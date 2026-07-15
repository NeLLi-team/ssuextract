#!/usr/bin/env python3
"""Run or validate an IMG centroid BLAST search completion contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Mapping, Sequence

from atomic_io import replace_and_fsync


SCHEMA_VERSION = 1
STATUS = "complete"
MAX_TARGETS = 501
MANIFEST_NAME = "manifest.json"
TAXONOMY_RELATIVE_PATH = Path("tables/preferred_taxonomy.parquet")
FILE_KEYS = (
    "curated_manifest",
    "preferred_taxonomy",
    "clusters",
    "source_fasta",
    "centroids_fasta",
    "blast_m8",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_marker_threads(marker: str, threads: int) -> None:
    if marker not in {"16S", "18S"}:
        raise ValueError(f"unsupported marker: {marker}")
    if not 1 <= threads <= 8:
        raise ValueError("threads must be an integer from 1 through 8")


def blastn_argv(
    marker: str,
    curated_profile: str | Path,
    centroids_fasta: str | Path,
    blast_m8: str | Path,
    threads: int,
) -> list[str]:
    _validate_marker_threads(marker, threads)
    profile = Path(curated_profile).resolve()
    centroids = Path(centroids_fasta).resolve()
    blast = Path(blast_m8).resolve()
    return [
        "blastn",
        "-task",
        "blastn",
        "-query",
        str(centroids),
        "-db",
        str(profile / "blast" / marker),
        "-evalue",
        "1e-5",
        "-max_target_seqs",
        str(MAX_TARGETS),
        "-max_hsps",
        "1",
        "-num_threads",
        str(threads),
        "-outfmt",
        "6 qseqid sseqid pident length qlen slen qcovs bitscore",
        "-out",
        str(blast),
    ]


def _file_paths(
    curated_profile: str | Path,
    clusters: str | Path,
    source_fasta: str | Path,
    centroids_fasta: str | Path,
    blast_m8: str | Path,
) -> dict[str, Path]:
    profile = Path(curated_profile).resolve()
    return {
        "curated_manifest": profile / MANIFEST_NAME,
        "preferred_taxonomy": profile / TAXONOMY_RELATIVE_PATH,
        "clusters": Path(clusters).resolve(),
        "source_fasta": Path(source_fasta).resolve(),
        "centroids_fasta": Path(centroids_fasta).resolve(),
        "blast_m8": Path(blast_m8).resolve(),
    }


def _file_records(paths: Mapping[str, Path]) -> dict[str, dict[str, str]]:
    missing = [
        f"{name}={path}"
        for name, path in paths.items()
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing:
        raise RuntimeError(
            "search contract file is missing or empty: " + ", ".join(missing)
        )
    return {
        name: {"path": str(path), "sha256": _sha256(path)}
        for name, path in paths.items()
    }


def _payload(
    marker: str,
    curated_profile: str | Path,
    clusters: str | Path,
    source_fasta: str | Path,
    centroids_fasta: str | Path,
    blast_m8: str | Path,
    threads: int,
) -> dict[str, object]:
    _validate_marker_threads(marker, threads)
    paths = _file_paths(
        curated_profile, clusters, source_fasta, centroids_fasta, blast_m8
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "status": STATUS,
        "marker": marker,
        "files": _file_records(paths),
        "blastn": {
            "argv": blastn_argv(
                marker, curated_profile, centroids_fasta, blast_m8, threads
            ),
            "max_target_seqs": MAX_TARGETS,
            "threads": threads,
        },
    }


def _write_json_durably(path: Path, value: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
        replace_and_fsync(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def run_search(
    marker: str,
    curated_profile: str | Path,
    clusters: str | Path,
    source_fasta: str | Path,
    centroids_fasta: str | Path,
    blast_m8: str | Path,
    sidecar: str | Path,
    threads: int,
) -> dict[str, object]:
    _validate_marker_threads(marker, threads)
    paths = _file_paths(
        curated_profile, clusters, source_fasta, centroids_fasta, blast_m8
    )
    _file_records({name: path for name, path in paths.items() if name != "blast_m8"})
    sidecar_path = Path(sidecar)
    blast_path = Path(blast_m8)
    sidecar_path.unlink(missing_ok=True)
    blast_path.parent.mkdir(parents=True, exist_ok=True)
    blast_path.unlink(missing_ok=True)
    command = blastn_argv(
        marker, curated_profile, centroids_fasta, blast_m8, threads
    )
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as error:
        raise RuntimeError("required executable not found: blastn") from error
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"blastn failed with exit code {error.returncode}"
        ) from error
    payload = _payload(
        marker,
        curated_profile,
        clusters,
        source_fasta,
        centroids_fasta,
        blast_m8,
        threads,
    )
    _write_json_durably(sidecar_path, payload)
    return payload


def validate_search(
    marker: str,
    curated_profile: str | Path,
    clusters: str | Path,
    source_fasta: str | Path,
    centroids_fasta: str | Path,
    blast_m8: str | Path,
    sidecar: str | Path,
    threads: int,
) -> dict[str, object]:
    sidecar_path = Path(sidecar)
    if not sidecar_path.is_file() or sidecar_path.stat().st_size == 0:
        raise RuntimeError(f"IMG search completion sidecar is missing or empty: {sidecar_path}")
    try:
        observed = json.loads(sidecar_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RuntimeError(f"invalid IMG search completion sidecar: {sidecar_path}") from error
    if not isinstance(observed, dict) or set(observed) != {
        "schema_version",
        "status",
        "marker",
        "files",
        "blastn",
    }:
        raise RuntimeError("invalid IMG search completion sidecar schema")
    if (
        observed["schema_version"] == 2
        and isinstance(observed.get("blastn"), dict)
        and observed["blastn"].get("mode") == "query_chunks"
    ):
        import img_chunked_search

        return img_chunked_search.validate_search_payload(
            observed,
            marker,
            curated_profile,
            clusters,
            source_fasta,
            centroids_fasta,
            blast_m8,
            sidecar_path,
            threads,
        )
    if observed["schema_version"] != SCHEMA_VERSION or observed["status"] != STATUS:
        raise RuntimeError("IMG search completion sidecar is not schema-1 complete")
    if observed["marker"] != marker:
        raise RuntimeError("IMG search completion sidecar marker mismatch")

    expected = _payload(
        marker,
        curated_profile,
        clusters,
        source_fasta,
        centroids_fasta,
        blast_m8,
        threads,
    )
    observed_files = observed["files"]
    expected_files = expected["files"]
    if not isinstance(observed_files, dict) or set(observed_files) != set(FILE_KEYS):
        raise RuntimeError("invalid IMG search completion file binding")
    for name in FILE_KEYS:
        if observed_files.get(name) != expected_files[name]:
            raise RuntimeError(f"IMG search completion file binding mismatch: {name}")
    if observed["blastn"] != expected["blastn"]:
        raise RuntimeError("IMG search completion BLAST command contract mismatch")
    return observed


def portable_provenance(
    payload: Mapping[str, object], sidecar_sha256: str
) -> dict[str, object]:
    if payload.get("schema_version") == 2:
        import img_chunked_search

        return img_chunked_search.portable_provenance(payload, sidecar_sha256)
    if len(sidecar_sha256) != 64 or any(
        character not in "0123456789abcdef" for character in sidecar_sha256
    ):
        raise ValueError("sidecar_sha256 must be a lowercase SHA256 digest")
    if (
        payload.get("schema_version") != SCHEMA_VERSION
        or payload.get("status") != STATUS
        or payload.get("marker") not in {"16S", "18S"}
    ):
        raise ValueError("cannot make portable provenance from an invalid search payload")
    files = payload.get("files")
    blastn = payload.get("blastn")
    if not isinstance(files, dict) or set(files) != set(FILE_KEYS):
        raise ValueError("search payload has an invalid file binding")
    if not isinstance(blastn, dict):
        raise ValueError("search payload has an invalid BLAST contract")
    hashes: dict[str, dict[str, str]] = {}
    for name in FILE_KEYS:
        record = files.get(name)
        if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
            raise ValueError(f"search payload has an invalid file record: {name}")
        sha256 = record.get("sha256")
        if (
            not isinstance(sha256, str)
            or len(sha256) != 64
            or any(character not in "0123456789abcdef" for character in sha256)
        ):
            raise ValueError(f"search payload has an invalid file SHA256: {name}")
        hashes[name] = {"sha256": sha256}
    threads = blastn.get("threads")
    max_targets = blastn.get("max_target_seqs")
    if type(threads) is not int:
        raise ValueError("search payload has an invalid BLAST thread contract")
    _validate_marker_threads(str(payload["marker"]), threads)
    if max_targets != MAX_TARGETS:
        raise ValueError("search payload has an invalid max_target_seqs contract")
    observed_argv = blastn.get("argv")
    if not isinstance(observed_argv, list):
        raise ValueError("search payload has an invalid BLAST argv contract")
    normalized_argv = list(observed_argv)
    portable_argv = blastn_argv(str(payload["marker"]), ".", ".", ".", threads)
    path_roles = {
        "-query": "@centroids_fasta",
        "-db": "@curated_blast_database",
        "-out": "@blast_m8",
    }
    for flag, role in path_roles.items():
        if normalized_argv.count(flag) != 1 or portable_argv.count(flag) != 1:
            raise ValueError("search payload has an invalid BLAST argv contract")
        normalized_argv[normalized_argv.index(flag) + 1] = role
        portable_argv[portable_argv.index(flag) + 1] = role
    if normalized_argv != portable_argv:
        raise ValueError("search payload has an invalid BLAST argv contract")
    return {
        "schema_version": SCHEMA_VERSION,
        "status": STATUS,
        "marker": payload["marker"],
        "sidecar_sha256": sidecar_sha256,
        "files": hashes,
        "blastn": {
            "argv": portable_argv,
            "max_target_seqs": MAX_TARGETS,
            "threads": threads,
        },
    }


def _add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--marker", required=True)
    parser.add_argument("--curated-profile", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--source-fasta", type=Path, required=True)
    parser.add_argument("--centroids-fasta", type=Path, required=True)
    parser.add_argument("--blast-m8", type=Path, required=True)
    parser.add_argument("--sidecar", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    for action in ("run", "validate"):
        child = subparsers.add_parser(action)
        _add_common_arguments(child)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    function = run_search if args.action == "run" else validate_search
    payload = function(
        args.marker,
        args.curated_profile,
        args.clusters,
        args.source_fasta,
        args.centroids_fasta,
        args.blast_m8,
        args.sidecar,
        args.threads,
    )
    print(
        json.dumps(
            {
                "marker": payload["marker"],
                "schema_version": payload["schema_version"],
                "status": payload["status"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
