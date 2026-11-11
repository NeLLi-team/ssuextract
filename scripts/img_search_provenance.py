#!/usr/bin/env python3
"""Validate the completion contract of a chunked IMG centroid BLAST search."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Mapping, Sequence

from atomic_io import replace_and_fsync


SCHEMA_VERSION = 2
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
    if observed["schema_version"] != SCHEMA_VERSION or observed["status"] != STATUS:
        raise RuntimeError("IMG search completion sidecar is not schema-2 complete")
    import img_chunked_search  # noqa: PLC0415  # Chunk helpers import this module.

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


def portable_provenance(
    payload: Mapping[str, object], sidecar_sha256: str
) -> dict[str, object]:
    import img_chunked_search  # noqa: PLC0415  # Chunk helpers import this module.

    return img_chunked_search.portable_provenance(payload, sidecar_sha256)


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
    child = subparsers.add_parser("validate")
    _add_common_arguments(child)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    payload = validate_search(
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
