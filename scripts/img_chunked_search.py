#!/usr/bin/env python3
"""Run deterministic, restartable IMG centroid BLAST query chunks."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Mapping, Sequence

from atomic_io import replace_and_fsync
import img_search_provenance as base


SCHEMA_VERSION = 2
PLAN_SCHEMA_VERSION = 1
RECEIPT_SCHEMA_VERSION = 1
STATUS = "complete"
MAX_CHUNKS = 64
PLAN_NAME = "plan.json"
QUERY_DIRECTORY = "queries"
RESULT_DIRECTORY = "results"
RECEIPT_DIRECTORY = "receipts"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_concatenation(paths: Sequence[Path]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def _file_record(path: Path, *, allow_empty: bool = False) -> dict[str, str]:
    resolved = path.resolve()
    if not resolved.is_file() or (not allow_empty and resolved.stat().st_size == 0):
        raise RuntimeError(f"chunked search file is missing or empty: {resolved}")
    return {"path": str(resolved), "sha256": _sha256(resolved)}


def _database_binding(curated_profile: Path, marker: str) -> dict[str, object]:
    profile = curated_profile.resolve()
    manifest_path = profile / base.MANIFEST_NAME
    manifest = _load_json(manifest_path, "curated profile manifest")
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise RuntimeError("curated profile manifest has no artifact inventory")
    prefix = f"blast/{marker}."
    selected: list[dict[str, object]] = []
    for record in artifacts:
        if not isinstance(record, dict) or set(record) != {"path", "bytes", "sha256"}:
            raise RuntimeError("invalid curated profile artifact record")
        relative = record.get("path")
        if not isinstance(relative, str) or not relative.startswith(prefix):
            continue
        path = (profile / relative).resolve()
        if not path.is_relative_to(profile) or not path.is_file():
            raise RuntimeError(f"invalid curated BLAST artifact path: {relative}")
        if (
            type(record.get("bytes")) is not int
            or path.stat().st_size != record["bytes"]
            or _sha256(path) != record.get("sha256")
        ):
            raise RuntimeError(f"curated BLAST artifact binding mismatch: {relative}")
        selected.append(dict(record))
    selected.sort(key=lambda item: str(item["path"]))
    suffixes = {Path(str(item["path"])).suffix for item in selected}
    if not {".nhr", ".nin", ".nsq"}.issubset(suffixes):
        raise RuntimeError(f"curated BLAST artifact inventory is incomplete: {marker}")
    return {
        "manifest": _file_record(manifest_path),
        "artifacts": selected,
    }


def _blastn_tool() -> dict[str, str]:
    executable_name = shutil.which("blastn")
    if executable_name is None:
        raise RuntimeError("required executable not found: blastn")
    executable = Path(executable_name).resolve()
    try:
        result = subprocess.run(
            [str(executable), "-version"],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as error:
        raise RuntimeError("could not record blastn version") from error
    version = "\n".join(
        part.strip() for part in (result.stdout, result.stderr) if part.strip()
    )
    if not version:
        raise RuntimeError("blastn returned an empty version string")
    return {"path": str(executable), "sha256": _sha256(executable), "version": version}


def _write_bytes_durably(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(value)
            handle.flush()
        replace_and_fsync(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _read_fasta_records(path: Path) -> list[tuple[str, bytes]]:
    records: list[tuple[str, bytes]] = []
    current = bytearray()
    current_id: str | None = None
    seen: set[str] = set()
    with path.open("rb") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(b">"):
                if current_id is not None:
                    records.append((current_id, bytes(current)))
                try:
                    current_id = line[1:].split(None, 1)[0].decode("ascii")
                except (IndexError, UnicodeDecodeError) as error:
                    raise ValueError(
                        f"invalid FASTA header on line {line_number}: {path}"
                    ) from error
                if not current_id or current_id in seen:
                    raise ValueError(f"empty or duplicate FASTA identifier: {current_id!r}")
                seen.add(current_id)
                current = bytearray(line)
            else:
                if current_id is None:
                    raise ValueError(f"FASTA sequence appears before a header: {path}")
                current.extend(line)
    if current_id is not None:
        records.append((current_id, bytes(current)))
    if not records:
        raise ValueError(f"FASTA contains no records: {path}")
    return records


def _validate_marker(marker: str) -> None:
    if marker not in {"16S", "18S"}:
        raise ValueError(f"unsupported marker: {marker}")


def _chunk_path(directory: Path, role: str, index: int) -> Path:
    suffix = "fna" if role == QUERY_DIRECTORY else "m8" if role == RESULT_DIRECTORY else "json"
    return directory / role / f"{index:03d}.{suffix}"


def prepare_chunks(
    marker: str,
    centroids_fasta: str | Path,
    chunk_directory: str | Path,
    chunk_count: int,
    *,
    force: bool = False,
) -> dict[str, object]:
    _validate_marker(marker)
    if not 1 <= chunk_count <= MAX_CHUNKS:
        raise ValueError(f"chunk_count must be from 1 through {MAX_CHUNKS}")
    centroids = Path(centroids_fasta).resolve()
    directory = Path(chunk_directory).resolve()
    plan_path = directory / PLAN_NAME
    if plan_path.exists() and not force:
        plan = validate_plan(marker, centroids, directory)
        if plan["chunk_count"] != chunk_count:
            raise RuntimeError(
                "requested chunk_count differs from the existing chunk plan"
            )
        return plan
    if directory.exists():
        if not force and any(directory.iterdir()):
            raise RuntimeError(f"chunk directory exists without a valid plan: {directory}")
        if force:
            shutil.rmtree(directory)
    directory.mkdir(parents=True, exist_ok=True)
    records = _read_fasta_records(centroids)
    if chunk_count > len(records):
        raise ValueError("chunk_count exceeds the number of centroid records")
    base_size, remainder = divmod(len(records), chunk_count)
    chunks: list[dict[str, object]] = []
    offset = 0
    for index in range(chunk_count):
        size = base_size + (1 if index < remainder else 0)
        selected = records[offset : offset + size]
        offset += size
        query_path = _chunk_path(directory, QUERY_DIRECTORY, index)
        _write_bytes_durably(query_path, b"".join(record for _, record in selected))
        chunks.append(
            {
                "index": index,
                "query_count": len(selected),
                "first_query_id": selected[0][0],
                "last_query_id": selected[-1][0],
                "query_fasta": _file_record(query_path),
            }
        )
    plan: dict[str, object] = {
        "schema_version": PLAN_SCHEMA_VERSION,
        "marker": marker,
        "centroids_fasta": {
            **_file_record(centroids),
            "record_count": len(records),
        },
        "chunk_count": chunk_count,
        "chunks": chunks,
    }
    base._write_json_durably(plan_path, plan)
    return validate_plan(marker, centroids, directory)


def _load_json(path: Path, label: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RuntimeError(f"invalid {label}: {path}") from error
    if not isinstance(value, dict):
        raise RuntimeError(f"invalid {label}: {path}")
    return value


def validate_plan(
    marker: str, centroids_fasta: str | Path, chunk_directory: str | Path
) -> dict[str, object]:
    _validate_marker(marker)
    centroids = Path(centroids_fasta).resolve()
    directory = Path(chunk_directory).resolve()
    plan = _load_json(directory / PLAN_NAME, "chunk plan")
    if set(plan) != {
        "schema_version",
        "marker",
        "centroids_fasta",
        "chunk_count",
        "chunks",
    }:
        raise RuntimeError("invalid chunk plan schema")
    if plan["schema_version"] != PLAN_SCHEMA_VERSION or plan["marker"] != marker:
        raise RuntimeError("chunk plan schema or marker mismatch")
    centroid_record = plan["centroids_fasta"]
    if not isinstance(centroid_record, dict) or set(centroid_record) != {
        "path",
        "sha256",
        "record_count",
    }:
        raise RuntimeError("invalid chunk plan centroid binding")
    expected_centroid = _file_record(centroids)
    if centroid_record.get("path") != expected_centroid["path"] or centroid_record.get(
        "sha256"
    ) != expected_centroid["sha256"]:
        raise RuntimeError("chunk plan centroid binding mismatch")
    chunk_count = plan["chunk_count"]
    chunks = plan["chunks"]
    if (
        type(chunk_count) is not int
        or not 1 <= chunk_count <= MAX_CHUNKS
        or not isinstance(chunks, list)
        or len(chunks) != chunk_count
    ):
        raise RuntimeError("invalid chunk plan count")
    paths: list[Path] = []
    total = 0
    for index, chunk in enumerate(chunks):
        if not isinstance(chunk, dict) or set(chunk) != {
            "index",
            "query_count",
            "first_query_id",
            "last_query_id",
            "query_fasta",
        }:
            raise RuntimeError("invalid chunk plan entry schema")
        if chunk["index"] != index or type(chunk["query_count"]) is not int or chunk[
            "query_count"
        ] < 1:
            raise RuntimeError("invalid chunk plan index or query count")
        query_path = _chunk_path(directory, QUERY_DIRECTORY, index).resolve()
        if chunk["query_fasta"] != _file_record(query_path):
            raise RuntimeError(f"chunk plan query binding mismatch: {index}")
        query_records = _read_fasta_records(query_path)
        if (
            len(query_records) != chunk["query_count"]
            or query_records[0][0] != chunk["first_query_id"]
            or query_records[-1][0] != chunk["last_query_id"]
        ):
            raise RuntimeError(f"chunk plan query content mismatch: {index}")
        paths.append(query_path)
        total += chunk["query_count"]
    if centroid_record.get("record_count") != total:
        raise RuntimeError("chunk plan record count mismatch")
    if _sha256_concatenation(paths) != expected_centroid["sha256"]:
        raise RuntimeError("chunk FASTA concatenation does not reproduce centroids")
    return plan


def _validate_chunk_receipt(
    marker: str,
    curated_profile: Path,
    directory: Path,
    plan: Mapping[str, object],
    index: int,
    threads: int,
    *,
    database_binding: Mapping[str, object] | None = None,
    blastn_tool: Mapping[str, str] | None = None,
) -> dict[str, object]:
    chunks = plan["chunks"]
    assert isinstance(chunks, list)
    chunk = chunks[index]
    assert isinstance(chunk, dict)
    query_path = _chunk_path(directory, QUERY_DIRECTORY, index).resolve()
    result_path = _chunk_path(directory, RESULT_DIRECTORY, index).resolve()
    receipt_path = _chunk_path(directory, RECEIPT_DIRECTORY, index).resolve()
    receipt = _load_json(receipt_path, "chunk receipt")
    database = dict(
        database_binding
        if database_binding is not None
        else _database_binding(curated_profile, marker)
    )
    tool = dict(blastn_tool if blastn_tool is not None else _blastn_tool())
    expected = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "status": STATUS,
        "marker": marker,
        "index": index,
        "query_count": chunk["query_count"],
        "database": database,
        "files": {
            "query_fasta": _file_record(query_path),
            "blast_m8": _file_record(result_path, allow_empty=True),
        },
        "blastn": {
            "argv": base.blastn_argv(
                marker, curated_profile, query_path, result_path, threads
            ),
            "max_target_seqs": base.MAX_TARGETS,
            "threads": threads,
            "tool": tool,
        },
    }
    if receipt != expected:
        raise RuntimeError(f"chunk receipt contract mismatch: {index}")
    return receipt


def run_chunk(
    marker: str,
    curated_profile: str | Path,
    centroids_fasta: str | Path,
    chunk_directory: str | Path,
    index: int,
    threads: int,
    *,
    force: bool = False,
) -> dict[str, object]:
    base._validate_marker_threads(marker, threads)
    profile = Path(curated_profile).resolve()
    directory = Path(chunk_directory).resolve()
    plan = validate_plan(marker, centroids_fasta, directory)
    if not 0 <= index < int(plan["chunk_count"]):
        raise ValueError("chunk index is outside the plan")
    result_path = _chunk_path(directory, RESULT_DIRECTORY, index).resolve()
    receipt_path = _chunk_path(directory, RECEIPT_DIRECTORY, index).resolve()
    if receipt_path.exists() and not force:
        return _validate_chunk_receipt(
            marker, profile, directory, plan, index, threads
        )
    receipt_path.unlink(missing_ok=True)
    result_path.unlink(missing_ok=True)
    query_path = _chunk_path(directory, QUERY_DIRECTORY, index).resolve()
    result_path.parent.mkdir(parents=True, exist_ok=True)
    database_before = _database_binding(profile, marker)
    tool = _blastn_tool()
    command = base.blastn_argv(marker, profile, query_path, result_path, threads)
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as error:
        raise RuntimeError("required executable not found: blastn") from error
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"blastn chunk {index} failed with exit code {error.returncode}"
        ) from error
    database_after = _database_binding(profile, marker)
    if database_after != database_before:
        raise RuntimeError(f"curated BLAST database changed during chunk {index}")
    chunks = plan["chunks"]
    assert isinstance(chunks, list) and isinstance(chunks[index], dict)
    receipt: dict[str, object] = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "status": STATUS,
        "marker": marker,
        "index": index,
        "query_count": chunks[index]["query_count"],
        "database": database_before,
        "files": {
            "query_fasta": _file_record(query_path),
            "blast_m8": _file_record(result_path, allow_empty=True),
        },
        "blastn": {
            "argv": command,
            "max_target_seqs": base.MAX_TARGETS,
            "threads": threads,
            "tool": tool,
        },
    }
    base._write_json_durably(receipt_path, receipt)
    return _validate_chunk_receipt(
        marker,
        profile,
        directory,
        plan,
        index,
        threads,
        database_binding=database_after,
        blastn_tool=tool,
    )


def _concatenate_results(paths: Sequence[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp"
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as destination:
            for path in paths:
                with path.open("rb") as source:
                    shutil.copyfileobj(source, destination, 1024 * 1024)
            destination.flush()
        replace_and_fsync(temporary, output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _chunked_payload(
    marker: str,
    curated_profile: str | Path,
    clusters: str | Path,
    source_fasta: str | Path,
    centroids_fasta: str | Path,
    blast_m8: str | Path,
    chunk_directory: str | Path,
    threads: int,
    *,
    validated_receipts: Sequence[Mapping[str, object]] | None = None,
) -> dict[str, object]:
    base._validate_marker_threads(marker, threads)
    profile = Path(curated_profile).resolve()
    directory = Path(chunk_directory).resolve()
    plan = validate_plan(marker, centroids_fasta, directory)
    chunk_count = int(plan["chunk_count"])
    if validated_receipts is None:
        database = _database_binding(profile, marker)
        tool = _blastn_tool()
        receipts = [
            _validate_chunk_receipt(
                marker,
                profile,
                directory,
                plan,
                index,
                threads,
                database_binding=database,
                blastn_tool=tool,
            )
            for index in range(chunk_count)
        ]
    else:
        receipts = [dict(receipt) for receipt in validated_receipts]
        if len(receipts) != chunk_count:
            raise RuntimeError("validated receipt count does not match the chunk plan")
    result_paths = [
        _chunk_path(directory, RESULT_DIRECTORY, index).resolve()
        for index in range(chunk_count)
    ]
    files = base._file_records(
        base._file_paths(
            profile, clusters, source_fasta, centroids_fasta, blast_m8
        )
    )
    if _sha256_concatenation(result_paths) != files["blast_m8"]["sha256"]:
        raise RuntimeError("combined BLAST output is not the ordered chunk concatenation")
    return {
        "schema_version": SCHEMA_VERSION,
        "status": STATUS,
        "marker": marker,
        "files": files,
        "blastn": {
            "mode": "query_chunks",
            "max_target_seqs": base.MAX_TARGETS,
            "threads": threads,
            "chunk_count": len(receipts),
            "chunks": receipts,
        },
    }


def finalize_search(
    marker: str,
    curated_profile: str | Path,
    clusters: str | Path,
    source_fasta: str | Path,
    centroids_fasta: str | Path,
    blast_m8: str | Path,
    sidecar: str | Path,
    chunk_directory: str | Path,
    threads: int,
) -> dict[str, object]:
    directory = Path(chunk_directory).resolve()
    sidecar_path = Path(sidecar).resolve()
    blast_path = Path(blast_m8).resolve()
    if directory != (sidecar_path.parent / "chunks").resolve():
        raise ValueError("chunk directory must be the sidecar directory's chunks path")
    sidecar_path.unlink(missing_ok=True)
    plan = validate_plan(marker, centroids_fasta, directory)
    profile = Path(curated_profile).resolve()
    database = _database_binding(profile, marker)
    tool = _blastn_tool()
    receipts = [
        _validate_chunk_receipt(
            marker,
            profile,
            directory,
            plan,
            index,
            threads,
            database_binding=database,
            blastn_tool=tool,
        )
        for index in range(int(plan["chunk_count"]))
    ]
    result_paths = [
        _chunk_path(directory, RESULT_DIRECTORY, index).resolve()
        for index in range(int(plan["chunk_count"]))
    ]
    _concatenate_results(result_paths, blast_path)
    payload = _chunked_payload(
        marker,
        curated_profile,
        clusters,
        source_fasta,
        centroids_fasta,
        blast_path,
        directory,
        threads,
        validated_receipts=receipts,
    )
    base._write_json_durably(sidecar_path, payload)
    return payload


def validate_search_payload(
    observed: Mapping[str, object],
    marker: str,
    curated_profile: str | Path,
    clusters: str | Path,
    source_fasta: str | Path,
    centroids_fasta: str | Path,
    blast_m8: str | Path,
    sidecar: str | Path,
    threads: int,
) -> dict[str, object]:
    directory = Path(sidecar).resolve().parent / "chunks"
    expected = _chunked_payload(
        marker,
        curated_profile,
        clusters,
        source_fasta,
        centroids_fasta,
        blast_m8,
        directory,
        threads,
    )
    if observed != expected:
        raise RuntimeError("chunked IMG search completion contract mismatch")
    return dict(observed)


def _portable_file(record: object, label: str) -> dict[str, str]:
    if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
        raise ValueError(f"invalid chunked search file record: {label}")
    path = record.get("path")
    if not isinstance(path, str) or not path:
        raise ValueError(f"invalid chunked search path: {label}")
    digest = record.get("sha256")
    if not isinstance(digest, str) or len(digest) != 64 or any(
        character not in "0123456789abcdef" for character in digest
    ):
        raise ValueError(f"invalid chunked search SHA256: {label}")
    return {"sha256": digest}


def portable_provenance(
    payload: Mapping[str, object], sidecar_sha256: str
) -> dict[str, object]:
    if len(sidecar_sha256) != 64 or any(
        character not in "0123456789abcdef" for character in sidecar_sha256
    ):
        raise ValueError("sidecar_sha256 must be a lowercase SHA256 digest")
    if (
        set(payload) != {"schema_version", "status", "marker", "files", "blastn"}
        or payload.get("schema_version") != SCHEMA_VERSION
        or payload.get("status") != STATUS
        or payload.get("marker") not in {"16S", "18S"}
    ):
        raise ValueError("invalid chunked search payload")
    files = payload.get("files")
    blastn = payload.get("blastn")
    if not isinstance(files, dict) or set(files) != set(base.FILE_KEYS):
        raise ValueError("invalid chunked search file binding")
    if not isinstance(blastn, dict) or set(blastn) != {
        "mode",
        "max_target_seqs",
        "threads",
        "chunk_count",
        "chunks",
    }:
        raise ValueError("invalid chunked BLAST contract")
    threads = blastn.get("threads")
    if type(threads) is not int:
        raise ValueError("invalid chunked BLAST thread contract")
    base._validate_marker_threads(str(payload["marker"]), threads)
    chunks = blastn.get("chunks")
    chunk_count = blastn.get("chunk_count")
    if (
        blastn.get("mode") != "query_chunks"
        or blastn.get("max_target_seqs") != base.MAX_TARGETS
        or type(chunk_count) is not int
        or not isinstance(chunks, list)
        or len(chunks) != chunk_count
    ):
        raise ValueError("invalid chunked BLAST execution contract")
    portable_chunks: list[dict[str, object]] = []
    for index, chunk in enumerate(chunks):
        if not isinstance(chunk, dict) or set(chunk) != {
            "schema_version",
            "status",
            "marker",
            "index",
            "query_count",
            "database",
            "files",
            "blastn",
        }:
            raise ValueError("invalid chunk receipt in search payload")
        chunk_files = chunk.get("files")
        database = chunk.get("database")
        command = chunk.get("blastn")
        if (
            chunk.get("schema_version") != RECEIPT_SCHEMA_VERSION
            or chunk.get("status") != STATUS
            or chunk.get("marker") != payload["marker"]
            or chunk.get("index") != index
            or type(chunk.get("query_count")) is not int
            or chunk["query_count"] < 1
            or not isinstance(chunk_files, dict)
            or set(chunk_files) != {"query_fasta", "blast_m8"}
            or not isinstance(database, dict)
            or set(database) != {"manifest", "artifacts"}
            or not isinstance(command, dict)
        ):
            raise ValueError("invalid chunk receipt contract")
        argv = command.get("argv")
        tool = command.get("tool")
        if (
            command.get("max_target_seqs") != base.MAX_TARGETS
            or command.get("threads") != threads
            or not isinstance(argv, list)
            or not isinstance(tool, dict)
            or set(tool) != {"path", "sha256", "version"}
        ):
            raise ValueError("invalid chunk BLAST command")
        portable_tool = {
            "sha256": _portable_file(
                {"path": tool.get("path"), "sha256": tool.get("sha256")},
                f"blastn tool {index}",
            )["sha256"],
            "version": tool.get("version"),
        }
        if not isinstance(portable_tool["version"], str) or not portable_tool[
            "version"
        ]:
            raise ValueError("invalid chunk BLAST version")
        normalized = list(argv)
        expected_argv = base.blastn_argv(
            str(payload["marker"]),
            ".",
            ".",
            ".",
            threads,
        )
        path_roles = {
            "-query": f"@chunk_query_{index:03d}",
            "-db": "@curated_blast_database",
            "-out": f"@chunk_output_{index:03d}",
        }
        for flag, role in path_roles.items():
            if normalized.count(flag) != 1 or expected_argv.count(flag) != 1:
                raise ValueError("invalid chunk BLAST path argument")
            normalized[normalized.index(flag) + 1] = role
            expected_argv[expected_argv.index(flag) + 1] = role
        if normalized != expected_argv:
            raise ValueError("invalid chunk BLAST argv")
        manifest = database.get("manifest")
        artifacts = database.get("artifacts")
        if not isinstance(artifacts, list) or not artifacts:
            raise ValueError("invalid chunk database artifact binding")
        portable_artifacts: list[dict[str, object]] = []
        artifact_paths: list[str] = []
        for artifact in artifacts:
            if not isinstance(artifact, dict) or set(artifact) != {
                "path",
                "bytes",
                "sha256",
            }:
                raise ValueError("invalid chunk database artifact record")
            artifact_path = str(artifact.get("path", ""))
            if not artifact_path.startswith(
                f"blast/{payload['marker']}."
            ):
                raise ValueError("wrong marker in chunk database artifact record")
            _portable_file(
                {"path": artifact["path"], "sha256": artifact["sha256"]},
                f"database artifact {index}",
            )
            if type(artifact.get("bytes")) is not int or artifact["bytes"] < 1:
                raise ValueError("invalid chunk database artifact size")
            artifact_paths.append(artifact_path)
            portable_artifacts.append(
                {
                    "name": artifact_path,
                    "bytes": artifact["bytes"],
                    "sha256": artifact["sha256"],
                }
            )
        if artifact_paths != sorted(set(artifact_paths)):
            raise ValueError("chunk database artifact records are not unique and sorted")
        if not {".nhr", ".nin", ".nsq"}.issubset(
            {Path(path).suffix for path in artifact_paths}
        ):
            raise ValueError("chunk database artifact records are incomplete")
        portable_database = {
            "manifest": _portable_file(manifest, f"database manifest {index}"),
            "artifacts": portable_artifacts,
        }
        portable_chunks.append(
            {
                "index": index,
                "query_count": chunk["query_count"],
                "database": portable_database,
                "files": {
                    "query_fasta": _portable_file(
                        chunk_files["query_fasta"], f"query_fasta {index}"
                    ),
                    "blast_m8": _portable_file(
                        chunk_files["blast_m8"], f"blast_m8 {index}"
                    ),
                },
                "blastn": {
                    "argv": normalized,
                    "max_target_seqs": base.MAX_TARGETS,
                    "threads": threads,
                    "tool": portable_tool,
                },
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "status": STATUS,
        "marker": payload["marker"],
        "sidecar_sha256": sidecar_sha256,
        "files": {
            name: _portable_file(files[name], name) for name in base.FILE_KEYS
        },
        "blastn": {
            "mode": "query_chunks",
            "max_target_seqs": base.MAX_TARGETS,
            "threads": threads,
            "chunk_count": chunk_count,
            "chunks": portable_chunks,
        },
    }


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--marker", required=True)
    parser.add_argument("--centroids-fasta", type=Path, required=True)
    parser.add_argument("--chunk-directory", type=Path, required=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    children = parser.add_subparsers(dest="action", required=True)
    prepare = children.add_parser("prepare")
    _add_common(prepare)
    prepare.add_argument("--chunks", type=int, required=True)
    prepare.add_argument("--force", action="store_true")
    run = children.add_parser("run-chunk")
    _add_common(run)
    run.add_argument("--curated-profile", type=Path, required=True)
    run.add_argument("--chunk-index", type=int, required=True)
    run.add_argument("--threads", type=int, required=True)
    run.add_argument("--force", action="store_true")
    finalize = children.add_parser("finalize")
    _add_common(finalize)
    finalize.add_argument("--curated-profile", type=Path, required=True)
    finalize.add_argument("--clusters", type=Path, required=True)
    finalize.add_argument("--source-fasta", type=Path, required=True)
    finalize.add_argument("--blast-m8", type=Path, required=True)
    finalize.add_argument("--sidecar", type=Path, required=True)
    finalize.add_argument("--threads", type=int, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.action == "prepare":
        payload = prepare_chunks(
            args.marker,
            args.centroids_fasta,
            args.chunk_directory,
            args.chunks,
            force=args.force,
        )
    elif args.action == "run-chunk":
        payload = run_chunk(
            args.marker,
            args.curated_profile,
            args.centroids_fasta,
            args.chunk_directory,
            args.chunk_index,
            args.threads,
            force=args.force,
        )
    else:
        payload = finalize_search(
            args.marker,
            args.curated_profile,
            args.clusters,
            args.source_fasta,
            args.centroids_fasta,
            args.blast_m8,
            args.sidecar,
            args.chunk_directory,
            args.threads,
        )
    print(
        json.dumps(
            {
                "action": args.action,
                "marker": args.marker,
                "schema_version": payload["schema_version"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
