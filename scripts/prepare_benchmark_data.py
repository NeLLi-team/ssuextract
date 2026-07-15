#!/usr/bin/env python3
"""Validate database benchmark artifacts and write documentation data tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Iterable, Mapping, Sequence


PROFILE_LABELS = {
    "curated-legacy": "Legacy 138.1/4.12",
    "curated": "Curated 138.2/5.1.1",
    "img": "IMG-enhanced 138.2/5.1.1",
}
BLAST_SUFFIXES = {
    ".ndb",
    ".nhr",
    ".nin",
    ".njs",
    ".not",
    ".nsq",
    ".ntf",
    ".nto",
}
FASTA_SUFFIXES = {".fa", ".fasta", ".fna"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", text=True
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            os.fchmod(handle.fileno(), 0o644)
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def load_timings(path: Path) -> list[dict[str, object]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"benchmark timing table is missing or empty: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"profile", "trial", "elapsed_seconds", "peak_rss_kb"}
        if not required <= set(reader.fieldnames or ()):
            raise ValueError("benchmark timing table has an invalid schema")
        observed: dict[tuple[str, int], dict[str, object]] = {}
        for row_number, row in enumerate(reader, 2):
            profile = row["profile"]
            try:
                trial = int(row["trial"])
                elapsed = float(row["elapsed_seconds"])
                peak_rss = int(row["peak_rss_kb"])
            except (TypeError, ValueError) as error:
                raise ValueError(f"invalid benchmark value on row {row_number}") from error
            if profile not in PROFILE_LABELS or trial not in range(4):
                raise ValueError(f"unexpected profile or trial on row {row_number}")
            if not math.isfinite(elapsed) or elapsed <= 0 or peak_rss <= 0:
                raise ValueError(f"invalid benchmark value on row {row_number}")
            key = (profile, trial)
            if key in observed:
                raise ValueError(f"duplicate benchmark profile/trial: {key}")
            observed[key] = {
                "profile": PROFILE_LABELS[profile],
                "trial": trial,
                "warmup": "true" if trial == 0 else "false",
                "elapsed_seconds": f"{elapsed:.2f}",
                "peak_rss_kb": peak_rss,
            }
    expected = {(profile, trial) for profile in PROFILE_LABELS for trial in range(4)}
    if set(observed) != expected:
        missing = sorted(expected - set(observed))
        raise ValueError(f"benchmark timing table is incomplete: {missing}")
    return [observed[key] for key in sorted(observed, key=lambda item: (list(PROFILE_LABELS).index(item[0]), item[1]))]


def _storage_components(profile: Path) -> dict[str, int]:
    if not profile.is_dir():
        raise ValueError(f"database profile directory is missing: {profile}")
    files = sorted(path for path in profile.rglob("*") if path.is_file())
    if not files:
        raise ValueError(f"database profile contains no runtime files: {profile}")
    components = {
        "source_fasta_bytes": 0,
        "blast_bytes": 0,
        "parquet_bytes": 0,
        "other_bytes": 0,
    }
    for path in files:
        relative = path.relative_to(profile)
        size = path.stat().st_size
        if path.suffix.lower() in FASTA_SUFFIXES:
            components["source_fasta_bytes"] += size
        elif path.suffix.lower() == ".parquet":
            components["parquet_bytes"] += size
        elif path.suffix.lower() in BLAST_SUFFIXES or relative.parts[0] == "blast":
            components["blast_bytes"] += size
        else:
            components["other_bytes"] += size
    components["installed_bytes"] = sum(components.values())
    return components


def storage_row(
    label: str,
    profile: Path,
    archive: Path | None,
) -> dict[str, object]:
    values: dict[str, object] = {"profile": label, **_storage_components(profile)}
    if archive is None:
        values.update({"archive_bytes": "", "archive_sha256": ""})
    else:
        if not archive.is_file() or archive.stat().st_size == 0:
            raise ValueError(f"database archive is missing or empty: {archive}")
        values.update({"archive_bytes": archive.stat().st_size, "archive_sha256": _sha256(archive)})
    return values


def _validate_profile_identity(profile: Path, expected_profile: str) -> None:
    manifest_path = profile / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid database manifest: {manifest_path}") from error
    if manifest.get("profile") != expected_profile:
        raise ValueError(f"database manifest profile mismatch: {manifest_path}")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-directory", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--curated-profile", type=Path, required=True)
    parser.add_argument("--curated-archive", type=Path, required=True)
    parser.add_argument("--img-profile", type=Path, required=True)
    parser.add_argument("--img-archive", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    timings = load_timings(args.benchmark_directory / "timing.tsv")
    _validate_profile_identity(args.curated_profile, "curated")
    _validate_profile_identity(args.img_profile, "img")
    storage = [
        storage_row(PROFILE_LABELS["curated-legacy"], args.legacy_root, None),
        storage_row(PROFILE_LABELS["curated"], args.curated_profile, args.curated_archive),
        storage_row(PROFILE_LABELS["img"], args.img_profile, args.img_archive),
    ]
    _atomic_tsv(
        args.output_directory / "database_profile_benchmark.tsv",
        ("profile", "trial", "warmup", "elapsed_seconds", "peak_rss_kb"),
        timings,
    )
    _atomic_tsv(
        args.output_directory / "database_storage.tsv",
        (
            "profile",
            "source_fasta_bytes",
            "blast_bytes",
            "parquet_bytes",
            "other_bytes",
            "installed_bytes",
            "archive_bytes",
            "archive_sha256",
        ),
        storage,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
