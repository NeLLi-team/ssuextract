#!/usr/bin/env python3
"""Assemble a validated SSUextract database profile and release archive."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping, Sequence

import build_database_release as builder
import database_manager as manager
from atomic_io import fsync_directory, fsync_file, replace_and_fsync


_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_TABLE_PATHS = {
    "sequences": "tables/sequences.parquet",
    "source_records": "tables/source_records.parquet",
    "taxonomy_assignments": "tables/taxonomy_assignments.parquet",
    "preferred": "tables/preferred_taxonomy.parquet",
    "img_location": "tables/img_location.parquet",
}


class AssemblyError(RuntimeError):
    """A profile or its archive could not be assembled safely."""


def _annotate_rollback_errors(
    publication_error: BaseException,
    rollback_errors: Sequence[BaseException],
    label: str,
) -> None:
    messages = tuple(
        f"{label}: {type(error).__name__}: {error}" for error in rollback_errors
    )
    setattr(publication_error, "rollback_errors", messages)
    add_note = getattr(publication_error, "add_note", None)
    if add_note is not None:
        for message in messages:
            add_note(message)


@dataclass(frozen=True)
class ArchiveMetadata:
    path: Path
    bytes: int
    sha256: str


@dataclass(frozen=True)
class AssemblyResult:
    profile_directory: Path
    manifest: dict[str, object]
    archive: ArchiveMetadata | None


def _require_identifier(value: str, label: str) -> str:
    if not isinstance(value, str) or not _IDENTIFIER.fullmatch(value):
        raise AssemblyError(f"{label} must be a safe, non-empty identifier")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _run(command: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            command,
            cwd=cwd,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise AssemblyError(f"Required executable not found: {command[0]}") from error
    except subprocess.CalledProcessError as error:
        detail = (error.stderr or error.stdout or "command returned an error").strip()
        raise AssemblyError(f"Command failed ({' '.join(command)}): {detail}") from error


def _makeblastdb_version(executable: str, *, cwd: Path) -> str:
    output = _run([executable, "-version"], cwd=cwd)
    lines = [line.strip() for line in (output.stdout + output.stderr).splitlines() if line.strip()]
    return " | ".join(lines)


def _marker_sequences(
    model: builder.DatabaseModel,
) -> tuple[dict[str, tuple[builder.SequenceRecord, ...]], dict[str, int]]:
    members: dict[str, list[builder.SequenceRecord]] = defaultdict(list)
    source_counts: Counter[str] = Counter()
    for sequence in model.sequences:
        for raw_marker in sequence.markers:
            marker = _require_identifier(raw_marker, "sequence marker")
            members[marker].append(sequence)
    for source in model.source_records:
        marker = _require_identifier(source.marker, "source marker")
        source_counts[marker] += 1
    if not members:
        raise AssemblyError("Database model has no marker memberships")
    marker_sequences = {
        marker: tuple(sorted(sequences, key=lambda record: record.sequence_id))
        for marker, sequences in sorted(members.items())
    }
    return marker_sequences, dict(source_counts)


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _artifact(path: Path, root: Path) -> dict[str, object]:
    relative = PurePosixPath(path.relative_to(root).as_posix())
    size = path.stat().st_size
    if size <= 0:
        raise AssemblyError(f"Generated artifact is empty: {relative}")
    return {"path": str(relative), "bytes": size, "sha256": _sha256(path)}


def _normalize_permissions(root: Path) -> None:
    root.chmod(0o755)
    for path in root.rglob("*"):
        if path.is_symlink():
            raise AssemblyError(f"Generated profile contains a symbolic link: {path}")
        path.chmod(0o755 if path.is_dir() else 0o644)


def _provenance(
    model: builder.DatabaseModel,
    profile: str,
    version: str,
    marker_sequences: dict[str, tuple[builder.SequenceRecord, ...]],
    source_counts: dict[str, int],
    makeblastdb_version: str,
    provenance_details: Mapping[str, object] | None,
) -> dict[str, object]:
    source_versions = Counter(
        (record.reference_source, record.source_version) for record in model.source_records
    )
    provenance = {
        "schema_version": 1,
        "profile": profile,
        "version": version,
        "counts": {
            "sequences": len(model.sequences),
            "source_records": len(model.source_records),
            "taxonomy_assignments": len(model.taxonomy_assignments),
            "preferred_taxonomy": len(model.preferred_taxonomy),
        },
        "markers": {
            marker: {
                "sequences": len(records),
                "source_records": source_counts[marker],
            }
            for marker, records in marker_sequences.items()
        },
        "sources": [
            {"name": source, "version": source_version, "source_records": count}
            for (source, source_version), count in sorted(source_versions.items())
        ],
        "tools": {
            "makeblastdb": {
                "version": makeblastdb_version,
                "command_template": [
                    "makeblastdb",
                    "-in",
                    "fasta/{marker}.fasta",
                    "-dbtype",
                    "nucl",
                    "-parse_seqids",
                    "-blastdb_version",
                    "5",
                    "-out",
                    "blast/{marker}",
                ],
            }
        },
    }
    if provenance_details:
        provenance["release_build"] = dict(provenance_details)
    return provenance


def _copy_release_files(staging: Path, files: Mapping[str, str | Path]) -> None:
    for relative_name, source_name in sorted(files.items()):
        relative = PurePosixPath(relative_name)
        if (
            relative.is_absolute()
            or not relative.parts
            or any(part in {"", ".", ".."} for part in relative.parts)
        ):
            raise AssemblyError(f"Unsafe supplemental release path: {relative_name!r}")
        source = Path(source_name)
        if not source.is_file() or source.is_symlink():
            raise AssemblyError(f"Supplemental release file is not a regular file: {source}")
        destination = staging / Path(relative)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)


def _assemble_staging_profile(
    staging: Path,
    model: builder.DatabaseModel,
    profile: str,
    version: str,
    img_locations: tuple[builder.ImgLocation, ...],
    makeblastdb: str,
    blastdbcmd: str,
    provenance_details: Mapping[str, object] | None,
    release_files: Mapping[str, str | Path],
) -> dict[str, object]:
    builder.validate_release(model, img_locations)
    marker_sequences, source_counts = _marker_sequences(model)

    builder.write_release_tables(staging / "tables", model, img_locations)
    blast_databases: dict[str, dict[str, str]] = {}
    for marker, sequences in marker_sequences.items():
        fasta = Path("fasta") / f"{marker}.fasta"
        prefix = Path("blast") / marker
        builder.write_marker_fasta(staging / fasta, sequences)
        (staging / prefix).parent.mkdir(parents=True, exist_ok=True)
        _run(
            [
                makeblastdb,
                "-in",
                fasta.as_posix(),
                "-dbtype",
                "nucl",
                "-parse_seqids",
                "-blastdb_version",
                "5",
                "-out",
                prefix.as_posix(),
                "-title",
                f"SSUextract {profile} {version} {marker}",
            ],
            cwd=staging,
        )
        blast_databases[marker] = {"prefix": prefix.as_posix()}
        (staging / fasta).unlink()
    fasta_directory = staging / "fasta"
    if fasta_directory.exists():
        fasta_directory.rmdir()

    _copy_release_files(staging, release_files)

    provenance_path = staging / "provenance.json"
    _write_json(
        provenance_path,
        _provenance(
            model,
            profile,
            version,
            marker_sequences,
            source_counts,
            _makeblastdb_version(makeblastdb, cwd=staging),
            provenance_details,
        ),
    )
    artifact_paths = sorted(
        path
        for path in staging.rglob("*")
        if path.is_file() and not path.is_symlink() and path.name != manager.MANIFEST_NAME
    )
    manifest: dict[str, object] = {
        "schema_version": manager.SCHEMA_VERSION,
        "profile": profile,
        "version": version,
        "artifacts": [_artifact(path, staging) for path in artifact_paths],
        "blast_databases": blast_databases,
        "taxonomy_database": dict(_TABLE_PATHS),
        "provenance": "provenance.json",
    }
    evidence_catalog = staging / "EVIDENCE" / "img_taxonomy_evidence.jsonl"
    if evidence_catalog.is_file():
        manifest["taxonomy_evidence_catalog"] = PurePosixPath(
            evidence_catalog.relative_to(staging).as_posix()
        ).as_posix()
    _write_json(staging / manager.MANIFEST_NAME, manifest)
    _normalize_permissions(staging)
    manager.validate_profile_directory(staging, profile, blastdbcmd)
    return manifest


def _write_deterministic_tar(profile_directory: Path, destination: Path) -> None:
    root_name = profile_directory.name
    entries = [
        profile_directory,
        *sorted(profile_directory.rglob("*"), key=lambda path: path.as_posix()),
    ]
    with tarfile.open(destination, mode="w", format=tarfile.PAX_FORMAT) as archive:
        for path in entries:
            if path.is_symlink() or not (path.is_dir() or path.is_file()):
                raise AssemblyError(f"Unsupported archive member: {path}")
            relative = path.relative_to(profile_directory)
            name = root_name if not relative.parts else f"{root_name}/{relative.as_posix()}"
            info = tarfile.TarInfo(name + ("/" if path.is_dir() else ""))
            info.uid = 0
            info.gid = 0
            info.uname = ""
            info.gname = ""
            info.mtime = 0
            if path.is_dir():
                info.type = tarfile.DIRTYPE
                info.mode = 0o755
                archive.addfile(info)
            else:
                info.size = path.stat().st_size
                info.mode = 0o644
                with path.open("rb") as handle:
                    archive.addfile(info, handle)


def package_profile(
    profile_directory: str | Path,
    archive_path: str | Path,
    *,
    zstd: str = "zstd",
    compression_level: int = 10,
) -> ArchiveMetadata:
    """Package exactly one profile directory into a deterministic ``tar.zst``."""

    profile = Path(profile_directory).resolve()
    if not profile.is_dir():
        raise AssemblyError(f"Profile directory does not exist: {profile}")
    destination = Path(archive_path).resolve()
    if destination == profile or profile in destination.parents:
        raise AssemblyError("Archive destination must be outside the profile directory")
    if destination.exists():
        raise AssemblyError(f"Archive destination already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{destination.name}.", dir=destination.parent
    ) as tmp:
        temporary_directory = Path(tmp)
        tar_path = temporary_directory / "profile.tar"
        compressed = temporary_directory / "profile.tar.zst"
        _write_deterministic_tar(profile, tar_path)
        _run(
            [
                zstd,
                f"-{compression_level}",
                "-T1",
                "-q",
                "-f",
                tar_path.name,
                "-o",
                compressed.name,
            ],
            cwd=temporary_directory,
        )
        linked = False
        try:
            # The temporary archive is created on the destination filesystem.
            # link(2) publishes it without the overwrite race of os.replace().
            fsync_file(compressed)
            os.link(compressed, destination)
            linked = True
            fsync_directory(destination.parent)
        except FileExistsError as error:
            raise AssemblyError(f"Archive destination already exists: {destination}") from error
        except BaseException as publication_error:
            rollback_errors: list[BaseException] = []
            if linked:
                try:
                    destination.unlink(missing_ok=True)
                except BaseException as error:
                    rollback_errors.append(error)
                try:
                    fsync_directory(destination.parent)
                except BaseException as error:
                    rollback_errors.append(error)
            _annotate_rollback_errors(
                publication_error,
                rollback_errors,
                "archive rollback also failed",
            )
            raise
    return ArchiveMetadata(destination, destination.stat().st_size, _sha256(destination))


def _stage_profile_archive(
    profile_directory: Path,
    archive_path: Path,
    *,
    zstd: str,
    compression_level: int = 10,
) -> tuple[Path, Path]:
    """Create but do not publish an archive on its destination filesystem."""

    archive_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_directory = Path(
        tempfile.mkdtemp(prefix=f".{archive_path.name}.staging-", dir=archive_path.parent)
    )
    tar_path = temporary_directory / "profile.tar"
    compressed = temporary_directory / archive_path.name
    try:
        _write_deterministic_tar(profile_directory, tar_path)
        _run(
            [
                zstd,
                f"-{compression_level}",
                "-T1",
                "-q",
                "-f",
                tar_path.name,
                "-o",
                compressed.name,
            ],
            cwd=temporary_directory,
        )
        return temporary_directory, compressed
    except BaseException:
        shutil.rmtree(temporary_directory, ignore_errors=True)
        raise


def assemble_database_profile(
    model: builder.DatabaseModel,
    profile_directory: str | Path,
    *,
    profile: str,
    version: str,
    img_locations: Iterable[builder.ImgLocation] = (),
    makeblastdb: str = "makeblastdb",
    blastdbcmd: str = "blastdbcmd",
    archive_path: str | Path | None = None,
    zstd: str = "zstd",
    provenance_details: Mapping[str, object] | None = None,
    release_files: Mapping[str, str | Path] | None = None,
) -> AssemblyResult:
    """Stage and validate both outputs, with rollback on reported publish failures."""

    profile = _require_identifier(profile, "profile")
    version = _require_identifier(version, "version")
    target = Path(profile_directory).resolve()
    if target.name != profile:
        raise AssemblyError(
            f"Profile directory name {target.name!r} must match profile {profile!r}"
        )
    if target.exists():
        raise AssemblyError(f"Profile directory already exists: {target}")
    archive_destination = Path(archive_path).resolve() if archive_path else None
    if archive_destination is not None:
        if archive_destination == target or target in archive_destination.parents:
            raise AssemblyError("Archive destination must be outside the profile directory")
        if archive_destination.exists():
            raise AssemblyError(
                f"Archive destination already exists: {archive_destination}"
            )
    target.parent.mkdir(parents=True, exist_ok=True)
    transaction = Path(
        tempfile.mkdtemp(prefix=f".{profile}.transaction-", dir=target.parent)
    )
    staging = transaction / profile
    staging.mkdir()
    archive_staging_directory: Path | None = None
    staged_archive: Path | None = None
    archive_bytes: int | None = None
    archive_sha256: str | None = None
    try:
        manifest = _assemble_staging_profile(
            staging,
            model,
            profile,
            version,
            tuple(img_locations),
            makeblastdb,
            blastdbcmd,
            provenance_details,
            release_files or {},
        )
        if archive_destination is not None:
            archive_staging_directory, staged_archive = _stage_profile_archive(
                staging, archive_destination, zstd=zstd
            )
            archive_bytes = staged_archive.stat().st_size
            archive_sha256 = _sha256(staged_archive)
            fsync_file(staged_archive)
        try:
            replace_and_fsync(staging, target)
        except BaseException as publication_error:
            rollback_errors: list[BaseException] = []
            if target.exists() and not staging.exists():
                try:
                    replace_and_fsync(target, staging)
                except BaseException as error:
                    rollback_errors.append(error)
            _annotate_rollback_errors(
                publication_error,
                rollback_errors,
                "profile publication rollback also failed",
            )
            raise
        if archive_destination is not None and staged_archive is not None:
            archive_linked = False
            try:
                # Publish without replacing an archive created by another process.
                os.link(staged_archive, archive_destination)
                archive_linked = True
                fsync_directory(archive_destination.parent)
            except BaseException as publication_error:
                # Roll back failures reported by the filesystem. A hard process
                # kill between the two operations can still leave a partial pair;
                # release publication therefore validates both paths separately.
                rollback_errors: list[BaseException] = []
                if archive_linked:
                    try:
                        archive_destination.unlink(missing_ok=True)
                    except BaseException as error:
                        rollback_errors.append(error)
                    try:
                        fsync_directory(archive_destination.parent)
                    except BaseException as error:
                        rollback_errors.append(error)
                try:
                    replace_and_fsync(target, staging)
                except BaseException as error:
                    rollback_errors.append(error)
                _annotate_rollback_errors(
                    publication_error,
                    rollback_errors,
                    "release rollback also failed",
                )
                raise
    except BaseException:
        shutil.rmtree(transaction, ignore_errors=True)
        if archive_staging_directory is not None:
            shutil.rmtree(archive_staging_directory, ignore_errors=True)
        raise
    shutil.rmtree(transaction, ignore_errors=True)
    if archive_staging_directory is not None:
        shutil.rmtree(archive_staging_directory, ignore_errors=True)
    archive = (
        ArchiveMetadata(
            archive_destination,
            archive_bytes,
            archive_sha256,
        )
        if archive_destination is not None
        and archive_bytes is not None
        and archive_sha256 is not None
        else None
    )
    return AssemblyResult(target, manifest, archive)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records-jsonl", type=Path, required=True)
    parser.add_argument("--img-metadata", type=Path)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--blastdbcmd", default="blastdbcmd")
    parser.add_argument("--zstd", default="zstd")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    model = builder.build_deduplicated_model(builder.read_prepared_jsonl(args.records_jsonl))
    locations = builder.read_img_metadata_tsv(args.img_metadata) if args.img_metadata else ()
    result = assemble_database_profile(
        model,
        args.output_root / args.profile,
        profile=args.profile,
        version=args.version,
        img_locations=locations,
        makeblastdb=args.makeblastdb,
        blastdbcmd=args.blastdbcmd,
        archive_path=args.archive,
        zstd=args.zstd,
    )
    summary = {
        "profile_directory": str(result.profile_directory),
        "archive": (
            {
                "path": str(result.archive.path),
                "bytes": result.archive.bytes,
                "sha256": result.archive.sha256,
            }
            if result.archive
            else None
        ),
    }
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
