#!/usr/bin/env python3

from __future__ import annotations

import argparse
import contextlib
import fcntl
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
import uuid
import warnings
from pathlib import Path, PurePosixPath
from typing import BinaryIO, Callable
from urllib.parse import urlparse

from atomic_io import fsync_directory, replace_and_fsync
from database_download import (
    DownloadError,
    DownloadIntegrityError,
    download_verified_archive,
)
from database_updates import (
    ReleaseDiscoveryError,
    compare_versions,
    discover_latest_catalog,
    validate_zenodo_config,
)


REPO = Path(__file__).resolve().parents[1]
DEFAULT_CATALOG = REPO / "config" / "database_catalog.json"
DEFAULT_MARKERS = REPO / "config" / "model_markers.json"
DEFAULT_ROOT = REPO / "resources" / "database"
MANIFEST_NAME = "manifest.json"
LEGACY_PREFIX = "silva-138-1_pr2-4-12"
SCHEMA_VERSION = 1

_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")

ProgressReporter = Callable[[str], None] | None


class DatabaseError(RuntimeError):
    """Base class for database catalog, profile, and installation failures."""


class CatalogError(DatabaseError):
    """The remote database catalog is malformed or incomplete."""


class ManifestError(DatabaseError):
    """A local profile manifest is malformed."""


class IntegrityError(DatabaseError):
    """An artifact does not match its recorded size or checksum."""


class InstallError(DatabaseError):
    """A database profile could not be installed safely."""


class LegacyDatabaseWarning(FutureWarning):
    """A pre-profile SSUextract database was selected."""


def _read_json(path: str | Path, label: str) -> object:
    source = Path(path)
    try:
        with source.open(encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise DatabaseError(f"Could not read {label} {source}: {error}") from error


def _require_object(value: object, label: str, error_type: type[DatabaseError]) -> dict:
    if not isinstance(value, dict):
        raise error_type(f"{label} must be a JSON object")
    return value


def _require_identifier(value: object, label: str, error_type: type[DatabaseError]) -> str:
    if not isinstance(value, str) or not _IDENTIFIER.fullmatch(value):
        raise error_type(f"{label} must be a safe, non-empty identifier")
    return value


def _safe_relative_path(
    value: object, label: str, error_type: type[DatabaseError]
) -> PurePosixPath:
    if not isinstance(value, str) or not value or "\\" in value:
        raise error_type(f"{label} must be a non-empty POSIX relative path")
    path = PurePosixPath(value)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise error_type(f"{label} escapes the profile directory: {value!r}")
    return path


def _require_size(value: object, label: str, error_type: type[DatabaseError]) -> int:
    if type(value) is not int or value <= 0:
        raise error_type(f"{label} must be a positive integer")
    return value


def _require_sha256(value: object, label: str, error_type: type[DatabaseError]) -> str:
    if not isinstance(value, str) or not _SHA256.fullmatch(value):
        raise error_type(f"{label} must be a lowercase SHA-256 digest")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _report(progress: ProgressReporter, message: str) -> None:
    if progress is not None:
        progress(message)


def _human_size(size: int) -> str:
    value = float(size)
    units = ("B", "KiB", "MiB", "GiB")
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024
    raise AssertionError("unreachable")


def load_catalog(path: str | Path = DEFAULT_CATALOG) -> dict:
    catalog = _require_object(_read_json(path, "database catalog"), "catalog", CatalogError)
    if catalog.get("schema_version") != SCHEMA_VERSION:
        raise CatalogError(f"catalog schema_version must be {SCHEMA_VERSION}")

    default_profile = _require_identifier(
        catalog.get("default_profile"), "catalog default_profile", CatalogError
    )
    profiles = _require_object(catalog.get("profiles"), "catalog profiles", CatalogError)
    if not profiles:
        raise CatalogError("catalog profiles must not be empty")
    if default_profile not in profiles:
        raise CatalogError(f"catalog default_profile {default_profile!r} is not listed")

    if "zenodo" in catalog:
        try:
            validate_zenodo_config(catalog["zenodo"])
        except ReleaseDiscoveryError as error:
            raise CatalogError(str(error)) from error

    for profile, raw_entry in profiles.items():
        _require_identifier(profile, "catalog profile name", CatalogError)
        entry = _require_object(raw_entry, f"catalog profile {profile!r}", CatalogError)
        _require_identifier(entry.get("version"), f"catalog profile {profile!r} version", CatalogError)
        description = entry.get("description")
        if description is not None and (
            not isinstance(description, str)
            or not description.strip()
            or any(character in description for character in "\r\n\t")
        ):
            raise CatalogError(
                f"catalog profile {profile!r} description must be non-empty, single-line text"
            )
        archive = _require_object(
            entry.get("archive"), f"catalog profile {profile!r} archive", CatalogError
        )
        values = (archive.get("url"), archive.get("bytes"), archive.get("sha256"))
        if values == (None, None, None):
            continue
        url, size, digest = values
        if not isinstance(url, str) or urlparse(url).scheme != "https" or not urlparse(url).netloc:
            raise CatalogError(f"catalog profile {profile!r} archive URL must use HTTPS")
        _require_size(size, f"catalog profile {profile!r} archive bytes", CatalogError)
        _require_sha256(digest, f"catalog profile {profile!r} archive sha256", CatalogError)
    return catalog


def _parse_marker_file(path: str | Path) -> dict[str, str]:
    source = Path(path)
    data = _require_object(_read_json(source, "model marker map"), "marker map", ManifestError)
    if "models" in data:
        if data.get("schema_version") != SCHEMA_VERSION:
            raise ManifestError(f"marker map schema_version must be {SCHEMA_VERSION}: {source}")
        raw_models = data.get("models")
    else:
        raw_models = data
    models = _require_object(raw_models, "marker map models", ManifestError)
    parsed: dict[str, str] = {}
    for model, marker in models.items():
        parsed[_require_identifier(model, "model name", ManifestError)] = _require_identifier(
            marker, f"marker for model {model!r}", ManifestError
        )
    return parsed


def load_marker_mapping(
    default_path: str | Path = DEFAULT_MARKERS,
    override_path: str | Path | None = None,
) -> dict[str, str]:
    models = _parse_marker_file(default_path)
    if override_path is not None:
        models.update(_parse_marker_file(override_path))
    return models


def load_manifest(profile_directory: str | Path, expected_profile: str | None = None) -> dict:
    directory = Path(profile_directory)
    manifest = _require_object(
        _read_json(directory / MANIFEST_NAME, "profile manifest"),
        "profile manifest",
        ManifestError,
    )
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ManifestError(f"manifest schema_version must be {SCHEMA_VERSION}")
    profile = _require_identifier(manifest.get("profile"), "manifest profile", ManifestError)
    if expected_profile is not None and profile != expected_profile:
        raise ManifestError(
            f"manifest profile {profile!r} does not match requested profile {expected_profile!r}"
        )
    _require_identifier(manifest.get("version"), "manifest version", ManifestError)

    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise ManifestError("manifest artifacts must be a non-empty array")
    artifact_paths: set[PurePosixPath] = set()
    for index, raw_artifact in enumerate(artifacts):
        artifact = _require_object(raw_artifact, f"artifact {index}", ManifestError)
        path = _safe_relative_path(artifact.get("path"), f"artifact {index} path", ManifestError)
        if path in artifact_paths:
            raise ManifestError(f"duplicate artifact path: {path}")
        artifact_paths.add(path)
        _require_size(artifact.get("bytes"), f"artifact {path} bytes", ManifestError)
        _require_sha256(artifact.get("sha256"), f"artifact {path} sha256", ManifestError)

    databases = _require_object(
        manifest.get("blast_databases"), "manifest blast_databases", ManifestError
    )
    if not databases:
        raise ManifestError("manifest blast_databases must not be empty")
    for marker, raw_database in databases.items():
        _require_identifier(marker, "BLAST marker", ManifestError)
        database = _require_object(raw_database, f"BLAST database {marker!r}", ManifestError)
        prefix = _safe_relative_path(
            database.get("prefix"), f"BLAST database {marker!r} prefix", ManifestError
        )
        if not any(
            path.parent == prefix.parent
            and (path.name == prefix.name or path.name.startswith(prefix.name + "."))
            for path in artifact_paths
        ):
            raise ManifestError(
                f"BLAST database {marker!r} prefix {prefix} has no manifest-listed artifacts"
            )

    taxonomy_database = _require_object(
        manifest.get("taxonomy_database"), "manifest taxonomy_database", ManifestError
    )
    preferred_taxonomy = _safe_relative_path(
        taxonomy_database.get("preferred"),
        "preferred taxonomy database",
        ManifestError,
    )
    if preferred_taxonomy not in artifact_paths:
        raise ManifestError(
            f"Preferred taxonomy database is not a manifest-listed artifact: "
            f"{preferred_taxonomy}"
        )
    return manifest


def _inside(root: Path, candidate: Path) -> bool:
    try:
        candidate.relative_to(root)
    except ValueError:
        return False
    return True


def _validate_blast_prefix(prefix: Path, blastdbcmd: str) -> None:
    try:
        subprocess.run(
            [blastdbcmd, "-db", str(prefix), "-info"],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise IntegrityError(f"BLAST validator not found: {blastdbcmd}") from error
    except subprocess.CalledProcessError as error:
        detail = error.stderr.strip() if error.stderr else "blastdbcmd returned an error"
        raise IntegrityError(f"Invalid BLAST database prefix {prefix}: {detail}") from error


def validate_profile_directory(
    profile_directory: str | Path,
    expected_profile: str | None = None,
    blastdbcmd: str = "blastdbcmd",
) -> dict:
    directory = Path(profile_directory).resolve()
    if not directory.is_dir():
        raise IntegrityError(f"Database profile directory does not exist: {directory}")
    manifest = load_manifest(directory, expected_profile=expected_profile)

    for artifact in manifest["artifacts"]:
        relative = _safe_relative_path(artifact["path"], "artifact path", ManifestError)
        path = (directory / Path(relative)).resolve()
        if not _inside(directory, path):
            raise IntegrityError(f"Artifact escapes profile directory: {relative}")
        if not path.is_file() or path.is_symlink():
            raise IntegrityError(f"Manifest artifact is not a regular file: {relative}")
        actual_size = path.stat().st_size
        if actual_size == 0:
            raise IntegrityError(f"Manifest artifact is empty: {relative}")
        if actual_size != artifact["bytes"]:
            raise IntegrityError(
                f"Artifact size mismatch for {relative}: expected {artifact['bytes']}, found {actual_size}"
            )
        actual_digest = _sha256(path)
        if actual_digest != artifact["sha256"]:
            raise IntegrityError(
                f"Artifact SHA-256 mismatch for {relative}: expected {artifact['sha256']}, "
                f"found {actual_digest}"
            )

    for database in manifest["blast_databases"].values():
        relative = _safe_relative_path(database["prefix"], "BLAST prefix", ManifestError)
        prefix = (directory / Path(relative)).resolve()
        if not _inside(directory, prefix):
            raise IntegrityError(f"BLAST prefix escapes profile directory: {relative}")
        _validate_blast_prefix(prefix, blastdbcmd)
    return manifest


def validate_profile(
    root: str | Path, profile: str, blastdbcmd: str = "blastdbcmd"
) -> dict:
    _require_identifier(profile, "profile", ManifestError)
    return validate_profile_directory(Path(root) / profile, profile, blastdbcmd)


def installed_version(root: str | Path, profile: str) -> str:
    _require_identifier(profile, "profile", ManifestError)
    manifest = load_manifest(Path(root) / profile, expected_profile=profile)
    return manifest["version"]


def check_profile_update(
    root: str | Path,
    profile: str,
    catalog_path: str | Path = DEFAULT_CATALOG,
    timeout: float = 5.0,
    opener: Callable[..., BinaryIO] = urllib.request.urlopen,
) -> dict[str, str]:
    installed = installed_version(root, profile)
    catalog = load_catalog(catalog_path)
    try:
        latest_catalog = discover_latest_catalog(catalog, timeout=timeout, opener=opener)
        latest = latest_catalog["profiles"][profile]["version"]
        comparison = compare_versions(installed, latest)
    except (KeyError, ReleaseDiscoveryError) as error:
        return {
            "status": "unavailable",
            "installed": installed,
            "latest": "",
            "reason": " ".join(str(error).split()),
        }
    if comparison < 0:
        status = "update_available"
    elif comparison > 0:
        status = "installed_newer"
    else:
        status = "current"
    return {
        "status": status,
        "installed": installed,
        "latest": latest,
        "reason": "",
    }


def detect_legacy_database(
    root: str | Path, blastdbcmd: str = "blastdbcmd"
) -> Path | None:
    prefix = Path(root).resolve() / LEGACY_PREFIX
    try:
        _validate_blast_prefix(prefix, blastdbcmd)
    except IntegrityError:
        return None
    return prefix


def resolve_database(
    root: str | Path,
    profile: str,
    model: str,
    marker_override: str | Path | None = None,
    marker_file: str | Path = DEFAULT_MARKERS,
    blastdbcmd: str = "blastdbcmd",
) -> Path:
    _require_identifier(profile, "profile", ManifestError)
    models = load_marker_mapping(marker_file, marker_override)
    try:
        marker = models[model]
    except KeyError as error:
        raise ManifestError(
            f"Unknown covariance model {model!r}; add it to a model-marker override file"
        ) from error

    directory = Path(root).resolve() / profile
    if directory.exists():
        manifest = validate_profile_directory(directory, profile, blastdbcmd)
        databases = manifest["blast_databases"]
        if marker not in databases:
            raise ManifestError(
                f"Profile {profile!r} has no BLAST database for marker {marker!r}"
            )
        relative = _safe_relative_path(databases[marker]["prefix"], "BLAST prefix", ManifestError)
        return directory / Path(relative)

    if profile == "curated":
        legacy = detect_legacy_database(root, blastdbcmd)
        if legacy is not None:
            warnings.warn(
                "Using the legacy unprofiled SSUextract database; reinstall the curated "
                "profile because legacy support will be removed in a future release",
                LegacyDatabaseWarning,
                stacklevel=2,
            )
            return legacy
    raise IntegrityError(f"Database profile is not installed: {directory}")


def download_archive(
    url: str,
    destination: str | Path,
    expected_size: int,
    expected_sha256: str,
    opener: Callable[..., BinaryIO] = urllib.request.urlopen,
    progress: ProgressReporter = None,
) -> None:
    parsed = urlparse(url)
    if parsed.scheme != "https" or not parsed.netloc:
        raise InstallError(f"Database archive URL must use HTTPS: {url}")
    _require_size(expected_size, "archive bytes", InstallError)
    _require_sha256(expected_sha256, "archive sha256", InstallError)

    try:
        download_verified_archive(
            url,
            destination,
            expected_size,
            expected_sha256,
            opener=opener,
            progress=progress,
        )
    except DownloadIntegrityError as error:
        raise IntegrityError(str(error)) from error
    except DownloadError as error:
        raise InstallError(str(error)) from error


def _catalog_for_install(
    catalog_path: str | Path,
    profile: str,
    latest: bool,
    timeout: float,
    opener: Callable[..., BinaryIO],
    progress: ProgressReporter,
) -> dict:
    catalog = load_catalog(catalog_path)
    if not latest:
        return catalog
    _report(progress, "Checking Zenodo for the latest database release...")
    try:
        discovered = discover_latest_catalog(catalog, timeout=timeout, opener=opener)
        pinned_version = catalog["profiles"][profile]["version"]
        discovered_version = discovered["profiles"][profile]["version"]
        if compare_versions(discovered_version, pinned_version) < 0:
            _report(
                progress,
                f"Zenodo reports v{discovered_version}, older than bundled "
                f"v{pinned_version}; using bundled release metadata.",
            )
            return catalog
        _report(progress, f"Latest verified Zenodo database release: v{discovered_version}")
        return discovered
    except (KeyError, ReleaseDiscoveryError) as error:
        bundled_version = catalog["profiles"].get(profile, {}).get("version", "unknown")
        _report(
            progress,
            f"Could not verify the latest Zenodo release ({' '.join(str(error).split())}); "
            f"using bundled release metadata v{bundled_version}.",
        )
        return catalog


def safe_extract_tar(archive: str | Path, destination: str | Path) -> None:
    destination_path = Path(destination).resolve()
    destination_path.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(archive, mode="r:*") as handle:
            for member in handle.getmembers():
                relative = _safe_relative_path(member.name, "archive member", InstallError)
                target = (destination_path / Path(relative)).resolve()
                if not _inside(destination_path, target):
                    raise InstallError(f"Archive member escapes extraction directory: {member.name}")
                if member.issym() or member.islnk():
                    raise InstallError(f"Archive links are not allowed: {member.name}")
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                if not member.isfile():
                    raise InstallError(f"Unsupported archive member type: {member.name}")
                target.parent.mkdir(parents=True, exist_ok=True)
                source = handle.extractfile(member)
                if source is None:
                    raise InstallError(f"Could not read archive member: {member.name}")
                with source, target.open("xb") as output:
                    shutil.copyfileobj(source, output)
    except DatabaseError:
        raise
    except (OSError, tarfile.TarError) as error:
        raise InstallError(f"Could not extract database archive {archive}: {error}") from error


def _find_extracted_profile(extracted: Path) -> Path:
    if (extracted / MANIFEST_NAME).is_file():
        return extracted
    children = [child for child in extracted.iterdir() if child.is_dir()]
    files = [child for child in extracted.iterdir() if child.is_file()]
    if not files and len(children) == 1 and (children[0] / MANIFEST_NAME).is_file():
        return children[0]
    raise InstallError(
        f"Archive must contain {MANIFEST_NAME} at its root or in one top-level directory"
    )


@contextlib.contextmanager
def _profile_install_lock(root: Path):
    """Serialize same-host installs; publication still defends cross-host races."""

    descriptor = os.open(root, os.O_RDONLY)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def _recover_interrupted_backup(
    root: Path, profile: str, blastdbcmd: str
) -> None:
    """Resolve an unambiguous backup left by an interrupted replacement."""

    target = root / profile
    backups = sorted(root.glob(f".{profile}.backup-*"))
    if not backups:
        return
    if len(backups) > 1:
        locations = ", ".join(str(path) for path in backups)
        raise InstallError(
            f"Interrupted replacement state for profile {profile!r}; inspect target "
            f"{target} and backup(s) {locations} before retrying"
        )
    if target.exists():
        try:
            validate_profile_directory(target, profile, blastdbcmd)
        except DatabaseError as error:
            raise InstallError(
                f"Interrupted replacement left both target {target} and backup "
                f"{backups[0]}, and the target does not validate: {error}"
            ) from error
        try:
            shutil.rmtree(backups[0])
            fsync_directory(root)
        except OSError as error:
            raise InstallError(
                f"Validated target {target}, but could not remove interrupted backup "
                f"{backups[0]}: {error}"
            ) from error
        warnings.warn(
            f"Retained validated profile {profile!r} and removed an interrupted "
            "replacement backup",
            RuntimeWarning,
            stacklevel=2,
        )
        return
    try:
        os.replace(backups[0], target)
        fsync_directory(root)
    except OSError as error:
        locations = [path for path in (target, backups[0]) if path.exists()]
        raise InstallError(
            f"Could not restore interrupted profile {profile!r}; surviving path(s): "
            f"{', '.join(str(path) for path in locations) or 'none'}: {error}"
        ) from error
    warnings.warn(
        f"Restored profile {profile!r} from an interrupted replacement backup",
        RuntimeWarning,
        stacklevel=2,
    )


def _remove_verified_download(path: Path, root: Path) -> None:
    try:
        path.unlink(missing_ok=True)
        fsync_directory(root)
    except OSError as error:
        warnings.warn(
            f"Could not remove verified download cache {path}: {error}",
            RuntimeWarning,
            stacklevel=2,
        )


def _install_profile_locked(
    root: str | Path,
    profile: str,
    catalog: dict,
    blastdbcmd: str = "blastdbcmd",
    replace: bool = False,
    opener: Callable[..., BinaryIO] = urllib.request.urlopen,
    progress: ProgressReporter = None,
) -> Path:
    _require_identifier(profile, "profile", CatalogError)
    try:
        entry = catalog["profiles"][profile]
    except KeyError as error:
        raise CatalogError(f"Unknown database profile {profile!r}") from error
    archive = entry["archive"]
    if archive["url"] is None:
        raise InstallError(
            f"Profile {profile!r} has no published archive URL yet; install cannot proceed"
        )

    root_path = Path(root).resolve()
    root_path.mkdir(parents=True, exist_ok=True)
    target = root_path / profile
    if target.exists() and not replace:
        raise InstallError(f"Database profile already exists: {target}; use --force to replace it")

    _report(
        progress,
        f"Installing profile {profile!r} v{entry['version']} "
        f"({_human_size(archive['bytes'])}) into {target}",
    )
    host = hashlib.sha256(socket.gethostname().encode()).hexdigest()[:12]
    archive_path = (
        root_path
        / f".{profile}.download-v{entry['version']}-{archive['sha256']}-{host}.part"
    )
    with contextlib.ExitStack() as cleanup:
        download_archive(
            archive["url"],
            archive_path,
            archive["bytes"],
            archive["sha256"],
            opener=opener,
            progress=progress,
        )
        cleanup.callback(_remove_verified_download, archive_path, root_path)
        temporary = cleanup.enter_context(
            tempfile.TemporaryDirectory(prefix=f".{profile}.staging-", dir=root_path)
        )
        staging = Path(temporary)
        extracted = staging / "extracted"
        _report(progress, "Extracting archive...")
        safe_extract_tar(archive_path, extracted)
        source = _find_extracted_profile(extracted)
        _report(progress, "Validating files and BLAST indexes...")
        manifest = validate_profile_directory(source, profile, blastdbcmd)
        if manifest["version"] != entry["version"]:
            raise IntegrityError(
                f"Profile version mismatch: catalog expects {entry['version']!r}, "
                f"archive contains {manifest['version']!r}"
            )

        incoming = staging / "ready"
        if target.exists() and not replace:
            raise InstallError(
                f"Database profile appeared during installation: {target}; "
                "use --force to replace it"
            )
        try:
            os.replace(source, incoming)
            if not target.exists():
                replace_and_fsync(incoming, target)
                _report(progress, f"Installed profile {profile!r} v{manifest['version']}.")
                return target
        except OSError as error:
            if target.exists() and not incoming.exists():
                raise InstallError(
                    f"Profile {profile!r} was moved into place, but its durability sync "
                    "failed; validate the installed profile before retrying: "
                    f"{error}"
                ) from error
            raise InstallError(f"Could not install database profile {profile!r}: {error}") from error

        # Advisory flock is host-local on some shared filesystems. Recheck here
        # so a cross-host winner is never converted into an implicit --force.
        if not replace:
            raise InstallError(
                f"Database profile appeared during publication: {target}; "
                "the existing profile was not replaced"
            )

        backup = root_path / f".{profile}.backup-{uuid.uuid4().hex}"
        try:
            os.replace(target, backup)
            replace_and_fsync(incoming, target)
        except OSError as error:
            if backup.exists():
                try:
                    if target.exists():
                        failed_incoming = staging / "failed-ready"
                        os.replace(target, failed_incoming)
                        fsync_directory(root_path)
                    os.replace(backup, target)
                    fsync_directory(root_path)
                except OSError as rollback_error:
                    locations = [path for path in (target, backup) if path.exists()]
                    raise InstallError(
                        f"Could not replace profile {profile!r}, and rollback failed; "
                        f"surviving path(s): "
                        f"{', '.join(str(path) for path in locations) or 'none'}: "
                        f"{rollback_error}"
                    ) from error
                raise InstallError(
                    f"Could not replace database profile {profile!r}; the previous profile "
                    f"was restored: {error}"
                ) from error
            raise InstallError(
                f"Could not replace database profile {profile!r}; the existing profile "
                f"remains in place: {error}"
            ) from error
        try:
            shutil.rmtree(backup)
            fsync_directory(root_path)
        except OSError as error:
            warnings.warn(
                f"Installed profile {profile!r}, but could not remove backup {backup}: {error}",
                RuntimeWarning,
                stacklevel=2,
            )
        _report(progress, f"Installed profile {profile!r} v{manifest['version']}.")
    return target


def install_profile(
    root: str | Path,
    profile: str,
    catalog_path: str | Path = DEFAULT_CATALOG,
    blastdbcmd: str = "blastdbcmd",
    replace: bool = False,
    opener: Callable[..., BinaryIO] = urllib.request.urlopen,
    latest: bool = False,
    timeout: float = 5.0,
    metadata_opener: Callable[..., BinaryIO] = urllib.request.urlopen,
    progress: ProgressReporter = None,
) -> Path:
    """Download, validate, and publish one profile under an install lock."""

    _require_identifier(profile, "profile", CatalogError)
    catalog = _catalog_for_install(
        catalog_path,
        profile,
        latest,
        timeout,
        metadata_opener,
        progress,
    )
    root_path = Path(root).resolve()
    root_path.mkdir(parents=True, exist_ok=True)
    with _profile_install_lock(root_path):
        _recover_interrupted_backup(root_path, profile, blastdbcmd)
        return _install_profile_locked(
            root_path,
            profile,
            catalog=catalog,
            blastdbcmd=blastdbcmd,
            replace=replace,
            opener=opener,
            progress=progress,
        )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Install and resolve SSUextract database profiles.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    profiles = subparsers.add_parser("profiles", help="list available database profiles")
    profiles.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)

    validate = subparsers.add_parser("validate", help="validate an installed profile")
    validate.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    validate.add_argument("--profile", required=True)
    validate.add_argument("--blastdbcmd", default="blastdbcmd")

    version = subparsers.add_parser("version", help="print an installed profile version")
    version.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    version.add_argument("--profile", required=True)

    update = subparsers.add_parser(
        "check-update", help="compare an installed profile with the latest Zenodo release"
    )
    update.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    update.add_argument("--profile", required=True)
    update.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    update.add_argument("--timeout", type=float, default=5.0)
    update.add_argument("--format", choices=("human", "tsv"), default="human")

    resolve = subparsers.add_parser("resolve", help="print the BLAST prefix for a model")
    resolve.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    resolve.add_argument("--profile", default="curated")
    resolve.add_argument("--model", required=True)
    resolve.add_argument("--markers", type=Path, default=DEFAULT_MARKERS)
    resolve.add_argument("--marker-override", type=Path)
    resolve.add_argument("--blastdbcmd", default="blastdbcmd")

    install = subparsers.add_parser("install", help="download and validate a profile installation")
    install.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    install.add_argument("--profile", default="curated")
    install.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    install.add_argument("--blastdbcmd", default="blastdbcmd")
    install.add_argument("--force", action="store_true")
    install.add_argument(
        "--latest",
        action="store_true",
        help="install the latest fully verified release from the Zenodo concept record",
    )
    install.add_argument("--timeout", type=float, default=5.0)
    return parser.parse_args(argv)


def _print_update_result(result: dict[str, str], profile: str, output_format: str) -> None:
    if output_format == "tsv":
        print(
            "\t".join(
                (
                    result["status"],
                    result["installed"],
                    result["latest"] or "-",
                    result["reason"] or "-",
                )
            )
        )
        return
    installed = result["installed"]
    latest = result["latest"]
    if result["status"] == "update_available":
        print(
            f"Database update available for {profile!r}: v{installed} -> v{latest}. "
            f"Run 'pixi run setup -- --database_profile {profile} --update'."
        )
    elif result["status"] == "current":
        print(f"Database profile {profile!r}: installed v{installed}; latest v{latest}.")
    elif result["status"] == "installed_newer":
        print(
            f"Database profile {profile!r}: installed v{installed}; "
            f"Zenodo latest is v{latest}."
        )
    else:
        print(
            f"Database profile {profile!r}: installed v{installed}. "
            f"Update check unavailable: {result['reason']}"
        )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.command == "profiles":
            catalog = load_catalog(args.catalog)
            for name, entry in catalog["profiles"].items():
                print(
                    "\t".join(
                        (
                            name,
                            entry["version"],
                            _human_size(entry["archive"]["bytes"]),
                            entry.get("description", ""),
                        )
                    )
                )
        elif args.command == "validate":
            validate_profile(args.root, args.profile, args.blastdbcmd)
            print((args.root / args.profile).resolve())
        elif args.command == "version":
            print(installed_version(args.root, args.profile))
        elif args.command == "check-update":
            result = check_profile_update(
                args.root,
                args.profile,
                catalog_path=args.catalog,
                timeout=args.timeout,
            )
            _print_update_result(result, args.profile, args.format)
        elif args.command == "resolve":
            print(
                resolve_database(
                    args.root,
                    args.profile,
                    args.model,
                    marker_override=args.marker_override,
                    marker_file=args.markers,
                    blastdbcmd=args.blastdbcmd,
                )
            )
        else:
            print(
                install_profile(
                    args.root,
                    args.profile,
                    catalog_path=args.catalog,
                    blastdbcmd=args.blastdbcmd,
                    replace=args.force,
                    latest=args.latest,
                    timeout=args.timeout,
                    progress=lambda message: print(
                        f"Database: {message}", file=sys.stderr, flush=True
                    ),
                )
            )
    except DatabaseError as error:
        print(f"database-manager: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
