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
import time
import urllib.request
import uuid
import warnings
from pathlib import Path, PurePosixPath
from typing import BinaryIO, Callable, TextIO
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
DEFAULT_ROOT = REPO / "resources" / "database"
MANIFEST_NAME = "manifest.json"
SCHEMA_VERSION = 1
EXIT_ERROR = 2
EXIT_NOT_INSTALLED = 3
EXIT_DAMAGED = 4

_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_DOWNLOAD_PROGRESS = re.compile(
    r"^Downloading archive: \[[=-]+\]\s+"
    r"(?P<percent>[0-9]+(?:\.[0-9]+)?)% "
    r"\((?P<written>.+?) of (?P<total>.+?)\) "
    r"(?P<rate>.+?/s) ETA (?P<eta>\S+)$"
)
_NON_TTY_PROGRESS_INTERVAL = 10.0

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


class NotInstalledError(DatabaseError):
    """The profile directory or its manifest is absent."""


class ArtifactDamagedError(IntegrityError):
    """A manifest artifact is missing, empty, or differs in size or SHA-256."""


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


class _ConsoleProgress:
    def __init__(
        self,
        stream: TextIO,
        clock: Callable[[], float] = time.monotonic,
        terminal_columns: Callable[[], int] | None = None,
        non_tty_interval: float = _NON_TTY_PROGRESS_INTERVAL,
    ) -> None:
        if non_tty_interval <= 0:
            raise ValueError("non-TTY progress interval must be positive")
        self.stream = stream
        self.clock = clock
        self.terminal_columns = terminal_columns or self._terminal_columns
        self.non_tty_interval = non_tty_interval
        self.last_non_tty_progress: float | None = None
        self.live = False

    def _terminal_columns(self) -> int:
        try:
            return os.get_terminal_size(self.stream.fileno()).columns
        except (AttributeError, OSError, ValueError):
            return shutil.get_terminal_size(fallback=(80, 24)).columns

    @staticmethod
    def _download_line(message: str, columns: int) -> str:
        limit = max(1, columns - 1)
        match = _DOWNLOAD_PROGRESS.fullmatch(message)
        if match is None:
            return f"Database: {message}"[:limit]

        percent = float(match.group("percent"))
        percent_text = f"{percent:5.1f}%"
        written = match.group("written")
        total = match.group("total")
        rate = match.group("rate")
        eta = match.group("eta")
        prefix = "Database: "
        tails = (
            f" {percent_text} {written}/{total} {rate} ETA {eta}",
            f" {percent_text} {written}/{total} {rate}",
            f" {percent_text} {written} {rate}",
            f" {percent_text} {written}",
        )
        for tail in tails:
            available = limit - len(prefix) - len(tail) - 2
            if available < 6:
                continue
            width = min(24, available)
            filled = min(width, int(percent * width / 100))
            bar = "=" * filled + "-" * (width - filled)
            return f"{prefix}[{bar}]{tail}"

        fallbacks = (
            f"{prefix}{percent_text} {written} {rate}",
            f"{prefix}{percent_text} {written}",
            f"{prefix}{percent_text}",
        )
        for rendered in fallbacks:
            if len(rendered) <= limit:
                return rendered
        return fallbacks[-1][:limit]

    def __call__(self, message: str) -> None:
        rendered = f"Database: {message}"
        is_download = message.startswith("Downloading archive: [")
        if self.stream.isatty() and is_download:
            rendered = self._download_line(message, self.terminal_columns())
            self.stream.write(f"\r{rendered}\x1b[K")
            self.stream.flush()
            self.live = True
            return
        if is_download:
            match = _DOWNLOAD_PROGRESS.fullmatch(message)
            percent = float(match.group("percent")) if match is not None else None
            now = self.clock()
            is_boundary = percent is not None and (percent <= 0 or percent >= 100)
            if (
                not is_boundary
                and self.last_non_tty_progress is not None
                and now - self.last_non_tty_progress < self.non_tty_interval
            ):
                return
            self.last_non_tty_progress = now
        self.finish()
        print(rendered, file=self.stream, flush=True)

    def finish(self) -> None:
        if self.live:
            self.stream.write("\n")
            self.stream.flush()
            self.live = False


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
        url, size, digest = (
            archive.get("url"),
            archive.get("bytes"),
            archive.get("sha256"),
        )
        if not isinstance(url, str) or urlparse(url).scheme != "https" or not urlparse(url).netloc:
            raise CatalogError(f"catalog profile {profile!r} archive URL must use HTTPS")
        _require_size(size, f"catalog profile {profile!r} archive bytes", CatalogError)
        _require_sha256(digest, f"catalog profile {profile!r} archive sha256", CatalogError)
    return catalog


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
    for key in ("preferred", "source_records"):
        taxonomy_path = _safe_relative_path(
            taxonomy_database.get(key),
            f"taxonomy_database.{key}",
            ManifestError,
        )
        if taxonomy_path not in artifact_paths:
            raise ManifestError(
                f"taxonomy_database.{key} is not a manifest-listed artifact: "
                f"{taxonomy_path}"
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
    *,
    quick: bool = False,
) -> dict:
    directory = Path(profile_directory).resolve()
    if not directory.is_dir():
        raise NotInstalledError(
            f"Database profile directory does not exist: {directory}"
        )
    if not (directory / MANIFEST_NAME).exists():
        raise NotInstalledError(f"Database profile {directory} has no {MANIFEST_NAME}")
    manifest = load_manifest(directory, expected_profile=expected_profile)

    for artifact in manifest["artifacts"]:
        relative = _safe_relative_path(artifact["path"], "artifact path", ManifestError)
        listed_path = directory / Path(relative)
        path = listed_path.resolve()
        if not _inside(directory, path):
            raise IntegrityError(f"Artifact escapes profile directory: {relative}")
        if listed_path.is_symlink():
            raise IntegrityError(f"Manifest artifact is a symbolic link: {relative}")
        if not path.is_file():
            raise ArtifactDamagedError(
                f"Manifest artifact is not a regular file: {relative}"
            )
        actual_size = path.stat().st_size
        if actual_size == 0:
            raise ArtifactDamagedError(f"Manifest artifact is empty: {relative}")
        if actual_size != artifact["bytes"]:
            raise ArtifactDamagedError(
                f"Artifact size mismatch for {relative}: expected {artifact['bytes']}, found {actual_size}"
            )
        if quick:
            continue
        actual_digest = _sha256(path)
        if actual_digest != artifact["sha256"]:
            raise ArtifactDamagedError(
                f"Artifact SHA-256 mismatch for {relative}: expected {artifact['sha256']}, "
                f"found {actual_digest}"
            )

    if quick:
        return manifest
    for database in manifest["blast_databases"].values():
        relative = _safe_relative_path(database["prefix"], "BLAST prefix", ManifestError)
        prefix = (directory / Path(relative)).resolve()
        if not _inside(directory, prefix):
            raise IntegrityError(f"BLAST prefix escapes profile directory: {relative}")
        _validate_blast_prefix(prefix, blastdbcmd)
    return manifest


def validate_profile(
    root: str | Path,
    profile: str,
    blastdbcmd: str = "blastdbcmd",
    *,
    quick: bool = False,
) -> dict:
    _require_identifier(profile, "profile", ManifestError)
    return validate_profile_directory(
        Path(root) / profile, profile, blastdbcmd, quick=quick
    )


def installed_version(root: str | Path, profile: str) -> str:
    _require_identifier(profile, "profile", ManifestError)
    manifest = load_manifest(Path(root) / profile, expected_profile=profile)
    return manifest["version"]


def _update_status(installed: str, latest: str) -> str:
    return {-1: "update_available", 0: "current", 1: "installed_newer"}[
        compare_versions(installed, latest)
    ]


def _read_update_cache(
    path: Path, installed: str, now: float, catalog_sha256: str
) -> dict[str, str] | None:
    try:
        cache = json.loads(path.read_text(encoding="utf-8"))
        checked_at = cache["checked_at"]
        result = cache["result"]
        if (
            cache["catalog_sha256"] == catalog_sha256
            and isinstance(checked_at, (int, float))
            and 0 <= now - checked_at < 24 * 60 * 60
            and isinstance(result, dict)
            and set(result) == {"status", "installed", "latest", "reason"}
            and all(isinstance(value, str) for value in result.values())
            and result["status"] in {"current", "installed_newer", "update_available"}
        ):
            result["installed"] = installed
            result["status"] = _update_status(installed, result["latest"])
            return result
    except (OSError, ValueError, KeyError, TypeError, ReleaseDiscoveryError):
        return None
    return None


def _write_update_cache(
    path: Path, result: dict[str, str], checked_at: float, catalog_sha256: str
) -> None:
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f"{path.name}.",
            delete=False,
        ) as handle:
            staged = Path(handle.name)
            json.dump(
                {
                    "checked_at": checked_at,
                    "catalog_sha256": catalog_sha256,
                    "result": result,
                },
                handle,
            )
        staged.replace(path)
    except OSError:
        # Read-only database roots remain usable for runs and setup checks.
        pass
    finally:
        if staged is not None:
            with contextlib.suppress(OSError):
                staged.unlink(missing_ok=True)


def check_profile_update(
    root: str | Path,
    profile: str,
    catalog_path: str | Path = DEFAULT_CATALOG,
    timeout: float = 5.0,
    opener: Callable[..., BinaryIO] = urllib.request.urlopen,
    *,
    use_cache: bool = False,
) -> dict[str, str]:
    installed = installed_version(root, profile)
    cache_path = Path(root) / f".{profile}.update-check.json"
    checked_at = time.time()
    catalog = load_catalog(catalog_path)
    catalog_sha256 = hashlib.sha256(
        json.dumps(catalog, sort_keys=True).encode("utf-8")
    ).hexdigest()
    if use_cache:
        cached = _read_update_cache(cache_path, installed, checked_at, catalog_sha256)
        if cached is not None:
            return cached
    print("Checking Zenodo for database updates...", file=sys.stderr)
    try:
        latest_catalog = discover_latest_catalog(catalog, timeout=timeout, opener=opener)
        latest = latest_catalog["profiles"][profile]["version"]
        status = _update_status(installed, latest)
    except (KeyError, ReleaseDiscoveryError) as error:
        result = {
            "status": "unavailable",
            "installed": installed,
            "latest": "",
            "reason": " ".join(str(error).split()),
        }
    else:
        result = {
            "status": status,
            "installed": installed,
            "latest": latest,
            "reason": "",
        }
        _write_update_cache(cache_path, result, checked_at, catalog_sha256)
    return result


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
        with tarfile.open(archive, mode="r|*") as handle:
            for member in handle:
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
    root_path = Path(root).resolve()
    root_path.mkdir(parents=True, exist_ok=True)
    target = root_path / profile
    if target.exists() and not replace:
        try:
            validate_profile_directory(target, profile, blastdbcmd)
        except DatabaseError:
            raise InstallError(
                f"Database profile already exists: {target}; use --force to replace it"
            ) from None
        _report(progress, f"Profile {profile!r} is already installed and valid.")
        return target

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
    require_version: str | None = None,
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
    if require_version is not None:
        expected = _require_identifier(
            require_version, "required database version", InstallError
        )
        selected = catalog["profiles"][profile]["version"]
        if selected != expected:
            raise InstallError(
                f"Database release changed during update: expected v{expected}, "
                f"selected v{selected}; run the command again"
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
    parser = argparse.ArgumentParser(
        description="Install and validate SSUextract database profiles."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    profiles = subparsers.add_parser("profiles", help="list available database profiles")
    profiles.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)

    validate = subparsers.add_parser("validate", help="validate an installed profile")
    validate.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    validate.add_argument("--profile", required=True)
    validate.add_argument("--blastdbcmd", default="blastdbcmd")
    validate.add_argument(
        "--quick", action="store_true", help="check manifest and artifact sizes only"
    )

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
    update.add_argument(
        "--cached",
        action="store_true",
        help="reuse update checks younger than 24 hours",
    )

    install = subparsers.add_parser("install", help="download and validate a profile installation")
    install.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    install.add_argument("--profile", default="curated")
    install.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    install.add_argument("--blastdbcmd", default="blastdbcmd")
    install.add_argument("--force", action="store_true")
    install.add_argument(
        "--require-version",
        help="install only when discovery selects this exact database version",
    )
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
            f"Run 'pixi run setup --database_profile {profile} --update'."
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


def _failure_status(command: str, error: Exception) -> int:
    """Return the exit status; `validate` separates absent and damaged profiles."""
    if command == "validate":
        if isinstance(error, NotInstalledError):
            return EXIT_NOT_INSTALLED
        if isinstance(error, ArtifactDamagedError):
            return EXIT_DAMAGED
    return EXIT_ERROR


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    progress = _ConsoleProgress(sys.stderr)
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
            validate_profile(args.root, args.profile, args.blastdbcmd, quick=args.quick)
            print((args.root / args.profile).resolve())
        elif args.command == "version":
            print(installed_version(args.root, args.profile))
        elif args.command == "check-update":
            result = check_profile_update(
                args.root,
                args.profile,
                catalog_path=args.catalog,
                timeout=args.timeout,
                use_cache=args.cached,
            )
            _print_update_result(result, args.profile, args.format)
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
                    progress=progress,
                    require_version=args.require_version,
                )
            )
    except (DatabaseError, OSError) as error:
        progress.finish()
        print(f"database-manager: {error}", file=sys.stderr)
        return _failure_status(args.command, error)
    finally:
        progress.finish()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
