#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
import queue
import re
import threading
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from typing import BinaryIO
from urllib.parse import urljoin, urlparse


MAX_JSON_BYTES = 2 * 1024 * 1024
MAX_CHECKSUM_BYTES = 64 * 1024
ZENODO_HOST = "zenodo.org"

_MD5 = re.compile(r"^[0-9a-f]{32}$")
_SEMANTIC_VERSION = re.compile(
    r"^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$"
)
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


class ReleaseDiscoveryError(RuntimeError):
    """Zenodo did not provide a complete, trusted database release contract."""


class _Deadline:
    def __init__(self, seconds: float):
        if seconds <= 0:
            raise ReleaseDiscoveryError("Zenodo timeout must be greater than zero")
        self._end = time.monotonic() + seconds

    def remaining(self) -> float:
        remaining = self._end - time.monotonic()
        if remaining <= 0:
            raise ReleaseDiscoveryError("Zenodo update check timed out")
        return remaining


def semantic_version(value: object, label: str = "version") -> tuple[int, int, int]:
    if not isinstance(value, str):
        raise ReleaseDiscoveryError(f"{label} must use MAJOR.MINOR.PATCH")
    match = _SEMANTIC_VERSION.fullmatch(value)
    if match is None:
        raise ReleaseDiscoveryError(f"{label} must use MAJOR.MINOR.PATCH: {value!r}")
    return tuple(int(part) for part in match.groups())


def compare_versions(installed: str, latest: str) -> int:
    installed_parts = semantic_version(installed, "installed database version")
    latest_parts = semantic_version(latest, "latest database version")
    return (installed_parts > latest_parts) - (installed_parts < latest_parts)


def _object_without_duplicate_keys(pairs: list[tuple[str, object]]) -> dict:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ReleaseDiscoveryError(f"Zenodo JSON contains a duplicate key: {key!r}")
        result[key] = value
    return result


def _require_object(value: object, label: str) -> dict:
    if not isinstance(value, dict):
        raise ReleaseDiscoveryError(f"{label} must be a JSON object")
    return value


def _require_positive_int(value: object, label: str) -> int:
    if type(value) is not int or value <= 0:
        raise ReleaseDiscoveryError(f"{label} must be a positive integer")
    return value


def _require_digest(value: object, pattern: re.Pattern[str], label: str) -> str:
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise ReleaseDiscoveryError(f"{label} is not a lowercase hexadecimal digest")
    return value


def _validate_zenodo_url(url: object, label: str) -> str:
    if not isinstance(url, str):
        raise ReleaseDiscoveryError(f"{label} must be a URL")
    parsed = urlparse(url)
    if (
        parsed.scheme != "https"
        or parsed.hostname != ZENODO_HOST
        or parsed.username is not None
        or parsed.password is not None
        or parsed.port is not None
        or parsed.query
        or parsed.fragment
        or not parsed.path.startswith("/api/records/")
    ):
        raise ReleaseDiscoveryError(
            f"{label} must be an unqualified https://{ZENODO_HOST}/api/records/ URL"
        )
    return url


def validate_zenodo_config(value: object) -> dict:
    config = _require_object(value, "catalog zenodo")
    record_id = _require_positive_int(config.get("record_id"), "Zenodo record_id")
    concept_record_id = config.get("concept_record_id")
    if not isinstance(concept_record_id, str) or not concept_record_id.isdigit():
        raise ReleaseDiscoveryError("Zenodo concept_record_id must contain only digits")
    api_url = _validate_zenodo_url(config.get("api_url"), "Zenodo api_url")
    if urlparse(api_url).path != f"/api/records/{record_id}":
        raise ReleaseDiscoveryError("Zenodo api_url does not match record_id")
    return config


def _read_url(
    url: str,
    label: str,
    limit: int,
    deadline: _Deadline,
    opener: Callable[..., BinaryIO],
) -> bytes:
    _validate_zenodo_url(url, label)
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "SSUextract database update check",
        },
    )
    chunks: list[bytes] = []
    size = 0
    try:
        with opener(request, timeout=deadline.remaining()) as response:
            geturl = getattr(response, "geturl", None)
            if callable(geturl):
                final_url = geturl()
                if final_url:
                    _validate_zenodo_url(final_url, f"redirected {label}")
            while True:
                chunk = response.read(min(64 * 1024, limit + 1 - size))
                if not chunk:
                    break
                chunks.append(chunk)
                size += len(chunk)
                if size > limit:
                    raise ReleaseDiscoveryError(f"{label} exceeds {limit} bytes")
                deadline.remaining()
    except ReleaseDiscoveryError:
        raise
    except (OSError, ValueError, urllib.error.URLError) as error:
        raise ReleaseDiscoveryError(f"Could not read {label}: {error}") from error
    return b"".join(chunks)


def _read_json_url(
    url: str,
    label: str,
    deadline: _Deadline,
    opener: Callable[..., BinaryIO],
) -> dict:
    data = _read_url(url, label, MAX_JSON_BYTES, deadline, opener)
    try:
        value = json.loads(data, object_pairs_hook=_object_without_duplicate_keys)
    except ReleaseDiscoveryError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ReleaseDiscoveryError(f"{label} is not valid UTF-8 JSON: {error}") from error
    return _require_object(value, label)


def _record_identity(record: dict, concept_record_id: str, expected_id: int | None) -> int:
    record_id = _require_positive_int(record.get("id"), "Zenodo release record id")
    if expected_id is not None and record_id != expected_id:
        raise ReleaseDiscoveryError(
            f"Zenodo record id mismatch: expected {expected_id}, found {record_id}"
        )
    if str(record.get("conceptrecid")) != concept_record_id:
        raise ReleaseDiscoveryError("Zenodo release belongs to a different concept record")
    return record_id


def _follow_body_location(
    value: dict,
    base_url: str,
    deadline: _Deadline,
    opener: Callable[..., BinaryIO],
) -> dict:
    if "id" in value:
        return value
    location = value.get("location")
    if not isinstance(location, str) or not location:
        raise ReleaseDiscoveryError("Zenodo latest-release response has no record")
    target = _validate_zenodo_url(urljoin(base_url, location), "Zenodo latest location")
    return _read_json_url(target, "Zenodo latest record", deadline, opener)


def _record_files(record: dict, record_id: int, version: str) -> dict[str, dict]:
    raw_files = record.get("files")
    if isinstance(raw_files, dict):
        raw_files = raw_files.get("entries")
    if not isinstance(raw_files, list):
        raise ReleaseDiscoveryError("Zenodo record files must be an array")

    files: dict[str, dict] = {}
    for raw_file in raw_files:
        file = _require_object(raw_file, "Zenodo file")
        name = file.get("key", file.get("filename"))
        if not isinstance(name, str) or not name or "/" in name or "\\" in name:
            raise ReleaseDiscoveryError("Zenodo file has an unsafe name")
        if name in files:
            raise ReleaseDiscoveryError(f"Zenodo record repeats file {name!r}")
        size = _require_positive_int(
            file.get("size", file.get("filesize")), f"Zenodo file {name!r} size"
        )
        checksum = file.get("checksum")
        if not isinstance(checksum, str) or not checksum.startswith("md5:"):
            raise ReleaseDiscoveryError(f"Zenodo file {name!r} has no MD5 checksum")
        md5 = _require_digest(checksum[4:], _MD5, f"Zenodo file {name!r} MD5")
        links = _require_object(file.get("links"), f"Zenodo file {name!r} links")
        content_url = links.get("content", links.get("self", links.get("download")))
        content_url = _validate_zenodo_url(content_url, f"Zenodo file {name!r} URL")
        expected_path = f"/api/records/{record_id}/files/{name}/content"
        if urlparse(content_url).path != expected_path:
            raise ReleaseDiscoveryError(
                f"Zenodo file {name!r} URL does not belong to release record {record_id}"
            )
        files[name] = {"bytes": size, "md5": md5, "url": content_url}

    # The publisher emits this exact inventory. Reject additions so an archive or
    # sidecar cannot enter the install path without a reviewed contract change.
    expected_names = {
        f"ssuextract-db-curated-v{version}.tar.zst",
        f"ssuextract-db-img-v{version}.tar.zst",
        f"ssuextract-db-release-v{version}.json",
        "SHA256SUMS",
    }
    if set(files) != expected_names:
        missing = sorted(expected_names - set(files))
        unexpected = sorted(set(files) - expected_names)
        raise ReleaseDiscoveryError(
            f"Zenodo file inventory mismatch; missing={missing}, unexpected={unexpected}"
        )
    return files


def _verified_file_content(
    file: dict,
    name: str,
    limit: int,
    deadline: _Deadline,
    opener: Callable[..., BinaryIO],
) -> bytes:
    if file["bytes"] > limit:
        raise ReleaseDiscoveryError(f"Zenodo file {name!r} exceeds {limit} bytes")
    content = _read_url(file["url"], f"Zenodo file {name!r}", limit, deadline, opener)
    if len(content) != file["bytes"]:
        raise ReleaseDiscoveryError(f"Zenodo file {name!r} size does not match its record")
    actual_md5 = hashlib.md5(content, usedforsecurity=False).hexdigest()
    if actual_md5 != file["md5"]:
        raise ReleaseDiscoveryError(f"Zenodo file {name!r} MD5 does not match its record")
    return content


def _parse_release_manifest(data: bytes, version: str) -> dict:
    try:
        value = json.loads(data, object_pairs_hook=_object_without_duplicate_keys)
    except ReleaseDiscoveryError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ReleaseDiscoveryError(f"Zenodo release manifest is invalid JSON: {error}") from error
    manifest = _require_object(value, "Zenodo release manifest")
    if manifest.get("schema_version") != 1:
        raise ReleaseDiscoveryError("Zenodo release manifest schema_version must be 1")
    release = _require_object(manifest.get("release"), "Zenodo release")
    if release.get("version") != version:
        raise ReleaseDiscoveryError("Zenodo record and release manifest versions differ")
    return manifest


def _validate_release_contract(
    manifest: dict,
    files: dict[str, dict],
    sums_data: bytes,
    version: str,
) -> dict[str, dict]:
    profiles = _require_object(manifest.get("profiles"), "Zenodo release profiles")
    if set(profiles) != {"curated", "img"}:
        raise ReleaseDiscoveryError("Zenodo release must contain curated and img profiles")

    archives: dict[str, dict] = {}
    for profile in ("curated", "img"):
        entry = _require_object(profiles[profile], f"Zenodo profile {profile!r}")
        if entry.get("version") != version:
            raise ReleaseDiscoveryError(f"Zenodo profile {profile!r} has the wrong version")
        archive = _require_object(entry.get("archive"), f"Zenodo profile {profile!r} archive")
        expected_name = f"ssuextract-db-{profile}-v{version}.tar.zst"
        if archive.get("filename") != expected_name:
            raise ReleaseDiscoveryError(f"Zenodo profile {profile!r} archive name is invalid")
        size = _require_positive_int(archive.get("bytes"), f"{profile} archive bytes")
        md5 = _require_digest(archive.get("md5"), _MD5, f"{profile} archive MD5")
        sha256 = _require_digest(
            archive.get("sha256"), _SHA256, f"{profile} archive SHA-256"
        )
        if files[expected_name]["bytes"] != size or files[expected_name]["md5"] != md5:
            raise ReleaseDiscoveryError(
                f"Zenodo profile {profile!r} archive metadata disagrees with the record"
            )
        archives[profile] = {
            "filename": expected_name,
            "url": files[expected_name]["url"],
            "bytes": size,
            "sha256": sha256,
        }

    checksum_contract = _require_object(manifest.get("checksums"), "checksum contract")
    if checksum_contract.get("filename") != "SHA256SUMS":
        raise ReleaseDiscoveryError("Checksum contract filename must be SHA256SUMS")
    checksum_bytes = _require_positive_int(
        checksum_contract.get("bytes"), "checksum contract bytes"
    )
    checksum_sha256 = _require_digest(
        checksum_contract.get("sha256"), _SHA256, "checksum contract SHA-256"
    )
    expected_names = {archive["filename"] for archive in archives.values()}
    covers = checksum_contract.get("covers")
    if not isinstance(covers, list) or set(covers) != expected_names or len(covers) != 2:
        raise ReleaseDiscoveryError("Checksum contract does not cover exactly both archives")
    if len(sums_data) != checksum_bytes or hashlib.sha256(sums_data).hexdigest() != checksum_sha256:
        raise ReleaseDiscoveryError("SHA256SUMS does not match the release manifest")

    try:
        sums_text = sums_data.decode("ascii")
    except UnicodeDecodeError as error:
        raise ReleaseDiscoveryError("SHA256SUMS is not ASCII") from error
    if not sums_text.endswith("\n"):
        raise ReleaseDiscoveryError("SHA256SUMS must end with a newline")
    sums: dict[str, str] = {}
    for line in sums_text.splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  ([A-Za-z0-9._-]+)", line)
        if match is None or match.group(2) in sums:
            raise ReleaseDiscoveryError("SHA256SUMS has an invalid or duplicate entry")
        sums[match.group(2)] = match.group(1)
    if set(sums) != expected_names:
        raise ReleaseDiscoveryError("SHA256SUMS does not list exactly both archives")
    for profile, archive in archives.items():
        if sums[archive["filename"]] != archive["sha256"]:
            raise ReleaseDiscoveryError(
                f"SHA256SUMS disagrees with the {profile!r} archive manifest"
            )
    return archives


def _discover_latest_catalog(
    catalog: dict,
    timeout: float,
    opener: Callable[..., BinaryIO],
) -> dict:
    zenodo = validate_zenodo_config(catalog.get("zenodo"))
    deadline = _Deadline(timeout)

    bootstrap = _read_json_url(
        zenodo["api_url"], "Zenodo database record", deadline, opener
    )
    _record_identity(bootstrap, zenodo["concept_record_id"], zenodo["record_id"])
    links = _require_object(bootstrap.get("links"), "Zenodo record links")
    latest_url = _validate_zenodo_url(links.get("latest"), "Zenodo latest-release URL")
    latest = _read_json_url(latest_url, "Zenodo latest-release record", deadline, opener)
    latest = _follow_body_location(latest, latest_url, deadline, opener)
    latest_record_id = _record_identity(latest, zenodo["concept_record_id"], None)

    metadata = _require_object(latest.get("metadata"), "Zenodo release metadata")
    version = metadata.get("version")
    semantic_version(version, "Zenodo release version")
    assert isinstance(version, str)
    files = _record_files(latest, latest_record_id, version)

    release_name = f"ssuextract-db-release-v{version}.json"
    manifest_data = _verified_file_content(
        files[release_name], release_name, MAX_JSON_BYTES, deadline, opener
    )
    sums_data = _verified_file_content(
        files["SHA256SUMS"], "SHA256SUMS", MAX_CHECKSUM_BYTES, deadline, opener
    )
    manifest = _parse_release_manifest(manifest_data, version)
    archives = _validate_release_contract(manifest, files, sums_data, version)

    descriptions = {
        name: entry.get("description", "")
        for name, entry in catalog.get("profiles", {}).items()
        if isinstance(entry, dict)
    }
    return {
        "schema_version": 1,
        "default_profile": catalog.get("default_profile", "curated"),
        "zenodo": {
            "record_id": latest_record_id,
            "concept_record_id": zenodo["concept_record_id"],
            "api_url": f"https://{ZENODO_HOST}/api/records/{latest_record_id}",
        },
        "profiles": {
            profile: {
                "version": version,
                "description": descriptions.get(profile, ""),
                "archive": {
                    "url": archives[profile]["url"],
                    "bytes": archives[profile]["bytes"],
                    "sha256": archives[profile]["sha256"],
                },
            }
            for profile in ("curated", "img")
        },
    }


def discover_latest_catalog(
    catalog: dict,
    timeout: float = 5.0,
    opener: Callable[..., BinaryIO] = urllib.request.urlopen,
) -> dict:
    """Return a verified latest-release catalog within a hard wall-time bound."""

    if timeout <= 0:
        raise ReleaseDiscoveryError("Zenodo timeout must be greater than zero")
    outcome: queue.Queue[tuple[bool, object]] = queue.Queue(maxsize=1)

    def discover() -> None:
        try:
            outcome.put((True, _discover_latest_catalog(catalog, timeout, opener)))
        except Exception as error:
            outcome.put((False, error))

    worker = threading.Thread(target=discover, name="zenodo-update-check", daemon=True)
    worker.start()
    worker.join(timeout)
    if worker.is_alive():
        # The daemon may finish its bounded read later, but it cannot delay this
        # process or mutate local database state.
        raise ReleaseDiscoveryError("Zenodo update check timed out")
    succeeded, value = outcome.get_nowait()
    if succeeded:
        assert isinstance(value, dict)
        return value
    if isinstance(value, ReleaseDiscoveryError):
        raise value
    assert isinstance(value, Exception)
    raise ReleaseDiscoveryError(f"Zenodo update check failed: {value}") from value
