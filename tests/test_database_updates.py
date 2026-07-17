import hashlib
import io
import json
import sys
import time
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import database_updates as updates


class FakeResponse(io.BytesIO):
    def __init__(self, data: bytes, url: str):
        super().__init__(data)
        self.url = url

    def geturl(self) -> str:
        return self.url

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class MappingOpener:
    def __init__(self, responses: dict[str, bytes]):
        self.responses = responses
        self.calls: list[tuple[str, float]] = []

    def __call__(self, request, timeout: float):
        url = request.full_url
        self.calls.append((url, timeout))
        try:
            data = self.responses[url]
        except KeyError as error:
            raise OSError(f"unexpected URL {url}") from error
        return FakeResponse(data, url)


def json_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def md5(data: bytes) -> str:
    return hashlib.md5(data, usedforsecurity=False).hexdigest()


def release_fixture() -> tuple[dict, dict[str, bytes]]:
    concept_id = "9"
    bootstrap_id = 10
    release_id = 20
    version = "1.2.0"
    base = "https://zenodo.org/api/records"
    names = {
        profile: f"ssuextract-db-{profile}-v{version}.tar.zst"
        for profile in ("curated", "img")
    }
    archive_data = {"curated": b"curated archive", "img": b"img archive"}
    archive_sha = {
        profile: hashlib.sha256(data).hexdigest()
        for profile, data in archive_data.items()
    }
    sums = "".join(
        f"{archive_sha[profile]}  {names[profile]}\n"
        for profile in ("curated", "img")
    ).encode()
    manifest = {
        "schema_version": 1,
        "release": {"version": version},
        "profiles": {
            profile: {
                "version": version,
                "archive": {
                    "filename": names[profile],
                    "bytes": len(archive_data[profile]),
                    "md5": md5(archive_data[profile]),
                    "sha256": archive_sha[profile],
                },
            }
            for profile in ("curated", "img")
        },
        "checksums": {
            "filename": "SHA256SUMS",
            "bytes": len(sums),
            "sha256": hashlib.sha256(sums).hexdigest(),
            "covers": [names["curated"], names["img"]],
        },
    }
    manifest_data = json_bytes(manifest)
    release_name = f"ssuextract-db-release-v{version}.json"

    def file_entry(name: str, data: bytes) -> dict:
        return {
            "key": name,
            "size": len(data),
            "checksum": f"md5:{md5(data)}",
            "links": {"self": f"{base}/{release_id}/files/{name}/content"},
        }

    latest = {
        "id": release_id,
        "conceptrecid": concept_id,
        "metadata": {"version": version},
        "files": [
            file_entry(names["curated"], archive_data["curated"]),
            file_entry(release_name, manifest_data),
            file_entry(names["img"], archive_data["img"]),
            file_entry("SHA256SUMS", sums),
        ],
    }
    bootstrap_url = f"{base}/{bootstrap_id}"
    latest_url = f"{base}/{bootstrap_id}/versions/latest"
    catalog = {
        "schema_version": 1,
        "default_profile": "curated",
        "zenodo": {
            "record_id": bootstrap_id,
            "concept_record_id": concept_id,
            "api_url": bootstrap_url,
        },
        "profiles": {
            "curated": {"description": "curated references"},
            "img": {"description": "IMG-enhanced references"},
        },
    }
    responses = {
        bootstrap_url: json_bytes(
            {
                "id": bootstrap_id,
                "conceptrecid": concept_id,
                "links": {"latest": latest_url},
            }
        ),
        latest_url: json_bytes(latest),
        f"{base}/{release_id}/files/{release_name}/content": manifest_data,
        f"{base}/{release_id}/files/SHA256SUMS/content": sums,
    }
    return catalog, responses


class DatabaseUpdateTests(unittest.TestCase):
    def test_complete_release_contract_produces_install_catalog(self) -> None:
        catalog, responses = release_fixture()
        result = updates.discover_latest_catalog(
            catalog, timeout=2, opener=MappingOpener(responses)
        )
        self.assertEqual(result["zenodo"]["record_id"], 20)
        self.assertEqual(result["profiles"]["curated"]["version"], "1.2.0")
        self.assertEqual(
            result["profiles"]["img"]["archive"]["sha256"],
            hashlib.sha256(b"img archive").hexdigest(),
        )
        self.assertEqual(
            result["profiles"]["curated"]["description"], "curated references"
        )

    def test_body_location_is_followed_and_revalidated(self) -> None:
        catalog, responses = release_fixture()
        latest_url = responses[catalog["zenodo"]["api_url"]]
        bootstrap = json.loads(latest_url)
        version_url = bootstrap["links"]["latest"]
        latest_record = responses[version_url]
        target = "https://zenodo.org/api/records/20"
        responses[version_url] = json_bytes({"location": target})
        responses[target] = latest_record

        result = updates.discover_latest_catalog(
            catalog, timeout=2, opener=MappingOpener(responses)
        )

        self.assertEqual(result["zenodo"]["record_id"], 20)

    def test_release_from_another_concept_is_rejected(self) -> None:
        catalog, responses = release_fixture()
        bootstrap = json.loads(responses[catalog["zenodo"]["api_url"]])
        latest_url = bootstrap["links"]["latest"]
        latest = json.loads(responses[latest_url])
        latest["conceptrecid"] = "other"
        responses[latest_url] = json_bytes(latest)

        with self.assertRaisesRegex(updates.ReleaseDiscoveryError, "different concept"):
            updates.discover_latest_catalog(
                catalog, timeout=2, opener=MappingOpener(responses)
            )

    def test_unexpected_file_inventory_is_rejected(self) -> None:
        catalog, responses = release_fixture()
        bootstrap = json.loads(responses[catalog["zenodo"]["api_url"]])
        latest_url = bootstrap["links"]["latest"]
        latest = json.loads(responses[latest_url])
        latest["files"].append(
            {
                "key": "unexpected.txt",
                "size": 1,
                "checksum": f"md5:{md5(b'x')}",
                "links": {
                    "self": "https://zenodo.org/api/records/20/files/unexpected.txt/content"
                },
            }
        )
        responses[latest_url] = json_bytes(latest)

        with self.assertRaisesRegex(updates.ReleaseDiscoveryError, "inventory mismatch"):
            updates.discover_latest_catalog(
                catalog, timeout=2, opener=MappingOpener(responses)
            )

    def test_release_manifest_md5_mismatch_is_rejected(self) -> None:
        catalog, responses = release_fixture()
        manifest_url = next(url for url in responses if "db-release" in url)
        responses[manifest_url] += b"corrupt"

        with self.assertRaisesRegex(updates.ReleaseDiscoveryError, "size does not match"):
            updates.discover_latest_catalog(
                catalog, timeout=2, opener=MappingOpener(responses)
            )

    def test_sha256sums_disagreement_is_rejected(self) -> None:
        catalog, responses = release_fixture()
        bootstrap = json.loads(responses[catalog["zenodo"]["api_url"]])
        latest_url = bootstrap["links"]["latest"]
        latest = json.loads(responses[latest_url])
        checksum_file = next(entry for entry in latest["files"] if entry["key"] == "SHA256SUMS")
        checksum_url = checksum_file["links"]["self"]
        bad_sums = responses[checksum_url].replace(b"a", b"b", 1)
        responses[checksum_url] = bad_sums
        checksum_file["checksum"] = f"md5:{md5(bad_sums)}"
        responses[latest_url] = json_bytes(latest)

        with self.assertRaisesRegex(updates.ReleaseDiscoveryError, "release manifest"):
            updates.discover_latest_catalog(
                catalog, timeout=2, opener=MappingOpener(responses)
            )

    def test_non_zenodo_content_url_is_rejected(self) -> None:
        catalog, responses = release_fixture()
        bootstrap = json.loads(responses[catalog["zenodo"]["api_url"]])
        latest_url = bootstrap["links"]["latest"]
        latest = json.loads(responses[latest_url])
        latest["files"][0]["links"]["self"] = "https://example.org/archive"
        responses[latest_url] = json_bytes(latest)

        with self.assertRaisesRegex(updates.ReleaseDiscoveryError, "zenodo.org"):
            updates.discover_latest_catalog(
                catalog, timeout=2, opener=MappingOpener(responses)
            )

    def test_semantic_versions_are_compared_numerically(self) -> None:
        self.assertEqual(updates.compare_versions("1.9.0", "1.10.0"), -1)
        self.assertEqual(updates.compare_versions("2.0.0", "1.10.0"), 1)
        self.assertEqual(updates.compare_versions("1.2.3", "1.2.3"), 0)

    def test_update_check_has_a_hard_wall_time_limit(self) -> None:
        catalog, _ = release_fixture()

        def slow_opener(request, timeout):
            time.sleep(1)
            raise OSError("still offline")

        started = time.monotonic()
        with self.assertRaisesRegex(updates.ReleaseDiscoveryError, "timed out"):
            updates.discover_latest_catalog(catalog, timeout=0.05, opener=slow_opener)
        self.assertLess(time.monotonic() - started, 0.5)


if __name__ == "__main__":
    unittest.main()
