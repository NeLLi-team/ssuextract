import contextlib
import hashlib
import io
import json
import re
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import database_download as downloader
import database_manager as manager


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")


def write_profile(
    root: Path,
    profile: str = "curated",
    marker: str = "16S",
    artifact_data: bytes = b"blast-index",
    artifact_bytes: int | None = None,
    artifact_sha256: str | None = None,
    artifact_path: str = "blast/ssu.unusual",
    version: str = "test-1",
) -> Path:
    directory = root / profile
    artifact = directory / artifact_path
    artifact.parent.mkdir(parents=True, exist_ok=True)
    artifact.write_bytes(artifact_data)
    write_json(
        directory / manager.MANIFEST_NAME,
        {
            "schema_version": 1,
            "profile": profile,
            "version": version,
            "artifacts": [
                {
                    "path": artifact_path,
                    "bytes": len(artifact_data)
                    if artifact_bytes is None
                    else artifact_bytes,
                    "sha256": sha256(artifact_data)
                    if artifact_sha256 is None
                    else artifact_sha256,
                }
            ],
            "blast_databases": {marker: {"prefix": "blast/ssu"}},
            "taxonomy_database": {
                "preferred": artifact_path,
                "source_records": artifact_path,
            },
        },
    )
    return directory


def blast_ok(*args, **kwargs) -> subprocess.CompletedProcess:
    return subprocess.CompletedProcess(args[0], 0, stdout="Database: test\n", stderr="")


def run_validate(root: Path, blastdbcmd: str = "blastdbcmd") -> tuple[int, str]:
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        status = manager.main(
            [
                "validate",
                "--root",
                str(root),
                "--profile",
                "curated",
                "--blastdbcmd",
                blastdbcmd,
            ]
        )
    return status, stderr.getvalue()


def write_published_catalog(path: Path) -> None:
    write_json(
        path,
        {
            "schema_version": 1,
            "default_profile": "curated",
            "profiles": {
                "curated": {
                    "version": "test-1",
                    "archive": {
                        "url": "https://example.org/curated.tar.gz",
                        "bytes": 1,
                        "sha256": sha256(b"x"),
                    },
                }
            },
        },
    )


class FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class FakeHTTPResponse(FakeResponse):
    def __init__(self, data: bytes, status: int, headers: dict[str, str]) -> None:
        super().__init__(data)
        self.status = status
        self.headers = headers

    def geturl(self) -> str:
        return "https://example.org/curated.tar.gz"


class InterruptingHTTPResponse(FakeHTTPResponse):
    def __init__(
        self,
        data: bytes,
        status: int,
        headers: dict[str, str],
        interrupt_after: int,
    ) -> None:
        super().__init__(data, status, headers)
        self.remaining = interrupt_after

    def read(self, size: int = -1) -> bytes:
        if self.remaining == 0:
            raise TimeoutError("simulated manager-level timeout")
        chunk = super().read(min(size, self.remaining))
        self.remaining -= len(chunk)
        return chunk

    def read1(self, size: int = -1) -> bytes:
        return self.read(size)


class TtyBuffer(io.StringIO):
    def isatty(self) -> bool:
        return True


def tar_bytes(files: dict[str, bytes]) -> bytes:
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name, data in files.items():
            member = tarfile.TarInfo(name)
            member.size = len(data)
            archive.addfile(member, io.BytesIO(data))
    return output.getvalue()


class DatabaseManagerTests(unittest.TestCase):
    def test_quick_validation_checks_sizes_without_hashing_or_blast(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            directory = write_profile(root, artifact_sha256="0" * 64)
            with (
                mock.patch.object(manager, "_sha256") as digest,
                mock.patch.object(manager.subprocess, "run") as blast,
                contextlib.redirect_stdout(io.StringIO()),
            ):
                status = manager.main(
                    ["validate", "--root", str(root), "--profile", "curated", "--quick"]
                )
            self.assertEqual(status, 0)
            digest.assert_not_called()
            blast.assert_not_called()
            for content in (b"wrong size", None):
                with self.subTest(content=content):
                    artifact = directory / "blast/ssu.unusual"
                    if content is None:
                        artifact.unlink()
                    else:
                        artifact.write_bytes(content)
                    with contextlib.redirect_stderr(io.StringIO()):
                        status = manager.main(
                            [
                                "validate",
                                "--root",
                                str(root),
                                "--profile",
                                "curated",
                                "--quick",
                            ]
                        )
                    self.assertEqual(status, 4)

    def test_manifest_requires_runtime_taxonomy_artifact_keys(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = write_profile(Path(tmp))
            manifest_path = directory / manager.MANIFEST_NAME
            original = json.loads(manifest_path.read_text())
            for key in ("preferred", "source_records"):
                for value in (None, "../outside", "unlisted.parquet"):
                    with self.subTest(key=key, value=value):
                        manifest = json.loads(json.dumps(original))
                        if value is None:
                            del manifest["taxonomy_database"][key]
                        else:
                            manifest["taxonomy_database"][key] = value
                        write_json(manifest_path, manifest)
                        with self.assertRaises(manager.ManifestError):
                            manager.load_manifest(directory)

    def test_tar_extraction_reads_members_in_one_streaming_pass(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root / "profile.tar.gz"
            archive.write_bytes(
                tar_bytes({"profile/first": b"one", "profile/second": b"two"})
            )
            real_open = manager.tarfile.open
            with (
                mock.patch.object(manager.tarfile, "open", wraps=real_open) as opened,
                mock.patch.object(
                    manager.tarfile.TarFile,
                    "getmembers",
                    side_effect=AssertionError("second pass"),
                ),
            ):
                manager.safe_extract_tar(archive, root / "extract")
            self.assertEqual(opened.call_args.kwargs["mode"], "r|*")
            self.assertEqual((root / "extract/profile/first").read_bytes(), b"one")
            self.assertEqual((root / "extract/profile/second").read_bytes(), b"two")

    @mock.patch.object(manager, "discover_latest_catalog")
    def test_cached_update_check_reuses_result_until_24_hours(self, discover) -> None:
        discover.return_value = {"profiles": {"curated": {"version": "1.1.0"}}}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_profile(root, version="1.0.0")
            with mock.patch.object(manager.time, "time", return_value=100_000):
                first = manager.check_profile_update(root, "curated", use_cache=True)
            cache = root / ".curated.update-check.json"
            self.assertEqual(json.loads(cache.read_text())["checked_at"], 100_000)
            with mock.patch.object(manager.time, "time", return_value=186_399):
                self.assertEqual(
                    manager.check_profile_update(root, "curated", use_cache=True), first
                )
            self.assertEqual(discover.call_count, 1)
            with mock.patch.object(manager.time, "time", return_value=186_400):
                manager.check_profile_update(root, "curated", use_cache=True)
            self.assertEqual(discover.call_count, 2)
            self.assertEqual(json.loads(cache.read_text())["checked_at"], 186_400)
            manager.check_profile_update(root, "curated")
            self.assertEqual(discover.call_count, 3)

    @mock.patch.object(manager, "discover_latest_catalog")
    def test_update_cache_is_scoped_to_root_and_profile_after_installation(
        self, discover
    ) -> None:
        discover.return_value = {
            "profiles": {"curated": {"version": "1.1.0"}, "img": {"version": "1.1.0"}}
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for database_root, profile, version in (
                (root, "curated", "1.0.0"),
                (root, "img", "1.0.0"),
                (root / "other", "curated", "1.0.0"),
                (root, "curated", "1.1.0"),
            ):
                write_profile(database_root, profile=profile, version=version)
                result = manager.check_profile_update(
                    database_root, profile, use_cache=True
                )
                self.assertEqual(result["installed"], version)
            self.assertEqual(discover.call_count, 3)
            self.assertEqual(result["status"], "current")

    @mock.patch.object(
        manager,
        "discover_latest_catalog",
        side_effect=manager.ReleaseDiscoveryError("offline"),
    )
    def test_unavailable_update_result_is_not_cached(self, discover) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_profile(root, version="1.0.0")
            for _ in range(2):
                result = manager.check_profile_update(root, "curated", use_cache=True)
                self.assertEqual(result["status"], "unavailable")
            self.assertEqual(discover.call_count, 2)
            self.assertFalse((root / ".curated.update-check.json").exists())

    @mock.patch.object(manager, "discover_latest_catalog")
    def test_update_cache_checks_catalog_contents(self, discover) -> None:
        discover.side_effect = [
            {"profiles": {"curated": {"version": version}}}
            for version in ("1.0.0", "1.1.0", "1.2.0")
        ]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_profile(root, version="1.0.0")
            first_catalog = root / "catalog.json"
            write_published_catalog(first_catalog)
            second_catalog = root / "other-catalog.json"
            catalog = json.loads(first_catalog.read_text(encoding="utf-8"))
            catalog["profiles"]["curated"]["version"] = "1.1.0"
            write_json(second_catalog, catalog)
            for path, expected in (
                (first_catalog, "1.0.0"),
                (second_catalog, "1.1.0"),
            ):
                with self.subTest(catalog=path):
                    result = manager.check_profile_update(
                        root, "curated", catalog_path=path, use_cache=True
                    )
                    self.assertEqual(result["latest"], expected)
            catalog["profiles"]["curated"]["version"] = "1.2.0"
            write_json(second_catalog, catalog)
            result = manager.check_profile_update(
                root, "curated", catalog_path=second_catalog, use_cache=True
            )
            self.assertEqual(result["latest"], "1.2.0")
            self.assertEqual(discover.call_count, 3)

    @mock.patch.object(manager, "discover_latest_catalog")
    def test_update_query_notice_is_emitted_only_for_discovery(self, discover) -> None:
        discover.return_value = {"profiles": {"curated": {"version": "1.0.0"}}}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_profile(root, version="1.0.0")
            for flags, queries in (([], 1), (["--cached"], 0), ([], 1)):
                stdout, stderr = io.StringIO(), io.StringIO()
                with (
                    contextlib.redirect_stdout(stdout),
                    contextlib.redirect_stderr(stderr),
                ):
                    status = manager.main(
                        [
                            "check-update",
                            "--root",
                            str(root),
                            "--profile",
                            "curated",
                            "--format",
                            "tsv",
                            *flags,
                        ]
                    )
                self.assertEqual(status, 0)
                self.assertEqual(stdout.getvalue(), "current\t1.0.0\t1.0.0\t-\n")
                self.assertEqual(
                    stderr.getvalue().count("Checking Zenodo for database updates..."),
                    queries,
                )
            self.assertEqual(discover.call_count, 2)

    @mock.patch.object(manager, "discover_latest_catalog")
    def test_bad_or_unwritable_update_cache_does_not_fail_check(self, discover) -> None:
        discover.return_value = {"profiles": {"curated": {"version": "1.0.0"}}}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_profile(root, version="1.0.0")
            cache = root / ".curated.update-check.json"
            for content in ("{", "[]", '{"checked_at":"bad","result":{}}'):
                cache.write_text(content)
                self.assertEqual(
                    manager.check_profile_update(root, "curated", use_cache=True)[
                        "status"
                    ],
                    "current",
                )
            with mock.patch.object(
                manager.tempfile,
                "NamedTemporaryFile",
                side_effect=PermissionError("read only"),
            ):
                self.assertEqual(
                    manager.check_profile_update(root, "curated")["status"], "current"
                )

    def test_repository_catalog_has_coherent_release_state(self) -> None:
        catalog = manager.load_catalog(REPO / "config" / "database_catalog.json")
        self.assertEqual(set(catalog["profiles"]), {"curated", "img"})
        self.assertEqual(catalog["zenodo"]["concept_record_id"], "21367798")
        for profile in ("curated", "img"):
            archive = catalog["profiles"][profile]["archive"]
            self.assertRegex(
                archive["url"],
                rf"^https://zenodo\.org/records/[0-9]+/files/"
                rf"ssuextract-db-{profile}-v[0-9]+\.[0-9]+\.[0-9]+"
                r"\.tar\.zst\?download=1$",
            )

    def test_tutorial_profile_versions_and_sizes_match_catalog(self) -> None:
        catalog = manager.load_catalog(REPO / "config" / "database_catalog.json")
        tutorial = (REPO / "docs" / "tutorials" / "first-run.md").read_text()
        for number, profile in enumerate(("curated", "img"), start=1):
            entry = catalog["profiles"][profile]
            expected = (
                f"{number}) {profile} v{entry['version']} "
                f"({manager._human_size(entry['archive']['bytes'])})"
            )
            self.assertIn(expected, tutorial)

    def test_human_size_formats_binary_unit_boundaries(self) -> None:
        self.assertEqual(manager._human_size(1023), "1023 B")
        self.assertEqual(manager._human_size(1024), "1.0 KiB")
        self.assertEqual(manager._human_size(1024**2), "1.0 MiB")
        self.assertEqual(manager._human_size(1024**3), "1.0 GiB")

    @mock.patch.object(manager.subprocess, "run", side_effect=blast_ok)
    def test_artifact_size_must_match_manifest(self, run) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = write_profile(Path(tmp), artifact_bytes=99)
            with self.assertRaisesRegex(manager.IntegrityError, "size mismatch"):
                manager.validate_profile_directory(directory, "curated")
        run.assert_not_called()

    @mock.patch.object(manager.subprocess, "run", side_effect=blast_ok)
    def test_artifact_hash_must_match_manifest(self, run) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = write_profile(Path(tmp), artifact_sha256="0" * 64)
            with self.assertRaisesRegex(manager.IntegrityError, "SHA-256 mismatch"):
                manager.validate_profile_directory(directory, "curated")
        run.assert_not_called()

    @mock.patch.object(manager.subprocess, "run", side_effect=blast_ok)
    def test_symlinked_artifact_inside_profile_is_rejected(self, run) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = write_profile(Path(tmp))
            artifact = directory / "blast" / "ssu.unusual"
            target = directory / "blast" / "real-file"
            artifact.rename(target)
            artifact.symlink_to(target.name)

            with self.assertRaisesRegex(
                manager.IntegrityError, "symbolic link"
            ) as caught:
                manager.validate_profile_directory(directory, "curated")
            status, _ = run_validate(Path(tmp))
        self.assertNotIsInstance(caught.exception, manager.ArtifactDamagedError)
        self.assertEqual(status, 2)
        run.assert_not_called()

    def test_manifest_rejects_artifact_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / "curated"
            write_json(
                directory / manager.MANIFEST_NAME,
                {
                    "schema_version": 1,
                    "profile": "curated",
                    "version": "test-1",
                    "artifacts": [
                        {"path": "../outside", "bytes": 1, "sha256": "0" * 64}
                    ],
                    "blast_databases": {"16S": {"prefix": "blast/ssu"}},
                    "taxonomy_database": {"preferred": "../outside"},
                },
            )
            with self.assertRaisesRegex(manager.ManifestError, "escapes"):
                manager.load_manifest(directory, "curated")

    def test_validate_reports_absent_profile_directory_as_not_installed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            status, stderr = run_validate(Path(tmp))
        self.assertEqual(status, 3)
        self.assertIn("does not exist", stderr)

    def test_validate_reports_absent_manifest_as_not_installed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "curated").mkdir()
            status, stderr = run_validate(Path(tmp))
        self.assertEqual(status, 3)
        self.assertIn(manager.MANIFEST_NAME, stderr)

    def test_validate_reports_size_mismatch_as_damaged(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_profile(Path(tmp), artifact_bytes=99)
            status, stderr = run_validate(Path(tmp))
        self.assertEqual(status, 4)
        self.assertIn("size mismatch", stderr)

    def test_validate_reports_hash_mismatch_as_damaged(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_profile(Path(tmp), artifact_sha256="0" * 64)
            status, stderr = run_validate(Path(tmp))
        self.assertEqual(status, 4)
        self.assertIn("SHA-256 mismatch", stderr)

    def test_validate_reports_missing_artifact_as_damaged(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = write_profile(Path(tmp))
            (directory / "blast" / "ssu.unusual").unlink()
            status, stderr = run_validate(Path(tmp))
        self.assertEqual(status, 4)
        self.assertIn("not a regular file", stderr)

    def test_validate_reports_invalid_manifest_with_generic_status(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = write_profile(Path(tmp))
            (directory / manager.MANIFEST_NAME).write_text("{", encoding="utf-8")
            status, stderr = run_validate(Path(tmp))
        self.assertEqual(status, 2)
        self.assertIn("database-manager:", stderr)

    def test_validate_reports_missing_blastdbcmd_with_generic_status(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_profile(Path(tmp))
            status, stderr = run_validate(
                Path(tmp), blastdbcmd="/nonexistent/blastdbcmd"
            )
        self.assertEqual(status, 2)
        self.assertIn("BLAST validator not found", stderr)

    def test_unusable_database_root_is_reported_without_traceback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            blocker = Path(tmp) / "file"
            blocker.write_text("not a directory", encoding="utf-8")
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                status = manager.main(
                    [
                        "install",
                        "--root",
                        str(blocker / "database"),
                        "--profile",
                        "curated",
                    ]
                )

        self.assertEqual(status, 2)
        self.assertIn("database-manager:", stderr.getvalue())

    def test_zstd_tar_archive_is_extracted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root / "profile.tar.zst"
            with tarfile.open(archive, mode="w:zst") as handle:
                member = tarfile.TarInfo("curated/file.txt")
                member.size = len(b"payload")
                handle.addfile(member, io.BytesIO(b"payload"))

            manager.safe_extract_tar(archive, root / "extract")

            self.assertEqual(
                (root / "extract" / "curated" / "file.txt").read_bytes(), b"payload"
            )

    def test_declared_python_range_can_read_zstd_archives(self) -> None:
        text = (REPO / "pixi.toml").read_text(encoding="utf-8")
        match = re.search(r'^python = ">=3\.(\d+),', text, flags=re.MULTILINE)
        self.assertIsNotNone(match)
        # tarfile reads .tar.zst archives from Python 3.14.
        self.assertGreaterEqual(int(match.group(1)), 14)

    def test_tar_extraction_rejects_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root / "bad.tar.gz"
            archive.write_bytes(tar_bytes({"../outside": b"payload"}))
            with self.assertRaisesRegex(manager.InstallError, "escapes"):
                manager.safe_extract_tar(archive, root / "extract")
            self.assertFalse((root / "outside").exists())

    def test_tar_extraction_rejects_links(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root / "link.tar"
            with tarfile.open(archive, mode="w") as handle:
                member = tarfile.TarInfo("profile/link")
                member.type = tarfile.SYMTYPE
                member.linkname = "../../outside"
                handle.addfile(member)
            with self.assertRaisesRegex(manager.InstallError, "links are not allowed"):
                manager.safe_extract_tar(archive, root / "extract")

    def test_download_requires_https(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(manager.InstallError, "must use HTTPS"):
                manager.download_archive(
                    "http://example.org/profile.tar.gz",
                    Path(tmp) / "profile.tar.gz",
                    1,
                    sha256(b"x"),
                )

    def test_download_reports_visible_progress(self) -> None:
        data = b"x" * 100
        messages: list[str] = []
        with tempfile.TemporaryDirectory() as tmp:
            manager.download_archive(
                "https://example.org/profile.tar.gz",
                Path(tmp) / "profile.tar.gz",
                len(data),
                sha256(data),
                opener=lambda url, timeout: FakeResponse(data),
                progress=messages.append,
            )
        self.assertTrue(any("[------------------------]" in message for message in messages))
        self.assertTrue(any("100.0%" in message for message in messages))
        self.assertTrue(any("ETA" in message for message in messages))

    def test_console_progress_redraws_terminal_bar_and_finishes_line(self) -> None:
        stream = TtyBuffer()
        progress = manager._ConsoleProgress(stream, terminal_columns=lambda: 60)
        progress(
            "Downloading archive: [------------------------]   0.0% "
            "(0 B of 828.8 MiB) -- B/s ETA --:--"
        )
        progress(
            "Downloading archive: [============------------]  50.0% "
            "(414.4 MiB of 828.8 MiB) 1.0 MiB/s ETA 06:54"
        )
        progress("Extracting archive...")
        rendered = stream.getvalue()
        updates = [
            part.split("\x1b[K", 1)[0]
            for part in rendered.split("\r")[1:]
        ]
        self.assertEqual(rendered.count("\rDatabase: ["), 2)
        self.assertTrue(all(len(update) <= 59 for update in updates))
        self.assertIn("\x1b[K\nDatabase: Extracting archive...\n", rendered)
        self.assertFalse(progress.live)

    def test_console_progress_adapts_detail_to_terminal_width(self) -> None:
        message = (
            "Downloading archive: [============------------]  50.0% "
            "(414.4 MiB of 828.8 MiB) 1.0 MiB/s ETA 06:54"
        )
        narrow = manager._ConsoleProgress._download_line(message, 40)
        wide = manager._ConsoleProgress._download_line(message, 120)
        self.assertLessEqual(len(narrow), 39)
        self.assertLessEqual(len(wide), 119)
        self.assertIn("[", narrow)
        self.assertIn("50.0%", narrow)
        self.assertIn("1.0 MiB/s", wide)
        self.assertIn("ETA 06:54", wide)

    def test_console_progress_rate_limits_nonterminal_updates(self) -> None:
        stream = io.StringIO()
        times = iter((0.0, 1.0, 11.0, 12.0))
        progress = manager._ConsoleProgress(stream, clock=lambda: next(times))
        progress(
            "Downloading archive: [------------------------]   0.0% "
            "(0 B of 828.8 MiB) -- B/s ETA --:--"
        )
        progress(
            "Downloading archive: [==----------------------]  10.0% "
            "(82.9 MiB of 828.8 MiB) 1.0 MiB/s ETA 12:26"
        )
        progress(
            "Downloading archive: [=====-------------------]  20.0% "
            "(165.8 MiB of 828.8 MiB) 1.0 MiB/s ETA 11:03"
        )
        progress(
            "Downloading archive: [========================] 100.0% "
            "(828.8 MiB of 828.8 MiB) 1.0 MiB/s ETA 00:00"
        )
        self.assertNotIn("\r", stream.getvalue())
        self.assertEqual(stream.getvalue().count("Database: Downloading archive:"), 3)
        self.assertNotIn(" 10.0%", stream.getvalue())
        self.assertIn(" 20.0%", stream.getvalue())
        self.assertIn("100.0%", stream.getvalue())

    @mock.patch.object(manager, "discover_latest_catalog")
    def test_update_check_reports_newer_verified_release(self, discover) -> None:
        discover.return_value = {
            "profiles": {"curated": {"version": "1.1.0"}}
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_profile(root, version="1.0.0")
            result = manager.check_profile_update(root, "curated")
        self.assertEqual(result["status"], "update_available")
        self.assertEqual(result["installed"], "1.0.0")
        self.assertEqual(result["latest"], "1.1.0")

    @mock.patch.object(manager, "discover_latest_catalog")
    def test_update_check_reports_current_and_installed_newer_states(self, discover) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_profile(root, version="1.1.0")
            for latest, expected in (("1.1.0", "current"), ("1.0.0", "installed_newer")):
                with self.subTest(latest=latest):
                    discover.return_value = {
                        "profiles": {"curated": {"version": latest}}
                    }
                    result = manager.check_profile_update(root, "curated")
                    self.assertEqual(result["status"], expected)

    @mock.patch.object(manager, "discover_latest_catalog")
    def test_latest_install_catalog_accepts_newer_verified_release(self, discover) -> None:
        discovered = {
            "profiles": {"curated": {"version": "1.1.0", "archive": {}}}
        }
        discover.return_value = discovered
        messages: list[str] = []
        selected = manager._catalog_for_install(
            manager.DEFAULT_CATALOG,
            "curated",
            latest=True,
            timeout=1,
            opener=mock.Mock(),
            progress=messages.append,
        )
        self.assertIs(selected, discovered)
        self.assertIn("Latest verified Zenodo database release: v1.1.0", messages)

    @mock.patch.object(manager, "discover_latest_catalog")
    def test_latest_install_catalog_never_downgrades_bundled_release(self, discover) -> None:
        discover.return_value = {
            "profiles": {"curated": {"version": "0.9.0", "archive": {}}}
        }
        bundled = manager.load_catalog(manager.DEFAULT_CATALOG)
        selected = manager._catalog_for_install(
            manager.DEFAULT_CATALOG,
            "curated",
            latest=True,
            timeout=1,
            opener=mock.Mock(),
            progress=None,
        )
        self.assertEqual(
            selected["profiles"]["curated"]["version"],
            bundled["profiles"]["curated"]["version"],
        )

    @mock.patch.object(
        manager,
        "discover_latest_catalog",
        side_effect=manager.ReleaseDiscoveryError("offline"),
    )
    def test_latest_install_catalog_falls_back_when_zenodo_is_unavailable(self, discover) -> None:
        messages: list[str] = []
        bundled = manager.load_catalog(manager.DEFAULT_CATALOG)
        selected = manager._catalog_for_install(
            manager.DEFAULT_CATALOG,
            "curated",
            latest=True,
            timeout=1,
            opener=mock.Mock(),
            progress=messages.append,
        )
        self.assertEqual(
            selected["profiles"]["curated"]["version"],
            bundled["profiles"]["curated"]["version"],
        )
        self.assertTrue(any("using bundled release metadata" in message for message in messages))

    @mock.patch.object(manager, "_install_profile_locked")
    @mock.patch.object(manager, "_catalog_for_install")
    def test_required_install_version_rejects_changed_release(
        self, select_catalog, install
    ) -> None:
        select_catalog.return_value = {
            "profiles": {"curated": {"version": "1.0.2"}}
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "database"
            with self.assertRaisesRegex(
                manager.InstallError,
                "expected v1.0.3, selected v1.0.2",
            ):
                manager.install_profile(
                    root,
                    "curated",
                    latest=True,
                    require_version="1.0.3",
                )
            self.assertFalse(root.exists())
        install.assert_not_called()

    @mock.patch.object(manager, "_install_profile_locked")
    @mock.patch.object(manager, "_catalog_for_install")
    def test_required_install_version_accepts_exact_release(
        self, select_catalog, install
    ) -> None:
        select_catalog.return_value = {
            "profiles": {"curated": {"version": "1.0.2"}}
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "database"
            expected = root / "curated"
            install.return_value = expected
            result = manager.install_profile(
                root,
                "curated",
                latest=True,
                require_version="1.0.2",
            )
        self.assertEqual(result, expected)
        install.assert_called_once()

    @mock.patch.object(
        manager,
        "discover_latest_catalog",
        side_effect=manager.ReleaseDiscoveryError("offline"),
    )
    def test_update_check_is_nonfatal_when_zenodo_is_unavailable(self, discover) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_profile(root, version="1.0.0")
            result = manager.check_profile_update(root, "curated")
        self.assertEqual(result["status"], "unavailable")
        self.assertEqual(result["installed"], "1.0.0")
        self.assertEqual(result["reason"], "offline")

    def test_interrupted_backup_is_restored_before_install_retry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp, self.assertWarns(RuntimeWarning):
            root = Path(tmp)
            backup = root / ".curated.backup-test"
            backup.mkdir()
            (backup / "sentinel").write_text("previous profile", encoding="utf-8")

            manager._recover_interrupted_backup(root, "curated", "blastdbcmd")

            self.assertFalse(backup.exists())
            self.assertEqual(
                (root / "curated" / "sentinel").read_text(encoding="utf-8"),
                "previous profile",
            )

    @mock.patch.object(manager.subprocess, "run", side_effect=blast_ok)
    def test_valid_target_wins_over_interrupted_backup(self, run) -> None:
        with tempfile.TemporaryDirectory() as tmp, self.assertWarns(RuntimeWarning):
            root = Path(tmp)
            target = write_profile(root)
            backup = root / ".curated.backup-test"
            backup.mkdir()
            (backup / "sentinel").write_text("previous profile", encoding="utf-8")

            manager._recover_interrupted_backup(root, "curated", "blastdbcmd")

            self.assertTrue(target.is_dir())
            self.assertFalse(backup.exists())

    @mock.patch.object(manager.subprocess, "run", side_effect=blast_ok)
    def test_install_keeps_profile_restored_from_interrupted_backup(self, run) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "databases"
            catalog = Path(tmp) / "catalog.json"
            write_published_catalog(catalog)
            backup = write_profile(root).rename(root / ".curated.backup-test")
            opener = mock.Mock(side_effect=AssertionError("must not download"))

            with self.assertWarns(RuntimeWarning):
                installed = manager.install_profile(
                    root, "curated", catalog_path=catalog, opener=opener
                )

            self.assertEqual(installed, root / "curated")
            self.assertFalse(backup.exists())
            opener.assert_not_called()

    @mock.patch.object(manager.subprocess, "run", side_effect=blast_ok)
    def test_install_without_replace_keeps_valid_existing_profile(self, run) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "databases"
            catalog = Path(tmp) / "catalog.json"
            write_published_catalog(catalog)
            write_profile(root)
            opener = mock.Mock(side_effect=AssertionError("must not download"))

            installed = manager.install_profile(
                root, "curated", catalog_path=catalog, opener=opener
            )

            self.assertEqual(installed, root / "curated")
            opener.assert_not_called()

    def test_install_without_replace_rejects_invalid_existing_profile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "databases"
            catalog = Path(tmp) / "catalog.json"
            write_published_catalog(catalog)
            (root / "curated").mkdir(parents=True)

            with self.assertRaisesRegex(manager.InstallError, "use --force"):
                manager.install_profile(root, "curated", catalog_path=catalog)

    def test_cross_host_winner_is_not_replaced_without_force(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "databases"
            catalog = Path(tmp) / "catalog.json"
            write_json(
                catalog,
                {
                    "schema_version": 1,
                    "default_profile": "curated",
                    "profiles": {
                        "curated": {
                            "version": "test-2",
                            "archive": {
                                "url": "https://example.org/curated.tar.gz",
                                "bytes": 1,
                                "sha256": sha256(b"x"),
                            },
                        }
                    },
                },
            )

            def extracted_profile(extracted):
                source = extracted / "curated"
                source.mkdir(parents=True)
                return source

            real_replace = manager.os.replace

            def competitor_after_staging(source, destination):
                result = real_replace(source, destination)
                if Path(destination).name == "ready":
                    competitor = root / "curated"
                    competitor.mkdir()
                    (competitor / "sentinel").write_text(
                        "concurrent install", encoding="utf-8"
                    )
                return result

            with mock.patch.object(manager, "download_archive"), \
                mock.patch.object(manager, "safe_extract_tar"), \
                mock.patch.object(
                    manager, "_find_extracted_profile", side_effect=extracted_profile
                ), \
                mock.patch.object(
                    manager,
                    "validate_profile_directory",
                    return_value={"version": "test-2"},
                ), \
                mock.patch.object(
                    manager.os, "replace", side_effect=competitor_after_staging
                ):
                with self.assertRaisesRegex(manager.InstallError, "appeared during publication"):
                    manager.install_profile(root, "curated", catalog_path=catalog)

            self.assertEqual(
                (root / "curated" / "sentinel").read_text(encoding="utf-8"),
                "concurrent install",
            )

    def test_install_failure_preserves_existing_profile_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "databases"
            old_profile = root / "curated"
            old_profile.mkdir(parents=True)
            sentinel = old_profile / "sentinel"
            sentinel.write_text("old installation", encoding="utf-8")

            artifact = b"corrupt"
            manifest = {
                "schema_version": 1,
                "profile": "curated",
                "version": "test-2",
                "artifacts": [
                    {
                        "path": "blast/ssu.fake",
                        "bytes": len(artifact),
                        "sha256": "0" * 64,
                    }
                ],
                "blast_databases": {"16S": {"prefix": "blast/ssu"}},
                "taxonomy_database": {
                    "preferred": "blast/ssu.fake",
                    "source_records": "blast/ssu.fake",
                },
            }
            archive = tar_bytes(
                {
                    "curated/manifest.json": json.dumps(manifest).encode(),
                    "curated/blast/ssu.fake": artifact,
                }
            )
            catalog = Path(tmp) / "catalog.json"
            write_json(
                catalog,
                {
                    "schema_version": 1,
                    "default_profile": "curated",
                    "profiles": {
                        "curated": {
                            "version": "test-2",
                            "archive": {
                                "url": "https://example.org/curated.tar.gz",
                                "bytes": len(archive),
                                "sha256": sha256(archive),
                            },
                        }
                    },
                },
            )

            def opener(url, timeout):
                self.assertEqual(url.full_url, "https://example.org/curated.tar.gz")
                self.assertEqual(timeout, 60)
                return FakeResponse(archive)

            with self.assertRaisesRegex(manager.IntegrityError, "SHA-256 mismatch"):
                manager.install_profile(
                    root,
                    "curated",
                    catalog_path=catalog,
                    replace=True,
                    opener=opener,
                )
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "old installation")
            self.assertEqual([path.name for path in root.iterdir()], ["curated"])

    @mock.patch.object(manager.subprocess, "run", side_effect=blast_ok)
    def test_successful_install_uses_durable_profile_publish(self, run) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "databases"
            artifact = b"valid-index"
            manifest = {
                "schema_version": 1,
                "profile": "curated",
                "version": "test-2",
                "artifacts": [
                    {
                        "path": "blast/ssu.fake",
                        "bytes": len(artifact),
                        "sha256": sha256(artifact),
                    }
                ],
                "blast_databases": {"16S": {"prefix": "blast/ssu"}},
                "taxonomy_database": {
                    "preferred": "blast/ssu.fake",
                    "source_records": "blast/ssu.fake",
                },
            }
            archive = tar_bytes(
                {
                    "curated/manifest.json": json.dumps(manifest).encode(),
                    "curated/blast/ssu.fake": artifact,
                }
            )
            catalog = Path(tmp) / "catalog.json"
            write_json(
                catalog,
                {
                    "schema_version": 1,
                    "default_profile": "curated",
                    "profiles": {
                        "curated": {
                            "version": "test-2",
                            "archive": {
                                "url": "https://example.org/curated.tar.gz",
                                "bytes": len(archive),
                                "sha256": sha256(archive),
                            },
                        }
                    },
                },
            )

            durable_publish = manager.replace_and_fsync
            with mock.patch.object(
                manager, "replace_and_fsync", wraps=durable_publish
            ) as replace_and_sync:
                installed = manager.install_profile(
                    root,
                    "curated",
                    catalog_path=catalog,
                    opener=lambda url, timeout: FakeResponse(archive),
                )

            self.assertEqual(installed, root / "curated")
            self.assertEqual(
                (installed / "blast" / "ssu.fake").read_bytes(), artifact
            )
            replace_and_sync.assert_called_once()

    @mock.patch.object(manager.subprocess, "run", side_effect=blast_ok)
    def test_manager_retry_resumes_retained_partial_and_removes_cache(self, run) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "databases"
            artifact = b"valid-index"
            manifest = {
                "schema_version": 1,
                "profile": "curated",
                "version": "test-2",
                "artifacts": [
                    {
                        "path": "blast/ssu.fake",
                        "bytes": len(artifact),
                        "sha256": sha256(artifact),
                    }
                ],
                "blast_databases": {"16S": {"prefix": "blast/ssu"}},
                "taxonomy_database": {
                    "preferred": "blast/ssu.fake",
                    "source_records": "blast/ssu.fake",
                },
            }
            archive = tar_bytes(
                {
                    "curated/manifest.json": json.dumps(manifest).encode(),
                    "curated/blast/ssu.fake": artifact,
                }
            )
            catalog = Path(tmp) / "catalog.json"
            write_json(
                catalog,
                {
                    "schema_version": 1,
                    "default_profile": "curated",
                    "profiles": {
                        "curated": {
                            "version": "test-2",
                            "archive": {
                                "url": "https://example.org/curated.tar.gz",
                                "bytes": len(archive),
                                "sha256": sha256(archive),
                            },
                        }
                    },
                },
            )
            first_ranges: list[str | None] = []

            def response_for(request, interrupt: bool):
                requested_range = request.get_header("Range")
                offset = (
                    int(requested_range.split("=")[1].split("-")[0])
                    if requested_range
                    else 0
                )
                body = archive[offset:]
                status = 206 if offset else 200
                headers = {"Content-Length": str(len(body))}
                if offset:
                    headers["Content-Range"] = (
                        f"bytes {offset}-{len(archive) - 1}/{len(archive)}"
                    )
                response_type = InterruptingHTTPResponse if interrupt else FakeHTTPResponse
                if interrupt:
                    return response_type(body, status, headers, 10)
                return response_type(body, status, headers)

            def interrupted_opener(request, timeout):
                first_ranges.append(request.get_header("Range"))
                return response_for(request, interrupt=True)

            with mock.patch.object(downloader.time, "sleep") as sleep:
                with self.assertRaisesRegex(manager.InstallError, "stopped after 5 attempts"):
                    manager.install_profile(
                        root,
                        "curated",
                        catalog_path=catalog,
                        opener=interrupted_opener,
                    )
            self.assertEqual(sleep.call_count, 4)
            self.assertEqual(
                first_ranges,
                [None]
                + [f"bytes={offset}-{len(archive) - 1}" for offset in (10, 20, 30, 40)],
            )
            partials = list(root.glob(".curated.download-*.part"))
            self.assertEqual(len(partials), 1)
            self.assertEqual(partials[0].stat().st_size, 50)
            self.assertFalse((root / "curated").exists())

            resumed_ranges: list[str | None] = []

            def resumed_opener(request, timeout):
                resumed_ranges.append(request.get_header("Range"))
                return response_for(request, interrupt=False)

            installed = manager.install_profile(
                root,
                "curated",
                catalog_path=catalog,
                opener=resumed_opener,
            )
            self.assertEqual((installed / "blast" / "ssu.fake").read_bytes(), artifact)
            self.assertEqual(resumed_ranges, [f"bytes=50-{len(archive) - 1}"])
            self.assertFalse(partials[0].exists())

    @mock.patch.object(manager.subprocess, "run", side_effect=blast_ok)
    def test_replacement_rename_failure_restores_previous_profile(self, run) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "databases"
            old_profile = root / "curated"
            old_profile.mkdir(parents=True)
            sentinel = old_profile / "sentinel"
            sentinel.write_text("old installation", encoding="utf-8")

            artifact = b"valid-index"
            manifest = {
                "schema_version": 1,
                "profile": "curated",
                "version": "test-2",
                "artifacts": [
                    {
                        "path": "blast/ssu.fake",
                        "bytes": len(artifact),
                        "sha256": sha256(artifact),
                    }
                ],
                "blast_databases": {"16S": {"prefix": "blast/ssu"}},
                "taxonomy_database": {
                    "preferred": "blast/ssu.fake",
                    "source_records": "blast/ssu.fake",
                },
            }
            archive = tar_bytes(
                {
                    "curated/manifest.json": json.dumps(manifest).encode(),
                    "curated/blast/ssu.fake": artifact,
                }
            )
            catalog = Path(tmp) / "catalog.json"
            write_json(
                catalog,
                {
                    "schema_version": 1,
                    "default_profile": "curated",
                    "profiles": {
                        "curated": {
                            "version": "test-2",
                            "archive": {
                                "url": "https://example.org/curated.tar.gz",
                                "bytes": len(archive),
                                "sha256": sha256(archive),
                            },
                        }
                    },
                },
            )

            real_replace = manager.os.replace
            call_count = 0

            def fail_new_profile_rename(source, destination):
                nonlocal call_count
                call_count += 1
                if call_count == 3:
                    raise OSError("simulated rename failure")
                return real_replace(source, destination)

            with mock.patch.object(manager.os, "replace", side_effect=fail_new_profile_rename):
                with self.assertRaisesRegex(manager.InstallError, "previous profile was restored"):
                    manager.install_profile(
                        root,
                        "curated",
                        catalog_path=catalog,
                        replace=True,
                        opener=lambda url, timeout: FakeResponse(archive),
                    )

            self.assertEqual(sentinel.read_text(encoding="utf-8"), "old installation")
            self.assertEqual([path.name for path in root.iterdir()], ["curated"])
            self.assertEqual(call_count, 4)

    @mock.patch.object(manager.subprocess, "run", side_effect=blast_ok)
    def test_post_rename_sync_failure_restores_previous_profile(self, run) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "databases"
            old_profile = root / "curated"
            old_profile.mkdir(parents=True)
            sentinel = old_profile / "sentinel"
            sentinel.write_text("old installation", encoding="utf-8")

            artifact = b"valid-index"
            manifest = {
                "schema_version": 1,
                "profile": "curated",
                "version": "test-2",
                "artifacts": [
                    {
                        "path": "blast/ssu.fake",
                        "bytes": len(artifact),
                        "sha256": sha256(artifact),
                    }
                ],
                "blast_databases": {"16S": {"prefix": "blast/ssu"}},
                "taxonomy_database": {
                    "preferred": "blast/ssu.fake",
                    "source_records": "blast/ssu.fake",
                },
            }
            archive = tar_bytes(
                {
                    "curated/manifest.json": json.dumps(manifest).encode(),
                    "curated/blast/ssu.fake": artifact,
                }
            )
            catalog = Path(tmp) / "catalog.json"
            write_json(
                catalog,
                {
                    "schema_version": 1,
                    "default_profile": "curated",
                    "profiles": {
                        "curated": {
                            "version": "test-2",
                            "archive": {
                                "url": "https://example.org/curated.tar.gz",
                                "bytes": len(archive),
                                "sha256": sha256(archive),
                            },
                        }
                    },
                },
            )

            durable_publish = manager.replace_and_fsync

            def publish_then_fail(source, destination):
                durable_publish(source, destination)
                raise OSError("simulated post-rename sync failure")

            with mock.patch.object(
                manager, "replace_and_fsync", side_effect=publish_then_fail
            ):
                with self.assertRaisesRegex(manager.InstallError, "previous profile was restored"):
                    manager.install_profile(
                        root,
                        "curated",
                        catalog_path=catalog,
                        replace=True,
                        opener=lambda url, timeout: FakeResponse(archive),
                    )

            self.assertEqual(sentinel.read_text(encoding="utf-8"), "old installation")
            self.assertEqual([path.name for path in root.iterdir()], ["curated"])



if __name__ == "__main__":
    unittest.main()
