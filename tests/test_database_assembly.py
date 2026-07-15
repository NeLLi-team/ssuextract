import hashlib
import json
import stat
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import assemble_database_profile as assembler
import build_database_release as builder
import database_manager as manager


def tiny_model() -> builder.DatabaseModel:
    records = [
        builder.PreparedSourceRecord(
            "SILVA",
            "138.2",
            "silva-16",
            "silva-16",
            "AAAA",
            "16S",
            ("Bacteria", "Firmicutes"),
            "SILVA",
        ),
        builder.PreparedSourceRecord(
            "PR2",
            "5.1.1",
            "pr2-plastid",
            "pr2-plastid",
            "ACGTACGT",
            "16S",
            ("Eukaryota", "Archaeplastida"),
            "PR2",
            "plastid",
        ),
        builder.PreparedSourceRecord(
            "PR2",
            "5.1.1",
            "pr2-nucleus",
            "pr2-nucleus",
            "ACGTACGT",
            "18S",
            ("Eukaryota", "Archaeplastida"),
            "PR2",
            "nucleus",
        ),
        builder.PreparedSourceRecord(
            "PR2",
            "5.1.1",
            "pr2-18",
            "pr2-18",
            "CCCC",
            "18S",
            ("Eukaryota", "TSAR"),
            "PR2",
            "nucleus",
        ),
    ]
    return builder.build_deduplicated_model(records)


class DatabaseAssemblyTests(unittest.TestCase):
    def test_assembles_marker_indexes_tables_manifest_and_validates(self) -> None:
        model = tiny_model()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "database"
            profile = root / "curated"
            result = assembler.assemble_database_profile(
                model,
                profile,
                profile="curated",
                version="2026.1",
                provenance_details={"taxonomy_policy_version": "1.0.0"},
            )

            manifest = manager.validate_profile(root, "curated")
            self.assertEqual(result.manifest, manifest)
            self.assertEqual(stat.S_IMODE(profile.stat().st_mode), 0o755)
            self.assertEqual(
                stat.S_IMODE((profile / "manifest.json").stat().st_mode), 0o644
            )
            self.assertEqual(set(manifest["blast_databases"]), {"16S", "18S"})
            self.assertEqual(
                manifest["taxonomy_database"]["preferred"],
                "tables/preferred_taxonomy.parquet",
            )

            source_members = {
                marker: {
                    source.sequence_id
                    for source in model.source_records
                    if source.marker == marker
                }
                for marker in ("16S", "18S")
            }
            for marker, expected in source_members.items():
                output = subprocess.run(
                    [
                        "blastdbcmd",
                        "-db",
                        str(profile / "blast" / marker),
                        "-entry",
                        "all",
                        "-outfmt",
                        "%a",
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(set(output.stdout.splitlines()), expected)
            self.assertFalse((profile / "fasta").exists())

            artifact_paths = [artifact["path"] for artifact in manifest["artifacts"]]
            self.assertEqual(artifact_paths, sorted(artifact_paths))
            self.assertIn("provenance.json", artifact_paths)
            self.assertIn("tables/preferred_taxonomy.parquet", artifact_paths)
            for artifact in manifest["artifacts"]:
                path = profile / artifact["path"]
                self.assertEqual(path.stat().st_size, artifact["bytes"])
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), artifact["sha256"])

            provenance = json.loads((profile / "provenance.json").read_text(encoding="utf-8"))
            self.assertEqual(provenance["counts"]["sequences"], 3)
            self.assertEqual(provenance["markers"]["16S"]["sequences"], 2)
            self.assertEqual(provenance["markers"]["18S"]["sequences"], 2)
            self.assertEqual(
                provenance["release_build"]["taxonomy_policy_version"], "1.0.0"
            )
            self.assertEqual(
                {(row["name"], row["version"]) for row in provenance["sources"]},
                {("PR2", "5.1.1"), ("SILVA", "138.2")},
            )

    def test_package_is_byte_deterministic_and_has_one_top_level_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "curated"
            assembler.assemble_database_profile(
                tiny_model(), profile, profile="curated", version="2026.1"
            )
            first = assembler.package_profile(profile, root / "first.tar.zst")
            second = assembler.package_profile(profile, root / "second.tar.zst")
            self.assertEqual(first.bytes, second.bytes)
            self.assertEqual(first.sha256, second.sha256)

            extracted = root / "extracted"
            manager.safe_extract_tar(first.path, extracted)
            self.assertEqual([path.name for path in extracted.iterdir()], ["curated"])
            manager.validate_profile(extracted, "curated")
            with tarfile.open(first.path, mode="r:*") as archive:
                roots = {Path(member.name).parts[0] for member in archive.getmembers()}
            self.assertEqual(roots, {"curated"})

    def test_rejects_mismatched_profile_directory_without_writing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "wrong-name"
            with self.assertRaisesRegex(assembler.AssemblyError, "must match profile"):
                assembler.assemble_database_profile(
                    tiny_model(), target, profile="curated", version="2026.1"
                )
            self.assertFalse(target.exists())

    def test_rejects_archive_inside_profile_tree(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            profile = Path(tmp) / "curated"
            profile.mkdir()
            with self.assertRaisesRegex(assembler.AssemblyError, "outside"):
                assembler.package_profile(profile, profile / "curated.tar.zst")

    def test_archive_failure_publishes_neither_profile_nor_archive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "profiles" / "curated"
            archive = root / "archives" / "curated.tar.zst"
            with mock.patch.object(
                assembler,
                "_stage_profile_archive",
                side_effect=assembler.AssemblyError("injected compression failure"),
            ):
                with self.assertRaisesRegex(assembler.AssemblyError, "injected"):
                    assembler.assemble_database_profile(
                        tiny_model(),
                        profile,
                        profile="curated",
                        version="2026.1",
                        archive_path=archive,
                    )
            self.assertFalse(profile.exists())
            self.assertFalse(archive.exists())

    def test_interrupt_during_archive_staging_publishes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "profiles" / "curated"
            archive = root / "archives" / "curated.tar.zst"
            with mock.patch.object(
                assembler, "_stage_profile_archive", side_effect=KeyboardInterrupt()
            ):
                with self.assertRaises(KeyboardInterrupt):
                    assembler.assemble_database_profile(
                        tiny_model(),
                        profile,
                        profile="curated",
                        version="2026.1",
                        archive_path=archive,
                    )
            self.assertFalse(profile.exists())
            self.assertFalse(archive.exists())

    def test_archive_publish_failure_rolls_back_profile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "profiles" / "curated"
            archive = root / "archives" / "curated.tar.zst"
            with mock.patch.object(
                assembler.os, "link", side_effect=OSError("injected publish failure")
            ):
                with self.assertRaisesRegex(OSError, "injected"):
                    assembler.assemble_database_profile(
                        tiny_model(),
                        profile,
                        profile="curated",
                        version="2026.1",
                        archive_path=archive,
                    )
            self.assertFalse(profile.exists())
            self.assertFalse(archive.exists())

    def test_profile_post_rename_sync_failure_rolls_back_for_retry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "profiles" / "curated"
            archive = root / "archives" / "curated.tar.zst"
            durable_publish = assembler.replace_and_fsync

            def publish_then_fail(source: Path, destination: Path) -> None:
                durable_publish(source, destination)
                if Path(destination).resolve() == profile.resolve():
                    raise OSError("injected post-rename profile sync failure")

            with (
                mock.patch.object(
                    assembler,
                    "replace_and_fsync",
                    side_effect=publish_then_fail,
                ),
                self.assertRaisesRegex(OSError, "post-rename profile sync failure"),
            ):
                assembler.assemble_database_profile(
                    tiny_model(),
                    profile,
                    profile="curated",
                    version="2026.1",
                    archive_path=archive,
                )

            self.assertFalse(profile.exists())
            self.assertFalse(archive.exists())

            result = assembler.assemble_database_profile(
                tiny_model(),
                profile,
                profile="curated",
                version="2026.1",
                archive_path=archive,
            )
            self.assertEqual(result.profile_directory, profile)
            self.assertEqual(result.archive.path, archive)
            manager.validate_profile(profile.parent, "curated")

    def test_post_link_fsync_failure_rolls_back_profile_and_archive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "profiles" / "curated"
            archive = root / "archives" / "curated.tar.zst"
            real_fsync_directory = assembler.fsync_directory

            def fail_once_link_is_visible(path: Path) -> None:
                if Path(path).resolve() == archive.parent.resolve() and archive.exists():
                    raise OSError("injected post-link sync failure")
                real_fsync_directory(path)

            with (
                mock.patch.object(
                    assembler,
                    "fsync_directory",
                    side_effect=fail_once_link_is_visible,
                ),
                self.assertRaisesRegex(OSError, "post-link sync failure"),
            ):
                assembler.assemble_database_profile(
                    tiny_model(),
                    profile,
                    profile="curated",
                    version="2026.1",
                    archive_path=archive,
                )

            self.assertFalse(profile.exists())
            self.assertFalse(archive.exists())

    def test_repeated_archive_fsync_failure_still_attempts_profile_rollback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "profiles" / "curated"
            archive = root / "archives" / "curated.tar.zst"
            real_fsync_directory = assembler.fsync_directory

            def fail_for_archive_directory(path: Path) -> None:
                if Path(path).resolve() == archive.parent.resolve():
                    raise OSError("injected repeated archive sync failure")
                real_fsync_directory(path)

            with (
                mock.patch.object(
                    assembler,
                    "fsync_directory",
                    side_effect=fail_for_archive_directory,
                ),
                self.assertRaisesRegex(OSError, "repeated archive sync failure") as raised,
            ):
                assembler.assemble_database_profile(
                    tiny_model(),
                    profile,
                    profile="curated",
                    version="2026.1",
                    archive_path=archive,
                )

            self.assertFalse(profile.exists())
            self.assertFalse(archive.exists())
            self.assertTrue(
                any(
                    "release rollback also failed" in message
                    for message in raised.exception.rollback_errors
                )
            )

    def test_archive_unlink_failure_still_attempts_profile_rollback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "profiles" / "curated"
            archive = root / "archives" / "curated.tar.zst"
            real_fsync_directory = assembler.fsync_directory
            real_unlink = Path.unlink

            def fail_once_link_is_visible(path: Path) -> None:
                if Path(path).resolve() == archive.parent.resolve() and archive.exists():
                    raise OSError("injected archive publication sync failure")
                real_fsync_directory(path)

            def fail_archive_unlink(path: Path, *args: object, **kwargs: object) -> None:
                if path.resolve() == archive.resolve():
                    raise OSError("injected archive unlink failure")
                real_unlink(path, *args, **kwargs)

            with (
                mock.patch.object(
                    assembler,
                    "fsync_directory",
                    side_effect=fail_once_link_is_visible,
                ),
                mock.patch.object(Path, "unlink", new=fail_archive_unlink),
                self.assertRaisesRegex(
                    OSError, "archive publication sync failure"
                ) as raised,
            ):
                assembler.assemble_database_profile(
                    tiny_model(),
                    profile,
                    profile="curated",
                    version="2026.1",
                    archive_path=archive,
                )

            self.assertFalse(profile.exists())
            self.assertTrue(archive.exists())
            self.assertTrue(
                any(
                    "archive unlink failure" in message
                    for message in raised.exception.rollback_errors
                )
            )

    def test_package_profile_post_link_failure_removes_archive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "curated"
            profile.mkdir()
            (profile / "manifest.json").write_text("{}\n", encoding="utf-8")
            archive = root / "curated.tar.zst"
            real_fsync_directory = assembler.fsync_directory

            def fail_once_link_is_visible(path: Path) -> None:
                if Path(path).resolve() == root.resolve() and archive.exists():
                    raise OSError("injected package sync failure")
                real_fsync_directory(path)

            with (
                mock.patch.object(
                    assembler,
                    "fsync_directory",
                    side_effect=fail_once_link_is_visible,
                ),
                self.assertRaisesRegex(OSError, "package sync failure"),
            ):
                assembler.package_profile(profile, archive)

            self.assertFalse(archive.exists())

    def test_existing_archive_is_never_replaced(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "profiles" / "curated"
            archive = root / "archives" / "curated.tar.zst"
            archive.parent.mkdir()
            archive.write_bytes(b"existing")
            with self.assertRaisesRegex(assembler.AssemblyError, "already exists"):
                assembler.assemble_database_profile(
                    tiny_model(),
                    profile,
                    profile="curated",
                    version="2026.1",
                    archive_path=archive,
                )
            self.assertFalse(profile.exists())
            self.assertEqual(archive.read_bytes(), b"existing")

    def test_copies_checksum_covered_release_files_and_rejects_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            notice = root / "NOTICE.txt"
            notice.write_text("source notice\n", encoding="utf-8")
            evidence = root / "evidence.jsonl"
            evidence.write_text('{"record_type":"catalog"}\n', encoding="utf-8")
            result = assembler.assemble_database_profile(
                tiny_model(),
                root / "curated",
                profile="curated",
                version="2026.1",
                release_files={
                    "LICENSES/NOTICE.txt": notice,
                    "EVIDENCE/img_taxonomy_evidence.jsonl": evidence,
                },
            )
            self.assertEqual(
                (root / "curated" / "LICENSES" / "NOTICE.txt").read_text(),
                "source notice\n",
            )
            self.assertIn(
                "LICENSES/NOTICE.txt",
                {row["path"] for row in result.manifest["artifacts"]},
            )
            self.assertEqual(
                result.manifest["taxonomy_evidence_catalog"],
                "EVIDENCE/img_taxonomy_evidence.jsonl",
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            notice = root / "NOTICE.txt"
            notice.write_text("source notice\n", encoding="utf-8")
            with self.assertRaisesRegex(assembler.AssemblyError, "Unsafe supplemental"):
                assembler.assemble_database_profile(
                    tiny_model(),
                    root / "curated",
                    profile="curated",
                    version="2026.1",
                    release_files={"../NOTICE.txt": notice},
                )


if __name__ == "__main__":
    unittest.main()
