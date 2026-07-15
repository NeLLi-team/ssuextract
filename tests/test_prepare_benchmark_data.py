import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import prepare_benchmark_data as benchmark


class PrepareBenchmarkDataTests(unittest.TestCase):
    def write_complete_timings(self, path: Path) -> None:
        path.write_text(
            "profile\ttrial\telapsed_seconds\tpeak_rss_kb\n"
            + "".join(
                f"{profile}\t{trial}\t{10 + trial}.5\t{1000 + trial}\n"
                for profile in benchmark.PROFILE_LABELS
                for trial in range(4)
            ),
            encoding="utf-8",
        )

    def test_timing_contract_requires_all_profiles_and_trials(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            timing = Path(temporary) / "timing.tsv"
            self.write_complete_timings(timing)
            rows = benchmark.load_timings(timing)
            self.assertEqual(len(rows), 12)
            self.assertEqual(sum(row["warmup"] == "true" for row in rows), 3)

            timing.write_text(
                "profile\ttrial\telapsed_seconds\tpeak_rss_kb\n"
                "curated\t0\t1\t1\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "incomplete"):
                benchmark.load_timings(timing)

    def test_timing_contract_rejects_invalid_rows(self) -> None:
        cases = [
            ("invalid schema", "profile\ttrial\telapsed_seconds\ncurated\t0\t1\n"),
            ("unexpected profile or trial", "profile\ttrial\telapsed_seconds\tpeak_rss_kb\nunknown\t0\t1\t1\n"),
            ("unexpected profile or trial", "profile\ttrial\telapsed_seconds\tpeak_rss_kb\ncurated\t4\t1\t1\n"),
            ("invalid benchmark value", "profile\ttrial\telapsed_seconds\tpeak_rss_kb\ncurated\t0\t0\t1\n"),
            ("invalid benchmark value", "profile\ttrial\telapsed_seconds\tpeak_rss_kb\ncurated\t0\tnan\t1\n"),
            ("invalid benchmark value", "profile\ttrial\telapsed_seconds\tpeak_rss_kb\ncurated\t0\tinf\t1\n"),
            ("invalid benchmark value", "profile\ttrial\telapsed_seconds\tpeak_rss_kb\ncurated\t0\t-inf\t1\n"),
            ("duplicate", "profile\ttrial\telapsed_seconds\tpeak_rss_kb\ncurated\t0\t1\t1\ncurated\t0\t1\t1\n"),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            timing = Path(temporary) / "timing.tsv"
            for expected_error, contents in cases:
                with self.subTest(expected_error=expected_error):
                    timing.write_text(contents, encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, expected_error):
                        benchmark.load_timings(timing)

    def test_manifest_and_archive_validation_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profile = root / "profile"
            profile.mkdir()
            with self.assertRaisesRegex(ValueError, "invalid database manifest"):
                benchmark._validate_profile_identity(profile, "curated")
            (profile / "manifest.json").write_text(
                json.dumps({"profile": "img"}), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "profile mismatch"):
                benchmark._validate_profile_identity(profile, "curated")
            with self.assertRaisesRegex(ValueError, "archive is missing or empty"):
                benchmark.storage_row("Curated", profile, root / "missing.tar.zst")

    def test_storage_components_use_one_all_files_rule(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            profile = Path(temporary)
            (profile / "blast").mkdir()
            (profile / "tables").mkdir()
            (profile / "sequence.fna").write_bytes(b"fasta")
            (profile / "blast" / "16S.nsq").write_bytes(b"blast")
            (profile / "blast" / "taxonomy.parquet").write_bytes(b"PAR1")
            (profile / "manifest.json").write_bytes(b"{}")
            components = benchmark._storage_components(profile)
            self.assertEqual(components["source_fasta_bytes"], 5)
            self.assertEqual(components["blast_bytes"], 5)
            self.assertEqual(components["parquet_bytes"], 4)
            self.assertEqual(components["other_bytes"], 2)
            self.assertEqual(components["installed_bytes"], 16)

    def test_main_writes_validated_three_profile_tables(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            benchmark_directory = root / "benchmark"
            benchmark_directory.mkdir()
            (benchmark_directory / "timing.tsv").write_text(
                "profile\ttrial\telapsed_seconds\tpeak_rss_kb\n"
                + "".join(
                    f"{profile}\t{trial}\t{20 + trial}.25\t{2000 + trial}\n"
                    for profile in benchmark.PROFILE_LABELS
                    for trial in range(4)
                ),
                encoding="utf-8",
            )
            legacy = root / "legacy"
            legacy.mkdir()
            (legacy / "old.fasta").write_text(">x\nACGT\n", encoding="ascii")
            (legacy / "old.nsq").write_bytes(b"blast")
            (legacy / "index.html").write_text("not a database file", encoding="utf-8")

            profiles: dict[str, Path] = {}
            archives: dict[str, Path] = {}
            for profile in ("curated", "img"):
                directory = root / profile
                (directory / "blast").mkdir(parents=True)
                (directory / "tables").mkdir()
                (directory / "manifest.json").write_text(
                    json.dumps({"profile": profile}), encoding="utf-8"
                )
                (directory / "blast" / "16S.nsq").write_bytes(b"blast-index")
                (directory / "tables" / "taxonomy.parquet").write_bytes(b"PAR1")
                archive = root / f"{profile}.tar.zst"
                archive.write_bytes(profile.encode("ascii"))
                profiles[profile] = directory
                archives[profile] = archive

            output = root / "docs-data"
            result = benchmark.main(
                [
                    "--benchmark-directory",
                    str(benchmark_directory),
                    "--legacy-root",
                    str(legacy),
                    "--curated-profile",
                    str(profiles["curated"]),
                    "--curated-archive",
                    str(archives["curated"]),
                    "--img-profile",
                    str(profiles["img"]),
                    "--img-archive",
                    str(archives["img"]),
                    "--output-directory",
                    str(output),
                ]
            )
            self.assertEqual(result, 0)
            with (output / "database_profile_benchmark.tsv").open(newline="") as handle:
                timing_rows = list(csv.DictReader(handle, delimiter="\t"))
            with (output / "database_storage.tsv").open(newline="") as handle:
                storage_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(timing_rows), 12)
            self.assertEqual(len(storage_rows), 3)
            legacy_row = storage_rows[0]
            self.assertEqual(int(legacy_row["source_fasta_bytes"]), 8)
            self.assertEqual(int(legacy_row["blast_bytes"]), 5)
            self.assertEqual(int(legacy_row["other_bytes"]), 19)
            self.assertEqual(int(legacy_row["installed_bytes"]), 32)
            self.assertEqual(legacy_row["archive_sha256"], "")
            self.assertEqual(len(storage_rows[1]["archive_sha256"]), 64)
            self.assertEqual(storage_rows[1]["parquet_bytes"], "4")


if __name__ == "__main__":
    unittest.main()
