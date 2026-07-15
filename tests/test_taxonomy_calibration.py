import json
import os
import sys
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path
from unittest import mock

import duckdb


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import calibrate_taxonomy as calibration
import classify_img_clusters as classifier


class CalibrationMathTests(unittest.TestCase):
    def test_wilson_lower_bound_matches_known_cases(self) -> None:
        self.assertAlmostEqual(calibration.wilson_lower_bound(100, 100), 0.9630065, places=6)
        self.assertAlmostEqual(calibration.wilson_lower_bound(98, 100), 0.9299882, places=6)
        self.assertEqual(calibration.wilson_lower_bound(0, 0), 0.0)

    def test_self_hit_is_additional_to_overflow_sentinel(self) -> None:
        def hit(subject: str) -> classifier.BlastHit:
            return classifier.BlastHit(
                query="query",
                subject=subject,
                percent_identity=Decimal("100"),
                alignment_length=100,
                query_length=100,
                subject_length=100,
                query_coverage=Decimal("100"),
                bit_score=Decimal("100"),
            )

        taxonomy = {
            f"ref-{index}": classifier.TaxonomyRecord(
                ("Bacteria", "Firmicutes"), "SILVA", "Bacteria", ""
            )
            for index in range(classifier.BLAST_MAX_TARGETS + 1)
        }
        exactly_at_policy = [hit("query")] + [
            hit(f"ref-{index}") for index in range(classifier.BLAST_MAX_TARGETS)
        ]
        overflow = exactly_at_policy + [hit(f"ref-{classifier.BLAST_MAX_TARGETS}")]

        self.assertEqual(
            calibration._prediction("query", exactly_at_policy, taxonomy),
            ("Bacteria", "Firmicutes"),
        )
        self.assertEqual(
            calibration._prediction("query", overflow, taxonomy),
            ("Bacteria",),
        )
        self.assertEqual(calibration.CALIBRATION_FETCH_TARGETS, 502)

    def test_low_coverage_raw_overflow_sentinel_still_backs_off(self) -> None:
        def hit(
            subject: str, query_coverage: Decimal = Decimal("100")
        ) -> classifier.BlastHit:
            return classifier.BlastHit(
                query="query",
                subject=subject,
                percent_identity=Decimal("100"),
                alignment_length=100,
                query_length=100,
                subject_length=100,
                query_coverage=query_coverage,
                bit_score=Decimal("100"),
            )

        taxonomy = {
            f"ref-{index}": classifier.TaxonomyRecord(
                ("Bacteria", "Firmicutes"), "SILVA", "Bacteria", ""
            )
            for index in range(classifier.BLAST_MAX_TARGETS)
        }
        raw_hits = [hit("query")] + [
            hit(f"ref-{index}") for index in range(classifier.BLAST_MAX_TARGETS)
        ] + [hit("low-coverage-sentinel", Decimal("79.99"))]

        self.assertEqual(
            calibration._prediction("query", raw_hits, taxonomy),
            ("Bacteria",),
        )

    def test_prediction_reports_missing_subject_taxonomy(self) -> None:
        hit = classifier.BlastHit(
            query="query",
            subject="missing-reference",
            percent_identity=Decimal("99"),
            alignment_length=100,
            query_length=100,
            subject_length=100,
            query_coverage=Decimal("100"),
            bit_score=Decimal("100"),
        )
        with self.assertRaisesRegex(
            RuntimeError, "taxonomy is missing.*missing-reference"
        ):
            calibration._prediction("query", [hit], {})


class CalibrationSelectionTests(unittest.TestCase):
    def test_parquet_selection_filters_and_samples_marker_strata(self) -> None:
        preferred_rows = [
            ("a", "Bacteria;A", "SILVA", "Bacteria", False),
            ("b", "Bacteria;B", "SILVA", "Bacteria", False),
            ("multi", "Eukaryota;SAR", "PR2", "Eukaryota", False),
            ("tokens", "Archaea;A", "SILVA", "Archaea", False),
            ("ambiguous", "Bacteria;X", "SILVA", "Bacteria", True),
            ("blank", "   ", "SILVA", "Bacteria", False),
            ("empty-markers", "Bacteria;Y", "SILVA", "Bacteria", False),
        ]
        sequence_rows = [
            ("a", "16S"),
            ("b", "16S"),
            ("multi", "16S;18S"),
            ("tokens", "16S;;junk;18S;16S"),
            ("ambiguous", "16S"),
            ("blank", "18S"),
            ("empty-markers", ";;junk;"),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            profile = Path(temporary)
            tables = profile / "tables"
            tables.mkdir()
            connection = duckdb.connect(":memory:")
            try:
                connection.execute(
                    """
                    CREATE TABLE preferred_taxonomy (
                        sequence_id VARCHAR,
                        taxonomy VARCHAR,
                        taxonomy_source VARCHAR,
                        domain VARCHAR,
                        cross_domain_conflict BOOLEAN
                    )
                    """
                )
                connection.executemany(
                    "INSERT INTO preferred_taxonomy VALUES (?, ?, ?, ?, ?)",
                    preferred_rows,
                )
                connection.execute(
                    "CREATE TABLE sequences (sequence_id VARCHAR, markers VARCHAR)"
                )
                connection.executemany(
                    "INSERT INTO sequences VALUES (?, ?)", sequence_rows
                )
                connection.execute(
                    "COPY preferred_taxonomy TO ? (FORMAT PARQUET)",
                    [str(tables / "preferred_taxonomy.parquet")],
                )
                connection.execute(
                    "COPY sequences TO ? (FORMAT PARQUET)",
                    [str(tables / "sequences.parquet")],
                )
            finally:
                connection.close()

            with (
                mock.patch.object(calibration, "MIN_CALIBRATION_CALLS", 1),
                mock.patch.dict(
                    os.environ,
                    {"SLURM_CPUS_PER_TASK": "", "OMP_NUM_THREADS": ""},
                ),
            ):
                rows = calibration.select_calibration_rows(profile, 1)

        observed = [
            (
                row["sequence_id"],
                row["marker"],
                row["taxonomy_source"],
                row["domain"],
            )
            for row in rows
        ]
        self.assertEqual(
            observed,
            [
                ("multi", "16S", "PR2", "Eukaryota"),
                ("tokens", "16S", "SILVA", "Archaea"),
                ("a", "16S", "SILVA", "Bacteria"),
                ("multi", "18S", "PR2", "Eukaryota"),
                ("tokens", "18S", "SILVA", "Archaea"),
            ],
        )


class CalibrationSchemaTests(unittest.TestCase):
    @staticmethod
    def _hit(query: str, subject: str) -> classifier.BlastHit:
        return classifier.BlastHit(
            query=query,
            subject=subject,
            percent_identity=Decimal("99"),
            alignment_length=100,
            query_length=100,
            subject_length=100,
            query_coverage=Decimal("100"),
            bit_score=Decimal("100"),
        )

    def test_failed_stratum_is_recorded_while_passing_stratum_sets_cap(self) -> None:
        passing_rows = [
            {
                "sequence_id": f"pass-{index}",
                "taxonomy": "Bacteria;Firmicutes",
                "taxonomy_source": "SILVA",
                "domain": "Bacteria",
                "marker": "16S",
            }
            for index in range(100)
        ]
        failing_rows = [
            {
                "sequence_id": f"fail-{index}",
                "taxonomy": "Eukaryota;SAR",
                "taxonomy_source": "PR2",
                "domain": "Eukaryota",
                "marker": "16S",
            }
            for index in range(100)
        ]
        hits = {
            row["sequence_id"]: (self._hit(row["sequence_id"], "pass-reference"),)
            for row in passing_rows
        }
        hits.update(
            {
                row["sequence_id"]: (self._hit(row["sequence_id"], "bad-reference"),)
                for row in failing_rows
            }
        )
        taxonomy = {
            "pass-reference": classifier.TaxonomyRecord(
                ("Bacteria", "Firmicutes"), "SILVA", "Bacteria", ""
            ),
            "bad-reference": classifier.TaxonomyRecord(
                ("Bacteria", "Firmicutes"), "SILVA", "Bacteria", ""
            ),
        }

        with tempfile.TemporaryDirectory() as temporary:
            blast_output = Path(temporary) / "leave_one_out.m8"
            blast_output.touch()
            with (
                mock.patch.object(classifier, "parse_blast_hits", return_value=hits),
                mock.patch.object(classifier, "load_taxonomy", return_value=taxonomy),
            ):
                result = calibration.evaluate_calibration(
                    temporary,
                    passing_rows + failing_rows,
                    {"16S": blast_output},
                    samples_per_stratum_requested=777,
                )

        passing_key = "16S|SILVA|Bacteria"
        failing_key = "16S|PR2|Eukaryota"
        self.assertEqual(result["schema_version"], 2)
        self.assertEqual(result["samples_per_stratum_requested"], 777)
        self.assertEqual(result["rank_caps"], {passing_key: 1})
        self.assertEqual(result["strata"][passing_key]["status"], "calibrated")
        self.assertEqual(result["strata"][passing_key]["rank_cap"], 1)
        self.assertEqual(result["strata"][passing_key]["reason"], "")
        self.assertEqual(result["strata"][failing_key]["status"], "failed")
        self.assertIsNone(result["strata"][failing_key]["rank_cap"])
        self.assertEqual(
            result["strata"][failing_key]["reason"],
            "domain_precision_below_threshold",
        )
        self.assertEqual(result["strata"][failing_key]["metrics"][0]["called"], 100)

    def test_no_rows_and_no_passing_strata_are_fatal(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "selected no rows"):
            calibration.evaluate_calibration("unused", [], {})

        row = {
            "sequence_id": "fail",
            "taxonomy": "Eukaryota;SAR",
            "taxonomy_source": "PR2",
            "domain": "Eukaryota",
            "marker": "18S",
        }
        with tempfile.TemporaryDirectory() as temporary:
            blast_output = Path(temporary) / "leave_one_out.m8"
            blast_output.touch()
            with (
                mock.patch.object(classifier, "parse_blast_hits", return_value={}),
                mock.patch.object(classifier, "load_taxonomy", return_value={}),
                self.assertRaisesRegex(RuntimeError, "no strata passed"),
            ):
                calibration.evaluate_calibration(
                    temporary, [row], {"18S": blast_output}
                )


class CalibrationReuseTests(unittest.TestCase):
    @staticmethod
    def _row() -> dict[str, str]:
        return {
            "sequence_id": "query-1",
            "taxonomy": "Bacteria;Firmicutes",
            "taxonomy_source": "SILVA",
            "domain": "Bacteria",
            "marker": "16S",
        }

    @staticmethod
    def _m8(query: str = "query-1", subject: str = "query-1") -> str:
        return f"{query}\t{subject}\t100\t4\t4\t4\t100\t8\n"

    @staticmethod
    def _write_manifest(profile: Path, version: str = "test-v1") -> None:
        profile.mkdir(parents=True, exist_ok=True)
        (profile / calibration.manager.MANIFEST_NAME).write_text(
            json.dumps({"version": version}) + "\n", encoding="utf-8"
        )

    def _prepare_retained(
        self,
        profile: Path,
        output: Path,
        row: dict[str, str],
        *,
        provenance: bool = True,
    ) -> Path:
        self._write_manifest(profile)
        marker_directory = output / "16S"
        marker_directory.mkdir(parents=True, exist_ok=True)
        ids = marker_directory / "query_ids.txt"
        ids.write_text(
            calibration._query_ids_text([row]), encoding="ascii"
        )
        truth = marker_directory / "truth.tsv"
        truth.write_text(
            calibration._truth_tsv_text([row]), encoding="utf-8"
        )
        fasta = marker_directory / "queries.fna"
        fasta.write_text(">query-1\nACGT\n", encoding="ascii")
        blast_output = marker_directory / "leave_one_out.m8"
        blast_output.write_text(self._m8(), encoding="ascii")
        if provenance:
            calibration._write_search_provenance(
                profile,
                "16S",
                marker_directory,
                ids,
                truth,
                fasta,
                blast_output,
                8,
            )
        else:
            (marker_directory / calibration.SEARCH_PROVENANCE_NAME).unlink(
                missing_ok=True
            )
        return blast_output

    def test_fresh_search_writes_provenance_then_reuse_skips_blastn(self) -> None:
        row = self._row()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profile = root / "profile"
            output = root / "output"
            self._write_manifest(profile)

            def search(command: list[str]) -> None:
                destination = Path(command[command.index("-out") + 1])
                if command[0] == "blastdbcmd":
                    destination.write_text(">query-1\nACGT\n", encoding="ascii")
                else:
                    self.assertEqual(command[0], "blastn")
                    destination.write_text(self._m8(), encoding="ascii")

            with mock.patch.object(calibration, "_run", side_effect=search) as runner:
                outputs = calibration._write_queries(
                    profile, output, [row], threads=8
                )

            retained = output / "16S" / "leave_one_out.m8"
            self.assertEqual(outputs, {"16S": retained})
            self.assertEqual(runner.call_count, 2)
            provenance_path = output / "16S" / calibration.SEARCH_PROVENANCE_NAME
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            self.assertEqual(provenance["schema_version"], 1)
            self.assertEqual(provenance["status"], "complete")
            self.assertEqual(
                set(provenance["files"]),
                {"query_ids.txt", "truth.tsv", "queries.fna", "leave_one_out.m8"},
            )
            self.assertEqual(set(provenance["commands"]), {"blastdbcmd", "blastn"})

            def regenerate(command: list[str]) -> None:
                self.assertEqual(command[0], "blastdbcmd")
                destination = Path(command[command.index("-out") + 1])
                destination.write_text(">query-1\nACGT\n", encoding="ascii")

            with mock.patch.object(calibration, "_run", side_effect=regenerate) as runner:
                outputs = calibration._write_queries(
                    profile,
                    output,
                    [row],
                    threads=8,
                    reuse_existing_blast=True,
                )

            self.assertEqual(outputs, {"16S": retained})
            runner.assert_called_once()

    def test_fresh_search_rejects_incomplete_results_before_provenance(self) -> None:
        rows = [
            self._row(),
            {
                **self._row(),
                "sequence_id": "query-2",
            },
        ]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profile = root / "profile"
            output = root / "output"
            self._write_manifest(profile)

            def incomplete_search(command: list[str]) -> None:
                destination = Path(command[command.index("-out") + 1])
                if command[0] == "blastdbcmd":
                    destination.write_text(
                        ">query-1\nACGT\n>query-2\nACGT\n", encoding="ascii"
                    )
                else:
                    destination.write_text(self._m8(), encoding="ascii")

            with (
                mock.patch.object(
                    calibration, "_run", side_effect=incomplete_search
                ),
                self.assertRaisesRegex(RuntimeError, "query coverage differs"),
            ):
                calibration._write_queries(profile, output, rows, threads=8)

            self.assertFalse(
                (output / "16S" / calibration.SEARCH_PROVENANCE_NAME).exists()
            )

    def test_failed_fresh_search_leaves_no_complete_provenance(self) -> None:
        row = self._row()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profile = root / "profile"
            output = root / "output"
            self._prepare_retained(profile, output, row)

            def fail_blastn(command: list[str]) -> None:
                if command[0] == "blastdbcmd":
                    destination = Path(command[command.index("-out") + 1])
                    destination.write_text(">query-1\nACGT\n", encoding="ascii")
                    return
                raise RuntimeError("synthetic blastn failure")

            with (
                mock.patch.object(calibration, "_run", side_effect=fail_blastn),
                self.assertRaisesRegex(RuntimeError, "synthetic blastn failure"),
            ):
                calibration._write_queries(profile, output, [row], threads=8)
            self.assertFalse(
                (output / "16S" / calibration.SEARCH_PROVENANCE_NAME).exists()
            )

    def test_reuse_rejects_missing_provenance_profile_drift_and_tampering(self) -> None:
        row = self._row()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profile = root / "profile"
            output = root / "output"
            self._prepare_retained(
                profile, output, row, provenance=False
            )
            with self.assertRaisesRegex(RuntimeError, "search_provenance.json"):
                calibration._write_queries(
                    profile,
                    output,
                    [row],
                    threads=8,
                    reuse_existing_blast=True,
                )

            self._prepare_retained(profile, output, row)
            self._write_manifest(profile, "test-v2")
            with self.assertRaisesRegex(RuntimeError, "profile manifest has changed"):
                calibration._write_queries(
                    profile,
                    output,
                    [row],
                    threads=8,
                    reuse_existing_blast=True,
                )

            self._prepare_retained(profile, output, row)
            blast_output = output / "16S" / "leave_one_out.m8"
            blast_output.write_text(
                blast_output.read_text(encoding="ascii") + self._m8("query-1", "other"),
                encoding="ascii",
            )
            with (
                mock.patch.object(classifier, "parse_blast_hits") as parser,
                self.assertRaisesRegex(RuntimeError, "hash mismatch"),
            ):
                calibration._write_queries(
                    profile,
                    output,
                    [row],
                    threads=8,
                    reuse_existing_blast=True,
                )
            parser.assert_not_called()

    def test_reuse_rechecks_truth_and_exact_query_coverage(self) -> None:
        row = self._row()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profile = root / "profile"
            output = root / "output"
            self._prepare_retained(profile, output, row)
            changed_truth = dict(row, taxonomy="Bacteria;Actinobacteriota")
            with self.assertRaisesRegex(RuntimeError, "truth.tsv"):
                calibration._write_queries(
                    profile,
                    output,
                    [changed_truth],
                    threads=8,
                    reuse_existing_blast=True,
                )

            blast_output = self._prepare_retained(profile, output, row)
            blast_output.write_text(
                self._m8(query="other", subject="other"), encoding="ascii"
            )
            marker_directory = output / "16S"
            calibration._write_search_provenance(
                profile,
                "16S",
                marker_directory,
                marker_directory / "query_ids.txt",
                marker_directory / "truth.tsv",
                marker_directory / "queries.fna",
                blast_output,
                8,
            )

            def regenerate(command: list[str]) -> None:
                destination = Path(command[command.index("-out") + 1])
                destination.write_text(">query-1\nACGT\n", encoding="ascii")

            with (
                mock.patch.object(calibration, "_run", side_effect=regenerate),
                self.assertRaisesRegex(RuntimeError, "query coverage differs"),
            ):
                calibration._write_queries(
                    profile,
                    output,
                    [row],
                    threads=8,
                    reuse_existing_blast=True,
                )


if __name__ == "__main__":
    unittest.main()
