import ast
import csv
import importlib.util
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from annotate_hits import SUMMARY_FIELDS
from tree_schema import SUMMARY_TREE_FIELDS

SCRIPT = REPO / "scripts" / "validate_example_output.py"
EXPECTATIONS = REPO / "data" / "example" / "expected_annotations.tsv"
TUTORIAL = REPO / "docs" / "tutorials" / "first-run.md"
RUN_ASSEMBLIES = REPO / "docs" / "how-to" / "run-assemblies.md"
CLI_REFERENCE = REPO / "docs" / "reference" / "cli.md"
OUTPUT_REFERENCE = REPO / "docs" / "reference" / "outputs.md"
TOP_HIT_REPORTING = REPO / "scripts" / "top_hit_reporting.py"
CONDA_ENVIRONMENT = REPO / "config" / "environment.yml"
DATABASE_VERSION = "1.0.2"
SPEC = importlib.util.spec_from_file_location("validate_example_output", SCRIPT)
validator = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(validator)


class ExampleOutputTests(unittest.TestCase):
    def expected_rows(self, profile: str) -> list[dict[str, str]]:
        rows = [
            row
            for row in validator.read_tsv(EXPECTATIONS)
            if row["profile"] == profile
            and row["database_version"] == DATABASE_VERSION
        ]
        for row in rows:
            for field in validator.CENTROID_FIELDS:
                row.setdefault(field, "")
        return rows

    def write_summary(self, path: Path, rows: list[dict[str, str]]) -> None:
        fields = [*validator.KEY_FIELDS, *validator.ANNOTATION_FIELDS]
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows({field: row[field] for field in fields} for row in rows)

    def write_expectations(self, path: Path, rows: list[dict[str, str]]) -> None:
        fields = [
            "profile",
            "database_version",
            *validator.KEY_FIELDS,
            *validator.ANNOTATION_FIELDS,
        ]
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows({field: row[field] for field in fields} for row in rows)

    def test_each_profile_defines_ten_exact_annotations(self) -> None:
        for profile in ("curated", "img"):
            with self.subTest(profile=profile):
                rows = self.expected_rows(profile)
                self.assertEqual(len(rows), 10)
                self.assertEqual(
                    Counter(row["model"] for row in rows),
                    {"RF00177": 9, "RF01960": 1},
                )
                self.assertEqual(len({row["contig_name"] for row in rows}), 10)

    def test_validator_accepts_each_exact_profile_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            for profile in ("curated", "img"):
                with self.subTest(profile=profile):
                    summary = Path(temporary) / f"{profile}.tsv"
                    self.write_summary(summary, self.expected_rows(profile))
                    self.assertEqual(
                        validator.validate_example(
                            summary,
                            EXPECTATIONS,
                            profile,
                            DATABASE_VERSION,
                        ),
                        {"RF00177": 9, "RF01960": 1},
                    )

    def test_validator_accepts_centroid_contracts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for profile in ("curated", "img"):
                with self.subTest(profile=profile):
                    rows = self.expected_rows(profile)
                    expectations = root / f"{profile}-expectations.tsv"
                    summary = root / f"{profile}-summary.tsv"
                    self.write_expectations(expectations, rows)
                    self.write_summary(summary, rows)
                    self.assertEqual(
                        validator.validate_example(
                            summary,
                            expectations,
                            profile,
                            DATABASE_VERSION,
                        ),
                        {"RF00177": 9, "RF01960": 1},
                    )

    def test_validator_rejects_empty_v102_img_centroid_evidence(self) -> None:
        rows = self.expected_rows("img")
        derived = next(
            row
            for row in rows
            if row["taxonomy_assignment_method"] == "updated_reference_cluster"
        )
        derived["centroid_names"] = ""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expectations = root / "expectations.tsv"
            summary = root / "summary.tsv"
            self.write_expectations(expectations, rows)
            self.write_summary(summary, rows)
            with self.assertRaisesRegex(ValueError, "nonempty centroid"):
                validator.validate_example(
                    summary,
                    expectations,
                    "img",
                    DATABASE_VERSION,
                )

    def test_validator_rejects_v102_curated_centroid_evidence(self) -> None:
        rows = self.expected_rows("curated")
        rows[0]["centroid_names"] = "unexpected_centroid"
        rows[0]["centroid_taxonomy"] = rows[0]["taxonomy"]
        rows[0]["centroid_taxonomy_source"] = rows[0]["taxonomy_source"]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expectations = root / "expectations.tsv"
            summary = root / "summary.tsv"
            self.write_expectations(expectations, rows)
            self.write_summary(summary, rows)
            with self.assertRaisesRegex(ValueError, "curated.*centroid evidence"):
                validator.validate_example(
                    summary,
                    expectations,
                    "curated",
                    DATABASE_VERSION,
                )

    def test_validator_rejects_uncapped_v102_img_member_taxonomy(self) -> None:
        rows = self.expected_rows("img")
        derived = next(
            row
            for row in rows
            if row["taxonomy_assignment_method"] == "updated_reference_cluster"
        )
        derived["taxonomy"] += ";Unsupported propagated lineage"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expectations = root / "expectations.tsv"
            summary = root / "summary.tsv"
            self.write_expectations(expectations, rows)
            self.write_summary(summary, rows)
            with self.assertRaisesRegex(ValueError, "domain-capped"):
                validator.validate_example(
                    summary,
                    expectations,
                    "img",
                    DATABASE_VERSION,
                )

    def test_validator_rejects_shallow_v102_prokaryotic_centroid_taxonomy(self) -> None:
        rows = self.expected_rows("img")
        prokaryotic = next(
            row
            for row in rows
            if row["taxonomy_assignment_method"] == "updated_reference_cluster"
            and row["taxonomy_domain"] in {"Bacteria", "Archaea"}
        )
        prokaryotic["centroid_taxonomy"] = prokaryotic["taxonomy_domain"]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expectations = root / "expectations.tsv"
            summary = root / "summary.tsv"
            self.write_expectations(expectations, rows)
            self.write_summary(summary, rows)
            with self.assertRaisesRegex(ValueError, "deeper centroid taxonomy"):
                validator.validate_example(
                    summary,
                    expectations,
                    "img",
                    DATABASE_VERSION,
                )

    def test_validator_reports_changed_taxonomy(self) -> None:
        rows = self.expected_rows("curated")
        rows[0]["taxonomy"] = "Bacteria"
        with tempfile.TemporaryDirectory() as temporary:
            summary = Path(temporary) / "changed.tsv"
            self.write_summary(summary, rows)
            with self.assertRaisesRegex(ValueError, "taxonomy.*expected.*observed"):
                validator.validate_example(
                    summary,
                    EXPECTATIONS,
                    "curated",
                    DATABASE_VERSION,
                )

    def test_validator_rejects_missing_and_extra_rows(self) -> None:
        expected = self.expected_rows("curated")
        extra = {
            **expected[0],
            "sample": "unexpected_sample",
            "contig_name": "unexpected_contig",
        }
        cases = {
            "missing": expected[:-1],
            "extra": [*expected, extra],
        }
        with tempfile.TemporaryDirectory() as temporary:
            for label, rows in cases.items():
                with self.subTest(label=label):
                    summary = Path(temporary) / f"{label}.tsv"
                    self.write_summary(summary, rows)
                    with self.assertRaisesRegex(ValueError, label):
                        validator.validate_example(
                            summary,
                            EXPECTATIONS,
                            "curated",
                            DATABASE_VERSION,
                        )

    def test_validator_rejects_unrecognized_database_version(self) -> None:
        with self.assertRaisesRegex(ValueError, "no annotation contract"):
            validator.validate_example(
                Path("unused.tsv"),
                EXPECTATIONS,
                "curated",
                "9.9.9",
            )

    def test_validator_rejects_a_truncated_expectation_contract(self) -> None:
        rows = self.expected_rows("curated")[:4] + self.expected_rows("curated")[-2:]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expectations = root / "expectations.tsv"
            summary = root / "summary.tsv"
            with expectations.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=list(rows[0]),
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerows(rows)
            self.write_summary(summary, rows)
            with self.assertRaisesRegex(ValueError, "must contain 10 native rows"):
                validator.validate_example(
                    summary,
                    expectations,
                    "curated",
                    DATABASE_VERSION,
                )

    def test_tutorial_cut_fields_match_the_summary_schema(self) -> None:
        summary_fields = SUMMARY_FIELDS
        one_based_columns = (2, 3, 5, 14, 15, 23, 24, 25, 26)
        self.assertEqual(
            [summary_fields[index - 1] for index in one_based_columns],
            [
                "sample",
                "model",
                "coordinates",
                "reference_source",
                "taxonomy",
                "centroid_names",
                "centroid_taxonomy",
                "centroid_taxonomy_source",
                "reference_identifiers",
            ],
        )
        self.assertIn(
            "cut -f2,3,5,14,15,23-26 results/smoke/cmsearch_summary.tsv",
            TUTORIAL.read_text(),
        )
        self.assertEqual(
            summary_fields[25:28],
            ["reference_identifiers", "reference_versions", "query_sequence"],
        )
        self.assertEqual(summary_fields[28:], SUMMARY_TREE_FIELDS)

    def test_ranked_hit_documentation_matches_the_output_schema(self) -> None:
        module = ast.parse(TOP_HIT_REPORTING.read_text())
        assignment = next(
            node
            for node in module.body
            if isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == "TOP_HIT_FIELDS"
                for target in node.targets
            )
        )
        fields = ast.literal_eval(assignment.value)
        one_based_columns = (1, 4, 5, 7, 9, 10, 17, 19, 20, 28)
        self.assertEqual(
            [fields[index - 1] for index in one_based_columns],
            [
                "name",
                "hit_rank",
                "selection_reason",
                "reference_identifiers",
                "reference_source",
                "taxonomy",
                "centroid_taxonomy",
                "blast_pident",
                "blast_length",
                "blast_bitscore",
            ],
        )
        command = "cut -f1,4-5,7,9-10,17,19-20,28"
        self.assertIn(command, TUTORIAL.read_text())
        self.assertIn(command, RUN_ASSEMBLIES.read_text())
        self.assertIn("`--top_hits` | `5`", CLI_REFERENCE.read_text())
        self.assertIn("`blast_top_hits.tsv`", OUTPUT_REFERENCE.read_text())
        self.assertIn("`--tree_classification` | off", CLI_REFERENCE.read_text())
        self.assertIn("`tree_nearest_neighbors.tsv`", OUTPUT_REFERENCE.read_text())

    def test_tree_dependencies_are_declared_in_the_conda_environment(self) -> None:
        environment = CONDA_ENVIRONMENT.read_text()
        self.assertIn("  - iqtree >=3.0,<4\n", environment)
        self.assertIn("  - ete4 >=4.4,<5\n", environment)


if __name__ == "__main__":
    unittest.main()
