import ast
import csv
import importlib.util
import tempfile
import unittest
from collections import Counter
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "validate_example_output.py"
EXPECTATIONS = REPO / "data" / "example" / "expected_annotations.tsv"
TUTORIAL = REPO / "docs" / "tutorials" / "first-run.md"
SPEC = importlib.util.spec_from_file_location("validate_example_output", SCRIPT)
validator = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(validator)


class ExampleOutputTests(unittest.TestCase):
    def expected_rows(self, profile: str) -> list[dict[str, str]]:
        return [
            row
            for row in validator.read_tsv(EXPECTATIONS)
            if row["profile"] == profile and row["database_version"] == "1.0.1"
        ]

    def write_summary(self, path: Path, rows: list[dict[str, str]]) -> None:
        fields = [*validator.KEY_FIELDS, *validator.ANNOTATION_FIELDS]
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
                            "1.0.1",
                        ),
                        {"RF00177": 9, "RF01960": 1},
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
                    "1.0.1",
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
                            "1.0.1",
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
            with self.assertRaisesRegex(ValueError, "exactly 9 RF00177"):
                validator.validate_example(
                    summary,
                    expectations,
                    "curated",
                    "1.0.1",
                )

    def test_tutorial_cut_fields_match_the_summary_schema(self) -> None:
        module = ast.parse((REPO / "scripts" / "annotate_hits.py").read_text())
        assignment = next(
            node
            for node in module.body
            if isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == "SUMMARY_FIELDS"
                for target in node.targets
            )
        )
        summary_fields = ast.literal_eval(assignment.value)
        one_based_columns = (2, 3, 5, 14, 15)
        self.assertEqual(
            [summary_fields[index - 1] for index in one_based_columns],
            ["sample", "model", "coordinates", "reference_source", "taxonomy"],
        )
        self.assertIn(
            "cut -f2,3,5,14,15 results/smoke/cmsearch_summary.tsv",
            TUTORIAL.read_text(),
        )


if __name__ == "__main__":
    unittest.main()
