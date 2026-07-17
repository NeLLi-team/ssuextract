from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import duckdb


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "notebooks"))

import prepare_database_composition as composition


class TaxonomyRankRowsTests(unittest.TestCase):
    def test_pr2_suffix_disclosure_uses_the_build_contract(self) -> None:
        expected = tuple(composition.PR2_TAXONOMY_SUFFIX_TO_COMPARTMENT)
        self.assertEqual(composition.PR2_COMPARTMENT_SUFFIXES, expected)
        documentation = (
            REPO / "docs" / "reference" / "database-composition.md"
        ).read_text(encoding="utf-8")
        notebook = json.loads(
            (REPO / "notebooks" / "database_composition.ipynb").read_text(
                encoding="utf-8"
            )
        )
        notebook_markdown = "\n".join(
            "".join(cell["source"])
            for cell in notebook["cells"]
            if cell["cell_type"] == "markdown"
        )
        for suffix in expected:
            with self.subTest(suffix=suffix):
                self.assertIn(f"`:{suffix}`", documentation)
                self.assertIn(f"`:{suffix}`", notebook_markdown)

    def test_pr2_placeholder_remains_a_reported_supergroup(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parquet = Path(temporary) / "preferred.parquet"
            connection = duckdb.connect(":memory:")
            connection.execute(
                "CREATE TABLE preferred (domain VARCHAR, taxonomy VARCHAR)"
            )
            connection.executemany(
                "INSERT INTO preferred VALUES (?, ?)",
                [
                    ("Bacteria", "Bacteria;Pseudomonadota"),
                    ("Archaea", "Archaea;Thermoproteota"),
                    ("Eukaryota", "Eukaryota;TSAR"),
                    ("Eukaryota", "Eukaryota;Eukaryota_X"),
                    ("Eukaryota", "Eukaryota"),
                ],
            )
            connection.execute(
                "COPY preferred TO ? (FORMAT PARQUET)", [str(parquet)]
            )

            resolution, lineages = composition.taxonomy_rank_rows(
                connection, "curated", "1.0.1", parquet
            )

        eukaryota = next(row for row in resolution if row["domain"] == "Eukaryota")
        self.assertEqual(eukaryota["total_sequences"], 3)
        self.assertEqual(eukaryota["resolved_sequences"], 2)
        self.assertEqual(eukaryota["unresolved_sequences"], 1)
        self.assertEqual(eukaryota["named_lineages"], 2)
        self.assertEqual(
            [row["lineage"] for row in lineages if row["domain"] == "Eukaryota"],
            ["Eukaryota_X", "TSAR"],
        )


if __name__ == "__main__":
    unittest.main()
