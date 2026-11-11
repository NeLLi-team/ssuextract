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


class CompositionRowsTests(unittest.TestCase):
    def test_coordinate_counts_include_taxa_without_a_location_row(self) -> None:
        """Count each IMG taxon, including incomplete and absent locations."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = {
                name: root / f"{name}.parquet"
                for name in (
                    "sequences",
                    "source_records",
                    "preferred_taxonomy",
                    "img_location",
                )
            }
            with duckdb.connect(":memory:") as connection:
                connection.execute(
                    "COPY (SELECT '16S' AS markers) TO ? (FORMAT PARQUET)",
                    [str(paths["sequences"])],
                )
                connection.execute(
                    """COPY (SELECT * FROM (VALUES
                           ('IMG', '16S', '1'), ('IMG', '18S', '1'),
                           ('IMG', '16S', '2'), ('IMG', '16S', '3'),
                           ('SILVA', '16S', '4')
                       ) AS sources(reference_source, marker, taxon_oid))
                       TO ? (FORMAT PARQUET)""",
                    [str(paths["source_records"])],
                )
                connection.execute(
                    """COPY (SELECT 'Bacteria' AS domain,
                           false AS cross_domain_conflict)
                       TO ? (FORMAT PARQUET)""",
                    [str(paths["preferred_taxonomy"])],
                )
                connection.execute(
                    """COPY (SELECT * FROM (VALUES
                           ('1', 12.0, 34.0), ('2', NULL, 10.0)
                       ) AS locations(taxon_oid, latitude, longitude))
                       TO ? (FORMAT PARQUET)""",
                    [str(paths["img_location"])],
                )
                qa = {
                    "counts": {"sequences": 1},
                    "img_location": {
                        "rows": 2,
                        "latitude_rows": 1,
                        "longitude_rows": 2,
                    },
                }

                rows = composition.composition_rows(
                    connection, "img", "1.0.2", paths, qa
                )

        coordinate_counts = {
            row["category"]: row["count"]
            for row in rows
            if row["dimension"] == "img_taxon_coordinates"
        }
        self.assertEqual(
            coordinate_counts, {"Located": 1, "Coordinates unavailable": 2}
        )
        self.assertEqual(sum(coordinate_counts.values()), 3)


if __name__ == "__main__":
    unittest.main()
