import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

import duckdb


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from annotate_hits import SUMMARY_FIELDS, annotate_hits
from finalize_summaries import (
    load_metadata,
    load_summary_rows,
    merge_m8_files,
    write_category_summary,
    write_detailed_summary,
)
from hit_processing import HIT_FIELDS, META_FIELDS


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_taxonomy_parquet(path: Path, rows: list[tuple[str, ...]]) -> None:
    connection = duckdb.connect(":memory:")
    connection.execute(
        """
        CREATE TABLE taxonomy (
            sequence_id VARCHAR,
            reference_source VARCHAR,
            taxonomy VARCHAR,
            taxonomy_source VARCHAR,
            domain VARCHAR,
            compartment VARCHAR,
            assignment_method VARCHAR,
            cross_domain_conflict BOOLEAN,
            taxonomy_alternatives VARCHAR
        )
        """
    )
    if rows:
        connection.executemany(
            "INSERT INTO taxonomy VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [(*row, False, "") if len(row) == 7 else row for row in rows],
        )
    connection.execute("COPY taxonomy TO ? (FORMAT PARQUET, COMPRESSION ZSTD)", [str(path)])
    connection.close()


class AnnotationTests(unittest.TestCase):
    def test_lca_does_not_collapse_missing_internal_ranks(self) -> None:
        from annotate_hits import lowest_common_taxonomy

        self.assertEqual(
            lowest_common_taxonomy(
                ["Bacteria;P;;OrderA", "Bacteria;P;ClassB;OrderB"]
            ),
            "Bacteria;P",
        )

    def test_best_bitscore_is_selected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "sample.hits.tsv"
            m8 = root / "sample.m8"
            output = root / "sample.summary.tsv"
            write_tsv(
                hits,
                HIT_FIELDS,
                [
                    {
                        "name": "query1",
                        "sample": "sample",
                        "model": "RFTEST",
                        "length": 4,
                        "coordinates": "1-4",
                        "strand": "+",
                        "sequence_type": "simple",
                        "contig_name": "contig1",
                        "is_assembled": "False",
                    }
                ],
            )
            m8.write_text(
                "query1\tlow\t90.0\t4\t0\t0\t1\t4\t1\t4\t1e-5\t10\n"
                "query1\thigh;Bacteria\t99.0\t4\t0\t0\t1\t4\t1\t4\t1e-9\t20\n"
            )
            annotate_hits(hits, m8, output)
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(rows[0]["blast_sseqid"], "high;Bacteria")
        self.assertEqual(rows[0]["blast_bitscore"], "20")

    def test_no_hits_still_writes_header(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "sample.hits.tsv"
            m8 = root / "sample.m8"
            output = root / "sample.summary.tsv"
            write_tsv(hits, HIT_FIELDS, [])
            m8.touch()
            annotate_hits(hits, m8, output)

            self.assertEqual(output.read_text().rstrip("\n").split("\t"), SUMMARY_FIELDS)
            self.assertNotIn("\r", output.read_text())

    def test_equal_bitscore_taxonomy_is_lowest_common_ancestor(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "sample.hits.tsv"
            m8 = root / "sample.m8"
            taxonomy = root / "taxonomy.parquet"
            output = root / "sample.summary.tsv"
            write_tsv(
                hits,
                HIT_FIELDS,
                [
                    {
                        "name": "query1",
                        "sample": "sample",
                        "model": "RF01960",
                        "length": 4,
                        "coordinates": "1-4",
                        "strand": "+",
                        "sequence_type": "simple",
                        "contig_name": "contig1",
                        "is_assembled": "False",
                    }
                ],
            )
            m8.write_text(
                "query1	SSU_b	98.0	4	0	0	1	4	1	4	1e-9	20\n"
                "query1	SSU_a	99.0	4	0	0	1	4	1	4	1e-9	20\n"
            )
            write_taxonomy_parquet(
                taxonomy,
                [
                    (
                        "SSU_a",
                        "PR2",
                        "Eukaryota;TSAR;Alveolata",
                        "PR2",
                        "Eukaryota",
                        "nucleus",
                        "native",
                    ),
                    (
                        "SSU_b",
                        "PR2",
                        "Eukaryota;TSAR;Rhizaria",
                        "PR2",
                        "Eukaryota",
                        "nucleus",
                        "native",
                    ),
                ],
            )
            annotate_hits(hits, m8, output, taxonomy)
            with output.open(newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(row["blast_sseqid"], "SSU_a")
        self.assertEqual(row["taxonomy"], "Eukaryota;TSAR")
        self.assertEqual(row["taxonomy_source"], "PR2")
        self.assertEqual(row["taxonomy_domain"], "Eukaryota")
        self.assertEqual(row["blast_tied_subjects"], "2")

    def test_manifest_database_requires_taxonomy_for_every_subject(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "sample.hits.tsv"
            m8 = root / "sample.m8"
            taxonomy = root / "taxonomy.parquet"
            output = root / "sample.summary.tsv"
            write_tsv(
                hits,
                HIT_FIELDS,
                [
                    {
                        "name": "query1",
                        "sample": "sample",
                        "model": "RF00177",
                        "length": 4,
                        "coordinates": "1-4",
                        "strand": "+",
                        "sequence_type": "simple",
                        "contig_name": "contig1",
                        "is_assembled": "False",
                    }
                ],
            )
            m8.write_text(
                "query1	missing	99.0	4	0	0	1	4	1	4	1e-9	20\n"
            )
            write_taxonomy_parquet(taxonomy, [])

            with self.assertRaisesRegex(ValueError, "Missing taxonomy metadata"):
                annotate_hits(hits, m8, output, taxonomy)

    def test_truncated_equal_best_set_backs_taxonomy_off_to_domain(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "sample.hits.tsv"
            m8 = root / "sample.m8"
            taxonomy = root / "taxonomy.parquet"
            output = root / "sample.summary.tsv"
            row = dict.fromkeys(HIT_FIELDS, "")
            row.update(
                name="query1",
                sample="sample",
                model="RF01960",
                length=4,
                coordinates="1-4",
                strand="+",
                sequence_type="simple",
                contig_name="contig1",
                is_assembled="False",
            )
            write_tsv(hits, HIT_FIELDS, [row])
            m8.write_text(
                "query1\tSSU_a\t99\t4\t0\t0\t1\t4\t1\t4\t1e-9\t20\n"
                "query1\tSSU_b\t99\t4\t0\t0\t1\t4\t1\t4\t1e-9\t20\n"
                "query1\tSSU_c\t99\t4\t0\t0\t1\t4\t1\t4\t1e-9\t20\n"
            )
            write_taxonomy_parquet(
                taxonomy,
                [
                    (
                        "SSU_a", "PR2", "Eukaryota;TSAR;Alveolata", "PR2",
                        "Eukaryota", "nucleus", "native",
                    ),
                    (
                        "SSU_b", "PR2", "Eukaryota;TSAR;Rhizaria", "PR2",
                        "Eukaryota", "nucleus", "native",
                    ),
                    (
                        "SSU_c", "PR2", "Eukaryota;TSAR;Stramenopiles", "PR2",
                        "Eukaryota", "nucleus", "native",
                    ),
                ],
            )
            annotate_hits(hits, m8, output, taxonomy, max_targets=2)
            with output.open(newline="") as handle:
                result = next(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(result["taxonomy"], "Eukaryota")
        self.assertEqual(result["compartment"], "")
        self.assertEqual(result["taxonomy_assignment_method"], "truncated_equal_best_lca")
        self.assertEqual(result["blast_ties_truncated"], "true")

    def test_cross_domain_subject_ambiguity_is_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "sample.hits.tsv"
            m8 = root / "sample.m8"
            taxonomy = root / "taxonomy.parquet"
            output = root / "sample.summary.tsv"
            row = dict.fromkeys(HIT_FIELDS, "")
            row.update(
                name="query1",
                sample="sample",
                model="RF00177",
                length=4,
                coordinates="1-4",
                strand="+",
                sequence_type="simple",
                contig_name="contig1",
                is_assembled="False",
            )
            write_tsv(hits, HIT_FIELDS, [row])
            m8.write_text(
                "query1\tSSU_ambiguous\t100\t4\t0\t0\t1\t4\t1\t4\t0\t20\n"
                "query1\tSSU_bacteria\t100\t4\t0\t0\t1\t4\t1\t4\t0\t20\n"
            )
            alternatives = json.dumps(
                [
                    {
                        "taxonomy_source": "SILVA",
                        "taxonomy": "Bacteria;P",
                        "domain": "Bacteria",
                    },
                    {
                        "taxonomy_source": "PR2",
                        "taxonomy": "Eukaryota;T",
                        "domain": "Eukaryota",
                    },
                ]
            )
            write_taxonomy_parquet(
                taxonomy,
                [
                    (
                        "SSU_ambiguous",
                        "PR2+SILVA",
                        "",
                        "PR2+SILVA",
                        "ambiguous",
                        "mixed",
                        "cross_domain_ambiguous_exact_sequence",
                        True,
                        alternatives,
                    ),
                    (
                        "SSU_bacteria",
                        "SILVA",
                        "Bacteria;P",
                        "SILVA",
                        "Bacteria",
                        "",
                        "native",
                        False,
                        "",
                    ),
                ],
            )
            annotate_hits(hits, m8, output, taxonomy)
            with output.open(newline="") as handle:
                result = next(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(result["taxonomy_domain"], "ambiguous")
        self.assertEqual(result["taxonomy"], "")
        self.assertEqual(result["compartment"], "mixed")
        self.assertEqual(
            result["taxonomy_assignment_method"],
            "cross_domain_ambiguous_equal_best",
        )
        self.assertEqual(
            {row["domain"] for row in json.loads(result["taxonomy_alternatives"])},
            {"Bacteria", "Eukaryota"},
        )


class FinalSummaryTests(unittest.TestCase):
    def test_no_hit_sample_is_present_in_category_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / "sample.meta.tsv"
            summary = root / "sample.summary.tsv"
            detailed_output = root / "detailed.tsv"
            category_output = root / "categories.tsv"
            write_tsv(metadata, META_FIELDS, [{"sample": "sample", "model": "RFTEST"}])
            write_tsv(summary, SUMMARY_FIELDS, [])

            metadata_rows = load_metadata(str(root / "*.meta.tsv"))
            summary_rows = load_summary_rows(str(root / "*.summary.tsv"))
            write_detailed_summary(summary_rows, detailed_output)
            write_category_summary(summary_rows, metadata_rows, category_output)

            self.assertEqual(detailed_output.read_text().splitlines()[0].split("\t"), SUMMARY_FIELDS)
            self.assertEqual(category_output.read_text(), "\nsample\n")

    def test_categories_are_deduplicated_per_sample_contig(self) -> None:
        rows = [
            {
                "name": "q1",
                "sample": "sample",
                "model": "RF1",
                "length": "4",
                "coordinates": "1-4",
                "strand": "+",
                "sequence_type": "simple",
                "contig_name": "contig1",
                "blast_sseqid": "ref;Bacteria;Patescibacteria",
                "blast_pident": "99.0",
                "blast_length": "4",
                "blast_bitscore": "20.0",
                "is_assembled": "False",
            },
            {
                "name": "q2",
                "sample": "sample",
                "model": "RF2",
                "length": "4",
                "coordinates": "2-5",
                "strand": "+",
                "sequence_type": "simple",
                "contig_name": "contig1",
                "blast_sseqid": "ref;Bacteria;Patescibacteria",
                "blast_pident": "99.0",
                "blast_length": "4",
                "blast_bitscore": "20.0",
                "is_assembled": "False",
            },
        ]
        metadata = [{"sample": "sample", "model": "RF1"}]
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "categories.tsv"
            write_category_summary(rows, metadata, output)
            lines = output.read_text().splitlines()

        self.assertEqual(lines[0], "\tBacteriaSSU\tPatescibacteriaSSU")
        self.assertEqual(lines[1], "sample\t1\t1")

    def test_categories_use_normalized_taxonomy_and_compartment(self) -> None:
        row = dict.fromkeys(SUMMARY_FIELDS, "")
        row.update(
            name="q1",
            sample="sample",
            model="RF00177",
            contig_name="contig1",
            coordinates="1-4",
            strand="+",
            blast_sseqid="SSU_stable_identifier",
            taxonomy="Eukaryota;Archaeplastida",
            taxonomy_domain="Eukaryota",
            compartment="plastid",
        )
        metadata = [{"sample": "sample", "model": "RF00177"}]
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "categories.tsv"
            write_category_summary([row], metadata, output)
            lines = output.read_text().splitlines()

        self.assertEqual(lines[0], "\tPlastidSSU\tEukaryotaSSU")
        self.assertEqual(lines[1], "sample\t1\t1")

    def test_categories_union_divergent_model_assignments_per_contig(self) -> None:
        rows = []
        for model, domain in (("RF00177", "Bacteria"), ("RF01960", "Eukaryota")):
            row = dict.fromkeys(SUMMARY_FIELDS, "")
            row.update(
                name=f"q-{model}",
                sample="sample",
                model=model,
                contig_name="contig1",
                coordinates="1-4",
                strand="+",
                blast_sseqid=f"ref-{model}",
                taxonomy=domain,
                taxonomy_domain=domain,
            )
            rows.append(row)
        metadata = [{"sample": "sample", "model": "RF00177"}]
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "categories.tsv"
            write_category_summary(rows, metadata, output)
            lines = output.read_text().splitlines()

        self.assertEqual(lines[0], "\tBacteriaSSU\tEukaryotaSSU")
        self.assertEqual(lines[1], "sample\t1\t1")

    def test_legacy_pipe_delimiters_are_tokenized(self) -> None:
        row = dict.fromkeys(SUMMARY_FIELDS, "")
        row.update(
            name="q1",
            sample="sample",
            model="RF00177",
            contig_name="contig1",
            coordinates="1-4",
            strand="+",
            blast_sseqid="legacy|Bacteria|Patescibacteria",
        )
        metadata = [{"sample": "sample", "model": "RF00177"}]
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "categories.tsv"
            write_category_summary([row], metadata, output)
            lines = output.read_text().splitlines()

        self.assertEqual(lines[0], "\tBacteriaSSU\tPatescibacteriaSSU")
        self.assertEqual(lines[1], "sample\t1\t1")

    def test_detailed_rows_are_sorted_deterministically(self) -> None:
        row_a = dict.fromkeys(SUMMARY_FIELDS, "")
        row_a.update(
            sample="b", model="RF2", contig_name="contig2", coordinates="10-20", strand="+"
        )
        row_b = dict.fromkeys(SUMMARY_FIELDS, "")
        row_b.update(
            sample="a", model="RF1", contig_name="contig1", coordinates="1-4", strand="+"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_tsv(root / "z.summary.tsv", SUMMARY_FIELDS, [row_a])
            write_tsv(root / "a.summary.tsv", SUMMARY_FIELDS, [row_b])
            rows = load_summary_rows(str(root / "*.summary.tsv"))

        self.assertEqual([row["sample"] for row in rows], ["a", "b"])

    def test_m8_merge_uses_filename_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "b.m8").write_text("b\n")
            (root / "a.m8").write_text("a\n")
            output = root / "merged.m8"
            output.write_text("stale\n")
            merge_m8_files(str(root / "*.m8"), output)
            self.assertEqual(output.read_text(), "a\nb\n")


if __name__ == "__main__":
    unittest.main()
