import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

import duckdb


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from annotate_hits import (
    SUMMARY_FIELDS,
    annotate_hits,
    load_taxonomy_records,
)
from finalize_summaries import (
    apply_tree_assignments,
    load_metadata,
    load_summary_rows,
    load_top_hit_rows,
    merge_m8_files,
    write_category_summary,
    write_detailed_summary,
    write_top_hit_summary,
)
from hit_processing import HIT_FIELDS, META_FIELDS
from top_hit_reporting import TOP_HIT_FIELDS, load_reference_records
from tree_schema import TREE_ASSIGNMENT_FIELDS


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_taxonomy_parquet(
    path: Path,
    rows: list[tuple[str, ...]],
    *,
    include_centroids: bool = False,
) -> None:
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
            {centroid_columns}
        )
        """.format(
            centroid_columns=(
                ", centroid_names VARCHAR, centroid_taxonomy VARCHAR, "
                "centroid_taxonomy_source VARCHAR"
                if include_centroids
                else ""
            )
        )
    )
    if rows:
        normalized = [(*row, False, "") if len(row) == 7 else row for row in rows]
        if include_centroids:
            normalized = [
                (*row, "", "", "") if len(row) == 9 else row for row in normalized
            ]
        connection.executemany(
            "INSERT INTO taxonomy VALUES ("
            + ", ".join("?" for _ in range(12 if include_centroids else 9))
            + ")",
            normalized,
        )
    connection.execute("COPY taxonomy TO ? (FORMAT PARQUET, COMPRESSION ZSTD)", [str(path)])
    connection.close()


def write_source_records_parquet(
    path: Path,
    rows: list[tuple[str, str, str, str]],
) -> None:
    connection = duckdb.connect(":memory:")
    connection.execute(
        """
        CREATE TABLE source_records (
            sequence_id VARCHAR,
            reference_source VARCHAR,
            source_version VARCHAR,
            source_identifier VARCHAR
        )
        """
    )
    if rows:
        connection.executemany(
            "INSERT INTO source_records VALUES (?, ?, ?, ?)", rows
        )
    connection.execute(
        "COPY source_records TO ? (FORMAT PARQUET, COMPRESSION ZSTD)", [str(path)]
    )
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
            top_hits = root / "sample.top_hits.tsv"
            write_tsv(hits, HIT_FIELDS, [])
            m8.touch()
            annotate_hits(hits, m8, output, top_hits_output=top_hits)

            self.assertEqual(output.read_text().rstrip("\n").split("\t"), SUMMARY_FIELDS)
            self.assertEqual(top_hits.read_text().rstrip("\n").split("\t"), TOP_HIT_FIELDS)
            self.assertNotIn("\r", output.read_text())

    def test_top_hits_join_metadata_and_retain_best_source_hit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "sample.hits.tsv"
            m8 = root / "sample.m8"
            taxonomy = root / "preferred_taxonomy.parquet"
            sources = root / "source_records.parquet"
            fasta = root / "sample.fna"
            summary = root / "sample.summary.tsv"
            top_hits = root / "sample.top_hits.tsv"
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
            fasta.write_text(">query1\nacGt\n", encoding="ascii")
            m8.write_text(
                "query1\tSSU_img\t99\t4\t0\t0\t1\t4\t1\t4\t1e-20\t100\n"
                "query1\tSSU_img_tie\t98\t4\t1\t0\t1\t4\t1\t4\t1e-19\t100\n"
                "query1\tSSU_pr2\t98\t4\t1\t0\t1\t4\t1\t4\t1e-10\t90\n"
                "query1\tSSU_third\t97\t4\t1\t0\t1\t4\t1\t4\t1e-5\t80\n",
                encoding="ascii",
            )
            write_taxonomy_parquet(
                taxonomy,
                [
                    (
                        "SSU_img", "IMG", "Unclassified", "SILVA+PR2",
                        "Unclassified", "", "updated_reference_unclassified",
                        False, "",
                        '["IMG_3300000001.a:contig_1"]',
                        "Eukaryota;Amoebozoa;Discosea;Echinamoebida", "PR2",
                    ),
                    (
                        "SSU_img_tie", "IMG", "Unclassified", "SILVA+PR2",
                        "Unclassified", "", "updated_reference_unclassified",
                        False, "",
                        '["IMG_3300000002.a:contig_2"]',
                        "Eukaryota;Amoebozoa;Discosea;Echinamoebida", "PR2",
                    ),
                    (
                        "SSU_pr2", "PR2",
                        "Eukaryota;Amoebozoa;Discosea;Echinamoebida;Echinamoeba",
                        "PR2", "Eukaryota", "nucleus", "native", False, "", "", "", "",
                    ),
                    (
                        "SSU_third", "PR2", "Eukaryota;Amoebozoa", "PR2",
                        "Eukaryota", "nucleus", "native", False, "", "", "", "",
                    ),
                ],
                include_centroids=True,
            )
            write_source_records_parquet(
                sources,
                [
                    ("SSU_img", "IMG", "2025", "IMG_3300000001.a:contig_9"),
                    ("SSU_img_tie", "IMG", "2025", "IMG_3300000002.a:contig_8"),
                    ("SSU_pr2", "PR2", "5.1.1", "AB123456.1.1800_U"),
                    ("SSU_third", "PR2", "5.1.1", "AB999999.1.1800_U"),
                ],
            )
            annotate_hits(
                hits,
                m8,
                summary,
                taxonomy,
                query_fasta=fasta,
                source_records_file=sources,
                top_hits_output=top_hits,
                top_hits=1,
            )
            with summary.open(newline="") as handle:
                summary_row = next(csv.DictReader(handle, delimiter="\t"))
            with top_hits.open(newline="") as handle:
                top_rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(summary_row["query_sequence"], "acGt")
        self.assertEqual(
            summary_row["reference_identifiers"],
            "IMG:IMG_3300000001.a:contig_9",
        )
        self.assertEqual(summary_row["reference_versions"], "IMG:2025")
        self.assertEqual(summary_row["taxonomy"], "Unclassified")
        self.assertIn("Echinamoebida", summary_row["centroid_taxonomy"])
        self.assertEqual(len(top_rows), 3)
        self.assertEqual([row["hit_rank"] for row in top_rows], ["1", "2", "3"])
        self.assertEqual(
            [row["blast_sseqid"] for row in top_rows],
            ["SSU_img", "SSU_img_tie", "SSU_pr2"],
        )
        self.assertEqual(top_rows[0]["taxonomy"], "Unclassified")
        self.assertEqual(
            top_rows[0]["selection_reason"],
            "overall_top_n|equal_best_assignment|best_IMG",
        )
        self.assertIn("Echinamoebida", top_rows[0]["centroid_taxonomy"])
        self.assertEqual(top_rows[1]["selection_reason"], "equal_best_assignment")
        self.assertIn("Echinamoeba", top_rows[2]["taxonomy"])
        self.assertEqual(top_rows[2]["selection_reason"], "best_PR2")
        self.assertEqual(
            top_rows[2]["reference_identifiers"], "PR2:AB123456.1.1800_U"
        )
        self.assertEqual({row["query_sequence"] for row in top_rows}, {"acGt"})

    def test_source_metadata_is_required_for_every_reported_subject(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sources = Path(tmp) / "source_records.parquet"
            write_source_records_parquet(
                sources, [("SSU_present", "PR2", "5.1.1", "AB1.1.100_U")]
            )
            with self.assertRaisesRegex(ValueError, "Missing source metadata"):
                load_reference_records(sources, {"SSU_present", "SSU_missing"})

    def test_exact_sequence_retains_all_public_source_identifiers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sources = Path(tmp) / "source_records.parquet"
            write_source_records_parquet(
                sources,
                [
                    ("SSU_shared", "PR2", "5.1.1", "AB1.1.100_U"),
                    ("SSU_shared", "PR2", "5.1.1", "AB2.1.100_U"),
                    ("SSU_shared", "SILVA", "138.2", "AB3.1.100"),
                ],
            )
            record = load_reference_records(sources, {"SSU_shared"})["SSU_shared"]

        self.assertEqual(
            record.identifiers,
            "PR2:AB1.1.100_U|PR2:AB2.1.100_U|SILVA:AB3.1.100",
        )
        self.assertEqual(record.versions, "PR2:5.1.1|SILVA:138.2")
        self.assertEqual(record.sources, ("PR2", "SILVA"))

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

    def test_selected_best_hit_reports_its_centroid_evidence(self) -> None:
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
                "query1\tSSU_img\t99\t4\t0\t0\t1\t4\t1\t4\t0\t20\n",
                encoding="ascii",
            )
            write_taxonomy_parquet(
                taxonomy,
                [
                    (
                        "SSU_img",
                        "IMG",
                        "Bacteria",
                        "SILVA",
                        "Bacteria",
                        "",
                        "updated_reference_cluster",
                        False,
                        "",
                        '["IMG_centroid_1","REF_SILVA_AB123.1.1500"]',
                        "Bacteria;Bacteroidota;Bacteroidia;Flavobacteriales;Flavobacteriaceae",
                        "SILVA",
                    )
                ],
                include_centroids=True,
            )
            annotate_hits(hits, m8, output, taxonomy)
            with output.open(newline="") as handle:
                result = next(csv.DictReader(handle, delimiter="\t"))
            raw_output = output.read_text()

        self.assertEqual(result["blast_sseqid"], "SSU_img")
        self.assertEqual(
            result["centroid_names"],
            "IMG_centroid_1|REF_SILVA_AB123.1.1500",
        )
        self.assertNotIn('""', raw_output)
        self.assertEqual(
            result["centroid_taxonomy"],
            "Bacteria;Bacteroidota;Bacteroidia;Flavobacteriales;Flavobacteriaceae",
        )
        self.assertEqual(result["centroid_taxonomy_source"], "SILVA")

    def test_invalid_centroid_names_json_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            taxonomy = Path(tmp) / "taxonomy.parquet"
            write_taxonomy_parquet(
                taxonomy,
                [
                    (
                        "SSU_img",
                        "IMG",
                        "Bacteria",
                        "SILVA",
                        "Bacteria",
                        "",
                        "updated_reference_cluster",
                        False,
                        "",
                        "not-json",
                        "Bacteria;Bacteroidota",
                        "SILVA",
                    )
                ],
                include_centroids=True,
            )
            with self.assertRaisesRegex(ValueError, "Invalid centroid_names JSON"):
                load_taxonomy_records(taxonomy, {"SSU_img"})

    def test_incomplete_centroid_schema_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            taxonomy = Path(tmp) / "taxonomy.parquet"
            connection = duckdb.connect(":memory:")
            connection.execute(
                "CREATE TABLE taxonomy (sequence_id VARCHAR, centroid_names VARCHAR)"
            )
            connection.execute(
                "COPY taxonomy TO ? (FORMAT PARQUET, COMPRESSION ZSTD)",
                [str(taxonomy)],
            )
            connection.close()
            with self.assertRaisesRegex(ValueError, "incomplete centroid schema"):
                load_taxonomy_records(taxonomy, {"SSU_img"})

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

    def test_top_hit_rows_are_merged_in_query_rank_order(self) -> None:
        row_a = dict.fromkeys(TOP_HIT_FIELDS, "")
        row_a.update(
            name="q2", sample="sample", model="RF01960", hit_rank="2",
            blast_sseqid="subject-b",
        )
        row_b = dict.fromkeys(TOP_HIT_FIELDS, "")
        row_b.update(
            name="q1", sample="sample", model="RF01960", hit_rank="1",
            blast_sseqid="subject-a",
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_tsv(root / "z.top_hits.tsv", TOP_HIT_FIELDS, [row_a])
            write_tsv(root / "a.top_hits.tsv", TOP_HIT_FIELDS, [row_b])
            rows = load_top_hit_rows(str(root / "*.top_hits.tsv"))
            output = root / "blast_top_hits.tsv"
            write_top_hit_summary(rows, output)
            with output.open(newline="") as handle:
                merged = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual([row["name"] for row in merged], ["q1", "q2"])
        self.assertEqual(merged[0]["blast_sseqid"], "subject-a")

    def test_tree_mode_replaces_selected_taxonomy_but_retains_blast_taxonomy(self) -> None:
        row = dict.fromkeys(SUMMARY_FIELDS, "")
        row.update(
            name="query1",
            sample="sample",
            model="RF01960",
            taxonomy="Eukaryota",
            taxonomy_source="PR2",
            taxonomy_domain="Eukaryota",
            taxonomy_assignment_method="lowest_common_ancestor",
            taxonomy_mode="blast",
            blast_taxonomy="Eukaryota",
            blast_taxonomy_source="PR2",
            blast_taxonomy_domain="Eukaryota",
            blast_compartment="nucleus",
            blast_taxonomy_assignment_method="lowest_common_ancestor",
            blast_taxonomy_alternatives='[{"taxonomy":"Eukaryota"}]',
        )
        assignment = dict.fromkeys(TREE_ASSIGNMENT_FIELDS, "")
        assignment.update(
            name="query1",
            sample="sample",
            model="RF01960",
            tree_model="RF01960",
            tree_marker="18S",
            tree_taxonomy=(
                "Eukaryota;Amoebozoa;Discosea;Echinamoebida;Echinamoeba"
            ),
            tree_taxonomy_source="PR2",
            tree_taxonomy_domain="Eukaryota",
            tree_compartment="nucleus",
            tree_assignment_method="tree_nearest_named_lca",
            tree_basis_neighbors="5",
        )
        merged = apply_tree_assignments([row], [assignment], "tree")[0]
        self.assertEqual(merged["taxonomy_mode"], "tree")
        self.assertEqual(
            merged["taxonomy"],
            "Eukaryota;Amoebozoa;Discosea;Echinamoebida;Echinamoeba",
        )
        self.assertEqual(merged["blast_taxonomy"], "Eukaryota")
        self.assertEqual(merged["blast_compartment"], "nucleus")
        self.assertEqual(
            merged["blast_taxonomy_alternatives"],
            '[{"taxonomy":"Eukaryota"}]',
        )
        self.assertEqual(
            merged["taxonomy_assignment_method"], "tree_nearest_named_lca"
        )

    def test_sparse_tree_query_retains_blast_assignment_and_skip_reason(self) -> None:
        row = dict.fromkeys(SUMMARY_FIELDS, "")
        row.update(
            name="query1",
            sample="sample",
            model="RF01960",
            taxonomy="Eukaryota;Amoebozoa",
            taxonomy_source="PR2",
            taxonomy_domain="Eukaryota",
            taxonomy_assignment_method="native",
            taxonomy_mode="blast",
            blast_taxonomy="Eukaryota;Amoebozoa",
        )
        assignment = dict.fromkeys(TREE_ASSIGNMENT_FIELDS, "")
        assignment.update(
            name="query1",
            sample="sample",
            model="RF01960",
            tree_model="RF01960",
            tree_marker="18S",
            tree_route_decision="majority_global_top_hits",
            tree_route_18s_votes="2",
            tree_assignment_method="tree_skipped_insufficient_references",
            tree_basis_neighbors="0",
        )

        merged = apply_tree_assignments([row], [assignment], "tree")[0]

        self.assertEqual(merged["taxonomy_mode"], "blast")
        self.assertEqual(merged["taxonomy"], "Eukaryota;Amoebozoa")
        self.assertEqual(merged["taxonomy_assignment_method"], "native")
        self.assertEqual(
            merged["tree_assignment_method"],
            "tree_skipped_insufficient_references",
        )


if __name__ == "__main__":
    unittest.main()
