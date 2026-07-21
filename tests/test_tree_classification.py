import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from annotate_hits import BlastHit
from tree_phylogeny import classify_tree, trim_alignment
from tree_reference_selection import (
    build_alignment_input,
    choose_marker,
    prepare_tree_tasks,
)
from tree_schema import REFERENCE_FIELDS, TREE_ASSIGNMENT_FIELDS, TREE_NEIGHBOR_FIELDS


def blast_hit(subject: str, score: float) -> BlastHit:
    return BlastHit(
        subject=subject,
        percent_identity=99.0,
        alignment_length=100,
        mismatches=1,
        gap_opens=0,
        query_start=1,
        query_end=100,
        subject_start=1,
        subject_end=100,
        evalue=1e-20,
        bit_score=score,
    )


def reference_row(
    leaf: str,
    subject: str,
    rank: int,
    *,
    taxonomy: str = "",
    centroid_taxonomy: str = "",
) -> dict[str, object]:
    return {
        "leaf_id": leaf,
        "blast_sseqid": subject,
        "hit_rank": rank,
        "reference_identifiers": f"PR2:{subject}",
        "reference_versions": "PR2:5.1.1",
        "reference_source": "PR2",
        "taxonomy": taxonomy,
        "taxonomy_source": "PR2" if taxonomy else "",
        "taxonomy_domain": "Eukaryota",
        "compartment": "nucleus",
        "taxonomy_assignment_method": "native" if taxonomy else "unclassified",
        "centroid_names": "centroid" if centroid_taxonomy else "",
        "centroid_taxonomy": centroid_taxonomy,
        "centroid_taxonomy_source": "PR2" if centroid_taxonomy else "",
        "blast_pident": "99",
        "blast_length": "100",
        "blast_evalue": "1e-20",
        "blast_bitscore": str(200 - rank),
    }


class TreeSchemaTests(unittest.TestCase):
    def test_empty_sentinels_match_code_owned_schemas(self) -> None:
        sentinels = {
            "empty.tree_assignment.tsv": TREE_ASSIGNMENT_FIELDS,
            "empty.tree_neighbors.tsv": TREE_NEIGHBOR_FIELDS,
        }
        for filename, expected_fields in sentinels.items():
            with self.subTest(filename=filename):
                with (REPO / "config" / filename).open(newline="") as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    self.assertEqual(reader.fieldnames, expected_fields)
                    self.assertEqual(list(reader), [])


class TreeRoutingTests(unittest.TestCase):
    def test_no_blast_hits_retains_detected_marker(self) -> None:
        marker, decision, votes, scores = choose_marker(
            {"16S": [], "18S": []},
            "18S",
            100,
        )
        self.assertEqual(marker, "18S")
        self.assertEqual(decision, "detected_model_no_blast_hits")
        self.assertEqual(votes, {"16S": 0, "18S": 0})
        self.assertEqual(scores, {"16S": float("-inf"), "18S": float("-inf")})

    def test_global_hit_majority_selects_marker(self) -> None:
        marker, decision, votes, _scores = choose_marker(
            {
                "16S": [blast_hit("b1", 100), blast_hit("b2", 80)],
                "18S": [
                    blast_hit("e1", 120),
                    blast_hit("e2", 110),
                    blast_hit("e3", 105),
                ],
            },
            "16S",
            4,
        )
        self.assertEqual(marker, "18S")
        self.assertEqual(decision, "majority_global_top_hits")
        self.assertEqual(votes, {"16S": 1, "18S": 3})

    def test_equal_route_evidence_uses_detected_model(self) -> None:
        marker, decision, votes, scores = choose_marker(
            {
                "16S": [blast_hit("bacterial", 100)],
                "18S": [blast_hit("eukaryotic", 100)],
            },
            "18S",
            100,
        )
        self.assertEqual(marker, "18S")
        self.assertEqual(decision, "detected_model_tiebreak")
        self.assertEqual(votes, {"16S": 1, "18S": 1})
        self.assertEqual(scores["16S"], scores["18S"])

    def test_equal_vote_count_uses_best_bitscore(self) -> None:
        marker, decision, votes, scores = choose_marker(
            {
                "16S": [blast_hit("bacterial", 150)],
                "18S": [blast_hit("eukaryotic", 100)],
            },
            "18S",
            100,
        )
        self.assertEqual(marker, "16S")
        self.assertEqual(decision, "best_bitscore_tiebreak")
        self.assertEqual(votes, {"16S": 1, "18S": 1})
        self.assertGreater(scores["16S"], scores["18S"])

    def test_alignment_input_relabels_query_and_references(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "query.fna").write_text(">query name\nACGT\n")
            rows = [
                reference_row("REF0001", "subject1", 1),
                reference_row("REF0002", "subject2", 2),
                reference_row("REF0003", "subject3", 3),
            ]
            with (root / "references.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=REFERENCE_FIELDS, delimiter="\t"
                )
                writer.writeheader()
                writer.writerows(rows)
            references = root / "references.fna"
            references.write_text(
                ">subject3 third\nACGA\n"
                ">subject1 first\nACGC\n"
                ">subject2 second\nACGG\n"
            )
            output = root / "input.fna"
            build_alignment_input(root, references, output)
            headers = [
                line[1:] for line in output.read_text().splitlines() if line.startswith(">")
            ]
        self.assertEqual(headers, ["QUERY", "REF0001", "REF0002", "REF0003"])

    def test_sparse_query_writes_blast_fallback_instead_of_aborting(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            skipped = root / "skipped.tsv"
            with (
                patch(
                    "tree_reference_selection.load_query_sequences",
                    return_value={"query1": "ACGT"},
                ),
                patch(
                    "tree_reference_selection.load_blast_hits",
                    side_effect=[
                        {"query1": []},
                        {"query1": [blast_hit("e1", 100), blast_hit("e2", 90)]},
                    ],
                ),
                patch(
                    "tree_reference_selection.load_taxonomy_records",
                    return_value={},
                ),
                patch(
                    "tree_reference_selection.load_reference_records",
                    return_value={},
                ),
            ):
                directories = prepare_tree_tasks(
                    query_fasta=root / "query.fna",
                    blast_files={"16S": root / "16S.m8", "18S": root / "18S.m8"},
                    taxonomy_file=root / "taxonomy.parquet",
                    source_records_file=root / "sources.parquet",
                    sample="sample",
                    detected_model="RF01960",
                    detected_marker="18S",
                    marker_models={"16S": "RF00177", "18S": "RF01960"},
                    output_directory=root / "tasks",
                    skipped_assignments_file=skipped,
                )
            with skipped.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(directories, [])
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["tree_marker"], "18S")
        self.assertEqual(rows[0]["tree_route_18s_votes"], "2")
        self.assertEqual(
            rows[0]["tree_assignment_method"],
            "tree_skipped_insufficient_references",
        )


class TreePhylogenyTests(unittest.TestCase):
    def test_mild_gap_trimming_retains_shared_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            alignment = root / "aligned.fna"
            alignment.write_text(
                ">QUERY\nAUgGT\n"
                ">REF0001\nAT.GT\n"
                ">REF0002\nAT.GT\n"
                ">REF0003\nAT.GT\n"
            )
            output = root / "trimmed.fna"
            qc_file = root / "qc.json"
            qc = trim_alignment(
                alignment, output, qc_file, maximum_gap_fraction=0.8
            )
            sequences = [
                line for line in output.read_text().splitlines() if not line.startswith(">")
            ]
        self.assertEqual(qc["input_columns"], 5)
        self.assertEqual(qc["retained_columns"], 4)
        self.assertEqual(qc["removed_insert_columns"], 1)
        self.assertEqual(qc["removed_high_gap_columns"], 0)
        self.assertEqual(sequences, ["ATGT", "ATGT", "ATGT", "ATGT"])

    def test_tree_taxonomy_uses_nearest_named_lineages_and_centroids(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            references = root / "references.tsv"
            lineage = "Eukaryota;Amoebozoa;Discosea;Echinamoebida;Echinamoeba"
            rows = [
                reference_row(
                    "REF0001", "img_exact", 1, centroid_taxonomy=lineage
                ),
                reference_row("REF0002", "pr2_named", 2, taxonomy=lineage),
                reference_row(
                    "REF0003",
                    "pr2_other",
                    3,
                    taxonomy="Eukaryota;Amoebozoa;Tubulinea",
                ),
            ]
            with references.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=REFERENCE_FIELDS, delimiter="\t"
                )
                writer.writeheader()
                writer.writerows(rows)
            task = root / "task.json"
            task.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "name": "query1",
                        "sample": "sample",
                        "detected_model": "RF01960",
                        "tree_model": "RF01960",
                        "tree_marker": "18S",
                        "tree_route_decision": "majority_global_top_hits",
                        "tree_route_16s_votes": 0,
                        "tree_route_18s_votes": 100,
                        "tree_route_16s_best_bitscore": None,
                        "tree_route_18s_best_bitscore": 5710.0,
                    }
                )
            )
            tree = root / "tree.nwk"
            tree.write_text(
                "(((QUERY:0.01,REF0001:0.01)95:0.02,REF0002:0.04)90:0.02,"
                "REF0003:0.30);\n"
            )
            assignment_file = root / "assignment.tsv"
            neighbors_file = root / "neighbors.tsv"
            assignment = classify_tree(
                tree_file=tree,
                references_file=references,
                task_file=task,
                assignment_output=assignment_file,
                neighbors_output=neighbors_file,
                assignment_neighbors=2,
            )
            with assignment_file.open(newline="") as handle:
                assignment_reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(assignment_reader.fieldnames, TREE_ASSIGNMENT_FIELDS)
            with neighbors_file.open(newline="") as handle:
                neighbor_reader = csv.DictReader(handle, delimiter="\t")
                neighbor_rows = list(neighbor_reader)
                self.assertEqual(neighbor_reader.fieldnames, TREE_NEIGHBOR_FIELDS)

        self.assertEqual(assignment["tree_taxonomy"], lineage)
        self.assertEqual(assignment["tree_taxonomy_source"], "PR2")
        self.assertEqual(assignment["tree_basis_neighbors"], "2")
        self.assertEqual(assignment["tree_route_16s_best_bitscore"], "")
        self.assertEqual(assignment["tree_route_18s_best_bitscore"], "5710")
        self.assertEqual(assignment["tree_query_edge_support"], "95")
        self.assertEqual(neighbor_rows[0]["blast_sseqid"], "img_exact")
        self.assertEqual(neighbor_rows[0]["tree_lineage_basis"], "centroid_taxonomy")
        self.assertEqual(neighbor_rows[0]["used_for_assignment"], "true")

    def test_tree_taxonomy_includes_equal_distance_boundary_neighbors(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            references = root / "references.tsv"
            rows = [
                reference_row(
                    "REF0001",
                    "amoeba_1",
                    1,
                    taxonomy="Eukaryota;Amoebozoa;Discosea;Echinamoebida",
                ),
                reference_row(
                    "REF0002",
                    "amoeba_2",
                    2,
                    taxonomy="Eukaryota;Amoebozoa;Discosea;Echinamoebida",
                ),
                reference_row(
                    "REF0003",
                    "amoeba_3",
                    3,
                    taxonomy="Eukaryota;Amoebozoa;Tubulinea",
                ),
            ]
            with references.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=REFERENCE_FIELDS, delimiter="\t"
                )
                writer.writeheader()
                writer.writerows(rows)
            task = root / "task.json"
            task.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "name": "query1",
                        "sample": "sample",
                        "detected_model": "RF01960",
                        "tree_model": "RF01960",
                        "tree_marker": "18S",
                        "tree_route_decision": "majority_global_top_hits",
                        "tree_route_16s_votes": 0,
                        "tree_route_18s_votes": 100,
                    }
                )
            )
            tree = root / "tree.nwk"
            tree.write_text(
                "(QUERY:0.01,REF0001:0.01,REF0002:0.01,REF0003:0.01);\n"
            )
            assignment = classify_tree(
                tree_file=tree,
                references_file=references,
                task_file=task,
                assignment_output=root / "assignment.tsv",
                neighbors_output=root / "neighbors.tsv",
                assignment_neighbors=2,
            )

        self.assertEqual(assignment["tree_taxonomy"], "Eukaryota;Amoebozoa")
        self.assertEqual(assignment["tree_basis_neighbors"], "3")


if __name__ == "__main__":
    unittest.main()
