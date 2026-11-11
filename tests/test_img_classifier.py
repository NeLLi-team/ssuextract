import csv
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import build_database_release as builder
import classify_img_clusters as classifier  # noqa: E402  # Add the scripts directory before local imports.
import database_contracts  # noqa: E402  # Add the scripts directory before local imports.
import database_release_io  # noqa: E402  # Add the scripts directory before local imports.
import img_chunked_search as chunked  # noqa: E402  # Add the scripts directory before local imports.
import img_classification_data as classification_data  # noqa: E402  # Add the scripts directory before local imports.
import img_search_provenance as search_provenance

EUK_PREFIX = (
    "Eukaryota",
    "TSAR",
    "Alveolata",
    "Dinoflagellata",
    "Dinophyceae",
    "Gymnodiniales",
    "Gymnodiniaceae",
    "Gymnodinium",
)


def clusters(
    *rows: tuple[str, str, list[str]],
) -> tuple[classification_data.Cluster, ...]:
    text = io.StringIO()
    writer = csv.writer(text, delimiter="\t", lineterminator="\n")
    writer.writerow(["cluster_id", "centroid", "legacy_taxonomy", "sequences"])
    for cluster_id, centroid, members in rows:
        writer.writerow([cluster_id, centroid, "Wrong;Legacy;Taxonomy", repr(members)])
    text.seek(0)
    return classification_data.parse_clusters(text)


def hits(
    *rows: tuple[object, ...],
) -> tuple[tuple[str, tuple[classification_data.BlastHit, ...]], ...]:
    expanded = []
    for query, subject, pident, qcov, bitscore in rows:
        coverage = float(qcov)
        query_length = 100 if coverage.is_integer() else 10_000
        alignment_length = round(query_length * coverage / 100)
        expanded.append(
            (
                query,
                subject,
                pident,
                alignment_length,
                query_length,
                query_length,
                qcov,
                bitscore,
            )
        )
    return tuple(
        classification_data.parse_blast_hits(
            "\t".join(map(str, row)) + "\n" for row in expanded
        )
    )


def taxonomy(
    path: tuple[str, ...], source: str = "PR2", compartment: str = "nucleus"
) -> classification_data.TaxonomyRecord:
    return classification_data.TaxonomyRecord(path, source, path[0], compartment)


def calibrated(rank_cap: int) -> classification_data.CalibrationStratum:
    return classification_data.CalibrationStratum("calibrated", rank_cap, "")


def failed(
    reason: str = "insufficient_domain_calls",
) -> classification_data.CalibrationStratum:
    return classification_data.CalibrationStratum("failed", None, reason)


def calibration_document(
    rank_caps: dict[str, int],
    strata: dict[str, dict[str, object]],
    rank_rules: dict[str, list[dict[str, object]]] | None = None,
) -> dict[str, object]:
    return {
        "schema_version": classification_data.CALIBRATION_SCHEMA_VERSION,
        "status": "complete" if rank_caps else "failed",
        "calibration_use_case": "img_centroid",
        "stratum_basis": "predicted_candidate_source_and_lca_domain",
        "classification_policy": classification_data.classification_policy(),
        "rank_caps": rank_caps,
        "rank_rules": rank_rules or {},
        "strata": strata,
    }


class BlastParsingTests(unittest.TestCase):
    def test_parser_streams_query_groups_and_rejects_reappearing_query(self) -> None:
        rows = [
            "C1\tref\t100\t100\t100\t100\t100\t100\n",
            "C2\tref\t100\t100\t100\t100\t100\t100\n",
            "C1\tother\t100\t100\t100\t100\t100\t100\n",
        ]
        parsed = classification_data.parse_blast_hits(rows)
        self.assertEqual(next(parsed)[0], "C1")
        with self.assertRaisesRegex(ValueError, "reappears"):
            list(parsed)

    def test_grouped_stream_matches_baseline_multi_query_outcomes(
        self,
    ) -> None:
        rows = [
            "C2\tref_a\t99\t100\t100\t100\t100\t100\n",
            "C2\tref_b\t100\t100\t100\t100\t79\t101\n",
            "C1\tref_b\t100\t100\t100\t100\t100\t100\n",
            "C1\tref_a\t99\t100\t100\t100\t100\t98\n",
        ]
        cluster_rows = clusters(
            ("C1", "first", ["IMG_1"]),
            ("C2", "second", ["IMG_2"]),
            ("C3", "third", []),
        )
        records = {
            "ref_a": taxonomy(("Bacteria", "Firmicutes"), "SILVA", ""),
            "ref_b": taxonomy(("Bacteria", "Proteobacteria"), "SILVA", ""),
        }
        streamed = classifier.classify_clusters(
            cluster_rows, classification_data.parse_blast_hits(rows), records
        )
        # SHA256 of the two classified outcomes produced at commit 58ccc07.
        expected_digest = (
            "dbe8795d3a0f73f4182c8cb14fd15f321deb6b518e41275b19fd21ffd70f6fce"
        )
        rendered = json.dumps(streamed[1][:2], sort_keys=True, separators=(",", ":"))
        self.assertEqual(hashlib.sha256(rendered.encode()).hexdigest(), expected_digest)
        self.assertEqual(streamed[1][2]["reason"], "no_img_members")

    def test_one_hsp_per_subject_is_enforced(self) -> None:
        with self.assertRaisesRegex(ValueError, "Multiple HSPs"):
            hits(
                ("C1", "subject", 100, 100, 50),
                ("C1", "subject", 99, 100, 49),
            )


class ClassificationPolicyTests(unittest.TestCase):
    def test_all_outcomes_have_evidence_ids_and_assignments_reuse_them(self) -> None:
        assignments, outcomes, _ = classifier.classify_clusters(
            clusters(
                ("C1", "first", ["IMG_1"]),
                ("C2", "second", []),
            ),
            (),
            {},
            marker="16S",
        )
        self.assertTrue(
            all(row["evidence_id"].startswith("IMGEV_") for row in outcomes)
        )
        self.assertEqual(assignments[0]["evidence_id"], outcomes[0]["evidence_id"])
        self.assertEqual(outcomes[1]["reason"], "no_img_members")

    def test_centroid_name_is_not_a_blast_query_alias(self) -> None:
        with self.assertRaisesRegex(ValueError, "cluster ID"):
            classifier.classify_clusters(
                clusters(("C1", "centroid", ["IMG_1"])),
                hits(("centroid", "ref", 100, 100, 100)),
                {"ref": taxonomy(("Bacteria",), "SILVA", "")},
            )

    def test_centroid_display_name_removes_only_silva_lineage_text(self) -> None:
        self.assertEqual(
            classifier.centroid_name(
                "REF_SILVA_AB123.1.1500;Bacteria;Bacteroidota"
            ),
            "REF_SILVA_AB123.1.1500",
        )
        spr = "REF_SPR_KC486120.1.1650_Eukaryota_SAR_Alveolata"
        self.assertEqual(classifier.centroid_name(spr), spr)
        with self.assertRaisesRegex(ValueError, "invalid display name"):
            classifier.centroid_name("IMG_centroid|second")

    def test_tied_candidates_use_lca_and_ignore_input_order(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_2", "IMG_1"]))
        records = {
            "ref_a": taxonomy(EUK_PREFIX[:-2] + ("FamilyA", "GenusA", "Species_a")),
            "ref_b": taxonomy(EUK_PREFIX[:-2] + ("FamilyB", "GenusB", "Species_b")),
        }
        forward = hits(
            ("C1", "ref_b", 99, 100, 200),
            ("C1", "ref_a", 98, 100, 200),
        )
        reverse = hits(
            ("C1", "ref_a", 98, 100, 200),
            ("C1", "ref_b", 99, 100, 200),
        )

        first = classifier.classify_clusters(cluster_rows, forward, records)
        second = classifier.classify_clusters(cluster_rows, reverse, records)

        self.assertEqual(first, second)
        assignments, outcomes, _ = first
        self.assertEqual([row["source_identifier"] for row in assignments], ["IMG_1", "IMG_2"])
        self.assertEqual(outcomes[0]["taxonomy"], ";".join(EUK_PREFIX[:-2]))
        self.assertEqual(outcomes[0]["candidate_count"], 2)

    def test_coverage_filter_is_applied_before_best_score_window(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        records = {
            "low_coverage": taxonomy(("Eukaryota", "Discarded")),
            "eligible": taxonomy(("Bacteria", "Firmicutes"), "SILVA", ""),
        }
        blast_hits = hits(
            ("C1", "low_coverage", 100, 79.99, 1000),
            ("C1", "eligible", 97, 80, 100),
        )

        _, outcomes, qc = classifier.classify_clusters(cluster_rows, blast_hits, records)

        self.assertEqual(outcomes[0]["taxonomy"], "Bacteria;Firmicutes")
        self.assertEqual(outcomes[0]["candidate_count"], 1)
        self.assertEqual(qc["coverage_filtered_hits"], 1)

    def test_candidate_window_includes_98_percent_and_excludes_lower_hit(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        records = {
            "best": taxonomy(("Bacteria", "Firmicutes"), "SILVA", ""),
            "edge": taxonomy(("Bacteria", "Proteobacteria"), "SILVA", ""),
            "below": taxonomy(("Archaea", "Thermoproteota"), "SILVA", ""),
        }
        blast_hits = hits(
            ("C1", "below", 100, 100, 97.99),
            ("C1", "edge", 100, 100, 98),
            ("C1", "best", 100, 100, 100),
        )

        _, outcomes, _ = classifier.classify_clusters(cluster_rows, blast_hits, records)

        self.assertEqual(outcomes[0]["candidate_count"], 2)
        self.assertEqual(outcomes[0]["taxonomy"], "Bacteria")

    def test_overflow_candidate_marks_truncation_and_backs_off_to_domain(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        records = {
            "ref_a": taxonomy(EUK_PREFIX + ("Species_a",)),
            "ref_b": taxonomy(EUK_PREFIX[:-2] + ("OtherFamily", "OtherGenus", "Species_b")),
            "ref_c": taxonomy(("Eukaryota", "Other")),
        }
        blast_hits = hits(
            ("C1", "ref_a", 100, 100, 100),
            ("C1", "ref_b", 100, 100, 99),
            ("C1", "ref_c", 100, 100, 98.5),
        )

        assignments, outcomes, qc = classifier.classify_clusters(
            cluster_rows, blast_hits, records, max_targets=2
        )

        self.assertEqual(outcomes[0]["taxonomy"], "Eukaryota")
        self.assertTrue(outcomes[0]["truncated"])
        self.assertEqual(assignments[0]["compartment"], "")
        self.assertEqual(qc["clusters_truncated"], 1)

    def test_last_returned_hit_outside_window_is_not_truncated(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        species_path = EUK_PREFIX + ("Species_a",)
        records = {
            "best": taxonomy(species_path),
            "outside": taxonomy(("Eukaryota", "Other")),
        }
        blast_hits = hits(
            ("C1", "best", 100, 100, 100),
            ("C1", "outside", 100, 100, 97.99),
        )

        _, outcomes, qc = classifier.classify_clusters(
            cluster_rows, blast_hits, records, max_targets=2
        )

        self.assertEqual(outcomes[0]["taxonomy"], ";".join(species_path))
        self.assertFalse(outcomes[0]["truncated"])
        self.assertNotIn("clusters_truncated", qc)

    def test_low_coverage_raw_boundary_still_marks_truncation(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        records = {
            "best": taxonomy(("Bacteria", "Firmicutes"), "SILVA", ""),
            "second": taxonomy(("Bacteria", "Firmicutes"), "SILVA", ""),
            "filtered": taxonomy(("Bacteria", "Proteobacteria"), "SILVA", ""),
        }
        blast_hits = classification_data.parse_blast_hits(
            [
                "C1\tbest\t100\t100\t100\t100\t100\t100\n",
                "C1\tfiltered\t100\t100\t100\t100\t79\t99\n",
                "C1\tsecond\t100\t100\t100\t100\t100\t99\n",
            ]
        )

        _, outcomes, qc = classifier.classify_clusters(
            cluster_rows, blast_hits, records, max_targets=2
        )

        self.assertEqual(outcomes[0]["taxonomy"], "Bacteria")
        self.assertTrue(outcomes[0]["truncated"])
        self.assertEqual(qc["clusters_truncated"], 1)

    def test_species_requires_every_candidate_to_be_exact_and_agree(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        species_path = EUK_PREFIX + ("Species_a",)
        records = {
            "ref_a": taxonomy(species_path),
            "ref_b": taxonomy(species_path),
        }
        exact_hits = hits(
            ("C1", "ref_a", 100, 100, 100),
            ("C1", "ref_b", 100, 100, 99),
        )
        non_exact_hits = hits(
            ("C1", "ref_a", 100, 100, 100),
            ("C1", "ref_b", 99.9, 100, 99),
        )
        longer_subject_hits = classification_data.parse_blast_hits(
            ["C1\tref_a\t100\t100\t100\t101\t100\t100\n"]
        )

        _, exact_outcomes, _ = classifier.classify_clusters(
            cluster_rows, exact_hits, records
        )
        _, guarded_outcomes, _ = classifier.classify_clusters(
            cluster_rows, non_exact_hits, records
        )
        _, longer_subject_outcomes, _ = classifier.classify_clusters(
            cluster_rows, longer_subject_hits, records
        )

        self.assertEqual(exact_outcomes[0]["taxonomy"], ";".join(species_path))
        self.assertTrue(exact_outcomes[0]["species_called"])
        self.assertEqual(guarded_outcomes[0]["taxonomy"], ";".join(EUK_PREFIX))
        self.assertTrue(guarded_outcomes[0]["species_guard_applied"])
        self.assertFalse(guarded_outcomes[0]["species_called"])
        self.assertFalse(longer_subject_outcomes[0]["species_called"])

    def test_calibrated_rank_cap_limits_non_exact_assignment(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        path = EUK_PREFIX + ("Species_a",)
        records = {"ref": taxonomy(path)}
        blast_hits = hits(("C1", "ref", 99, 100, 100))
        assignments, outcomes, _ = classifier.classify_clusters(
            cluster_rows,
            blast_hits,
            records,
            marker="18S",
            calibration_strata={"18S|PR2|Eukaryota": calibrated(6)},
        )
        self.assertEqual(outcomes[0]["taxonomy"], ";".join(path[:7]))
        self.assertEqual(outcomes[0]["calibration_rank_cap"], 6)
        self.assertEqual(assignments[0]["taxonomy"], ";".join(path[:7]))

    def test_identity_rules_require_every_candidate_source(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        path = ("Eukaryota", "TSAR", "Alveolata")
        records = {
            "pr2": taxonomy(path, "PR2"),
            "silva": taxonomy(path, "SILVA", ""),
        }
        strata = {
            "16S|PR2|Eukaryota": calibrated(1),
            "16S|SILVA|Eukaryota": calibrated(0),
        }
        rules = {
            "16S|PR2|Eukaryota": (
                classification_data.CalibrationRankRule(0, Decimal(95)),
                classification_data.CalibrationRankRule(1, Decimal(97)),
            ),
            "16S|SILVA|Eukaryota": (
                classification_data.CalibrationRankRule(0, Decimal(90)),
            ),
        }
        assignments, outcomes, _ = classifier.classify_clusters(
            cluster_rows,
            hits(
                ("C1", "pr2", 96, 100, 100),
                ("C1", "silva", 96, 100, 100),
            ),
            records,
            marker="16S",
            calibration_strata=strata,
            calibration_rank_rules=rules,
        )
        self.assertEqual(assignments[0]["taxonomy"], "Eukaryota")
        self.assertEqual(outcomes[0]["calibration_rank_cap"], 0)
        self.assertEqual(outcomes[0]["calibration_identity_threshold"], "95")

        assignments, outcomes, _ = classifier.classify_clusters(
            cluster_rows,
            hits(
                ("C1", "pr2", 94, 100, 100),
                ("C1", "silva", 94, 100, 100),
            ),
            records,
            marker="16S",
            calibration_strata=strata,
            calibration_rank_rules=rules,
        )
        self.assertEqual(assignments[0]["taxonomy"], "Unclassified")
        self.assertEqual(outcomes[0]["reason"], "calibration_identity_below_domain")

    def test_cluster_propagation_cap_overrides_exact_species_call(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        species_path = EUK_PREFIX + ("Species_a",)
        assignments, outcomes, _ = classifier.classify_clusters(
            cluster_rows,
            hits(("C1", "ref", 100, 100, 100)),
            {"ref": taxonomy(species_path)},
            marker="18S",
            calibration_strata={"18S|PR2|Eukaryota": calibrated(7)},
            propagation_rank_cap=0,
        )
        self.assertEqual(outcomes[0]["taxonomy"], "Eukaryota")
        self.assertEqual(outcomes[0]["centroid_taxonomy"], ";".join(species_path))
        self.assertFalse(outcomes[0]["species_called"])
        self.assertEqual(assignments[0]["taxonomy"], "Eukaryota")
        self.assertEqual(assignments[0]["centroid"], "centroid")
        self.assertEqual(assignments[0]["centroid_name"], "centroid")
        self.assertEqual(
            assignments[0]["centroid_taxonomy"], ";".join(species_path)
        )
        self.assertEqual(assignments[0]["centroid_taxonomy_source"], "PR2")

    def test_failed_pr2_16s_stratum_unclassifies_exact_and_nonexact_hits(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        species_path = EUK_PREFIX + ("Species_a",)
        strata = {"16S|PR2|Eukaryota": failed()}
        source_records = (
            database_contracts.SourceRecord(
                "SRC_1", "SSU_1", "IMG", "2025", "IMG_1", "IMG_1", "16S"
            ),
        )
        for percent_identity in (100, 99):
            with self.subTest(percent_identity=percent_identity):
                assignments, outcomes, qc = classifier.classify_clusters(
                    cluster_rows,
                    hits(("C1", "ref", percent_identity, 100, 100)),
                    {"ref": taxonomy(species_path)},
                    marker="16S",
                    calibration_strata=strata,
                )
                self.assertEqual(outcomes[0]["classification_status"], "unclassified")
                self.assertEqual(outcomes[0]["reason"], "calibration_stratum_failed")
                self.assertEqual(assignments[0]["taxonomy"], "Unclassified")
                self.assertEqual(assignments[0]["marker"], "16S")
                self.assertEqual(
                    qc["unclassified_clusters_by_reason"],
                    {"calibration_stratum_failed": 1},
                )
                derived = builder.ingest_derived_cluster_assignments(
                    assignments, source_records
                )
                self.assertEqual(derived[0].taxonomy, ("Unclassified",))
        model = builder.build_deduplicated_model(
            [
                database_contracts.PreparedSourceRecord(
                    "PR2",
                    "5.1.1",
                    "pr2-ref",
                    "pr2-ref",
                    "AAAA",
                    "16S",
                    species_path,
                    "PR2",
                    "nucleus",
                ),
                database_contracts.PreparedSourceRecord(
                    "IMG",
                    "2025",
                    "IMG_1",
                    "IMG_1",
                    "AAAA",
                    "16S",
                    taxon_oid="1",
                ),
            ]
        )
        derived = builder.ingest_derived_cluster_assignments(
            assignments, model.source_records
        )
        updated = builder.add_taxonomy_assignments(model, derived)
        self.assertEqual(updated.preferred_taxonomy[0].taxonomy, species_path)
        self.assertEqual(updated.preferred_taxonomy[0].assignment_method, "native")

    def test_propagation_cap_cannot_override_failed_calibration(self) -> None:
        assignments, outcomes, _ = classifier.classify_clusters(
            clusters(("C1", "centroid", ["IMG_1"])),
            hits(("C1", "ref", 100, 100, 100)),
            {"ref": taxonomy(EUK_PREFIX + ("Species_a",))},
            marker="16S",
            calibration_strata={"16S|PR2|Eukaryota": failed()},
            propagation_rank_cap=0,
        )
        self.assertEqual(outcomes[0]["reason"], "calibration_stratum_failed")
        self.assertEqual(assignments[0]["taxonomy"], "Unclassified")

    def test_mixed_sources_fail_if_any_required_stratum_failed(self) -> None:
        path = ("Eukaryota", "TSAR")
        assignments, outcomes, _ = classifier.classify_clusters(
            clusters(("C1", "centroid", ["IMG_1"])),
            hits(
                ("C1", "pr2", 99, 100, 100),
                ("C1", "silva", 99, 100, 100),
            ),
            {
                "pr2": taxonomy(path, "PR2"),
                "silva": taxonomy(path, "SILVA", ""),
            },
            marker="16S",
            calibration_strata={
                "16S|PR2|Eukaryota": failed("domain_precision_below_threshold"),
                "16S|SILVA|Eukaryota": calibrated(0),
            },
        )
        self.assertEqual(assignments[0]["taxonomy"], "Unclassified")
        self.assertEqual(
            outcomes[0]["failed_calibration_strata"],
            [
                {
                    "key": "16S|PR2|Eukaryota",
                    "reason": "domain_precision_below_threshold",
                }
            ],
        )

    def test_missing_calibration_stratum_abstains(self) -> None:
        assignments, outcomes, _ = classifier.classify_clusters(
            clusters(("C1", "centroid", ["IMG_1"])),
            hits(("C1", "ref", 99, 100, 100)),
            {"ref": taxonomy(("Eukaryota", "TSAR"))},
            marker="16S",
            calibration_strata={"16S|SILVA|Bacteria": calibrated(0)},
        )
        self.assertEqual(assignments[0]["taxonomy"], "Unclassified")
        self.assertEqual(outcomes[0]["reason"], "calibration_stratum_missing")

    def test_no_hit_and_filtered_hit_are_explicitly_unclassified(self) -> None:
        cluster_rows = clusters(
            ("C1", "centroid_1", ["IMG_1"]),
            ("C2", "centroid_2", ["IMG_2"]),
        )
        blast_hits = hits(("C2", "ref", 100, 79, 100))

        assignments, outcomes, qc = classifier.classify_clusters(
            cluster_rows, blast_hits, {"ref": taxonomy(("Eukaryota", "TSAR"))}
        )

        self.assertEqual(len(assignments), 2)
        self.assertEqual(
            {row["assignment_method"] for row in assignments},
            {"updated_reference_unclassified"},
        )
        self.assertEqual({row["taxonomy"] for row in assignments}, {"Unclassified"})
        self.assertEqual(
            [outcome["classification_status"] for outcome in outcomes],
            ["unclassified", "unclassified"],
        )
        self.assertEqual(
            qc["unclassified_clusters_by_reason"],
            {"no_hits": 1, "no_hits_at_minimum_query_coverage": 1},
        )

    def test_cross_domain_candidate_set_is_unclassified(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        blast_hits = hits(
            ("C1", "bacteria", 100, 100, 100),
            ("C1", "eukaryote", 100, 100, 100),
        )
        records = {
            "bacteria": taxonomy(("Bacteria", "Firmicutes"), "SILVA", ""),
            "eukaryote": taxonomy(("Eukaryota", "TSAR")),
        }

        assignments, outcomes, _ = classifier.classify_clusters(
            cluster_rows, blast_hits, records
        )

        self.assertEqual(assignments[0]["assignment_method"], "updated_reference_unclassified")
        self.assertEqual(outcomes[0]["reason"], "no_common_domain")

    def test_ambiguous_reference_subject_is_explicitly_unclassified(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        blast_hits = hits(("C1", "ambiguous", 100, 100, 100))
        alternatives = json.dumps(
            [
                {"domain": "Bacteria", "taxonomy": "Bacteria;Firmicutes"},
                {"domain": "Eukaryota", "taxonomy": "Eukaryota;TSAR"},
            ]
        )
        records = {
            "ambiguous": classification_data.TaxonomyRecord(
                (), "PR2+SILVA", "ambiguous", "mixed", True, alternatives
            )
        }

        assignments, outcomes, qc = classifier.classify_clusters(
            cluster_rows, blast_hits, records, marker="16S"
        )

        self.assertEqual(assignments[0]["assignment_method"], "updated_reference_unclassified")
        self.assertEqual(outcomes[0]["reason"], "ambiguous_reference_taxonomy")
        self.assertEqual(
            outcomes[0]["candidates"][0]["taxonomy_alternatives"][0]["domain"],
            "Bacteria",
        )
        self.assertEqual(
            qc["unclassified_clusters_by_reason"],
            {"ambiguous_reference_taxonomy": 1},
        )

    def test_missing_reference_taxonomy_is_explicitly_unclassified(self) -> None:
        cluster_rows = clusters(("C1", "centroid", ["IMG_1"]))
        blast_hits = hits(("C1", "missing", 100, 100, 100))

        assignments, outcomes, _ = classifier.classify_clusters(
            cluster_rows, blast_hits, {}
        )

        self.assertEqual(assignments[0]["assignment_method"], "updated_reference_unclassified")
        self.assertEqual(outcomes[0]["reason"], "missing_reference_taxonomy")


class InputValidationTests(unittest.TestCase):
    def test_calibration_loader_accepts_an_explicit_failed_stratum(self) -> None:
        key = "16S|PR2|Eukaryota"
        document = calibration_document(
            {},
            {
                key: {
                    "status": "failed",
                    "rank_cap": None,
                    "reason": "insufficient_domain_calls",
                    "metrics": [{"rank": "domain", "called": 42}],
                }
            },
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calibration.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            loaded = classification_data.load_calibration(path)
        self.assertEqual(loaded.strata[key], failed())
        self.assertEqual(loaded.rank_rules, {})

    def test_calibration_loader_accepts_identity_rules(self) -> None:
        key = "16S|SILVA|Bacteria"
        document = calibration_document(
            {key: 1},
            {key: {"status": "calibrated", "rank_cap": 1, "reason": ""}},
            {
                key: [
                    {"rank_index": 0, "min_candidate_identity": 90},
                    {"rank_index": 1, "min_candidate_identity": 97.0},
                ]
            },
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calibration.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            loaded = classification_data.load_calibration(path)
        self.assertEqual(
            loaded.rank_rules[key],
            (
                classification_data.CalibrationRankRule(0, Decimal(90)),
                classification_data.CalibrationRankRule(1, Decimal("97.0")),
            ),
        )

    def test_calibration_loader_rejects_invalid_schema_v3_contracts(self) -> None:
        key = "16S|PR2|Eukaryota"
        failed_document = calibration_document(
            {}, {key: {"status": "failed", "rank_cap": None, "reason": "failed"}}
        )
        invalid_documents = {
            "old schema": {**failed_document, "schema_version": 2},
            "wrong use case": {
                **failed_document,
                "calibration_use_case": "runtime_query",
            },
            "wrong stratum basis": {
                **failed_document,
                "stratum_basis": "true_query_source_and_domain",
            },
            "wrong policy": {
                **failed_document,
                "classification_policy": {
                    **classification_data.classification_policy(),
                    "min_query_coverage_percent": "79",
                },
            },
            "unknown status": {
                **failed_document,
                "strata": {
                    key: {"status": "unknown", "rank_cap": None, "reason": "failed"}
                },
            },
            "failed with cap": {
                **failed_document,
                "status": "complete",
                "rank_caps": {key: 0},
                "strata": {
                    key: {"status": "failed", "rank_cap": 0, "reason": "failed"}
                },
            },
            "rank caps mismatch": {
                **failed_document,
                "rank_caps": {},
                "strata": {key: {"status": "calibrated", "rank_cap": 0, "reason": ""}},
            },
            "wrong root status": {**failed_document, "status": "complete"},
        }
        calibrated_document = calibration_document(
            {key: 1},
            {key: {"status": "calibrated", "rank_cap": 1, "reason": ""}},
        )
        invalid_documents.update(
            {
                "rule for failed stratum": {
                    **failed_document,
                    "rank_rules": {
                        key: [{"rank_index": 0, "min_candidate_identity": 90}]
                    },
                },
                "noncontiguous rules": {
                    **calibrated_document,
                    "rank_rules": {
                        key: [
                            {"rank_index": 0, "min_candidate_identity": 90},
                            {"rank_index": 2, "min_candidate_identity": 95},
                        ]
                    },
                },
                "decreasing rules": {
                    **calibrated_document,
                    "rank_rules": {
                        key: [
                            {"rank_index": 0, "min_candidate_identity": 95},
                            {"rank_index": 1, "min_candidate_identity": 90},
                        ]
                    },
                },
                "boolean threshold": {
                    **calibrated_document,
                    "rank_rules": {
                        key: [
                            {"rank_index": 0, "min_candidate_identity": True},
                            {"rank_index": 1, "min_candidate_identity": 95},
                        ]
                    },
                },
                "rules do not reach cap": {
                    **calibrated_document,
                    "rank_rules": {
                        key: [{"rank_index": 0, "min_candidate_identity": 90}]
                    },
                },
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calibration.json"
            for label, document in invalid_documents.items():
                with self.subTest(label=label):
                    path.write_text(json.dumps(document), encoding="utf-8")
                    with self.assertRaises(ValueError):
                        classification_data.load_calibration(path)

    def test_short_cluster_row_has_structured_error(self) -> None:
        table = io.StringIO("cluster_id\tcentroid\tsequences\nC1\n")

        with self.assertRaisesRegex(ValueError, "lacks cluster_id or centroid"):
            classification_data.parse_clusters(table)

    def test_internal_empty_taxonomy_rank_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "empty rank"):
            classification_data._taxonomy_record(
                {
                    "taxonomy": "Bacteria;;Firmicutes",
                    "taxonomy_source": "SILVA",
                    "domain": "Bacteria",
                },
                "ref",
            )


class OutputContractTests(unittest.TestCase):
    def test_cli_requires_release_classification_inputs(self) -> None:
        arguments = [
            "--blast",
            "hits.m8",
            "--blast-fetch-targets",
            "501",
            "--taxonomy",
            "taxonomy.parquet",
            "--clusters",
            "clusters.tsv",
            "--assignments-tsv",
            "assignments.tsv",
            "--outcomes-jsonl",
            "outcomes.jsonl",
            "--qc-json",
            "qc.json",
            "--marker",
            "16S",
            "--calibration",
            "calibration.json",
            "--search-provenance",
            "search.json",
            "--propagation-rank-cap",
            "0",
        ]
        for option in (
            "--marker",
            "--calibration",
            "--search-provenance",
            "--propagation-rank-cap",
        ):
            with self.subTest(option=option):
                missing = list(arguments)
                index = missing.index(option)
                del missing[index : index + 2]
                with self.assertRaises(SystemExit) as raised:
                    classifier._parser().parse_args(missing)
                self.assertEqual(raised.exception.code, 2)

    def test_shell_rejects_header_only_assignments(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stub = root / "python"
            stub.write_text(
                "#!/usr/bin/env bash\n"
                "set -eu\n"
                'if [[ "$1" == scripts/img_search_provenance.py ]]; then exit 0; fi\n'
                "while (($#)); do\n"
                '  case "$1" in\n'
                "    --assignments-tsv)\n"
                '      printf "marker\\tsource_identifier\\n" > "$2";;\n'
                '    --outcomes-jsonl|--qc-json) printf "{}\\n" > "$2";;\n'
                "  esac\n"
                "  shift\n"
                "done\n",
                encoding="utf-8",
            )
            stub.chmod(0o755)
            environment = dict(os.environ, PATH=f"{root}:{os.environ['PATH']}")
            completed = subprocess.run(
                [
                    "bash",
                    "scripts/classify_img_marker.sh",
                    "16S",
                    "clusters.tsv",
                    "source.fna",
                    "profile",
                    "calibration.json",
                    str(root / "output"),
                    "1",
                ],
                cwd=REPO,
                env=environment,
                check=False,
                capture_output=True,
                text=True,
            )
        self.assertEqual(completed.returncode, 1, completed.stderr)

    def test_cli_rejects_missing_overflow_target_contract(self) -> None:
        with self.assertRaisesRegex(ValueError, "must be 501"):
            classifier.main(
                [
                    "--blast",
                    "missing.m8",
                    "--blast-fetch-targets",
                    "500",
                    "--taxonomy",
                    "missing.parquet",
                    "--clusters",
                    "missing.tsv",
                    "--assignments-tsv",
                    "assignments.tsv",
                    "--outcomes-jsonl",
                    "outcomes.jsonl",
                    "--qc-json",
                    "qc.json",
                    "--marker",
                    "16S",
                    "--calibration",
                    "calibration.json",
                    "--search-provenance",
                    "search.json",
                    "--propagation-rank-cap",
                    "0",
                ]
            )

    def test_cli_writes_assignments_outcomes_and_qc(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            blast = root / "centroids.m8"
            cluster_tsv = root / "clusters.tsv"
            outputs = [
                root / "assignments.tsv",
                root / "outcomes.jsonl",
                root / "qc.json",
            ]
            calibration = root / "calibration.json"
            calibration.write_text(
                json.dumps(
                    calibration_document(
                        {"16S|SILVA|Bacteria": 0},
                        {
                            "16S|SILVA|Bacteria": {
                                "status": "calibrated",
                                "rank_cap": 0,
                                "reason": "",
                            }
                        },
                    )
                ),
                encoding="utf-8",
            )
            blast.write_text("C1\tref\t100\t4\t4\t4\t100\t100\n")
            cluster_tsv.write_text(
                "cluster_id\tcentroid\tsequences\n"
                "C1\tcentroid\t['IMG_1']\n"
            )
            profile = root / "curated"
            (profile / "tables").mkdir(parents=True)
            (profile / "blast").mkdir()
            (profile / "manifest.json").write_text("{}\n", encoding="utf-8")
            preferred = profile / "tables/preferred_taxonomy.parquet"
            database_release_io.write_preferred_taxonomy_parquet(
                preferred,
                [
                    database_contracts.PreferredTaxonomy(
                        sequence_id="ref",
                        reference_source="SILVA",
                        taxonomy=("Bacteria", "Firmicutes"),
                        taxonomy_source="SILVA",
                        domain="Bacteria",
                        compartment="",
                        assignment_method="native_reference",
                        cross_domain_conflict=False,
                    )
                ],
            )
            source_fasta = root / "source.fna"
            centroids_fasta = root / "centroids.fna"
            source_fasta.write_text(">IMG_1\nACGT\n", encoding="ascii")
            centroids_fasta.write_text(">C1\nACGT\n", encoding="ascii")
            search_sidecar = root / "search_provenance.json"
            for suffix in ("nhr", "nin", "nsq"):
                (profile / "blast" / f"16S.{suffix}").write_bytes(b"database")
            artifacts = [
                {
                    "path": f"blast/16S.{suffix}",
                    "bytes": 8,
                    "sha256": hashlib.sha256(b"database").hexdigest(),
                }
                for suffix in ("nhr", "nin", "nsq")
            ]
            (profile / "manifest.json").write_text(
                json.dumps({"artifacts": artifacts}), encoding="utf-8"
            )
            chunk_directory = root / "chunks"
            chunked.prepare_chunks("16S", centroids_fasta, chunk_directory, 1)
            tool = {"path": "/fake/blastn", "sha256": "a" * 64, "version": "test"}

            def fake_blast(command: list[str], check: bool) -> None:
                Path(command[command.index("-out") + 1]).write_bytes(blast.read_bytes())

            with (
                mock.patch.object(chunked, "_blastn_tool", return_value=tool),
                mock.patch.object(chunked.subprocess, "run", side_effect=fake_blast),
            ):
                chunked.run_chunk(
                    "16S", profile, centroids_fasta, chunk_directory, 0, 8
                )
                raw_search = chunked.finalize_search(
                    "16S",
                    profile,
                    cluster_tsv,
                    source_fasta,
                    centroids_fasta,
                    blast,
                    search_sidecar,
                    chunk_directory,
                    8,
                )
            portable_search = search_provenance.portable_provenance(
                raw_search, classifier._file_sha256(search_sidecar)
            )

            arguments = [
                "--blast",
                str(blast),
                "--blast-fetch-targets",
                str(classifier.BLAST_FETCH_TARGETS),
                "--taxonomy",
                str(preferred),
                "--clusters",
                str(cluster_tsv),
                "--marker",
                "16S",
                "--calibration",
                str(calibration),
                "--search-provenance",
                str(search_sidecar),
                "--propagation-rank-cap",
                "0",
                "--assignments-tsv",
                str(outputs[0]),
                "--outcomes-jsonl",
                str(outputs[1]),
                "--qc-json",
                str(outputs[2]),
            ]
            result = classifier.main(arguments)

            self.assertEqual(result, 0)
            self.assertTrue(all(path.is_file() and path.stat().st_size for path in outputs))
            qc = json.loads(outputs[2].read_text())
            self.assertEqual(qc["assignments_written"], 1)
            self.assertEqual(
                qc["classification_binding"],
                {
                    "schema_version": 1,
                    "calibration": {
                        "sha256": classifier._file_sha256(calibration),
                        "schema_version": 3,
                    },
                    "marker": "16S",
                    "propagation_rank_cap": 0,
                    "policy": classification_data.classification_policy(),
                    "outputs": {
                        "assignments_tsv": {
                            "sha256": classifier._file_sha256(outputs[0])
                        },
                        "outcomes_jsonl": {
                            "sha256": classifier._file_sha256(outputs[1])
                        },
                    },
                    "search_provenance": portable_search,
                },
            )

            mismatched_blast = root / "mismatched.m8"
            mismatched_blast.write_text(
                "C1\tref\t100\t4\t4\t4\t100\t99\n", encoding="ascii"
            )
            mismatched_arguments = list(arguments)
            mismatched_arguments[1] = str(mismatched_blast)
            with self.assertRaisesRegex(ValueError, "blast_m8"):
                classifier.main(mismatched_arguments)

    def test_current_builder_parquet_preserves_cross_domain_ambiguity(self) -> None:
        alternatives = json.dumps(
            [
                {"domain": "Bacteria", "taxonomy": "Bacteria;Firmicutes"},
                {"domain": "Eukaryota", "taxonomy": "Eukaryota;TSAR"},
            ],
            separators=(",", ":"),
            sort_keys=True,
        )
        preferred = database_contracts.PreferredTaxonomy(
            sequence_id="SSU_ambiguous",
            reference_source="PR2+SILVA",
            taxonomy=(),
            taxonomy_source="PR2+SILVA",
            domain="ambiguous",
            compartment="mixed",
            assignment_method="cross_domain_ambiguous_exact_sequence",
            cross_domain_conflict=True,
            taxonomy_alternatives=alternatives,
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "preferred_taxonomy.parquet"
            database_release_io.write_preferred_taxonomy_parquet(path, [preferred])
            with mock.patch.dict(
                os.environ,
                {"SLURM_CPUS_PER_TASK": "", "OMP_NUM_THREADS": ""},
            ):
                loaded = classification_data.load_taxonomy(path, {"SSU_ambiguous"})

        self.assertTrue(loaded["SSU_ambiguous"].cross_domain_conflict)
        self.assertEqual(
            json.loads(loaded["SSU_ambiguous"].taxonomy_alternatives)[1]["domain"],
            "Eukaryota",
        )

    def test_assignments_are_builder_ingestible_and_outputs_are_stable(self) -> None:
        cluster_rows = clusters(
            ("C1", "centroid", ["REF_legacy", "IMG_2", "IMG_1"])
        )
        records = {"ref": taxonomy(("Bacteria", "Firmicutes"), "SILVA", "")}
        blast_hits = hits(("C1", "ref", 100, 100, 100))
        assignments, outcomes, qc = classifier.classify_clusters(
            cluster_rows, blast_hits, records, marker="16S"
        )
        source_records = (
            database_contracts.SourceRecord(
                "SRC_1", "SSU_1", "IMG", "2025", "IMG_1", "IMG_1", "16S"
            ),
            database_contracts.SourceRecord(
                "SRC_2", "SSU_2", "IMG", "2025", "IMG_2", "IMG_2", "16S"
            ),
        )

        derived = builder.ingest_derived_cluster_assignments(assignments, source_records)

        self.assertEqual(len(derived), 2)
        self.assertEqual({item.taxonomy for item in derived}, {("Bacteria", "Firmicutes")})
        self.assertNotIn("Wrong", assignments[0]["taxonomy"])
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = [
                root / name for name in ("assignments.tsv", "outcomes.jsonl", "qc.json")
            ]
            with mock.patch.object(classifier.os, "fsync") as fsync:
                classifier._write_outputs(*paths, assignments, outcomes, qc)
            # Each output flushes its open handle, syncs the staged file after
            # close, then syncs the published directory entry.
            self.assertEqual(fsync.call_count, 3 * len(paths))
            with paths[0].open(newline="") as handle:
                written = list(csv.DictReader(handle, delimiter="\t"))
            outcome_json = [
                json.loads(line) for line in paths[1].read_text().splitlines()
            ]
            qc_json = json.loads(paths[2].read_text())

        self.assertEqual(written, assignments)
        self.assertEqual(outcome_json, outcomes)
        self.assertEqual(qc_json, qc)


if __name__ == "__main__":
    unittest.main()
