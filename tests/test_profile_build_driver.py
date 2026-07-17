import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import build_database_profiles as profiles


class ProfileBuildDriverTests(unittest.TestCase):
    def test_default_profile_version_tracks_taxonomy_policy_repair(self) -> None:
        args = profiles._parser().parse_args(
            [
                "curated",
                "--source-directory",
                "sources",
                "--output-root",
                "profiles",
                "--archive-directory",
                "archives",
            ]
        )
        self.assertEqual(args.version, "1.0.1")

    def test_search_and_classification_qc_bind_the_same_portable_provenance(self) -> None:
        calibration = {"sha256": "a" * 64, "schema_version": 2}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "curated"
            (profile / "tables").mkdir(parents=True)
            (profile / "blast").mkdir()
            (profile / "manifest.json").write_text("{}\n", encoding="utf-8")
            (profile / "tables/preferred_taxonomy.parquet").write_bytes(b"PAR1tax")
            clusters = root / "clusters.tsv"
            source_fasta = root / "source.fna"
            clusters.write_text("cluster_id\tcentroid\tsequences\n0\tIMG_1\t['IMG_1']\n")
            source_fasta.write_text(">IMG_1\nACGT\n", encoding="ascii")
            marker = root / "16S"
            marker.mkdir()
            centroids = marker / "centroids.fna"
            blast = marker / "centroids.m8"
            sidecar = marker / "search_provenance.json"
            assignments = marker / "assignments.tsv"
            outcomes = marker / "outcomes.jsonl"
            centroids.write_text(">0\nACGT\n", encoding="ascii")
            blast.write_text("0\tref\t100\t4\t4\t4\t100\t8\n")
            assignments.write_text("source_identifier\n", encoding="utf-8")
            outcomes.write_text("{}\n", encoding="utf-8")
            raw_search = profiles.img_search_provenance._payload(
                "16S", profile, clusters, source_fasta, centroids, blast, 8
            )
            sidecar.write_text(
                json.dumps(raw_search, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            portable = profiles.validate_classification_search(
                assignments,
                "16S",
                profile,
                clusters,
                source_fasta,
            )
            binding = {
                "schema_version": 1,
                "calibration": {"sha256": "a" * 64, "schema_version": 2},
                "marker": "16S",
                "propagation_rank_cap": 0,
                "policy": profiles.classifier.classification_policy(),
                "outputs": {
                    "assignments_tsv": {"sha256": profiles._sha256(assignments)},
                    "outcomes_jsonl": {"sha256": profiles._sha256(outcomes)},
                },
                "search_provenance": portable,
            }
            (marker / "qc.json").write_text(
                json.dumps({"classification_binding": binding}), encoding="utf-8"
            )
            self.assertEqual(
                profiles.validate_classification_qc(
                    assignments, "16S", calibration, portable
                ),
                marker / "qc.json",
            )
            self.assertNotIn("path", json.dumps(portable, sort_keys=True))
            blast.write_text("tampered\n", encoding="utf-8")
            with self.assertRaisesRegex(Exception, "blast_m8"):
                profiles.validate_classification_search(
                    assignments,
                    "16S",
                    profile,
                    clusters,
                    source_fasta,
                )

    def test_calibration_provenance_keeps_rank_caps_and_full_strata(self) -> None:
        calibrated_key = "18S|PR2|Eukaryota"
        failed_key = "16S|PR2|Eukaryota"
        strata = {
            calibrated_key: {
                "status": "calibrated",
                "rank_cap": 6,
                "reason": "",
                "sampled": 1000,
                "metrics": [{"rank": "domain", "accepted": True}],
            },
            failed_key: {
                "status": "failed",
                "rank_cap": None,
                "reason": "insufficient_domain_calls",
                "sampled": 42,
                "metrics": [{"rank": "domain", "accepted": False}],
            },
        }
        document = {
            "schema_version": 2,
            "rank_caps": {calibrated_key: 6},
            "strata": strata,
            "curated_profile": {
                "version": "1.0.0",
                "manifest_sha256": "a" * 64,
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calibration.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            provenance = profiles.calibration_provenance(path)

        self.assertEqual(provenance["schema_version"], 2)
        self.assertEqual(provenance["rank_caps"], {calibrated_key: 6})
        self.assertEqual(provenance["strata"], strata)
        self.assertEqual(provenance["curated_profile"], document["curated_profile"])
        self.assertEqual(len(provenance["sha256"]), 64)
        payload = profiles._safe_evidence_payload(
            "16S",
            {
                "classification_status": "unclassified",
                "reason": "calibration_stratum_failed",
                "failed_calibration_strata": [
                    {"key": failed_key, "reason": "insufficient_domain_calls"}
                ],
            },
        )
        self.assertEqual(
            payload["failed_calibration_strata"][0]["key"], failed_key
        )

    def test_curated_manifest_must_match_calibration_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.json"
            manifest.write_text('{"profile":"curated"}\n', encoding="utf-8")
            digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
            calibration = {
                "curated_profile": {
                    "version": "1.0.0",
                    "manifest_sha256": digest,
                }
            }
            self.assertEqual(
                profiles.validate_curated_manifest(manifest, calibration),
                {"sha256": digest},
            )
            calibration["curated_profile"]["manifest_sha256"] = "0" * 64
            with self.assertRaisesRegex(Exception, "does not match"):
                profiles.validate_curated_manifest(manifest, calibration)

    def test_img_cli_requires_curated_manifest(self) -> None:
        with self.assertRaisesRegex(SystemExit, "curated-manifest"):
            profiles.main(
                [
                    "img",
                    "--source-directory",
                    "sources",
                    "--output-root",
                    "profiles",
                    "--archive-directory",
                    "archives",
                    "--assignments",
                    "16S/assignments.tsv",
                    "--calibration",
                    "calibration.json",
                ]
            )

    def test_classification_qc_binding_rejects_ungated_or_mismatched_runs(self) -> None:
        calibration = {
            "sha256": "a" * 64,
            "schema_version": 2,
        }
        expected = {
            "schema_version": 1,
            "calibration": {"sha256": "a" * 64, "schema_version": 2},
            "marker": "16S",
            "propagation_rank_cap": 0,
            "policy": profiles.classifier.classification_policy(),
        }
        with tempfile.TemporaryDirectory() as tmp:
            marker = Path(tmp) / "16S"
            marker.mkdir()
            assignment = marker / "assignments.tsv"
            assignment.write_text("source_identifier\n", encoding="utf-8")
            outcomes = marker / "outcomes.jsonl"
            outcomes.write_text("{}\n", encoding="utf-8")
            expected["outputs"] = {
                "assignments_tsv": {
                    "sha256": hashlib.sha256(assignment.read_bytes()).hexdigest()
                },
                "outcomes_jsonl": {
                    "sha256": hashlib.sha256(outcomes.read_bytes()).hexdigest()
                },
            }
            with self.assertRaisesRegex(Exception, "classification QC"):
                profiles.validate_classification_qc(assignment, "16S", calibration)
            qc_path = marker / "qc.json"
            qc_path.write_text(
                json.dumps({"classification_binding": expected}), encoding="utf-8"
            )
            self.assertEqual(
                profiles.validate_classification_qc(assignment, "16S", calibration),
                qc_path,
            )
            assignment.write_text("substituted\n", encoding="utf-8")
            with self.assertRaisesRegex(Exception, "does not match"):
                profiles.validate_classification_qc(assignment, "16S", calibration)
            assignment.write_text("source_identifier\n", encoding="utf-8")
            outcomes.write_text('{"substituted":true}\n', encoding="utf-8")
            with self.assertRaisesRegex(Exception, "does not match"):
                profiles.validate_classification_qc(assignment, "16S", calibration)
            outcomes.write_text("{}\n", encoding="utf-8")
            mutations = (
                (
                    "wrong calibration",
                    {
                        **expected,
                        "calibration": {
                            "sha256": "b" * 64,
                            "schema_version": 2,
                        },
                    },
                ),
                ("wrong marker", {**expected, "marker": "18S"}),
                ("no propagation cap", {**expected, "propagation_rank_cap": None}),
                (
                    "wrong fetch count",
                    {
                        **expected,
                        "policy": {**expected["policy"], "blast_fetch_targets": 500},
                    },
                ),
                (
                    "wrong classification count",
                    {**expected, "policy": {**expected["policy"], "max_targets": 499}},
                ),
            )
            for label, binding in mutations:
                with self.subTest(label=label):
                    qc_path.write_text(
                        json.dumps({"classification_binding": binding}),
                        encoding="utf-8",
                    )
                    with self.assertRaisesRegex(Exception, "does not match"):
                        profiles.validate_classification_qc(
                            assignment, "16S", calibration
                        )

    def test_evidence_catalog_rejects_classified_failed_or_wrong_cap_outcomes(self) -> None:
        base_outcome = {
            "classification_status": "classified",
            "reason": "",
            "candidates": [
                {
                    "subject": "SSU_ref",
                    "taxonomy": "Bacteria;Firmicutes",
                    "taxonomy_source": "SILVA",
                }
            ],
            "taxonomy": "Bacteria",
            "taxonomy_source": "SILVA",
            "domain": "Bacteria",
            "assignment_method": "updated_reference_cluster",
            "calibration_rank_cap": 0,
            "propagation_rank_cap": 0,
            "evidence_id": "IMGEV_" + "a" * 64,
        }
        key = "16S|SILVA|Bacteria"
        cases = (
            (
                "failed stratum",
                {
                    "strata": {
                        key: {
                            "status": "failed",
                            "rank_cap": None,
                            "reason": "domain_precision_below_threshold",
                        }
                    }
                },
                base_outcome,
                "failed stratum",
            ),
            (
                "wrong cap",
                {
                    "strata": {
                        key: {"status": "calibrated", "rank_cap": 1, "reason": ""}
                    }
                },
                base_outcome,
                "wrong calibration rank cap",
            ),
        )
        with tempfile.TemporaryDirectory() as tmp:
            marker = Path(tmp) / "16S"
            marker.mkdir()
            assignments = marker / "assignments.tsv"
            assignments.write_text("source_identifier\n", encoding="utf-8")
            outcomes = marker / "outcomes.jsonl"
            for label, calibration, outcome, message in cases:
                with self.subTest(label=label):
                    outcomes.write_text(json.dumps(outcome) + "\n", encoding="utf-8")
                    with self.assertRaisesRegex(Exception, message):
                        profiles.create_img_evidence_catalog(
                            [assignments], marker / "evidence.jsonl", calibration
                        )

    def test_evidence_catalog_rejects_forged_failed_stratum_evidence(self) -> None:
        key = "16S|PR2|Eukaryota"
        outcome = {
            "classification_status": "unclassified",
            "reason": "calibration_stratum_failed",
            "candidates": [
                {
                    "subject": "SSU_ref",
                    "taxonomy": "Eukaryota;TSAR",
                    "taxonomy_source": "PR2",
                }
            ],
            "failed_calibration_strata": [
                {"key": key, "reason": "forged_reason"}
            ],
            "evidence_id": "IMGEV_" + "a" * 64,
        }
        calibration = {
            "strata": {
                key: {
                    "status": "failed",
                    "rank_cap": None,
                    "reason": "insufficient_domain_calls",
                }
            }
        }
        with tempfile.TemporaryDirectory() as tmp:
            marker = Path(tmp) / "16S"
            marker.mkdir()
            assignments = marker / "assignments.tsv"
            assignments.write_text("source_identifier\n", encoding="utf-8")
            (marker / "outcomes.jsonl").write_text(
                json.dumps(outcome) + "\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(Exception, "inconsistent evidence"):
                profiles.create_img_evidence_catalog(
                    [assignments], marker / "evidence.jsonl", calibration
                )
            outcome["failed_calibration_strata"][0]["reason"] = (
                "insufficient_domain_calls"
            )
            (marker / "outcomes.jsonl").write_text(
                json.dumps(outcome) + "\n", encoding="utf-8"
            )
            mappings, metadata = profiles.create_img_evidence_catalog(
                [assignments], marker / "evidence.jsonl", calibration
            )
            self.assertEqual(len(mappings[assignments.resolve()]), 1)
            self.assertEqual(metadata["records"], 1)

    def test_source_path_resolution_uses_source_integrity_validator(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source_directory = Path(tmp)
            source = source_directory / "source.dat"
            source.write_bytes(b"pinned-source")
            catalog = {
                "sources": {
                    "test_source": {
                        "filename": source.name,
                        "bytes": source.stat().st_size,
                        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                    }
                }
            }

            resolved = profiles._source_paths(
                catalog, source_directory, ("test_source",)
            )

            self.assertEqual(resolved, {"test_source": source})

    def test_assignment_reader_requires_builder_contract_and_streams_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "assignments.tsv"
            fields = [
                "source_identifier",
                "taxonomy",
                "taxonomy_source",
                "assignment_method",
                "evidence",
                "compartment",
            ]
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow(
                    {
                        "source_identifier": "IMG_1",
                        "taxonomy": "Unclassified",
                        "taxonomy_source": "SILVA+PR2",
                        "assignment_method": "updated_reference_unclassified",
                        "evidence": "no hit",
                        "compartment": "",
                    }
                )
            rows = list(profiles.read_assignment_rows([path]))
        self.assertEqual(rows[0]["source_identifier"], "IMG_1")

    def test_provenance_marks_cluster_inputs_as_build_only(self) -> None:
        catalog = {
            "taxonomy_policy_version": "1.0.0",
            "sources": {
                name: {
                    "filename": name,
                    "version": (
                        "138.2"
                        if name.startswith("silva_")
                        else "5.1.1"
                        if name == "pr2_ssu_fasta"
                        else "2025"
                    ),
                }
                for name in profiles.IMG_SOURCE_NAMES
            },
        }
        details = profiles._provenance_details(
            catalog,
            profiles.IMG_SOURCE_NAMES,
            "img",
            [{"action": "dropped_invalid_pair"}],
        )
        self.assertFalse(details["cluster_tables_packaged"])
        self.assertFalse(details["embedded_reference_records_retained"])
        self.assertEqual(
            details["metadata_corrections"],
            [{"action": "dropped_invalid_pair"}],
        )
        self.assertEqual(details["taxonomy_policy"]["Bacteria"], "SILVA 138.2")
        self.assertEqual(details["taxonomy_policy"]["Eukaryota"], "PR2 5.1.1")
        integrity = details["build_integrity"]
        self.assertEqual(integrity["command"]["profile_argument"], "img")
        self.assertEqual(len(integrity["repository"]["source_tree_sha256"]), 64)
        self.assertEqual(len(integrity["environment"]["pixi_lock_sha256"]), 64)
        self.assertIn("scripts/database_sources.py", profiles.SOURCE_TREE_FILES)
        self.assertIn("scripts/database_release_io.py", profiles.SOURCE_TREE_FILES)
        self.assertIn("scripts/img_classification_data.py", profiles.SOURCE_TREE_FILES)
        self.assertIn("scripts/taxonomy_utils.py", profiles.SOURCE_TREE_FILES)

    def test_img_evidence_catalog_is_content_addressed_and_privacy_safe(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            marker = root / "16S"
            marker.mkdir()
            old_id = "IMGEV_" + "a" * 64
            outcome = {
                "cluster_id": "private-cluster",
                "centroid": "private-centroid",
                "img_member_count": 17,
                "classification_status": "classified",
                "reason": "",
                "returned_hit_count": 1,
                "eligible_hit_count": 1,
                "coverage_filtered_hit_count": 0,
                "candidate_count": 1,
                "candidates": [
                    {
                        "subject": "SSU_public_reference",
                        "percent_identity": "100",
                        "taxonomy": "Bacteria;Firmicutes",
                        "taxonomy_source": "SILVA",
                    }
                ],
                "taxonomy": "Bacteria",
                "taxonomy_source": "SILVA",
                "domain": "Bacteria",
                "assignment_method": "updated_reference_cluster",
                "truncated": False,
                "calibration_rank_cap": 0,
                "propagation_rank_cap": 0,
                "evidence_id": old_id,
            }
            (marker / "outcomes.jsonl").write_text(
                json.dumps(outcome, sort_keys=True) + "\n", encoding="utf-8"
            )
            assignments = marker / "assignments.tsv"
            fields = [
                "source_identifier",
                "taxonomy",
                "taxonomy_source",
                "assignment_method",
                "evidence",
                "compartment",
                "evidence_id",
            ]
            with assignments.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow(
                    {
                        "source_identifier": "IMG_1",
                        "taxonomy": "Bacteria",
                        "taxonomy_source": "SILVA",
                        "assignment_method": "updated_reference_cluster",
                        "evidence": f"centroid_blast_lca evidence_id={old_id}",
                        "compartment": "",
                        "evidence_id": old_id,
                    }
                )
            catalog = root / "catalog.jsonl"
            mappings, metadata = profiles.create_img_evidence_catalog(
                [assignments],
                catalog,
                {
                    "strata": {
                        "16S|SILVA|Bacteria": {
                            "status": "calibrated",
                            "rank_cap": 0,
                            "reason": "",
                        }
                    }
                },
            )
            text = catalog.read_text(encoding="utf-8")
            records = [json.loads(line) for line in text.splitlines()]
            evidence = records[1]
            payload_bytes = profiles._canonical_json(evidence["decision"]).encode("utf-8")
            expected_hash = hashlib.sha256(payload_bytes).hexdigest()
            self.assertEqual(evidence["content_sha256"], expected_hash)
            self.assertEqual(evidence["evidence_id"], f"IMGEV_{expected_hash}")
            self.assertNotIn("private-cluster", text)
            self.assertNotIn("private-centroid", text)
            self.assertNotIn("IMG_1", text)
            self.assertNotIn("cluster_id", evidence["decision"])
            self.assertEqual(metadata["sha256"], hashlib.sha256(catalog.read_bytes()).hexdigest())
            rows = list(profiles.read_assignment_rows([assignments], mappings))
            self.assertIn(evidence["evidence_id"], rows[0]["evidence"])
            self.assertEqual(rows[0]["evidence_id"], evidence["evidence_id"])
            assignment_text = assignments.read_text(encoding="utf-8")
            assignments.write_text(
                assignment_text.replace("\tBacteria\t", "\tArchaea\t", 1),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(Exception, "does not match evidence"):
                list(profiles.read_assignment_rows([assignments], mappings))

    def test_img_evidence_rejects_nested_private_identifiers(self) -> None:
        with self.assertRaisesRegex(Exception, "forbidden"):
            profiles._safe_evidence_payload(
                "18S",
                {
                    "classification_status": "classified",
                    "candidates": [{"source_identifier": "IMG_123"}],
                },
            )

    def test_load_source_catalog_rejects_wrong_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sources.json"
            path.write_text(json.dumps({"schema_version": 2, "sources": {}}))
            with self.assertRaisesRegex(Exception, "schema_version"):
                profiles.load_source_catalog(path)


if __name__ == "__main__":
    unittest.main()
