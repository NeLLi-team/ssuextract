import ast
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import build_database_profiles as profiles  # noqa: E402  # Standalone scripts are loaded from the repository path.
import build_database_release as builder  # noqa: E402  # Standalone scripts are loaded from the repository path.
import classify_img_clusters as classifier  # noqa: E402  # Load repository scripts.
import database_contracts as contracts  # noqa: E402  # Load repository scripts.
import database_release_io as release_io  # noqa: E402  # Load repository scripts.
import img_chunked_search as chunked  # noqa: E402  # Standalone scripts are loaded from the repository path.
import img_classification_data as classification_data  # noqa: E402  # Standalone scripts are loaded from the repository path.


def calibration_document(
    rank_caps: dict[str, int], strata: dict[str, dict[str, object]]
) -> dict[str, object]:
    return {
        "schema_version": classification_data.CALIBRATION_SCHEMA_VERSION,
        "status": "complete" if rank_caps else "failed",
        "calibration_use_case": "img_centroid",
        "stratum_basis": "predicted_candidate_source_and_lca_domain",
        "classification_policy": classification_data.classification_policy(),
        "rank_caps": rank_caps,
        "strata": strata,
    }


class ProfileBuildDriverTests(unittest.TestCase):
    def test_source_fingerprint_covers_driver_import_closure(self) -> None:
        pending = [REPO / "scripts/build_database_profiles.py"]
        closure = set()
        while pending:
            path = pending.pop()
            if path in closure:
                continue
            closure.add(path)
            for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
                if isinstance(node, ast.Import):
                    modules = [alias.name for alias in node.names]
                elif isinstance(node, ast.ImportFrom) and node.module:
                    modules = [node.module]
                else:
                    continue
                pending.extend(
                    candidate
                    for module in modules
                    if (candidate := REPO / "scripts" / f"{module}.py").is_file()
                )
        self.assertLessEqual(
            {str(path.relative_to(REPO)) for path in closure},
            set(profiles.SOURCE_TREE_FILES),
        )

    def test_git_provenance_rejects_missing_unknown_or_dirty_state(self) -> None:
        clean_commit = subprocess.CompletedProcess([], 0, "a" * 40 + "\n", "")
        clean_status = subprocess.CompletedProcess([], 0, "", "")
        cases = (
            [FileNotFoundError("git")],
            [subprocess.CompletedProcess([], 0, "\n", ""), clean_status],
            [clean_commit, subprocess.CalledProcessError(1, ["git", "status"])],
            *(
                [clean_commit, subprocess.CompletedProcess([], 0, status, "")]
                for status in (
                    " M tracked.txt\n",
                    "M  tracked.txt\n",
                    " D tracked.txt\n",
                )
            ),
        )
        for responses in cases:
            with (
                self.subTest(responses=responses),
                mock.patch.object(profiles.subprocess, "run", side_effect=responses),
                mock.patch.object(profiles, "_tool_version", return_value="test"),
                self.assertRaisesRegex(Exception, "git|Git|dirty"),
            ):
                profiles._source_tree_provenance()
        with (
            mock.patch.object(
                profiles.subprocess,
                "run",
                side_effect=[clean_commit, clean_status],
            ) as run,
            mock.patch.object(profiles, "_tool_version", return_value="test"),
        ):
            provenance = profiles._source_tree_provenance()
        self.assertEqual(provenance["repository"]["commit"], "a" * 40)
        self.assertEqual(
            run.call_args_list[1].args[0],
            ["git", "status", "--porcelain", "--untracked-files=no"],
        )
        self.assertFalse(provenance["repository"]["tracked_worktree_dirty"])

    def test_build_checks_destinations_and_version_before_source_loading(self) -> None:
        for profile in ("curated", "img"):
            for invalid in ("profile", "archive", "version"):
                with (
                    self.subTest(profile=profile, invalid=invalid),
                    tempfile.TemporaryDirectory() as temporary,
                ):
                    root = Path(temporary)
                    args = profiles._parser().parse_args(
                        [
                            profile,
                            "--source-directory",
                            str(root / "sources"),
                            "--output-root",
                            str(root / "profiles"),
                            "--archive-directory",
                            str(root / "archives"),
                        ]
                    )
                    if invalid == "profile":
                        (args.output_root / profile).mkdir(parents=True)
                    elif invalid == "archive":
                        args.archive_directory.mkdir()
                        (
                            args.archive_directory
                            / f"ssuextract-db-{profile}-v{args.version}.tar.zst"
                        ).touch()
                    else:
                        args.version = "bad/version"
                    build = (
                        profiles.build_curated
                        if profile == "curated"
                        else profiles.build_img
                    )
                    with (
                        mock.patch.object(profiles, "load_source_catalog") as load,
                        self.assertRaisesRegex(
                            Exception, "already exists|safe, non-empty"
                        ),
                    ):
                        build(args)
                    load.assert_not_called()

    def test_location_corrections_only_include_released_taxa(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "metadata.tsv"
            path.write_text(
                "taxon_oid\tlatitude\tlongitude\n1\t999\t10\n2\t999\t10\n",
                encoding="utf-8",
            )
            locations, corrections = profiles.read_img_locations(path, {"1"})
        self.assertEqual([location.taxon_oid for location in locations], ["1"])
        self.assertEqual([correction["taxon_oid"] for correction in corrections], ["1"])

    def test_outcome_evidence_requires_classifier_identifier(self) -> None:
        with self.assertRaisesRegex(Exception, "evidence_id"):
            profiles._outcome_evidence_id({"classification_status": "unclassified"})

    def test_unclassified_classifier_evidence_roundtrips_through_builder(self) -> None:
        clusters = classification_data.parse_clusters(
            [
                "cluster_id\tcentroid\tsequences\n",
                "C1\tIMG_1\t['IMG_1']\n",
                "C2\tREF_A\t['REF_A']\n",
            ]
        )
        assignments, outcomes, _ = classifier.classify_clusters(
            clusters, (), {}, marker="16S"
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            marker = root / "16S"
            marker.mkdir()
            path = marker / "assignments.tsv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=classifier.ASSIGNMENT_FIELDS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(assignments)
            (marker / "outcomes.jsonl").write_text(
                "".join(json.dumps(outcome) + "\n" for outcome in outcomes),
                encoding="utf-8",
            )
            mappings, metadata = profiles.create_img_evidence_catalog(
                [path], root / "evidence.jsonl", {}
            )
            rows = list(profiles.read_assignment_rows([path], mappings))
        self.assertEqual(metadata["records"], 2)
        self.assertEqual(rows[0]["taxonomy"], "Unclassified")
        self.assertEqual(rows[0]["marker"], "16S")

    def test_curated_source_rebuild_matches_classification_reference(self) -> None:
        reference = builder.build_deduplicated_model(
            [
                contracts.PreparedSourceRecord(
                    "SILVA", "test", "ref", "ref", "ACGT", "16S", ("Bacteria",), "SILVA"
                )
            ]
        )
        other = builder.build_deduplicated_model(
            [
                contracts.PreparedSourceRecord(
                    "SILVA", "test", "ref", "ref", "ACGA", "16S", ("Bacteria",), "SILVA"
                )
            ]
        )
        with tempfile.TemporaryDirectory() as temporary:
            profile = Path(temporary)
            (profile / "tables").mkdir()
            release_io.write_sequences_parquet(
                profile / "tables/sequences.parquet", reference.sequences
            )
            profiles.validate_curated_sequences(reference, profile)
            with self.assertRaisesRegex(Exception, "curated.*sequences"):
                profiles.validate_curated_sequences(other, profile)

    def test_curated_source_filter_counts_reach_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = {
                "silva_ssu_nr99_fasta": root / "silva.fasta",
                "silva_ssu_taxonomy": root / "silva.tsv",
                "pr2_ssu_fasta": root / "pr2.fasta",
            }
            paths["silva_ssu_nr99_fasta"].write_text(
                ">S1.1.4 Bacteria;Firmicutes;\nACGT\n>E1.1.4 Eukaryota;TSAR;\nACGA\n",
                encoding="utf-8",
            )
            paths["silva_ssu_taxonomy"].write_text(
                "Bacteria;\t2\tdomain\n", encoding="utf-8"
            )
            paths["pr2_ssu_fasta"].write_text(
                ">P1.1.4_U|18S_rRNA|nucleus|9|Eukaryota|TSAR|D|S|C|O|F|G|Species\nACGG\n",
                encoding="utf-8",
            )
            catalog = {"sources": {name: {"version": "test"} for name in paths}}
            skipped = {}
            list(profiles.iter_curated_records(paths, catalog, skipped))
        self.assertEqual(skipped, {"silva_ssu_nr99_fasta": 1, "pr2_ssu_fasta": 0})

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
        self.assertEqual(args.version, "1.0.2")

    def test_search_and_classification_qc_bind_the_same_portable_provenance(self) -> None:
        calibration = {"sha256": "a" * 64, "schema_version": 3}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "curated"
            (profile / "tables").mkdir(parents=True)
            (profile / "blast").mkdir()
            artifacts = []
            for suffix in ("nhr", "nin", "nsq"):
                artifact = profile / "blast" / f"16S.{suffix}"
                artifact.write_bytes(b"blast fixture")
                artifacts.append(
                    {
                        "path": f"blast/16S.{suffix}",
                        "bytes": artifact.stat().st_size,
                        "sha256": profiles._sha256(artifact),
                    }
                )
            (profile / "manifest.json").write_text(
                json.dumps({"artifacts": artifacts}), encoding="utf-8"
            )
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
            chunk_directory = marker / "chunks"
            chunked.prepare_chunks("16S", centroids, chunk_directory, 1)
            tool = {"path": "/fake/blastn", "version": "test", "sha256": "a" * 64}

            def fake_blastn(command: list[str], *, check: bool) -> None:
                self.assertTrue(check)
                Path(command[command.index("-out") + 1]).write_text(
                    "0\tref\t100\t4\t4\t4\t100\t8\n", encoding="ascii"
                )

            with (
                mock.patch.object(chunked, "_blastn_tool", return_value=tool),
                mock.patch.object(chunked.subprocess, "run", side_effect=fake_blastn),
            ):
                chunked.run_chunk("16S", profile, centroids, chunk_directory, 0, 8)
                chunked.finalize_search(
                    "16S",
                    profile,
                    clusters,
                    source_fasta,
                    centroids,
                    blast,
                    sidecar,
                    chunk_directory,
                    8,
                )
            with mock.patch.object(chunked, "_blastn_tool", return_value=tool):
                portable = profiles.validate_classification_search(
                    assignments, "16S", profile, clusters, source_fasta
                )
            binding = {
                "schema_version": 1,
                "calibration": {"sha256": "a" * 64, "schema_version": 3},
                "marker": "16S",
                "propagation_rank_cap": 0,
                "policy": classification_data.classification_policy(),
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
            with (
                mock.patch.object(chunked, "_blastn_tool", return_value=tool),
                self.assertRaisesRegex(Exception, "ordered chunk concatenation"),
            ):
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
            **calibration_document({calibrated_key: 6}, strata),
            "curated_profile": {
                "version": "1.0.0",
                "manifest_sha256": "a" * 64,
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calibration.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            provenance = profiles.calibration_provenance(path)

        self.assertEqual(provenance["schema_version"], 3)
        self.assertEqual(provenance["status"], "complete")
        self.assertEqual(provenance["calibration_use_case"], "img_centroid")
        self.assertEqual(
            provenance["stratum_basis"],
            "predicted_candidate_source_and_lca_domain",
        )
        self.assertEqual(
            provenance["classification_policy"],
            classification_data.classification_policy(),
        )
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

    def test_build_provenance_rejects_schema_2_calibration(self) -> None:
        document = {
            **calibration_document(
                {},
                {
                    "16S|SILVA|Bacteria": {
                        "status": "failed",
                        "rank_cap": None,
                        "reason": "insufficient_domain_calls",
                    }
                },
            ),
            "schema_version": 2,
            "curated_profile": {
                "version": "1.0.2",
                "manifest_sha256": "a" * 64,
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calibration.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaisesRegex(Exception, "schema_version must be 3"):
                profiles.calibration_provenance(path)

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
            "schema_version": 3,
        }
        search = {"schema_version": 2}
        expected = {
            "search_provenance": search,
            "schema_version": 1,
            "calibration": {"sha256": "a" * 64, "schema_version": 3},
            "marker": "16S",
            "propagation_rank_cap": 0,
            "policy": classification_data.classification_policy(),
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
                profiles.validate_classification_qc(
                    assignment, "16S", calibration, search
                )
            qc_path = marker / "qc.json"
            qc_path.write_text(
                json.dumps({"classification_binding": expected}), encoding="utf-8"
            )
            self.assertEqual(
                profiles.validate_classification_qc(
                    assignment, "16S", calibration, search
                ),
                qc_path,
            )
            assignment.write_text("substituted\n", encoding="utf-8")
            with self.assertRaisesRegex(Exception, "does not match"):
                profiles.validate_classification_qc(
                    assignment, "16S", calibration, search
                )
            assignment.write_text("source_identifier\n", encoding="utf-8")
            outcomes.write_text('{"substituted":true}\n', encoding="utf-8")
            with self.assertRaisesRegex(Exception, "does not match"):
                profiles.validate_classification_qc(
                    assignment, "16S", calibration, search
                )
            outcomes.write_text("{}\n", encoding="utf-8")
            mutations = (
                (
                    "wrong calibration",
                    {
                        **expected,
                        "calibration": {
                            "sha256": "b" * 64,
                            "schema_version": 3,
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
                            assignment, "16S", calibration, search
                        )

    def test_evidence_catalog_rejects_classified_failed_or_wrong_cap_outcomes(self) -> None:
        base_outcome = {
            "cluster_id": "cluster-1",
            "centroid": "centroid-1",
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
            "centroid_taxonomy": "Bacteria",
            "centroid_taxonomy_source": "SILVA",
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
            "cluster_id": "cluster-1",
            "centroid": "centroid-1",
            "classification_status": "unclassified",
            "reason": "calibration_stratum_failed",
            "candidates": [
                {
                    "subject": "SSU_ref",
                    "taxonomy": "Eukaryota;TSAR",
                    "taxonomy_source": "PR2",
                }
            ],
            "failed_calibration_strata": [{"key": key, "reason": "forged_reason"}],
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

    def test_identity_gated_outcome_matches_calibration_rule(self) -> None:
        key = "16S|SILVA|Bacteria"
        outcome = {
            "cluster_id": "cluster-1",
            "centroid": "centroid-1",
            "classification_status": "unclassified",
            "reason": "calibration_identity_below_domain",
            "candidates": [
                {
                    "subject": "SSU_ref",
                    "taxonomy": "Bacteria;Firmicutes",
                    "taxonomy_source": "SILVA",
                    "percent_identity": "94",
                }
            ],
            "required_calibration_strata": [key],
            "minimum_candidate_identity": "94",
        }
        calibration = {
            "strata": {key: {"status": "calibrated", "rank_cap": 0, "reason": ""}},
            "rank_rules": {
                key: [{"rank_index": 0, "min_candidate_identity": 95}]
            },
        }
        profiles.validate_outcome_calibration("16S", outcome, calibration)
        with self.assertRaisesRegex(Exception, "violates calibration identity rules"):
            profiles.validate_outcome_calibration(
                "16S",
                {
                    **outcome,
                    "classification_status": "classified",
                    "reason": "",
                    "calibration_rank_cap": 0,
                },
                calibration,
            )

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
            path = Path(tmp) / "16S" / "assignments.tsv"
            path.parent.mkdir()
            fields = [
                "marker",
                "evidence_id",
                "source_identifier",
                "taxonomy",
                "taxonomy_source",
                "assignment_method",
                "evidence",
                "compartment",
                "centroid",
                "centroid_name",
                "centroid_taxonomy",
                "centroid_taxonomy_source",
            ]
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow(
                    {
                        "evidence_id": "old-id",
                        "marker": "16S",
                        "source_identifier": "IMG_1",
                        "taxonomy": "Unclassified",
                        "taxonomy_source": "SILVA+PR2",
                        "assignment_method": "updated_reference_unclassified",
                        "evidence": "no hit old-id",
                        "compartment": "",
                        "centroid": "centroid-1",
                        "centroid_name": "centroid-1",
                        "centroid_taxonomy": "",
                        "centroid_taxonomy_source": "",
                    }
                )
            resolution = {
                "evidence_id": "new-id",
                "taxonomy": "Unclassified",
                "taxonomy_source": "SILVA+PR2",
                "assignment_method": "updated_reference_unclassified",
                "compartment": "",
                "centroid": "centroid-1",
                "centroid_name": "centroid-1",
                "centroid_taxonomy": "",
                "centroid_taxonomy_source": "",
            }
            rows = list(
                profiles.read_assignment_rows(
                    [path], {path.resolve(): {"old-id": resolution}}
                )
            )
            path.write_text(
                path.read_text(encoding="utf-8").replace("16S", "18S"), encoding="utf-8"
            )
            with self.assertRaisesRegex(Exception, "marker"):
                list(
                    profiles.read_assignment_rows(
                        [path], {path.resolve(): {"old-id": resolution}}
                    )
                )
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
        with (
            mock.patch.object(
                profiles.subprocess,
                "run",
                side_effect=[
                    subprocess.CompletedProcess([], 0, "a" * 40, ""),
                    subprocess.CompletedProcess([], 0, "", ""),
                ],
            ),
            mock.patch.object(profiles, "_tool_version", return_value="test"),
        ):
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
                "centroid_taxonomy": "Bacteria",
                "centroid_taxonomy_source": "SILVA",
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
                "marker",
                "source_identifier",
                "taxonomy",
                "taxonomy_source",
                "assignment_method",
                "evidence",
                "compartment",
                "centroid",
                "centroid_name",
                "centroid_taxonomy",
                "centroid_taxonomy_source",
                "evidence_id",
            ]
            with assignments.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow(
                    {
                        "marker": "16S",
                        "source_identifier": "IMG_1",
                        "taxonomy": "Bacteria",
                        "taxonomy_source": "SILVA",
                        "assignment_method": "updated_reference_cluster",
                        "evidence": f"centroid_blast_lca evidence_id={old_id}",
                        "compartment": "",
                        "centroid": "private-centroid",
                        "centroid_name": "private-centroid",
                        "centroid_taxonomy": "Bacteria",
                        "centroid_taxonomy_source": "SILVA",
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
