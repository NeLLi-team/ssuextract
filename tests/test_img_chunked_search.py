import getpass
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

import img_chunked_search as chunked
import img_search_provenance as provenance
import build_database_profiles as profiles


class ImgChunkedSearchTests(unittest.TestCase):
    TOOL_RECORD = {
        "path": "/fake/blastn",
        "sha256": "a" * 64,
        "version": "blastn: test",
    }

    @staticmethod
    def _fixture(root: Path) -> dict[str, Path]:
        profile = root / "profile"
        (profile / "tables").mkdir(parents=True)
        (profile / "blast").mkdir()
        artifacts = []
        for marker in ("16S", "18S"):
            for suffix in ("nhr", "nin", "nsq"):
                path = profile / "blast" / f"{marker}.{suffix}"
                path.write_bytes(f"blast-{marker}-{suffix}".encode("ascii"))
                artifacts.append(
                    {
                        "path": f"blast/{marker}.{suffix}",
                        "bytes": path.stat().st_size,
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    }
                )
        (profile / provenance.MANIFEST_NAME).write_text(
            json.dumps({"version": "test", "artifacts": artifacts}) + "\n",
            encoding="utf-8",
        )
        (profile / provenance.TAXONOMY_RELATIVE_PATH).write_bytes(b"PAR1taxonomy")
        clusters = root / "clusters.tsv"
        source_fasta = root / "source.fna"
        centroids = root / "output" / "centroids.fna"
        centroids.parent.mkdir()
        clusters.write_text(
            "cluster_id\tcentroid\n"
            + "".join(f"C{index}\tIMG_{index}\n" for index in range(5)),
            encoding="utf-8",
        )
        source_fasta.write_text(
            "".join(f">IMG_{index}\nACGT{index}\n" for index in range(5)),
            encoding="ascii",
        )
        centroids.write_text(
            "".join(f">C{index}\nACGT{index}\n" for index in range(5)),
            encoding="ascii",
        )
        return {
            "profile": profile,
            "clusters": clusters,
            "source_fasta": source_fasta,
            "centroids": centroids,
            "chunk_directory": centroids.parent / "chunks",
            "blast_m8": centroids.parent / "centroids.m8",
            "sidecar": centroids.parent / "search_provenance.json",
        }

    @staticmethod
    def _fake_blastn(command: list[str], check: bool) -> None:
        if not check or command[0] != "blastn":
            raise AssertionError("unexpected subprocess contract")
        query = Path(command[command.index("-query") + 1])
        output = Path(command[command.index("-out") + 1])
        identifiers = [
            line[1:].split()[0]
            for line in query.read_text(encoding="ascii").splitlines()
            if line.startswith(">")
        ]
        output.write_text(
            "".join(
                f"{identifier}\tref\t100\t5\t5\t5\t100\t10\n"
                for identifier in identifiers
            ),
            encoding="ascii",
        )

    def _complete(
        self, paths: dict[str, Path], marker: str = "16S"
    ) -> dict[str, object]:
        plan = chunked.prepare_chunks(
            marker, paths["centroids"], paths["chunk_directory"], 2
        )
        self.assertEqual([item["query_count"] for item in plan["chunks"]], [3, 2])
        with (
            mock.patch.object(chunked, "_blastn_tool", return_value=self.TOOL_RECORD),
            mock.patch.object(chunked.subprocess, "run", side_effect=self._fake_blastn),
        ):
            for index in range(2):
                chunked.run_chunk(
                    marker,
                    paths["profile"],
                    paths["centroids"],
                    paths["chunk_directory"],
                    index,
                    8,
                )
            return chunked.finalize_search(
                marker,
                paths["profile"],
                paths["clusters"],
                paths["source_fasta"],
                paths["centroids"],
                paths["blast_m8"],
                paths["sidecar"],
                paths["chunk_directory"],
                8,
            )

    def test_chunks_preserve_query_order_and_validate_complete_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            payload = self._complete(paths)
            with mock.patch.object(
                chunked, "_blastn_tool", return_value=self.TOOL_RECORD
            ):
                observed = provenance.validate_search(
                    "16S",
                    paths["profile"],
                    paths["clusters"],
                    paths["source_fasta"],
                    paths["centroids"],
                    paths["blast_m8"],
                    paths["sidecar"],
                    8,
                )
                with self.assertRaisesRegex(RuntimeError, "receipt contract mismatch"):
                    provenance.validate_search(
                        "16S",
                        paths["profile"],
                        paths["clusters"],
                        paths["source_fasta"],
                        paths["centroids"],
                        paths["blast_m8"],
                        paths["sidecar"],
                        7,
                    )
            query_order = [
                line.split("\t", 1)[0]
                for line in paths["blast_m8"].read_text(encoding="ascii").splitlines()
            ]
            self.assertEqual(observed, payload)
            self.assertEqual(payload["schema_version"], 2)
            self.assertEqual(payload["blastn"]["mode"], "query_chunks")
            self.assertEqual(payload["blastn"]["chunk_count"], 2)
            self.assertEqual(query_order, ["C0", "C1", "C2", "C3", "C4"])

            sidecar_sha256 = hashlib.sha256(paths["sidecar"].read_bytes()).hexdigest()
            portable = provenance.portable_provenance(payload, sidecar_sha256)
            rendered = json.dumps(portable, sort_keys=True)
            self.assertEqual(portable["schema_version"], 2)
            self.assertEqual(portable["blastn"]["chunk_count"], 2)
            self.assertNotIn("path", rendered)
            self.assertNotIn(str(Path(temporary)), rendered)
            self.assertNotIn("/clusterfs", rendered)
            self.assertNotIn(getpass.getuser(), rendered)

    def test_schema_2_sidecar_is_consumed_by_profile_builder(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            payload = self._complete(paths)
            assignments = paths["sidecar"].parent / "assignments.tsv"
            assignments.write_text("source_identifier\n", encoding="utf-8")
            expected = provenance.portable_provenance(
                payload, hashlib.sha256(paths["sidecar"].read_bytes()).hexdigest()
            )
            with mock.patch.object(
                chunked, "_blastn_tool", return_value=self.TOOL_RECORD
            ):
                observed = profiles.validate_classification_search(
                    assignments,
                    "16S",
                    paths["profile"],
                    paths["clusters"],
                    paths["source_fasta"],
                )
            self.assertEqual(observed, expected)
            self.assertTrue(paths["chunk_directory"].is_dir())

    def test_validation_hashes_database_and_tool_once_for_all_chunks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            self._complete(paths)
            with (
                mock.patch.object(
                    chunked,
                    "_database_binding",
                    wraps=chunked._database_binding,
                ) as database_binding,
                mock.patch.object(
                    chunked, "_blastn_tool", return_value=self.TOOL_RECORD
                ) as blastn_tool,
            ):
                provenance.validate_search(
                    "16S",
                    paths["profile"],
                    paths["clusters"],
                    paths["source_fasta"],
                    paths["centroids"],
                    paths["blast_m8"],
                    paths["sidecar"],
                    8,
                )
            self.assertEqual(database_binding.call_count, 1)
            self.assertEqual(blastn_tool.call_count, 1)

    def test_18s_uses_its_own_database_binding(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            payload = self._complete(paths, marker="18S")
            self.assertEqual(payload["marker"], "18S")
            for receipt in payload["blastn"]["chunks"]:
                self.assertTrue(
                    all(
                        artifact["path"].startswith("blast/18S.")
                        for artifact in receipt["database"]["artifacts"]
                    )
                )

    def test_existing_plan_rejects_count_or_content_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            chunked.prepare_chunks(
                "16S", paths["centroids"], paths["chunk_directory"], 2
            )
            with self.assertRaisesRegex(RuntimeError, "chunk_count differs"):
                chunked.prepare_chunks(
                    "16S", paths["centroids"], paths["chunk_directory"], 3
                )
            plan_path = paths["chunk_directory"] / chunked.PLAN_NAME
            plan = json.loads(plan_path.read_text(encoding="utf-8"))
            plan["chunks"][0]["first_query_id"] = "tampered"
            plan_path.write_text(json.dumps(plan), encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "query content mismatch"):
                chunked.validate_plan(
                    "16S", paths["centroids"], paths["chunk_directory"]
                )

    def test_valid_receipt_skips_and_force_reruns_chunk(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            chunked.prepare_chunks(
                "16S", paths["centroids"], paths["chunk_directory"], 2
            )
            with (
                mock.patch.object(
                    chunked, "_blastn_tool", return_value=self.TOOL_RECORD
                ),
                mock.patch.object(
                    chunked.subprocess, "run", side_effect=self._fake_blastn
                ) as runner,
            ):
                chunked.run_chunk(
                    "16S",
                    paths["profile"],
                    paths["centroids"],
                    paths["chunk_directory"],
                    0,
                    8,
                )
                chunked.run_chunk(
                    "16S",
                    paths["profile"],
                    paths["centroids"],
                    paths["chunk_directory"],
                    0,
                    8,
                )
                self.assertEqual(runner.call_count, 1)
                chunked.run_chunk(
                    "16S",
                    paths["profile"],
                    paths["centroids"],
                    paths["chunk_directory"],
                    0,
                    8,
                    force=True,
                )
                self.assertEqual(runner.call_count, 2)

    def test_failed_or_tampered_chunk_cannot_finalize(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            chunked.prepare_chunks(
                "16S", paths["centroids"], paths["chunk_directory"], 2
            )
            failure = subprocess.CalledProcessError(7, ["blastn"])
            with (
                mock.patch.object(
                    chunked, "_blastn_tool", return_value=self.TOOL_RECORD
                ),
                mock.patch.object(chunked.subprocess, "run", side_effect=failure),
                self.assertRaisesRegex(RuntimeError, "chunk 0 failed"),
            ):
                chunked.run_chunk(
                    "16S",
                    paths["profile"],
                    paths["centroids"],
                    paths["chunk_directory"],
                    0,
                    8,
                )
            self.assertFalse(
                (paths["chunk_directory"] / "receipts" / "000.json").exists()
            )

        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            chunked.prepare_chunks(
                "16S", paths["centroids"], paths["chunk_directory"], 2
            )
            with (
                mock.patch.object(
                    chunked, "_blastn_tool", return_value=self.TOOL_RECORD
                ),
                mock.patch.object(
                    chunked.subprocess, "run", side_effect=self._fake_blastn
                ),
            ):
                chunked.run_chunk(
                    "16S",
                    paths["profile"],
                    paths["centroids"],
                    paths["chunk_directory"],
                    0,
                    8,
                )
                with self.assertRaisesRegex(RuntimeError, "invalid chunk receipt"):
                    chunked.finalize_search(
                        "16S",
                        paths["profile"],
                        paths["clusters"],
                        paths["source_fasta"],
                        paths["centroids"],
                        paths["blast_m8"],
                        paths["sidecar"],
                        paths["chunk_directory"],
                        8,
                    )

        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            self._complete(paths)
            result = paths["chunk_directory"] / "results" / "000.m8"
            result.write_text("tampered\n", encoding="ascii")
            with (
                mock.patch.object(
                    chunked, "_blastn_tool", return_value=self.TOOL_RECORD
                ),
                self.assertRaisesRegex(RuntimeError, "receipt contract mismatch"),
            ):
                provenance.validate_search(
                    "16S",
                    paths["profile"],
                    paths["clusters"],
                    paths["source_fasta"],
                    paths["centroids"],
                    paths["blast_m8"],
                    paths["sidecar"],
                    8,
                )

        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            self._complete(paths)
            (paths["profile"] / "blast" / "16S.nsq").write_bytes(b"tampered")
            with (
                mock.patch.object(
                    chunked, "_blastn_tool", return_value=self.TOOL_RECORD
                ),
                self.assertRaisesRegex(RuntimeError, "artifact binding mismatch"),
            ):
                provenance.validate_search(
                    "16S",
                    paths["profile"],
                    paths["clusters"],
                    paths["source_fasta"],
                    paths["centroids"],
                    paths["blast_m8"],
                    paths["sidecar"],
                    8,
                )


if __name__ == "__main__":
    unittest.main()
