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

import img_search_provenance as provenance


class ImgSearchProvenanceTests(unittest.TestCase):
    @staticmethod
    def _fixture(root: Path) -> dict[str, Path]:
        profile = root / "profile"
        (profile / "tables").mkdir(parents=True)
        (profile / "blast").mkdir()
        paths = {
            "profile": profile,
            "clusters": root / "clusters.tsv",
            "source_fasta": root / "source.fna",
            "centroids_fasta": root / "output" / "centroids.fna",
            "blast_m8": root / "output" / "centroids.m8",
            "sidecar": root / "output" / "search_provenance.json",
        }
        paths["centroids_fasta"].parent.mkdir()
        (profile / provenance.MANIFEST_NAME).write_text(
            '{"version":"test"}\n', encoding="utf-8"
        )
        (profile / provenance.TAXONOMY_RELATIVE_PATH).write_bytes(b"PAR1taxonomy")
        paths["clusters"].write_text(
            "cluster_id\tcentroid\nC1\tIMG_1\n", encoding="utf-8"
        )
        paths["source_fasta"].write_text(">IMG_1\nACGT\n", encoding="ascii")
        paths["centroids_fasta"].write_text(">C1\nACGT\n", encoding="ascii")
        return paths

    @staticmethod
    def _run_successfully(paths: dict[str, Path]) -> dict[str, object]:
        def fake_blastn(command: list[str], check: bool) -> None:
            if not check or command[0] != "blastn":
                raise AssertionError("unexpected subprocess contract")
            Path(command[command.index("-out") + 1]).write_text(
                "C1\tref\t100\t4\t4\t4\t100\t8\n", encoding="ascii"
            )

        with mock.patch.object(
            provenance.subprocess, "run", side_effect=fake_blastn
        ):
            return provenance.run_search(
                "16S",
                paths["profile"],
                paths["clusters"],
                paths["source_fasta"],
                paths["centroids_fasta"],
                paths["blast_m8"],
                paths["sidecar"],
                8,
            )

    @staticmethod
    def _validate(paths: dict[str, Path], *, marker: str = "16S", threads: int = 8):
        return provenance.validate_search(
            marker,
            paths["profile"],
            paths["clusters"],
            paths["source_fasta"],
            paths["centroids_fasta"],
            paths["blast_m8"],
            paths["sidecar"],
            threads,
        )

    def test_success_writes_exact_durable_contract_and_restart_validates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            payload = self._run_successfully(paths)
            observed = json.loads(paths["sidecar"].read_text(encoding="utf-8"))

            self.assertEqual(observed, payload)
            self.assertEqual(observed["schema_version"], 1)
            self.assertEqual(observed["status"], "complete")
            self.assertEqual(observed["marker"], "16S")
            self.assertEqual(observed["blastn"]["max_target_seqs"], 501)
            self.assertEqual(observed["blastn"]["threads"], 8)
            self.assertIn("501", observed["blastn"]["argv"])
            self.assertIn("8", observed["blastn"]["argv"])
            self.assertEqual(self._validate(paths), observed)
            self.assertEqual(self._validate(paths), observed)
            self.assertFalse(
                list(paths["sidecar"].parent.glob(".search_provenance.json.*.tmp"))
            )

    def test_failed_or_empty_search_never_leaves_completion_sidecar(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            paths["sidecar"].write_text('{"status":"stale"}\n')
            paths["blast_m8"].write_text("stale\n")
            failure = subprocess.CalledProcessError(9, ["blastn"])
            with (
                mock.patch.object(provenance.subprocess, "run", side_effect=failure),
                self.assertRaisesRegex(RuntimeError, "exit code 9"),
            ):
                provenance.run_search(
                    "16S",
                    paths["profile"],
                    paths["clusters"],
                    paths["source_fasta"],
                    paths["centroids_fasta"],
                    paths["blast_m8"],
                    paths["sidecar"],
                    8,
                )
            self.assertFalse(paths["sidecar"].exists())
            self.assertFalse(paths["blast_m8"].exists())

        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            with (
                mock.patch.object(provenance.subprocess, "run", return_value=None),
                self.assertRaisesRegex(RuntimeError, "blast_m8"),
            ):
                provenance.run_search(
                    "16S",
                    paths["profile"],
                    paths["clusters"],
                    paths["source_fasta"],
                    paths["centroids_fasta"],
                    paths["blast_m8"],
                    paths["sidecar"],
                    8,
                )
            self.assertFalse(paths["sidecar"].exists())

    def test_validator_rejects_every_bound_file_after_tampering(self) -> None:
        mutations = {
            "curated_manifest": lambda paths: paths["profile"]
            / provenance.MANIFEST_NAME,
            "preferred_taxonomy": lambda paths: paths["profile"]
            / provenance.TAXONOMY_RELATIVE_PATH,
            "clusters": lambda paths: paths["clusters"],
            "source_fasta": lambda paths: paths["source_fasta"],
            "centroids_fasta": lambda paths: paths["centroids_fasta"],
            "blast_m8": lambda paths: paths["blast_m8"],
        }
        for label, target in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                paths = self._fixture(Path(temporary))
                self._run_successfully(paths)
                with target(paths).open("ab") as handle:
                    handle.write(b"tampered")
                with self.assertRaisesRegex(RuntimeError, label):
                    self._validate(paths)

    def test_validator_rejects_marker_thread_command_and_schema_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            self._run_successfully(paths)
            with self.assertRaisesRegex(RuntimeError, "marker mismatch"):
                self._validate(paths, marker="18S")
            with self.assertRaisesRegex(RuntimeError, "command contract mismatch"):
                self._validate(paths, threads=7)

            sidecar = json.loads(paths["sidecar"].read_text(encoding="utf-8"))
            sidecar["blastn"]["argv"][10] = "500"
            paths["sidecar"].write_text(json.dumps(sidecar), encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "command contract mismatch"):
                self._validate(paths)

            sidecar["schema_version"] = 2
            paths["sidecar"].write_text(json.dumps(sidecar), encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "not schema-1 complete"):
                self._validate(paths)

    def test_portable_provenance_contains_hashes_but_no_local_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._fixture(Path(temporary))
            payload = self._run_successfully(paths)
            sidecar_sha256 = hashlib.sha256(paths["sidecar"].read_bytes()).hexdigest()
            portable = provenance.portable_provenance(payload, sidecar_sha256)
            rendered = json.dumps(portable, sort_keys=True)

        self.assertEqual(portable["sidecar_sha256"], sidecar_sha256)
        self.assertEqual(
            portable["blastn"]["argv"][4], "@centroids_fasta"
        )
        self.assertEqual(
            portable["blastn"]["argv"][6], "@curated_blast_database"
        )
        self.assertNotIn("path", rendered)
        self.assertNotIn(str(Path(temporary)), rendered)
        self.assertNotIn("/clusterfs", rendered)
        self.assertNotIn(getpass.getuser(), rendered)


if __name__ == "__main__":
    unittest.main()
