import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import extract_img_cluster_centroids as extractor


class ExtractCentroidTests(unittest.TestCase):
    def test_extracts_exact_headers_and_uses_cluster_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            clusters = root / "clusters.tsv"
            fasta = root / "source.fna"
            output = root / "centroids.fna"
            clusters.write_text(
                "cluster_id\tcentroid\tsequences\n"
                "C2\tREF_two;Taxonomy\t[]\n"
                "C1\tIMG_1\t[]\n"
            )
            fasta.write_text(
                ">unused\nAAAA\n>REF_two;Taxonomy\naugu\n>IMG_1\nCCCC\n"
            )
            report = extractor.extract_centroids(clusters, fasta, output)
            content = output.read_text()
        self.assertEqual(report, {"clusters": 2, "centroids_written": 2})
        self.assertEqual(content, ">C1\nCCCC\n>C2\nATGT\n")

    def test_missing_centroid_is_blocking(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            clusters = root / "clusters.tsv"
            fasta = root / "source.fna"
            clusters.write_text("cluster_id\tcentroid\nC1\tmissing\n")
            fasta.write_text(">other\nAAAA\n")
            with self.assertRaisesRegex(ValueError, "missing FASTA sequence"):
                extractor.extract_centroids(clusters, fasta, root / "out.fna")


if __name__ == "__main__":
    unittest.main()
