import csv
import statistics
import unittest
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data" / "example_performance.tsv"
PAGE = ROOT / "docs" / "reference" / "performance.md"


class PerformanceDocumentationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with DATA.open(encoding="utf-8", newline="") as handle:
            cls.rows = list(csv.DictReader(handle, delimiter="\t"))
        cls.page = PAGE.read_text(encoding="utf-8")

    def test_measurements_cover_every_profile_count_and_trial(self) -> None:
        conditions: dict[tuple[str, int], list[dict[str, str]]] = defaultdict(list)
        for row in self.rows:
            conditions[(row["profile"], int(row["query_sequences"]))].append(row)

        self.assertEqual(
            set(conditions),
            {
                (profile, count)
                for profile in ("curated", "img")
                for count in (1, 10, 100, 1000)
            },
        )
        for rows in conditions.values():
            rows.sort(key=lambda row: int(row["trial"]))
            self.assertEqual([int(row["trial"]) for row in rows], [0, 1, 2, 3])
            self.assertEqual(
                [row["warmup"] for row in rows],
                ["true", "false", "false", "false"],
            )
            self.assertTrue(all(float(row["elapsed_seconds"]) > 0 for row in rows))
            self.assertTrue(all(int(row["peak_rss_kb"]) > 0 for row in rows))
            self.assertTrue(all(int(row["hit_count"]) > 0 for row in rows))

    def test_page_table_matches_measured_medians(self) -> None:
        measured: dict[tuple[str, int], list[dict[str, str]]] = defaultdict(list)
        for row in self.rows:
            if row["warmup"] == "false":
                measured[(row["profile"], int(row["query_sequences"]))].append(row)

        for count in (1, 10, 100, 1000):
            curated = measured[("curated", count)]
            img = measured[("img", count)]
            self.assertEqual(len(curated), 3)
            self.assertEqual(len(img), 3)
            input_nucleotides = {int(row["total_nucleotides"]) for row in curated + img}
            accepted_loci = {int(row["hit_count"]) for row in curated + img}
            self.assertEqual(len(input_nucleotides), 1)
            self.assertEqual(len(accepted_loci), 1)

            curated_seconds = statistics.median(
                float(row["elapsed_seconds"]) for row in curated
            )
            curated_rss = statistics.median(
                int(row["peak_rss_kb"]) for row in curated
            ) / (1024 * 1024)
            img_seconds = statistics.median(float(row["elapsed_seconds"]) for row in img)
            img_rss = statistics.median(
                int(row["peak_rss_kb"]) for row in img
            ) / (1024 * 1024)
            expected = (
                f"| {count:,} | {input_nucleotides.pop():,} | {accepted_loci.pop():,} | "
                f"{curated_seconds:.2f} s | {curated_rss:.3f} GiB | "
                f"{img_seconds:.2f} s | {img_rss:.3f} GiB |"
            )
            self.assertIn(expected, self.page)

        self.assertNotIn("Measured elapsed times", self.page)
        self.assertIn("1,000 query sequences", self.page)


if __name__ == "__main__":
    unittest.main()
