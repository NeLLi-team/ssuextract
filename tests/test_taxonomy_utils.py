import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import taxonomy_utils


class TaxonomyPathTests(unittest.TestCase):
    def test_string_and_rank_sequence_have_one_normalization_contract(self) -> None:
        expected = ("Bacteria", "Firmicutes", "", "Bacillales")
        self.assertEqual(
            taxonomy_utils.taxonomy_path(" Bacteria ; Firmicutes ;; Bacillales ; "),
            expected,
        )
        self.assertEqual(
            taxonomy_utils.taxonomy_path([" Bacteria ", "Firmicutes", "", "Bacillales", ""]),
            expected,
        )

    def test_empty_domain_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-empty domain"):
            taxonomy_utils.taxonomy_path(";Firmicutes")


class LowestCommonAncestorTests(unittest.TestCase):
    def test_internal_missing_rank_ends_the_common_path(self) -> None:
        self.assertEqual(
            taxonomy_utils.lowest_common_ancestor(
                [
                    ("Bacteria", "Firmicutes", "", "OrderA"),
                    ("Bacteria", "Firmicutes", "ClassB", "OrderB"),
                ]
            ),
            ("Bacteria", "Firmicutes"),
        )

    def test_no_paths_has_no_common_ancestor(self) -> None:
        self.assertEqual(taxonomy_utils.lowest_common_ancestor([]), ())


class CommonValueTests(unittest.TestCase):
    def test_callers_choose_their_conflict_sentinel_explicitly(self) -> None:
        values = ["SILVA", "PR2", ""]
        self.assertEqual(taxonomy_utils.common_value(values), "")
        self.assertEqual(taxonomy_utils.common_value(values, conflict="mixed"), "mixed")

    def test_one_non_empty_value_is_returned(self) -> None:
        self.assertEqual(taxonomy_utils.common_value(["", "SILVA", "SILVA"]), "SILVA")


if __name__ == "__main__":
    unittest.main()
