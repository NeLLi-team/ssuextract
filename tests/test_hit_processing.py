import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from get_cmsequences import parse_seqmap
from hit_processing import (
    CmHit,
    ExtractionRegion,
    extract_regions,
    parse_cmsearch_tblout,
    resolve_extraction_regions,
)


def hit(
    *,
    subject: str = "contig1",
    model_from: int,
    model_to: int,
    sequence_from: int,
    sequence_to: int,
    strand: str = "+",
    included: bool = True,
) -> CmHit:
    return CmHit(
        subject=subject,
        model="RFTEST",
        model_from=model_from,
        model_to=model_to,
        sequence_from=sequence_from,
        sequence_to=sequence_to,
        strand=strand,
        included=included,
    )


class CmsearchParsingTests(unittest.TestCase):
    def test_description_with_spaces_stays_in_last_field(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tblout = Path(tmp) / "hits.out"
            tblout.write_text(
                "contig1 - model RFTEST cm 1 10 3 12 + no 1 0.50 0.0 "
                "42.0 1e-9 ! description with spaces\n"
            )
            parsed = parse_cmsearch_tblout(tblout)

        self.assertEqual(len(parsed), 1)
        self.assertEqual(parsed[0].subject, "contig1")
        self.assertTrue(parsed[0].included)

    def test_reported_but_not_included_hits_are_ignored(self) -> None:
        regions = resolve_extraction_regions(
            [
                hit(
                    model_from=1,
                    model_to=10,
                    sequence_from=1,
                    sequence_to=10,
                    included=False,
                )
            ],
            model_length=10,
        )
        self.assertEqual(regions, [])


class RegionResolutionTests(unittest.TestCase):
    def test_collinear_fragments_become_one_assembled_region(self) -> None:
        regions = resolve_extraction_regions(
            [
                hit(
                    model_from=1,
                    model_to=5,
                    sequence_from=10,
                    sequence_to=14,
                ),
                hit(
                    model_from=6,
                    model_to=9,
                    sequence_from=15,
                    sequence_to=18,
                ),
            ],
            model_length=10,
        )
        self.assertEqual(
            regions,
            [
                ExtractionRegion(
                    subject="contig1",
                    start=10,
                    end=18,
                    strand="+",
                    sequence_type="assembled",
                    is_assembled=True,
                )
            ],
        )

    def test_reverse_collinear_fragments_become_one_assembled_region(self) -> None:
        regions = resolve_extraction_regions(
            [
                hit(
                    model_from=1,
                    model_to=5,
                    sequence_from=30,
                    sequence_to=26,
                    strand="-",
                ),
                hit(
                    model_from=6,
                    model_to=9,
                    sequence_from=25,
                    sequence_to=22,
                    strand="-",
                ),
            ],
            model_length=10,
        )
        self.assertEqual(len(regions), 1)
        self.assertEqual((regions[0].start, regions[0].end), (22, 30))
        self.assertEqual(regions[0].strand, "-")
        self.assertTrue(regions[0].is_assembled)

    def test_out_of_order_fragments_remain_independent(self) -> None:
        regions = resolve_extraction_regions(
            [
                hit(
                    model_from=1,
                    model_to=5,
                    sequence_from=20,
                    sequence_to=24,
                ),
                hit(
                    model_from=6,
                    model_to=9,
                    sequence_from=10,
                    sequence_to=13,
                ),
            ],
            model_length=10,
        )
        self.assertEqual(len(regions), 2)
        self.assertTrue(all(region.sequence_type == "simple" for region in regions))

    def test_full_hit_is_not_reintroduced_into_fragment_assembly(self) -> None:
        regions = resolve_extraction_regions(
            [
                hit(
                    model_from=1,
                    model_to=10,
                    sequence_from=1,
                    sequence_to=10,
                ),
                hit(
                    model_from=2,
                    model_to=5,
                    sequence_from=20,
                    sequence_to=23,
                ),
            ],
            model_length=10,
        )
        self.assertEqual(len(regions), 2)
        self.assertEqual([(r.start, r.end) for r in regions], [(1, 10), (20, 23)])


class SequenceExtractionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.fasta = Path(self.tempdir.name) / "input.fna"
        self.fasta.write_text(">contig1 description retained\nAACCGGTT\n")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_one_based_inclusive_plus_interval(self) -> None:
        records = extract_regions(
            self.fasta,
            [ExtractionRegion("contig1", 1, 4, "+", "simple", False)],
            minimum_length=0,
        )
        self.assertEqual(records[0].sequence, "AACC")

    def test_one_based_inclusive_minus_interval(self) -> None:
        records = extract_regions(
            self.fasta,
            [ExtractionRegion("contig1", 1, 4, "-", "simple", False)],
            minimum_length=0,
        )
        self.assertEqual(records[0].sequence, "GGTT")

    def test_fasta_description_does_not_change_identifier(self) -> None:
        records = extract_regions(
            self.fasta,
            [ExtractionRegion("contig1", 2, 5, "+", "simple", False)],
            minimum_length=0,
        )
        self.assertEqual(len(records), 1)

    def test_missing_subject_fails_instead_of_silently_disappearing(self) -> None:
        with self.assertRaisesRegex(ValueError, "missing from FASTA"):
            extract_regions(
                self.fasta,
                [ExtractionRegion("missing", 1, 4, "+", "simple", False)],
                minimum_length=0,
            )

    def test_seqmap_without_final_newline_keeps_last_record(self) -> None:
        seqmap = Path(self.tempdir.name) / "input.seqmap"
        seqmap.write_text("contig1\t1\t4\t+\tsimple")
        regions = parse_seqmap(str(seqmap))
        records = extract_regions(self.fasta, regions, minimum_length=0)
        self.assertEqual(records[0].sequence, "AACC")


if __name__ == "__main__":
    unittest.main()

