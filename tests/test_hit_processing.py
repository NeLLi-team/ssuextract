import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from get_cmsequences import parse_seqmap
from hit_processing import (
    CmHit,
    CmModel,
    ExtractionRegion,
    extract_regions,
    parse_cmsearch_tblout,
    read_accepted_hits,
    read_covariance_model,
    resolve_competing_model_hits,
    resolve_extraction_regions,
    select_model_hits,
    write_accepted_hits,
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
    model_accession: str = "RF00177",
    model_name: str | None = None,
    bit_score: float = 100.0,
    e_value: float = 1e-20,
) -> CmHit:
    if model_name is None:
        model_name = {
            "RF00177": "SSU_rRNA_bacteria",
            "RF01960": "SSU_rRNA_eukarya",
        }.get(model_accession, "RFTEST")
    return CmHit(
        subject=subject,
        model=model_name,
        model_accession=model_accession,
        model_from=model_from,
        model_to=model_to,
        sequence_from=sequence_from,
        sequence_to=sequence_to,
        strand=strand,
        bit_score=bit_score,
        e_value=e_value,
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
        self.assertEqual(parsed[0].model_accession, "RFTEST")
        self.assertEqual(parsed[0].bit_score, 42.0)
        self.assertEqual(parsed[0].e_value, 1e-9)
        self.assertTrue(parsed[0].included)

    def test_known_bundled_model_without_accession_fails_loudly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tblout = Path(tmp) / "hits.out"
            tblout.write_text(
                "contig1 - SSU_rRNA_bacteria - cm 1 10 3 12 + no 1 0.50 0.0 "
                "42.0 1e-9 ! description\n"
            )
            with self.assertRaisesRegex(ValueError, "expected 'RF00177'"):
                parse_cmsearch_tblout(tblout)

    def test_unknown_custom_model_without_accession_is_retained(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tblout = Path(tmp) / "hits.out"
            tblout.write_text(
                "contig1 - custom_ssu - cm 1 10 3 12 + no 1 0.50 0.0 "
                "42.0 1e-9 ! description\n"
            )
            parsed = parse_cmsearch_tblout(tblout)

        self.assertEqual(parsed[0].model_accession, "custom_ssu")

    def test_bundled_covariance_model_identity_is_explicit(self) -> None:
        bacterial = read_covariance_model(REPO / "resources/models/RF00177.cm")
        eukaryotic = read_covariance_model(REPO / "resources/models/RF01960.cm")
        self.assertEqual(
            bacterial,
            CmModel("SSU_rRNA_bacteria", "RF00177", 1533),
        )
        self.assertEqual(
            eukaryotic,
            CmModel("SSU_rRNA_eukarya", "RF01960", 1851),
        )

    def test_covariance_model_file_must_contain_one_model(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            model = Path(tmp) / "multiple.cm"
            model.write_text(
                "INFERNAL1/a\nNAME first\nCLEN 10\n//\n"
                "INFERNAL1/a\nNAME second\nCLEN 20\n//\n"
            )
            with self.assertRaisesRegex(ValueError, "contains 2 models"):
                read_covariance_model(model)

    def test_bundled_covariance_model_filename_must_match_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            model = Path(tmp) / "renamed.cm"
            model.write_text(
                "INFERNAL1/a\nNAME SSU_rRNA_bacteria\n"
                "ACC RF00177\nCLEN 1533\n//\n"
            )
            with self.assertRaisesRegex(ValueError, "must use filename RF00177.cm"):
                read_covariance_model(model)

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


class CompetingModelTests(unittest.TestCase):
    def test_accepted_hit_table_round_trips_exactly(self) -> None:
        accepted = hit(
            model_accession="RF01960",
            model_from=2,
            model_to=1850,
            sequence_from=1734,
            sequence_to=3528,
            bit_score=1776.7,
            e_value=3.8e-216,
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "accepted.tsv"
            write_accepted_hits([accepted], path)
            observed = read_accepted_hits(path)
        self.assertEqual(observed, [accepted])

    def test_accepted_hit_table_rejects_malformed_rows(self) -> None:
        accepted = hit(
            model_accession="RF00177",
            model_from=1,
            model_to=10,
            sequence_from=1,
            sequence_to=10,
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "accepted.tsv"
            write_accepted_hits([accepted], path)
            text = path.read_text()
            path.write_text(text.replace("\t+\t", "\t?\t"))
            with self.assertRaisesRegex(ValueError, "strand is invalid"):
                read_accepted_hits(path)

    def test_bacterial_model_wins_user_example_overlap(self) -> None:
        bacterial = hit(
            model_accession="RF00177",
            model_from=1,
            model_to=1533,
            sequence_from=14154,
            sequence_to=15690,
            bit_score=1612.8,
            e_value=0.0,
        )
        eukaryotic = hit(
            model_accession="RF01960",
            model_from=1,
            model_to=1851,
            sequence_from=14159,
            sequence_to=15685,
            bit_score=719.0,
            e_value=3.8e-216,
        )
        self.assertEqual(resolve_competing_model_hits([eukaryotic, bacterial]), [bacterial])

    def test_eukaryotic_model_wins_when_it_has_better_evidence(self) -> None:
        bacterial = hit(
            model_accession="RF00177",
            model_from=1,
            model_to=1533,
            sequence_from=1729,
            sequence_to=3533,
            bit_score=573.7,
            e_value=2.7e-182,
        )
        eukaryotic = hit(
            model_accession="RF01960",
            model_from=1,
            model_to=1851,
            sequence_from=1734,
            sequence_to=3528,
            bit_score=1776.7,
            e_value=0.0,
        )
        self.assertEqual(resolve_competing_model_hits([bacterial, eukaryotic]), [eukaryotic])

    def test_equal_evalue_uses_higher_bit_score(self) -> None:
        lower = hit(
            model_accession="RF01960",
            model_from=1,
            model_to=10,
            sequence_from=1,
            sequence_to=10,
            bit_score=90.0,
            e_value=0.0,
        )
        higher = hit(
            model_accession="RF00177",
            model_from=1,
            model_to=10,
            sequence_from=1,
            sequence_to=10,
            bit_score=100.0,
            e_value=0.0,
        )
        self.assertEqual(resolve_competing_model_hits([lower, higher]), [higher])

    def test_nonoverlap_opposite_strand_and_same_model_are_preserved(self) -> None:
        hits = [
            hit(
                model_accession="RF00177",
                model_from=1,
                model_to=10,
                sequence_from=1,
                sequence_to=10,
            ),
            hit(
                model_accession="RF01960",
                model_from=1,
                model_to=10,
                sequence_from=20,
                sequence_to=30,
            ),
            hit(
                model_accession="RF01960",
                model_from=1,
                model_to=10,
                sequence_from=1,
                sequence_to=10,
                strand="-",
            ),
            hit(
                model_accession="RF00177",
                model_from=1,
                model_to=10,
                sequence_from=5,
                sequence_to=15,
                bit_score=99.0,
            ),
        ]
        self.assertCountEqual(resolve_competing_model_hits(hits), hits)

    def test_overlapping_unconfigured_model_fails_loudly(self) -> None:
        bundled = hit(
            model_accession="RF00177",
            model_from=1,
            model_to=10,
            sequence_from=1,
            sequence_to=10,
        )
        custom = hit(
            model_accession="RF01959",
            model_name="SSU_rRNA_archaea",
            model_from=1,
            model_to=10,
            sequence_from=2,
            sequence_to=9,
        )
        with self.assertRaisesRegex(ValueError, "without an explicit competition rule"):
            resolve_competing_model_hits([bundled, custom])

    def test_two_overlapping_unconfigured_models_fail_loudly(self) -> None:
        left = hit(
            model_accession="custom_left",
            model_name="custom_left",
            model_from=1,
            model_to=10,
            sequence_from=1,
            sequence_to=10,
        )
        right = hit(
            model_accession="custom_right",
            model_name="custom_right",
            model_from=1,
            model_to=10,
            sequence_from=2,
            sequence_to=9,
        )
        with self.assertRaisesRegex(ValueError, "without an explicit competition rule"):
            resolve_competing_model_hits([left, right])

    def test_exact_evidence_tie_fails_loudly(self) -> None:
        left = hit(
            model_accession="RF00177",
            model_from=1,
            model_to=10,
            sequence_from=1,
            sequence_to=10,
        )
        right = hit(
            model_accession="RF01960",
            model_from=1,
            model_to=10,
            sequence_from=2,
            sequence_to=9,
        )
        with self.assertRaisesRegex(ValueError, "Indistinguishable competing"):
            resolve_competing_model_hits([left, right])

    def test_select_model_hits_uses_name_and_accession(self) -> None:
        accepted = hit(
            model_accession="RF00177",
            model_from=1,
            model_to=10,
            sequence_from=1,
            sequence_to=10,
        )
        rejected = hit(
            model_accession="RF01960",
            model_from=1,
            model_to=10,
            sequence_from=1,
            sequence_to=10,
        )
        model = CmModel("SSU_rRNA_bacteria", "RF00177", 1533)
        self.assertEqual(select_model_hits([accepted, rejected], model), [accepted])

        conflicting = replace(accepted, model="SSU_rRNA_eukarya")
        with self.assertRaisesRegex(ValueError, "identity disagrees"):
            select_model_hits([conflicting], model)


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
