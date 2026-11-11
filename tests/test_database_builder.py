import gzip
import hashlib
import io
import json
import sys
import tempfile
import unittest
from dataclasses import astuple, replace
from pathlib import Path
from unittest import mock

import duckdb

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import build_database_release as builder
import database_contracts as contracts  # noqa: E402  # needs the scripts path set above
import database_release_io as release_io  # noqa: E402  # needs the scripts path set above
import database_sources as sources  # noqa: E402  # needs the scripts path set above
from database_contracts import (  # noqa: E402  # needs the scripts path set above
    ImgLocation,
    SourceRecord,
)
from database_release_io import (  # noqa: E402  # needs the scripts path set above
    write_img_location_parquet,
    write_source_records_parquet,
)

PR2_HEADER = (
    "AB353770.1.1740_U|18S_rRNA|nucleus||Eukaryota|TSAR|Alveolata|"
    "Dinoflagellata|Dinophyceae|Peridiniales|Kryptoperidiniaceae|"
    "Unruhdinium|Unruhdinium_kevei"
)


class FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def assignment(
    identifier: str,
    taxonomy: tuple[str, ...],
    source: str,
    *,
    method: str = "native",
    compartment: str = "",
    centroid_name: str = "",
    centroid_taxonomy: tuple[str, ...] = (),
    centroid_taxonomy_source: str = "",
) -> contracts.TaxonomyAssignment:
    return contracts.TaxonomyAssignment(
        taxonomy_assignment_id="TAX_" + identifier,
        source_record_id="SRC_" + identifier,
        sequence_id="SSU_test",
        taxonomy=taxonomy,
        taxonomy_source=source,
        domain=taxonomy[0],
        compartment=compartment,
        assignment_method=method,
        evidence="test evidence" if method != "native" else "",
        centroid_name=centroid_name,
        centroid_taxonomy=centroid_taxonomy,
        centroid_taxonomy_source=centroid_taxonomy_source,
    )


class HeaderParserTests(unittest.TestCase):
    def test_silva_rank_table_rejects_internal_empty_rank(self) -> None:
        with self.assertRaisesRegex(contracts.TaxonomyError, "empty rank"):
            sources.parse_silva_rank_table(["Bacteria;;Firmicutes;\t3\tphylum\n"])

    def test_pr2_511_header_parses_fixed_fields(self) -> None:
        parsed = sources.parse_pr2_header(PR2_HEADER)
        self.assertEqual(parsed.accession, "AB353770")
        self.assertEqual(parsed.coordinates, (1, 1740))
        self.assertEqual(parsed.orientation, "U")
        self.assertEqual(parsed.gene, "18S_rRNA")
        self.assertEqual(parsed.compartment, "nucleus")
        self.assertEqual(len(parsed.taxonomy), 9)
        self.assertEqual(parsed.taxonomy[:3], ("Eukaryota", "TSAR", "Alveolata"))
        self.assertTrue(sources.include_pr2_header(parsed))

    def test_pr2_rejects_wrong_taxonomy_width(self) -> None:
        with self.assertRaisesRegex(contracts.TaxonomyError, "13 pipe-delimited"):
            sources.parse_pr2_header("AB1.1.10_U|18S_rRNA|nucleus||Eukaryota")

    def test_pr2_accepts_reverse_source_coordinates(self) -> None:
        parsed = sources.parse_pr2_header(
            "AB512182.1072.1_U|18S_rRNA|nucleus||Eukaryota|TSAR|Alveolata|"
            "Dinoflagellata|Dinophyceae|Gymnodiniales|Gymnodiniaceae|Gymnodinium|"
            "Gymnodinium_sp."
        )
        self.assertEqual(parsed.accession, "AB512182")
        self.assertEqual(parsed.coordinates, (1072, 1))

    def test_pr2_empty_compartment_is_explicitly_unknown_for_included_eukaryote(self) -> None:
        header = PR2_HEADER.replace("|nucleus|", "||")
        parsed = sources.parse_pr2_header(header)
        self.assertEqual(parsed.compartment, "")
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "pr2.fasta"
            fasta.write_text(">" + header + "\nATGC\n", encoding="utf-8")
            records = list(
                sources.iter_curated_pr2_records(
                    fasta, source_version="5.1.1", skipped_counts={}
                )
            )
        self.assertEqual(records[0].compartment, "unknown")

    def test_pr2_preserves_zero_coordinate_sentinel(self) -> None:
        header = PR2_HEADER.replace("AB353770.1.1740", "AAYK01000203.1.0")
        parsed = sources.parse_pr2_header(header)
        self.assertEqual(parsed.coordinates, (1, 0))

    def test_pr2_organellar_domain_suffix_is_moved_to_compartment(self) -> None:
        taxonomy = (
            "Eukaryota:plas",
            *sources.parse_pr2_header(PR2_HEADER).taxonomy[1:],
        )
        self.assertEqual(
            sources.normalize_pr2_taxonomy(taxonomy, compartment="plastid")[0],
            "Eukaryota",
        )

    def test_pr2_compartment_suffix_is_removed_from_every_rank(self) -> None:
        taxonomy = tuple(
            f"{taxon}:plas"
            for taxon in (
                "Eukaryota",
                "Archaeplastida",
                "Streptophyta",
                "Streptophyta_X",
                "Embryophyceae",
                "Embryophyceae_X",
                "Embryophyceae_XX",
                "Silene",
                "Silene_latifolia",
            )
        )
        self.assertEqual(
            sources.normalize_pr2_taxonomy(taxonomy, compartment="plastid"),
            (
                "Eukaryota",
                "Archaeplastida",
                "Streptophyta",
                "Streptophyta_X",
                "Embryophyceae",
                "Embryophyceae_X",
                "Embryophyceae_XX",
                "Silene",
                "Silene_latifolia",
            ),
        )

    def test_pr2_suffix_must_match_compartment(self) -> None:
        with self.assertRaisesRegex(
            contracts.TaxonomyError, "requires compartment 'plastid'"
        ):
            sources.normalize_pr2_taxonomy(
                ("Eukaryota:plas", "Archaeplastida:plas"),
                compartment="nucleus",
            )

    def test_pr2_rejects_mixed_taxonomy_suffixes(self) -> None:
        with self.assertRaisesRegex(
            contracts.TaxonomyError, "mixed compartment suffixes"
        ):
            sources.normalize_pr2_taxonomy(
                ("Eukaryota:plas", "Archaeplastida:mito"),
                compartment="plastid",
            )

    def test_pr2_recognizes_every_compartment_suffix(self) -> None:
        cases = {
            "apic": "apicoplast",
            "chro": "chromatophore",
            "chrom": "chromatophore",
            "mito": "mitochondria",
            "nucl": "nucleomorph",
            "plas": "plastid",
        }
        for suffix, compartment in cases.items():
            with self.subTest(suffix=suffix):
                self.assertEqual(
                    sources.normalize_pr2_taxonomy(
                        (f"Eukaryota:{suffix}", f"lineage:{suffix}"),
                        compartment=compartment,
                    ),
                    ("Eukaryota", "lineage"),
                )

    def test_pr2_rejects_unknown_taxonomy_suffix(self) -> None:
        with self.assertRaisesRegex(
            contracts.TaxonomyError, "Unknown PR2 taxonomy suffix"
        ):
            sources.normalize_pr2_taxonomy(
                ("Eukaryota:future", "lineage:future"),
                compartment="plastid",
            )

    def test_pr2_rejects_repeated_taxonomy_suffix(self) -> None:
        with self.assertRaisesRegex(
            contracts.TaxonomyError, "Unknown PR2 taxonomy suffix"
        ):
            sources.normalize_pr2_taxonomy(
                ("Eukaryota:plas:plas",),
                compartment="plastid",
            )

    def test_pr2_accepts_clean_taxonomy_with_organelle_compartment(self) -> None:
        self.assertEqual(
            sources.normalize_pr2_taxonomy(
                ("Eukaryota", "Archaeplastida"), compartment="plastid"
            ),
            ("Eukaryota", "Archaeplastida"),
        )

    def test_silva_header_and_rank_aware_six_rank_path(self) -> None:
        parsed = sources.parse_silva_header(
            "AY846379.1.1791 Bacteria;Firmicutes;Bacilli;Bacillales;"
            "Bacillaceae;Bacillus;Bacillus subtilis"
        )
        rank_table = sources.parse_silva_rank_table(
            [
                "Bacteria;\t2\tdomain\t\t138.2\n",
                "Bacteria;Firmicutes;\t1239\tphylum\t\t138.2\n",
                "Bacteria;Firmicutes;Bacilli;\t91061\tclass\t\t138.2\n",
                "Bacteria;Firmicutes;Bacilli;Bacillales;\t1385\torder\t\t138.2\n",
                "Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;\t186817\tfamily\t\t138.2\n",
                "Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;"
                "\t1386\tgenus\t\t138.2\n",
            ]
        )
        self.assertEqual(parsed.accession, "AY846379")
        self.assertEqual(parsed.coordinates, (1, 1791))
        self.assertTrue(sources.include_silva_header(parsed))
        self.assertEqual(
            sources.rank_silva_prokaryotic_path(parsed.taxonomy, rank_table),
            (
                "Bacteria",
                "Firmicutes",
                "Bacilli",
                "Bacillales",
                "Bacillaceae",
                "Bacillus",
            ),
        )

    def test_source_filters_exclude_non_prok_silva_and_bacterial_pr2(self) -> None:
        silva = sources.parse_silva_header("X.1.10 Eukaryota;TSAR;Alveolata")
        bacterial_pr2 = sources.parse_pr2_header(
            "X.1.10_U|16S_rRNA|prokaryote||Bacteria|Firmicutes|Bacilli|"
            "Bacillales|Bacillaceae|Bacillus|Bacillus_subtilis|strain|isolate"
        )
        self.assertFalse(sources.include_silva_header(silva))
        self.assertFalse(sources.include_pr2_header(bacterial_pr2))

    def test_img_identifier_exposes_taxon_oid(self) -> None:
        parsed = sources.parse_img_identifier(
            ">IMG_3300000305.a_bgg_mtDRAFT_1001257 optional description"
        )
        self.assertEqual(parsed.taxon_oid, "3300000305")
        self.assertEqual(parsed.source_identifier, "IMG_3300000305.a_bgg_mtDRAFT_1001257")

    def test_curated_iterators_apply_source_and_marker_policy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            silva = root / "silva.fasta"
            ranks = root / "ranks.tsv"
            pr2 = root / "pr2.fasta"
            silva.write_text(
                ">A.1.4 Bacteria;Firmicutes;Bacilli\nATGC\n"
                ">E.1.4 Eukaryota;TSAR\nATGC\n",
                encoding="utf-8",
            )
            ranks.write_text(
                "Bacteria;\t2\tdomain\t\t138.2\n"
                "Bacteria;Firmicutes;\t1239\tphylum\t\t138.2\n"
                "Bacteria;Firmicutes;Bacilli;\t91061\tclass\t\t138.2\n",
                encoding="utf-8",
            )
            pr2.write_text(
                ">"
                + PR2_HEADER
                + "\nATGC\n>"
                + PR2_HEADER.replace("Eukaryota", "Bacteria").replace(
                    "nucleus", "unknown"
                )
                + "\nATGC\n",
                encoding="utf-8",
            )
            skipped = {}
            silva_records = list(
                sources.iter_curated_silva_records(
                    silva, ranks, source_version="138.2", skipped_counts=skipped
                )
            )
            pr2_records = list(
                sources.iter_curated_pr2_records(
                    pr2, source_version="5.1.1", skipped_counts=skipped
                )
            )
        self.assertEqual(len(silva_records), 1)
        self.assertEqual(silva_records[0].marker, "16S")
        self.assertEqual(silva_records[0].taxonomy[:3], ("Bacteria", "Firmicutes", "Bacilli"))
        self.assertEqual(len(pr2_records), 1)
        self.assertEqual(pr2_records[0].marker, "18S")
        self.assertEqual(pr2_records[0].compartment, "nucleus")
        self.assertEqual(skipped, {"SILVA": 1, "PR2": 1})

    def test_img_iterator_removes_embedded_legacy_references(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "mixed.fna"
            fasta.write_text(
                ">REF_SILVA_legacy\nATGC\n"
                ">REF_SPR_legacy\nATGC\n"
                ">IMG_3300000305.a_contig\nATGC\n",
                encoding="utf-8",
            )
            skipped = {}
            records = list(
                sources.iter_img_records(
                    fasta, "16S", source_version="2025", skipped_counts=skipped
                )
            )
        self.assertEqual([record.source_identifier for record in records], ["IMG_3300000305.a_contig"])
        self.assertEqual(records[0].taxon_oid, "3300000305")
        self.assertEqual(skipped, {"IMG": 2})


class SequenceModelTests(unittest.TestCase):
    def test_u_t_and_case_normalization_deduplicates_without_trimming(self) -> None:
        first = contracts.PreparedSourceRecord(
            "PR2",
            "5.1.1",
            "pr2-a",
            "pr2-a",
            "au gc",
            "18S",
            ("Eukaryota", "TSAR"),
            "PR2",
            "nucleus",
        )
        second = contracts.PreparedSourceRecord(
            "SILVA",
            "138.2",
            "silva-a",
            "silva-a",
            "ATGC",
            "18S",
            ("Eukaryota", "TSAR"),
            "SILVA",
            "",
        )
        model = builder.build_deduplicated_model([second, first])
        expected_id = sources.sequence_identifier("ATGC")
        self.assertEqual([record.sequence_id for record in model.sequences], [expected_id])
        self.assertEqual(len(model.source_records), 2)
        self.assertEqual(len(model.taxonomy_assignments), 2)
        self.assertEqual(model.sequences[0].sequence, "ATGC")
        self.assertEqual(model.sequences[0].markers, ("18S",))
        self.assertEqual(len(expected_id), 47)

    def test_hash_does_not_trim_or_reverse_complement(self) -> None:
        self.assertNotEqual(
            sources.sequence_identifier("ATGC"), sources.sequence_identifier("ATG")
        )
        self.assertNotEqual(
            sources.sequence_identifier("ATGC"), sources.sequence_identifier("GCAT")
        )

    def test_streams_plain_and_gzip_fasta_and_rejects_invalid_records(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            plain = root / "records.fasta"
            compressed = root / "records.fasta.gz"
            plain.write_text(">one\nAT\nGC\n", encoding="utf-8")
            with gzip.open(compressed, "wt") as handle:
                handle.write(">two\nAUGC\n")
            self.assertEqual(next(iter(sources.iter_fasta(plain))).sequence, "ATGC")
            self.assertEqual(
                next(iter(sources.iter_fasta(compressed))).sequence, "AUGC"
            )
        with self.assertRaisesRegex(contracts.FastaFormatError, "Empty nucleotide"):
            list(sources.iter_fasta_lines([">empty\n"]))

    def test_rejects_invalid_nucleotide_record(self) -> None:
        with self.assertRaisesRegex(contracts.FastaFormatError, "Invalid nucleotide"):
            list(sources.iter_fasta_lines([">protein\n", "MEEP\n"]))

    def test_marker_fasta_is_stable_and_sorted(self) -> None:
        records = [
            contracts.SequenceRecord(sources.sequence_identifier("TTTT"), "TTTT"),
            contracts.SequenceRecord(sources.sequence_identifier("AAAA"), "AAAA"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / "first.fasta"
            second = Path(tmp) / "second.fasta"
            release_io.write_marker_fasta(first, records, line_width=2)
            release_io.write_marker_fasta(second, reversed(records), line_width=2)
            content = first.read_bytes()
            self.assertEqual(content, second.read_bytes())
            headers = [line for line in first.read_text().splitlines() if line.startswith(">")]
            self.assertEqual(headers, sorted(headers))

    def test_exact_sequence_can_belong_to_both_marker_databases(self) -> None:
        records = [
            contracts.PreparedSourceRecord(
                "PR2",
                "5.1.1",
                "pr2-16",
                "pr2-16",
                "ATGC",
                "16S",
                ("Eukaryota", "Archaeplastida"),
                "PR2",
                "plastid",
            ),
            contracts.PreparedSourceRecord(
                "PR2",
                "5.1.1",
                "pr2-18",
                "pr2-18",
                "ATGC",
                "18S",
                ("Eukaryota", "Archaeplastida"),
                "PR2",
                "nucleus",
            ),
        ]
        model = builder.build_deduplicated_model(records)
        self.assertEqual(len(model.sequences), 1)
        self.assertEqual(model.sequences[0].markers, ("16S", "18S"))

    def test_real_pr2_dual_marker_headers_keep_full_host_taxonomy(self) -> None:
        plastid = (
            "AB189069.42332.43801_U|16S_rRNA|plastid||Eukaryota:plas|"
            "Archaeplastida:plas|Streptophyta:plas|Streptophyta_X:plas|"
            "Embryophyceae:plas|Embryophyceae_X:plas|Embryophyceae_XX:plas|"
            "Silene:plas|Silene_latifolia:plas"
        )
        nucleus = (
            "QBIE01005410.154239.155708_U|18S_rRNA|nucleus||Eukaryota|"
            "Archaeplastida|Streptophyta|Streptophyta_X|Embryophyceae|"
            "Embryophyceae_X|Embryophyceae_XX|Silene|Silene_latifolia"
        )
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "pr2.fasta"
            fasta.write_text(
                f">{plastid}\nATGC\n>{nucleus}\nATGC\n", encoding="utf-8"
            )
            records = list(
                sources.iter_curated_pr2_records(
                    fasta, source_version="5.1.1", skipped_counts={}
                )
            )
        model = builder.build_deduplicated_model(records)
        self.assertEqual(model.sequences[0].markers, ("16S", "18S"))
        self.assertEqual(
            model.preferred_taxonomy[0].taxonomy,
            (
                "Eukaryota",
                "Archaeplastida",
                "Streptophyta",
                "Streptophyta_X",
                "Embryophyceae",
                "Embryophyceae_X",
                "Embryophyceae_XX",
                "Silene",
                "Silene_latifolia",
            ),
        )
        self.assertEqual(model.preferred_taxonomy[0].compartment, "mixed")
        release_io.validate_release(model)


class TaxonomyPolicyTests(unittest.TestCase):
    def test_preferred_taxonomy_rejects_an_internal_empty_rank(self) -> None:
        with self.assertRaisesRegex(contracts.TaxonomyError, "empty rank"):
            builder.select_preferred_taxonomy(
                "SSU_test", [assignment("a", ("Bacteria", "", "Genus"), "SILVA")]
            )

    def test_single_assignment_keeps_its_full_lineage(self) -> None:
        lineage = ("Bacteria", "Firmicutes", "Bacilli")
        preferred = builder.select_preferred_taxonomy(
            "SSU_test", [assignment("a", lineage, "SILVA")]
        )
        self.assertEqual(preferred.taxonomy, lineage)

    def test_cross_domain_img_reference_source_is_img(self) -> None:
        preferred = builder.select_preferred_taxonomy(
            "SSU_test",
            [
                assignment(
                    "a", ("Bacteria",), "SILVA", method="updated_reference_derived"
                ),
                assignment(
                    "b", ("Eukaryota",), "PR2", method="updated_reference_derived"
                ),
            ],
        )
        self.assertEqual(preferred.reference_source, "IMG")

    def test_img_assignments_match_both_marker_and_identifier(self) -> None:
        records = [
            contracts.PreparedSourceRecord(
                "IMG", "1", "IMG_1.a", "IMG_1.a", "AAAA", "16S"
            ),
            contracts.PreparedSourceRecord(
                "IMG", "1", "IMG_1.a", "IMG_1.a", "CCCC", "18S"
            ),
        ]
        model = builder.build_deduplicated_model(records)
        derived = builder.ingest_derived_cluster_assignments(
            [
                {
                    "marker": "18S",
                    "source_identifier": "IMG_1.a",
                    "taxonomy": "Eukaryota;TSAR",
                    "taxonomy_source": "PR2",
                    "assignment_method": "updated_reference_derived",
                    "evidence": "evidence",
                }
            ],
            model.source_records,
        )
        self.assertEqual(len(derived), 1)
        self.assertEqual(derived[0].sequence_id, sources.sequence_identifier("CCCC"))

    def test_classified_method_cannot_assign_unclassified_taxonomy(self) -> None:
        record = contracts.PreparedSourceRecord(
            "SILVA", "1", "silva", "silva", "AAAA", "16S", ("Unclassified",), "SILVA"
        )
        with self.assertRaisesRegex(contracts.TaxonomyError, "Unclassified"):
            builder.build_deduplicated_model([record])

    def test_release_rejects_suffix_in_combined_source_assignment(self) -> None:
        record = contracts.PreparedSourceRecord(
            "PR2",
            "1",
            "pr2",
            "pr2",
            "AAAA",
            "18S",
            ("Eukaryota", "TSAR"),
            "PR2",
            "nucleus",
        )
        model = builder.build_deduplicated_model([record])
        malformed = replace(
            model,
            taxonomy_assignments=(
                replace(
                    model.taxonomy_assignments[0],
                    taxonomy=("Eukaryota", "TSAR:plas"),
                    taxonomy_source="PR2+SILVA",
                ),
            ),
        )
        with self.assertRaisesRegex(
            contracts.ReleaseValidationError, "normalized taxonomy assignments"
        ):
            release_io.validate_release(malformed)

    def test_missing_compartment_does_not_conflict_with_known_compartment(self) -> None:
        preferred = builder.select_preferred_taxonomy(
            "SSU_test",
            [
                assignment(
                    "a",
                    ("Eukaryota", "Archaeplastida"),
                    "PR2",
                    compartment="plastid",
                ),
                assignment(
                    "b",
                    ("Eukaryota", "Archaeplastida"),
                    "PR2",
                    compartment="",
                ),
            ],
        )
        self.assertEqual(preferred.compartment, "plastid")

    def test_same_priority_disagreement_uses_lca(self) -> None:
        preferred = builder.select_preferred_taxonomy(
            "SSU_test",
            [
                assignment("a", ("Eukaryota", "TSAR", "Alveolata"), "PR2", compartment="nucleus"),
                assignment("b", ("Eukaryota", "TSAR", "Rhizaria"), "PR2", compartment="nucleus"),
            ],
        )
        self.assertEqual(preferred.taxonomy, ("Eukaryota", "TSAR"))
        self.assertEqual(preferred.assignment_method, "lowest_common_ancestor")
        self.assertFalse(preferred.cross_domain_conflict)

    def test_native_source_priority_beats_derived_img_taxonomy(self) -> None:
        preferred = builder.select_preferred_taxonomy(
            "SSU_test",
            [
                assignment("native", ("Bacteria", "Firmicutes"), "SILVA"),
                assignment(
                    "img",
                    ("Bacteria", "Proteobacteria"),
                    "PR2",
                    method="updated_reference_cluster",
                ),
            ],
        )
        self.assertEqual(preferred.taxonomy, ("Bacteria", "Firmicutes"))
        self.assertEqual(preferred.taxonomy_source, "SILVA")
        self.assertEqual(preferred.reference_source, "SILVA")

    def test_native_source_priority_beats_cross_domain_derived_img_taxonomy(self) -> None:
        preferred = builder.select_preferred_taxonomy(
            "SSU_test",
            [
                assignment("native", ("Bacteria", "Firmicutes"), "SILVA"),
                assignment(
                    "img",
                    ("Eukaryota", "TSAR"),
                    "PR2",
                    method="updated_reference_cluster",
                ),
            ],
        )
        self.assertEqual(preferred.taxonomy, ("Bacteria", "Firmicutes"))
        self.assertEqual(preferred.taxonomy_source, "SILVA")
        self.assertFalse(preferred.cross_domain_conflict)

    def test_unclassified_img_assignment_does_not_conflict_with_exact_native_reference(self) -> None:
        preferred = builder.select_preferred_taxonomy(
            "SSU_test",
            [
                assignment("native", ("Bacteria", "Firmicutes"), "SILVA"),
                assignment(
                    "img",
                    ("Unclassified",),
                    "SILVA+PR2",
                    method="updated_reference_unclassified",
                ),
            ],
        )
        self.assertEqual(preferred.taxonomy, ("Bacteria", "Firmicutes"))
        self.assertFalse(preferred.cross_domain_conflict)

    def test_pr2_native_wins_for_eukaryota(self) -> None:
        preferred = builder.select_preferred_taxonomy(
            "SSU_test",
            [
                assignment("pr2", ("Eukaryota", "TSAR", "Alveolata"), "PR2", compartment="nucleus"),
                assignment("other", ("Eukaryota", "TSAR", "Rhizaria"), "SILVA"),
            ],
        )
        self.assertEqual(preferred.taxonomy, ("Eukaryota", "TSAR", "Alveolata"))
        self.assertEqual(preferred.taxonomy_source, "PR2")

    def test_cross_domain_exact_sequence_is_explicitly_ambiguous(self) -> None:
        records = [
            contracts.PreparedSourceRecord(
                "PR2",
                "5.1.1",
                "euk",
                "euk",
                "ATGC",
                "18S",
                ("Eukaryota", "TSAR"),
                "PR2",
                "nucleus",
            ),
            contracts.PreparedSourceRecord(
                "SILVA",
                "138.2",
                "prok",
                "prok",
                "ATGC",
                "16S",
                ("Bacteria", "Firmicutes"),
                "SILVA",
                "",
            ),
        ]
        model = builder.build_deduplicated_model(records)
        preferred = model.preferred_taxonomy[0]
        self.assertTrue(preferred.cross_domain_conflict)
        self.assertEqual(preferred.domain, "ambiguous")
        self.assertEqual(preferred.taxonomy, ())
        alternatives = json.loads(preferred.taxonomy_alternatives)
        self.assertEqual({row["domain"] for row in alternatives}, {"Bacteria", "Eukaryota"})
        release_io.validate_release(model)

    def test_img_taxonomy_must_be_explicit_and_evidence_bearing(self) -> None:
        record = contracts.PreparedSourceRecord(
            "IMG",
            "2025",
            "IMG_3300000305.a_contig",
            "IMG_3300000305.a_contig",
            "ATGC",
            "18S",
            taxon_oid="3300000305",
        )
        model = builder.build_deduplicated_model([record])
        with self.assertRaisesRegex(
            contracts.ReleaseValidationError, "IMG taxonomy is missing"
        ):
            release_io.validate_release(model)
        with self.assertRaisesRegex(contracts.TaxonomyError, "require taxonomy_source"):
            builder.ingest_derived_cluster_assignments(
                [
                    {
                        "source_identifier": record.source_identifier,
                        "marker": record.marker,
                        "taxonomy": "Eukaryota;TSAR",
                    }
                ],
                model.source_records,
            )
        derived = builder.ingest_derived_cluster_assignments(
            [
                {
                    "source_identifier": record.source_identifier,
                    "marker": record.marker,
                    "taxonomy": "Eukaryota;TSAR",
                    "taxonomy_source": "PR2",
                    "assignment_method": "updated_reference_cluster",
                    "evidence": "cluster C1 matched PR2 sequence X",
                    "compartment": "nucleus",
                    "centroid_name": "IMG_centroid_1",
                    "centroid_taxonomy": "Eukaryota;TSAR;Alveolata",
                    "centroid_taxonomy_source": "PR2",
                }
            ],
            model.source_records,
        )
        completed = builder.add_taxonomy_assignments(model, derived)
        release_io.validate_release(completed)
        with self.assertRaisesRegex(contracts.TaxonomyError, "prefix of the centroid"):
            builder.ingest_derived_cluster_assignments(
                [
                    {
                        "source_identifier": record.source_identifier,
                        "marker": record.marker,
                        "taxonomy": "Eukaryota;TSAR",
                        "taxonomy_source": "PR2",
                        "assignment_method": "updated_reference_cluster",
                        "evidence": "test evidence",
                        "centroid_name": "IMG_centroid_1",
                        "centroid_taxonomy": "Eukaryota;Archaeplastida",
                        "centroid_taxonomy_source": "PR2",
                    }
                ],
                model.source_records,
            )
        with self.assertRaisesRegex(contracts.TaxonomyError, "require a centroid name"):
            builder.ingest_derived_cluster_assignments(
                [
                    {
                        "source_identifier": record.source_identifier,
                        "marker": record.marker,
                        "taxonomy": "Unclassified",
                        "taxonomy_source": "SILVA+PR2",
                        "assignment_method": "updated_reference_unclassified",
                        "evidence": "test evidence",
                    }
                ],
                model.source_records,
            )

    def test_release_rejects_malformed_centroid_evidence(self) -> None:
        record = contracts.PreparedSourceRecord(
            "IMG",
            "2025",
            "IMG_3300000305.a_contig",
            "IMG_3300000305.a_contig",
            "ATGC",
            "18S",
            taxon_oid="3300000305",
        )
        model = builder.build_deduplicated_model([record])
        derived = builder.ingest_derived_cluster_assignments(
            [
                {
                    "source_identifier": record.source_identifier,
                    "marker": record.marker,
                    "taxonomy": "Eukaryota;TSAR",
                    "taxonomy_source": "PR2",
                    "assignment_method": "updated_reference_cluster",
                    "evidence": "test evidence",
                    "centroid_name": "IMG_centroid_1",
                    "centroid_taxonomy": "Eukaryota;TSAR;Alveolata",
                    "centroid_taxonomy_source": "PR2",
                }
            ],
            model.source_records,
        )
        completed = builder.add_taxonomy_assignments(model, derived)
        assignment_row = completed.taxonomy_assignments[0]
        preferred_row = completed.preferred_taxonomy[0]
        cases = {
            "invalid centroid name": replace(
                completed,
                taxonomy_assignments=(
                    replace(assignment_row, centroid_name="IMG_centroid|second"),
                ),
            ),
            "invalid cluster centroid evidence": replace(
                completed,
                taxonomy_assignments=(
                    replace(
                        assignment_row,
                        centroid_taxonomy=("Eukaryota", "Archaeplastida"),
                    ),
                ),
            ),
            "missing cluster centroid name": replace(
                completed,
                taxonomy_assignments=(replace(assignment_row, centroid_name=""),),
            ),
            "invalid centroid names": replace(
                completed,
                preferred_taxonomy=(
                    replace(preferred_row, centroid_names='["b","a","a"]'),
                ),
            ),
            "invalid centroid taxonomy": replace(
                completed,
                preferred_taxonomy=(
                    replace(preferred_row, centroid_taxonomy=("Bacteria",)),
                ),
            ),
            "orphan centroid taxonomy": replace(
                completed,
                preferred_taxonomy=(replace(preferred_row, centroid_names=""),),
            ),
        }
        for error, malformed in cases.items():
            with (
                self.subTest(error=error),
                self.assertRaisesRegex(contracts.ReleaseValidationError, error),
            ):
                release_io.validate_release(malformed)

    def test_multiple_img_centroids_are_retained_with_taxonomy_lca(self) -> None:
        preferred = builder.select_preferred_taxonomy(
            "SSU_test",
            [
                assignment(
                    "a",
                    ("Bacteria",),
                    "SILVA",
                    method="updated_reference_cluster",
                    centroid_name="centroid-b",
                    centroid_taxonomy=("Bacteria", "Patescibacteria", "Parcubacteria"),
                    centroid_taxonomy_source="SILVA",
                ),
                assignment(
                    "b",
                    ("Bacteria",),
                    "SILVA",
                    method="updated_reference_cluster",
                    centroid_name="centroid-a",
                    centroid_taxonomy=(
                        "Bacteria",
                        "Patescibacteria",
                        "Parcubacteria",
                        "Candidatus Campbellbacteria",
                    ),
                    centroid_taxonomy_source="SILVA",
                ),
            ],
        )
        self.assertEqual(
            json.loads(preferred.centroid_names),
            ["centroid-a", "centroid-b"],
        )
        self.assertEqual(
            preferred.centroid_taxonomy,
            ("Bacteria", "Patescibacteria", "Parcubacteria"),
        )
        self.assertEqual(preferred.centroid_taxonomy_source, "SILVA")

    def test_centroid_taxonomy_sources_are_tokenized_before_merging(self) -> None:
        preferred = builder.select_preferred_taxonomy(
            "SSU_test",
            [
                assignment(
                    "a",
                    ("Eukaryota",),
                    "PR2",
                    method="updated_reference_cluster",
                    centroid_name="centroid-a",
                    centroid_taxonomy=("Eukaryota", "TSAR"),
                    centroid_taxonomy_source="PR2+SILVA",
                ),
                assignment(
                    "b",
                    ("Eukaryota",),
                    "PR2",
                    method="updated_reference_cluster",
                    centroid_name="centroid-b",
                    centroid_taxonomy=("Eukaryota", "TSAR", "Alveolata"),
                    centroid_taxonomy_source="SILVA",
                ),
            ],
        )
        self.assertEqual(preferred.centroid_taxonomy_source, "PR2+SILVA")

    def test_cross_domain_centroids_do_not_report_a_source_without_taxonomy(self) -> None:
        preferred = builder.select_preferred_taxonomy(
            "SSU_test",
            [
                assignment(
                    "a",
                    ("Bacteria",),
                    "SILVA",
                    method="updated_reference_cluster",
                    centroid_name="centroid-a",
                    centroid_taxonomy=("Bacteria", "Patescibacteria"),
                    centroid_taxonomy_source="SILVA",
                ),
                assignment(
                    "b",
                    ("Eukaryota",),
                    "PR2",
                    method="updated_reference_cluster",
                    centroid_name="centroid-b",
                    centroid_taxonomy=("Eukaryota", "TSAR"),
                    centroid_taxonomy_source="PR2",
                ),
            ],
        )
        self.assertTrue(preferred.cross_domain_conflict)
        self.assertEqual(preferred.centroid_taxonomy, ())
        self.assertEqual(preferred.centroid_taxonomy_source, "")

    def test_release_rejects_pr2_compartment_suffixes(self) -> None:
        model = builder.build_deduplicated_model(
            [
                contracts.PreparedSourceRecord(
                    "PR2",
                    "5.1.1",
                    "pr2",
                    PR2_HEADER,
                    "ATGC",
                    "18S",
                    ("Eukaryota", "Archaeplastida:plas"),
                    "PR2",
                    "plastid",
                )
            ]
        )
        with self.assertRaisesRegex(
            contracts.ReleaseValidationError,
            "PR2 compartment suffix remains in 1 normalized taxonomy assignments",
        ):
            release_io.validate_release(model)

    def test_release_rejects_preferred_domain_mismatch(self) -> None:
        model = builder.build_deduplicated_model(
            [
                contracts.PreparedSourceRecord(
                    "PR2",
                    "5.1.1",
                    "pr2",
                    PR2_HEADER,
                    "ATGC",
                    "18S",
                    ("Eukaryota", "TSAR"),
                    "PR2",
                    "nucleus",
                )
            ]
        )
        preferred = model.preferred_taxonomy[0]
        malformed = contracts.DatabaseModel(
            model.sequences,
            model.source_records,
            model.taxonomy_assignments,
            (
                contracts.PreferredTaxonomy(
                    preferred.sequence_id,
                    preferred.reference_source,
                    preferred.taxonomy,
                    preferred.taxonomy_source,
                    "Bacteria",
                    preferred.compartment,
                    preferred.assignment_method,
                    preferred.cross_domain_conflict,
                    preferred.taxonomy_alternatives,
                ),
            ),
        )
        with self.assertRaisesRegex(
            contracts.ReleaseValidationError,
            "preferred taxonomy domain mismatch in 1 rows",
        ):
            release_io.validate_release(malformed)

    def test_release_rejects_unknown_preferred_compartment(self) -> None:
        model = builder.build_deduplicated_model(
            [
                contracts.PreparedSourceRecord(
                    "PR2",
                    "5.1.1",
                    "pr2",
                    PR2_HEADER,
                    "ATGC",
                    "18S",
                    ("Eukaryota", "TSAR"),
                    "PR2",
                    "nucleus",
                )
            ]
        )
        preferred = model.preferred_taxonomy[0]
        malformed = contracts.DatabaseModel(
            model.sequences,
            model.source_records,
            model.taxonomy_assignments,
            (
                contracts.PreferredTaxonomy(
                    preferred.sequence_id,
                    preferred.reference_source,
                    preferred.taxonomy,
                    preferred.taxonomy_source,
                    preferred.domain,
                    "chloroplast",
                    preferred.assignment_method,
                    preferred.cross_domain_conflict,
                    preferred.taxonomy_alternatives,
                ),
            ),
        )
        with self.assertRaisesRegex(
            contracts.ReleaseValidationError,
            "invalid preferred taxonomy compartment",
        ):
            release_io.validate_release(malformed)

    def test_release_rejects_suffix_in_combined_source_preferred_row(self) -> None:
        model = builder.build_deduplicated_model(
            [
                contracts.PreparedSourceRecord(
                    "PR2",
                    "5.1.1",
                    "pr2",
                    PR2_HEADER,
                    "ATGC",
                    "18S",
                    ("Eukaryota", "Archaeplastida"),
                    "PR2",
                    "plastid",
                )
            ]
        )
        preferred = model.preferred_taxonomy[0]
        malformed = contracts.DatabaseModel(
            model.sequences,
            model.source_records,
            model.taxonomy_assignments,
            (
                contracts.PreferredTaxonomy(
                    preferred.sequence_id,
                    "PR2+SILVA",
                    ("Eukaryota", "Archaeplastida:plas"),
                    "PR2+SILVA",
                    preferred.domain,
                    preferred.compartment,
                    preferred.assignment_method,
                    preferred.cross_domain_conflict,
                    preferred.taxonomy_alternatives,
                ),
            ),
        )
        with self.assertRaisesRegex(
            contracts.ReleaseValidationError,
            "PR2 compartment suffix remains in 1 preferred taxonomy rows",
        ):
            release_io.validate_release(malformed)


class PrivacyAndOutputTests(unittest.TestCase):
    def test_img_writer_checks_the_stored_metadata_header(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "metadata.parquet"
            with duckdb.connect(":memory:", config={"threads": 1}) as connection:
                connection.execute(
                    "COPY (SELECT '1' AS taxon_oid, 'private' AS contact_email) "
                    "TO ? (FORMAT PARQUET)",
                    [str(output)],
                )
            with (
                mock.patch.object(release_io, "_write_parquet", return_value=output),
                self.assertRaisesRegex(
                    contracts.ReleaseValidationError, "Forbidden PII"
                ),
            ):
                release_io.write_img_location_parquet(output, [])

    def test_img_metadata_requires_canonical_taxon_oid_column(self) -> None:
        with self.assertRaisesRegex(contracts.ReleaseValidationError, "taxon_oid"):
            sources.clean_img_metadata_row({"img genome id": "3300000305"})

    def test_img_metadata_allowlist_nulls_and_ranges(self) -> None:
        location = sources.clean_img_metadata_row(
            {
                "taxon_oid": "3300000305",
                "Latitude": "",
                "Longitude": "-73.91",
                "Contact Email": "must-not-propagate@example.org",
            }
        )
        self.assertEqual(location, contracts.ImgLocation("3300000305", None, -73.91))
        self.assertEqual(set(vars(location)), set(contracts.IMG_LOCATION_COLUMNS))
        with self.assertRaisesRegex(contracts.ReleaseValidationError, "outside"):
            sources.clean_img_metadata_row(
                {"taxon_oid": "1", "latitude": "190.1", "longitude": "0"}
            )
        with self.assertRaisesRegex(contracts.ReleaseValidationError, "outside"):
            sources.clean_img_metadata_row(
                {"taxon_oid": "1", "latitude": "nan", "longitude": "0"}
            )
        corrections = []
        cleaned = sources.clean_img_metadata_row(
            {"taxon_oid": "3300016796", "latitude": "120.06", "longitude": "-34.175"},
            corrections=corrections,
        )
        self.assertEqual(cleaned, contracts.ImgLocation("3300016796", None, None))
        self.assertEqual(corrections[0]["action"], "dropped_invalid_pair")
        with self.assertRaisesRegex(
            contracts.ReleaseValidationError, "latitude outside"
        ):
            sources.clean_img_metadata_row(
                {
                    "taxon_oid": "3300016796",
                    "latitude": "120.06",
                    "longitude": "-34.175",
                }
            )
        with self.assertRaisesRegex(contracts.ReleaseValidationError, "Forbidden PII"):
            sources.validate_privacy_columns(["taxon_oid", "contact_email"])

    def test_parquet_writers_emit_exact_tables_and_privacy_schema(self) -> None:
        record = contracts.PreparedSourceRecord(
            "PR2",
            "5.1.1",
            "euk",
            PR2_HEADER,
            "ATGC",
            "18S",
            ("Eukaryota", "TSAR"),
            "PR2",
            "nucleus",
        )
        model = builder.build_deduplicated_model([record])
        locations = [contracts.ImgLocation("3300000305", 1.5, None)]
        expected_rows = {
            "sequences": [
                (
                    item.sequence_id,
                    len(item.sequence),
                    hashlib.sha256(item.sequence.encode("ascii")).hexdigest(),
                    ";".join(item.markers),
                )
                for item in model.sequences
            ],
            "source_records": [astuple(item) for item in model.source_records],
            "taxonomy_assignments": [
                tuple(
                    ";".join(value) if isinstance(value, tuple) else value
                    for value in astuple(item)
                )
                for item in model.taxonomy_assignments
            ],
            "preferred_taxonomy": [
                tuple(
                    ";".join(value) if isinstance(value, tuple) else value
                    for value in astuple(item)
                )
                for item in model.preferred_taxonomy
            ],
            "img_location": [astuple(item) for item in locations],
        }
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            release_io.write_release_tables(output, model, locations)
            with duckdb.connect(":memory:", config={"threads": 1}) as connection:
                for table, expected in expected_rows.items():
                    with self.subTest(table=table):
                        stored = connection.execute(
                            "SELECT * FROM read_parquet(?)",
                            [str(output / f"{table}.parquet")],
                        ).fetchall()
                        self.assertEqual(stored, expected)
                columns = connection.execute(
                    "DESCRIBE SELECT * FROM read_parquet(?)",
                    [str(output / "img_location.parquet")],
                ).fetchall()
            self.assertEqual(
                {column[0] for column in columns}, set(contracts.IMG_LOCATION_COLUMNS)
            )

    def test_parquet_rows_round_trip_nulls_and_special_characters(self) -> None:
        headers = [
            "plain",
            "",
            'quote " inside',
            "comma, inside",
            "line\nbreak",
            "carriage\rreturn",
            "Ünïcode τ",
            "\\N",
            "NULL",
            " leading space",
        ]
        source_records = [
            SourceRecord(
                f"SRC_{index:02d}",
                f"SSU_{index:02d}",
                "IMG",
                "1",
                f"id_{index}",
                header,
                "16S",
                None if index % 2 else header,
            )
            for index, header in enumerate(headers)
        ]
        locations = [
            ImgLocation("3300000001", 1.5, None),
            ImgLocation("3300000002", None, -122.25),
            ImgLocation("3300000003", 0.1 + 0.2, 1e-05),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            write_source_records_parquet(output / "source.parquet", source_records)
            write_img_location_parquet(output / "location.parquet", locations)
            connection = duckdb.connect(":memory:")
            try:
                source_rows = connection.execute(
                    "SELECT * FROM read_parquet(?)", [str(output / "source.parquet")]
                ).fetchall()
                location_rows = connection.execute(
                    "SELECT * FROM read_parquet(?)", [str(output / "location.parquet")]
                ).fetchall()
            finally:
                connection.close()

        self.assertEqual(
            source_rows,
            [
                (
                    record.source_record_id,
                    record.sequence_id,
                    record.reference_source,
                    record.source_version,
                    record.source_identifier,
                    record.original_header,
                    record.marker,
                    record.taxon_oid,
                )
                for record in source_records
            ],
        )
        self.assertEqual(
            location_rows,
            [(item.taxon_oid, item.latitude, item.longitude) for item in locations],
        )

    def test_parquet_row_longer_than_the_default_csv_line_limit_is_written(
        self,
    ) -> None:
        header = "x" * 2_100_000
        record = SourceRecord("SRC_long", "SSU_long", "IMG", "1", "id", header, "16S")
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "source.parquet"
            write_source_records_parquet(target, [record])
            connection = duckdb.connect(":memory:")
            try:
                stored = connection.execute(
                    "SELECT original_header FROM read_parquet(?)", [str(target)]
                ).fetchone()[0]
            finally:
                connection.close()

        self.assertEqual(stored, header)


class FetchTests(unittest.TestCase):
    def test_fetch_rejects_non_https_source(self) -> None:
        with (
            tempfile.TemporaryDirectory() as tmp,
            self.assertRaisesRegex(contracts.SourceIntegrityError, "must use HTTPS"),
        ):
            sources.fetch_artifact(
                {
                    "url": "http://example.org/source.gz",
                    "bytes": 1,
                    "sha256": "0" * 64,
                },
                Path(tmp) / "source.gz",
            )

    def test_checksum_failure_does_not_install_partial_file(self) -> None:
        payload = b"downloaded bytes"
        source = {
            "url": "https://example.org/source.gz",
            "bytes": len(payload),
            "sha256": "0" * 64,
            "md5": hashlib.md5(payload).hexdigest(),
        }

        def opener(url, timeout):
            self.assertEqual(url, source["url"])
            self.assertEqual(timeout, 60)
            return FakeResponse(payload)

        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "source.gz"
            with self.assertRaisesRegex(
                contracts.SourceIntegrityError, "SHA-256 mismatch"
            ):
                sources.fetch_artifact(source, target, opener=opener)
            self.assertFalse(target.exists())
            self.assertEqual(list(Path(tmp).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
