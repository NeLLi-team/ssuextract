import gzip
import hashlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import build_database_release as builder


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
) -> builder.TaxonomyAssignment:
    return builder.TaxonomyAssignment(
        taxonomy_assignment_id="TAX_" + identifier,
        source_record_id="SRC_" + identifier,
        sequence_id="SSU_test",
        taxonomy=taxonomy,
        taxonomy_source=source,
        domain=taxonomy[0],
        compartment=compartment,
        assignment_method=method,
        evidence="test evidence" if method != "native" else "",
    )


class HeaderParserTests(unittest.TestCase):
    def test_pr2_511_header_parses_fixed_fields(self) -> None:
        parsed = builder.parse_pr2_header(PR2_HEADER)
        self.assertEqual(parsed.accession, "AB353770")
        self.assertEqual(parsed.coordinates, (1, 1740))
        self.assertEqual(parsed.orientation, "U")
        self.assertEqual(parsed.gene, "18S_rRNA")
        self.assertEqual(parsed.compartment, "nucleus")
        self.assertEqual(len(parsed.taxonomy), 9)
        self.assertEqual(parsed.taxonomy[:3], ("Eukaryota", "TSAR", "Alveolata"))
        self.assertTrue(builder.include_pr2_header(parsed))

    def test_pr2_rejects_wrong_taxonomy_width(self) -> None:
        with self.assertRaisesRegex(builder.TaxonomyError, "13 pipe-delimited"):
            builder.parse_pr2_header("AB1.1.10_U|18S_rRNA|nucleus||Eukaryota")

    def test_pr2_accepts_reverse_source_coordinates(self) -> None:
        parsed = builder.parse_pr2_header(
            "AB512182.1072.1_U|18S_rRNA|nucleus||Eukaryota|TSAR|Alveolata|"
            "Dinoflagellata|Dinophyceae|Gymnodiniales|Gymnodiniaceae|Gymnodinium|"
            "Gymnodinium_sp."
        )
        self.assertEqual(parsed.accession, "AB512182")
        self.assertEqual(parsed.coordinates, (1072, 1))

    def test_pr2_empty_compartment_is_explicitly_unknown_for_included_eukaryote(self) -> None:
        header = PR2_HEADER.replace("|nucleus|", "||")
        parsed = builder.parse_pr2_header(header)
        self.assertEqual(parsed.compartment, "")
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "pr2.fasta"
            fasta.write_text(">" + header + "\nATGC\n", encoding="utf-8")
            records = list(
                builder.iter_curated_pr2_records(fasta, source_version="5.1.1")
            )
        self.assertEqual(records[0].compartment, "unknown")

    def test_pr2_preserves_zero_coordinate_sentinel(self) -> None:
        header = PR2_HEADER.replace("AB353770.1.1740", "AAYK01000203.1.0")
        parsed = builder.parse_pr2_header(header)
        self.assertEqual(parsed.coordinates, (1, 0))

    def test_pr2_organellar_domain_suffix_is_moved_to_compartment(self) -> None:
        taxonomy = ("Eukaryota:plas", *builder.parse_pr2_header(PR2_HEADER).taxonomy[1:])
        self.assertEqual(builder.normalize_pr2_taxonomy(taxonomy)[0], "Eukaryota")

    def test_silva_header_and_rank_aware_six_rank_path(self) -> None:
        parsed = builder.parse_silva_header(
            "AY846379.1.1791 Bacteria;Firmicutes;Bacilli;Bacillales;"
            "Bacillaceae;Bacillus;Bacillus subtilis"
        )
        rank_table = builder.parse_silva_rank_table(
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
        self.assertTrue(builder.include_silva_header(parsed))
        self.assertEqual(
            builder.rank_silva_prokaryotic_path(parsed.taxonomy, rank_table),
            ("Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus"),
        )

    def test_source_filters_exclude_non_prok_silva_and_bacterial_pr2(self) -> None:
        silva = builder.parse_silva_header("X.1.10 Eukaryota;TSAR;Alveolata")
        bacterial_pr2 = builder.parse_pr2_header(
            "X.1.10_U|16S_rRNA|prokaryote||Bacteria|Firmicutes|Bacilli|"
            "Bacillales|Bacillaceae|Bacillus|Bacillus_subtilis|strain|isolate"
        )
        self.assertFalse(builder.include_silva_header(silva))
        self.assertFalse(builder.include_pr2_header(bacterial_pr2))

    def test_img_identifier_exposes_taxon_oid(self) -> None:
        parsed = builder.parse_img_identifier(
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
            pr2.write_text(">" + PR2_HEADER + "\nATGC\n", encoding="utf-8")
            silva_records = list(
                builder.iter_curated_silva_records(
                    silva, ranks, source_version="138.2"
                )
            )
            pr2_records = list(
                builder.iter_curated_pr2_records(pr2, source_version="5.1.1")
            )
        self.assertEqual(len(silva_records), 1)
        self.assertEqual(silva_records[0].marker, "16S")
        self.assertEqual(silva_records[0].taxonomy[:3], ("Bacteria", "Firmicutes", "Bacilli"))
        self.assertEqual(len(pr2_records), 1)
        self.assertEqual(pr2_records[0].marker, "18S")
        self.assertEqual(pr2_records[0].compartment, "nucleus")

    def test_img_iterator_removes_embedded_legacy_references(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "mixed.fna"
            fasta.write_text(
                ">REF_SILVA_legacy\nATGC\n"
                ">REF_SPR_legacy\nATGC\n"
                ">IMG_3300000305.a_contig\nATGC\n",
                encoding="utf-8",
            )
            records = list(
                builder.iter_img_records(fasta, "16S", source_version="2025")
            )
        self.assertEqual([record.source_identifier for record in records], ["IMG_3300000305.a_contig"])
        self.assertEqual(records[0].taxon_oid, "3300000305")


class SequenceModelTests(unittest.TestCase):
    def test_u_t_and_case_normalization_deduplicates_without_trimming(self) -> None:
        first = builder.PreparedSourceRecord(
            "PR2", "5.1.1", "pr2-a", "pr2-a", "au gc", "18S",
            ("Eukaryota", "TSAR"), "PR2", "nucleus",
        )
        second = builder.PreparedSourceRecord(
            "SILVA", "138.2", "silva-a", "silva-a", "ATGC", "18S",
            ("Eukaryota", "TSAR"), "SILVA", "",
        )
        model = builder.build_deduplicated_model([second, first])
        expected_id = builder.sequence_identifier("ATGC")
        self.assertEqual([record.sequence_id for record in model.sequences], [expected_id])
        self.assertEqual(len(model.source_records), 2)
        self.assertEqual(len(model.taxonomy_assignments), 2)
        self.assertEqual(model.sequences[0].sequence, "ATGC")
        self.assertEqual(model.sequences[0].markers, ("18S",))
        self.assertEqual(len(expected_id), 47)

    def test_hash_does_not_trim_or_reverse_complement(self) -> None:
        self.assertNotEqual(builder.sequence_identifier("ATGC"), builder.sequence_identifier("ATG"))
        self.assertNotEqual(builder.sequence_identifier("ATGC"), builder.sequence_identifier("GCAT"))

    def test_streams_plain_and_gzip_fasta_and_rejects_invalid_records(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            plain = root / "records.fasta"
            compressed = root / "records.fasta.gz"
            plain.write_text(">one\nAT\nGC\n", encoding="utf-8")
            with gzip.open(compressed, "wt") as handle:
                handle.write(">two\nAUGC\n")
            self.assertEqual(list(builder.iter_fasta(plain))[0].sequence, "ATGC")
            self.assertEqual(list(builder.iter_fasta(compressed))[0].sequence, "AUGC")
        with self.assertRaisesRegex(builder.FastaFormatError, "Empty nucleotide"):
            list(builder.iter_fasta_lines([">empty\n"]))

    def test_rejects_invalid_nucleotide_record(self) -> None:
        with self.assertRaisesRegex(builder.FastaFormatError, "Invalid nucleotide"):
            list(builder.iter_fasta_lines([">protein\n", "MEEP\n"]))

    def test_marker_fasta_is_stable_and_sorted(self) -> None:
        records = [
            builder.SequenceRecord(builder.sequence_identifier("TTTT"), "TTTT"),
            builder.SequenceRecord(builder.sequence_identifier("AAAA"), "AAAA"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / "first.fasta"
            second = Path(tmp) / "second.fasta"
            builder.write_marker_fasta(first, records, line_width=2)
            builder.write_marker_fasta(second, reversed(records), line_width=2)
            content = first.read_bytes()
            self.assertEqual(content, second.read_bytes())
            headers = [line for line in first.read_text().splitlines() if line.startswith(">")]
            self.assertEqual(headers, sorted(headers))

    def test_exact_sequence_can_belong_to_both_marker_databases(self) -> None:
        records = [
            builder.PreparedSourceRecord(
                "PR2", "5.1.1", "pr2-16", "pr2-16", "ATGC", "16S",
                ("Eukaryota", "Archaeplastida"), "PR2", "plastid",
            ),
            builder.PreparedSourceRecord(
                "PR2", "5.1.1", "pr2-18", "pr2-18", "ATGC", "18S",
                ("Eukaryota", "Archaeplastida"), "PR2", "nucleus",
            ),
        ]
        model = builder.build_deduplicated_model(records)
        self.assertEqual(len(model.sequences), 1)
        self.assertEqual(model.sequences[0].markers, ("16S", "18S"))


class TaxonomyPolicyTests(unittest.TestCase):
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
            builder.PreparedSourceRecord(
                "PR2", "5.1.1", "euk", "euk", "ATGC", "18S",
                ("Eukaryota", "TSAR"), "PR2", "nucleus",
            ),
            builder.PreparedSourceRecord(
                "SILVA", "138.2", "prok", "prok", "ATGC", "16S",
                ("Bacteria", "Firmicutes"), "SILVA", "",
            ),
        ]
        model = builder.build_deduplicated_model(records)
        preferred = model.preferred_taxonomy[0]
        self.assertTrue(preferred.cross_domain_conflict)
        self.assertEqual(preferred.domain, "ambiguous")
        self.assertEqual(preferred.taxonomy, ())
        alternatives = json.loads(preferred.taxonomy_alternatives)
        self.assertEqual({row["domain"] for row in alternatives}, {"Bacteria", "Eukaryota"})
        builder.validate_release(model)

    def test_img_taxonomy_must_be_explicit_and_evidence_bearing(self) -> None:
        record = builder.PreparedSourceRecord(
            "IMG", "2025", "IMG_3300000305.a_contig", "IMG_3300000305.a_contig",
            "ATGC", "18S", taxon_oid="3300000305",
        )
        model = builder.build_deduplicated_model([record])
        with self.assertRaisesRegex(builder.ReleaseValidationError, "IMG taxonomy is missing"):
            builder.validate_release(model)
        with self.assertRaisesRegex(builder.TaxonomyError, "require taxonomy_source"):
            builder.ingest_derived_cluster_assignments(
                [{"source_identifier": record.source_identifier, "taxonomy": "Eukaryota;TSAR"}],
                model.source_records,
            )
        derived = builder.ingest_derived_cluster_assignments(
            [
                {
                    "source_identifier": record.source_identifier,
                    "taxonomy": "Eukaryota;TSAR",
                    "taxonomy_source": "PR2",
                    "method": "updated_reference_cluster",
                    "evidence": "cluster C1 matched PR2 sequence X",
                    "compartment": "nucleus",
                }
            ],
            model.source_records,
        )
        completed = builder.add_taxonomy_assignments(model, derived)
        builder.validate_release(completed)


class PrivacyAndOutputTests(unittest.TestCase):
    def test_img_metadata_allowlist_nulls_and_ranges(self) -> None:
        location = builder.clean_img_metadata_row(
            {
                "taxon_oid": "3300000305",
                "Latitude": "",
                "Longitude": "-73.91",
                "Contact Email": "must-not-propagate@example.org",
            }
        )
        self.assertEqual(location, builder.ImgLocation("3300000305", None, -73.91))
        self.assertEqual(set(vars(location)), set(builder.IMG_LOCATION_COLUMNS))
        with self.assertRaisesRegex(builder.ReleaseValidationError, "outside"):
            builder.clean_img_metadata_row(
                {"taxon_oid": "1", "latitude": "190.1", "longitude": "0"}
            )
        with self.assertRaisesRegex(builder.ReleaseValidationError, "outside"):
            builder.clean_img_metadata_row(
                {"taxon_oid": "1", "latitude": "nan", "longitude": "0"}
            )
        corrections = []
        cleaned = builder.clean_img_metadata_row(
            {"taxon_oid": "3300016796", "latitude": "120.06", "longitude": "-34.175"},
            corrections=corrections,
        )
        self.assertEqual(cleaned, builder.ImgLocation("3300016796", None, None))
        self.assertEqual(corrections[0]["action"], "dropped_invalid_pair")
        with self.assertRaisesRegex(builder.ReleaseValidationError, "latitude outside"):
            builder.clean_img_metadata_row(
                {
                    "taxon_oid": "3300016796",
                    "latitude": "120.06",
                    "longitude": "-34.175",
                }
            )
        with self.assertRaisesRegex(builder.ReleaseValidationError, "Forbidden PII"):
            builder.validate_privacy_columns(["taxon_oid", "contact_email"])

    def test_parquet_writers_emit_exact_tables_and_privacy_schema(self) -> None:
        record = builder.PreparedSourceRecord(
            "PR2", "5.1.1", "euk", PR2_HEADER, "ATGC", "18S",
            ("Eukaryota", "TSAR"), "PR2", "nucleus",
        )
        model = builder.build_deduplicated_model([record])
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            builder.write_release_tables(
                output, model, [builder.ImgLocation("3300000305", 1.5, None)]
            )
            import duckdb

            connection = duckdb.connect(":memory:")
            try:
                self.assertEqual(
                    connection.execute(
                        "SELECT count(*) FROM read_parquet(?)", [str(output / "sequences.parquet")]
                    ).fetchone()[0],
                    1,
                )
                sequence_columns = {
                    row[0]
                    for row in connection.execute(
                        "DESCRIBE SELECT * FROM read_parquet(?)",
                        [str(output / "sequences.parquet")],
                    ).fetchall()
                }
                columns = {
                    row[0]
                    for row in connection.execute(
                        "DESCRIBE SELECT * FROM read_parquet(?)",
                        [str(output / "img_location.parquet")],
                    ).fetchall()
                }
            finally:
                connection.close()
            self.assertEqual(columns, set(builder.IMG_LOCATION_COLUMNS))
            self.assertEqual(sequence_columns, {"sequence_id", "length", "sha256", "markers"})


class FetchTests(unittest.TestCase):
    def test_fetch_rejects_non_https_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(builder.SourceIntegrityError, "must use HTTPS"):
                builder.fetch_artifact(
                    {"url": "http://example.org/source.gz", "bytes": 1, "sha256": "0" * 64},
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
            with self.assertRaisesRegex(builder.SourceIntegrityError, "SHA-256 mismatch"):
                builder.fetch_artifact(source, target, opener=opener)
            self.assertFalse(target.exists())
            self.assertEqual(list(Path(tmp).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
