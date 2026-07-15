#!/usr/bin/env python3
"""Deterministic, provenance-preserving SSU reference database builder.

Source acquisition and source-specific parsing live in :mod:`database_sources`.
Release validation and serialization live in :mod:`database_release_io`. This
module owns exact-sequence deduplication, taxonomy selection, and CLI orchestration.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

from database_contracts import (
    BuildError,
    DatabaseModel,
    FastaFormatError,
    FastaRecord,
    IMG_LOCATION_COLUMNS,
    ImgIdentifier,
    ImgLocation,
    PR2_COMPARTMENTS,
    PR2_RANKS,
    PROKARYOTIC_DOMAINS,
    PreferredTaxonomy,
    PreparedSourceRecord,
    Pr2Header,
    ReleaseValidationError,
    SILVA_PROKARYOTIC_RANKS,
    SequenceRecord,
    SilvaHeader,
    SourceIntegrityError,
    SourceRecord,
    TaxonomyAssignment,
    TaxonomyError,
)
from database_sources import (
    DEFAULT_SOURCES,
    clean_img_metadata,
    clean_img_metadata_row,
    fetch_artifact,
    fetch_configured_source,
    include_pr2_header,
    include_silva_header,
    iter_curated_pr2_records,
    iter_curated_silva_records,
    iter_fasta,
    iter_fasta_lines,
    iter_img_records,
    load_silva_rank_table,
    load_source_config,
    normalize_pr2_taxonomy,
    normalize_sequence_for_hashing,
    parse_img_identifier,
    parse_pr2_header,
    parse_silva_header,
    parse_silva_rank_table,
    rank_silva_prokaryotic_path,
    sequence_identifier,
    validate_nucleotide_sequence,
    validate_privacy_columns,
)
from database_release_io import (
    validate_release,
    write_img_location_parquet,
    write_marker_fasta,
    write_preferred_taxonomy_parquet,
    write_release_tables,
    write_sequences_parquet,
    write_source_records_parquet,
    write_taxonomy_assignments_parquet,
)
from taxonomy_utils import lowest_common_ancestor as _shared_lca
from taxonomy_utils import taxonomy_path


def _stable_id(prefix: str, *values: str) -> str:
    payload = "\0".join(values).encode("utf-8")
    return prefix + hashlib.sha256(payload).hexdigest()


def _taxonomy_tuple(value: Sequence[str] | str) -> tuple[str, ...]:
    try:
        return taxonomy_path(value)
    except ValueError as error:
        raise TaxonomyError(str(error)) from error


def lowest_common_ancestor(taxonomies: Iterable[Sequence[str] | str]) -> tuple[str, ...]:
    paths = [_taxonomy_tuple(taxonomy) for taxonomy in taxonomies]
    if not paths:
        raise TaxonomyError("Cannot compute an LCA without taxonomy paths")
    return _shared_lca(paths)


def taxonomy_priority(assignment: TaxonomyAssignment) -> int:
    """Return deterministic policy priority; lower values win."""

    method = assignment.assignment_method.lower()
    native = method == "native"
    source = assignment.taxonomy_source.upper()
    if assignment.domain in PROKARYOTIC_DOMAINS and native and source == "SILVA":
        return 0
    if (
        assignment.domain == "Eukaryota" or assignment.compartment.lower() in PR2_COMPARTMENTS
    ) and native and source == "PR2":
        return 0
    if method in {"updated_reference_cluster", "updated_reference_derived"}:
        return 10
    return 20


def select_preferred_taxonomy(
    sequence_id: str, assignments: Iterable[TaxonomyAssignment]
) -> PreferredTaxonomy:
    candidates = list(assignments)
    if not candidates:
        raise TaxonomyError(f"No taxonomy assignments for {sequence_id}")
    # Apply the source policy before testing for a domain conflict. Otherwise a
    # lower-priority, cluster-derived IMG assignment could turn an exact native
    # SILVA or PR2 match into a false cross-domain ambiguity.
    best_priority = min(taxonomy_priority(assignment) for assignment in candidates)
    best = [assignment for assignment in candidates if taxonomy_priority(assignment) == best_priority]
    domains = {
        assignment.domain
        for assignment in best
        if assignment.domain and assignment.domain != "Unclassified"
    }
    cross_domain_conflict = len(domains) > 1
    if cross_domain_conflict:
        alternatives = [
            {
                "taxonomy_source": assignment.taxonomy_source,
                "taxonomy": ";".join(assignment.taxonomy),
                "domain": assignment.domain,
                "compartment": assignment.compartment,
                "assignment_method": assignment.assignment_method,
            }
            for assignment in sorted(
                best,
                key=lambda item: (
                    item.taxonomy_source,
                    item.taxonomy,
                    item.compartment,
                    item.source_record_id,
                ),
            )
        ]
        return PreferredTaxonomy(
            sequence_id=sequence_id,
            reference_source="+".join(
                sorted({assignment.taxonomy_source for assignment in best})
            ),
            taxonomy=(),
            taxonomy_source="+".join(
                sorted({assignment.taxonomy_source for assignment in best})
            ),
            domain="ambiguous",
            compartment="mixed",
            assignment_method="cross_domain_ambiguous_exact_sequence",
            cross_domain_conflict=True,
            taxonomy_alternatives=json.dumps(
                alternatives, ensure_ascii=True, separators=(",", ":"), sort_keys=True
            ),
        )
    taxonomy = lowest_common_ancestor(assignment.taxonomy for assignment in best)
    if not taxonomy:
        raise TaxonomyError(f"No common taxonomy domain for {sequence_id}")
    disagreement = len({assignment.taxonomy for assignment in best}) > 1
    taxonomy_sources = sorted({assignment.taxonomy_source for assignment in best})
    compartments = sorted({assignment.compartment for assignment in best if assignment.compartment})
    reference_sources = sorted(
        {
            "IMG" if assignment.assignment_method.lower() != "native" else assignment.taxonomy_source
            for assignment in best
        }
    )
    return PreferredTaxonomy(
        sequence_id=sequence_id,
        reference_source="+".join(reference_sources),
        taxonomy=taxonomy,
        taxonomy_source="+".join(taxonomy_sources),
        domain=taxonomy[0],
        compartment=compartments[0] if len(compartments) == 1 else "",
        assignment_method="lowest_common_ancestor" if disagreement else best[0].assignment_method,
        cross_domain_conflict=cross_domain_conflict,
    )


def build_deduplicated_model(records: Iterable[PreparedSourceRecord]) -> DatabaseModel:
    """Deduplicate by normalized exact sequence while preserving every source row."""

    prepared = sorted(
        records,
        key=lambda record: (
            record.reference_source,
            record.source_version,
            record.source_identifier,
            record.original_header,
            record.marker,
        ),
    )
    sequences_by_id: dict[str, list[PreparedSourceRecord]] = defaultdict(list)
    source_records: list[SourceRecord] = []
    assignments: list[TaxonomyAssignment] = []
    seen_source_ids: set[str] = set()

    for record in prepared:
        sequence_id = sequence_identifier(record.sequence)
        sequences_by_id[sequence_id].append(record)
        source_record_id = _stable_id(
            "SRC_",
            record.reference_source,
            record.source_version,
            record.source_identifier,
            record.original_header,
            record.marker,
        )
        if source_record_id in seen_source_ids:
            raise ReleaseValidationError(
                f"Duplicate source record identity: {record.reference_source}:{record.source_identifier}"
            )
        seen_source_ids.add(source_record_id)
        source_record = SourceRecord(
            source_record_id,
            sequence_id,
            record.reference_source,
            record.source_version,
            record.source_identifier,
            record.original_header,
            record.marker,
            record.taxon_oid,
        )
        source_records.append(source_record)
        if record.taxonomy:
            taxonomy = _taxonomy_tuple(record.taxonomy)
            assignment_id = _stable_id(
                "TAX_",
                source_record_id,
                ";".join(taxonomy),
                record.taxonomy_source,
                record.compartment,
                record.assignment_method,
                record.evidence,
            )
            assignments.append(
                TaxonomyAssignment(
                    assignment_id,
                    source_record_id,
                    sequence_id,
                    taxonomy,
                    record.taxonomy_source,
                    taxonomy[0],
                    record.compartment,
                    record.assignment_method,
                    record.evidence,
                )
            )

    sequence_records = []
    for sequence_id, members in sorted(sequences_by_id.items()):
        normalized_sequences = {
            normalize_sequence_for_hashing(record.sequence) for record in members
        }
        if len(normalized_sequences) != 1:  # defensive: IDs already derive from this value
            raise ReleaseValidationError(f"hash collision for {sequence_id}")
        markers = tuple(sorted({record.marker for record in members}))
        sequence_records.append(
            SequenceRecord(sequence_id, normalized_sequences.pop(), markers)
        )

    by_sequence: dict[str, list[TaxonomyAssignment]] = defaultdict(list)
    for assignment in assignments:
        by_sequence[assignment.sequence_id].append(assignment)
    preferred = [
        select_preferred_taxonomy(sequence.sequence_id, by_sequence[sequence.sequence_id])
        for sequence in sequence_records
        if by_sequence[sequence.sequence_id]
    ]
    return DatabaseModel(
        tuple(sequence_records),
        tuple(sorted(source_records, key=lambda record: record.source_record_id)),
        tuple(sorted(assignments, key=lambda assignment: assignment.taxonomy_assignment_id)),
        tuple(sorted(preferred, key=lambda taxonomy: taxonomy.sequence_id)),
    )


def ingest_derived_cluster_assignments(
    rows: Iterable[Mapping[str, object]], source_records: Iterable[SourceRecord]
) -> tuple[TaxonomyAssignment, ...]:
    """Ingest explicit IMG cluster-derived assignments; never infer or invent them."""

    by_identifier: dict[str, list[SourceRecord]] = defaultdict(list)
    for source_record in source_records:
        if source_record.reference_source.upper() == "IMG":
            by_identifier[source_record.source_identifier].append(source_record)
    assignments: list[TaxonomyAssignment] = []
    for row_number, row in enumerate(rows, 1):
        source_identifier = str(row.get("source_identifier", "")).strip()
        taxonomy_source = str(row.get("taxonomy_source", "")).strip()
        method = str(row.get("method", row.get("assignment_method", ""))).strip()
        evidence = str(row.get("evidence", "")).strip()
        if source_identifier not in by_identifier:
            raise TaxonomyError(f"Unknown IMG source_identifier on assignment row {row_number}")
        if not taxonomy_source or method not in {
            "updated_reference_cluster",
            "updated_reference_derived",
            "updated_reference_unclassified",
        } or not evidence:
            raise TaxonomyError(
                "Derived IMG assignments require taxonomy_source, an updated-reference method, "
                "and evidence"
            )
        taxonomy = _taxonomy_tuple(row.get("taxonomy", ""))
        compartment = str(row.get("compartment", "")).strip()
        for source_record in by_identifier[source_identifier]:
            assignment_id = _stable_id(
                "TAX_",
                source_record.source_record_id,
                ";".join(taxonomy),
                taxonomy_source,
                compartment,
                method,
                evidence,
            )
            assignments.append(
                TaxonomyAssignment(
                    assignment_id,
                    source_record.source_record_id,
                    source_record.sequence_id,
                    taxonomy,
                    taxonomy_source,
                    taxonomy[0],
                    compartment,
                    method,
                    evidence,
                )
            )
    return tuple(sorted(assignments, key=lambda assignment: assignment.taxonomy_assignment_id))


def add_taxonomy_assignments(
    model: DatabaseModel, extra_assignments: Iterable[TaxonomyAssignment]
) -> DatabaseModel:
    """Return a model with explicit assignments merged and preferred rows recomputed."""

    assignments = tuple(
        sorted(
            (*model.taxonomy_assignments, *tuple(extra_assignments)),
            key=lambda assignment: assignment.taxonomy_assignment_id,
        )
    )
    by_sequence: dict[str, list[TaxonomyAssignment]] = defaultdict(list)
    for assignment in assignments:
        by_sequence[assignment.sequence_id].append(assignment)
    preferred = tuple(
        select_preferred_taxonomy(sequence.sequence_id, by_sequence[sequence.sequence_id])
        for sequence in model.sequences
        if by_sequence[sequence.sequence_id]
    )
    return DatabaseModel(model.sequences, model.source_records, assignments, preferred)



def read_prepared_jsonl(path: Path) -> Iterator[PreparedSourceRecord]:
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
                raw_taxonomy = row.pop("taxonomy", ())
                taxonomy = _taxonomy_tuple(raw_taxonomy) if raw_taxonomy else ()
                yield PreparedSourceRecord(taxonomy=taxonomy, **row)
            except (TypeError, ValueError, json.JSONDecodeError) as error:
                raise BuildError(f"Invalid prepared JSONL row {line_number}: {error}") from error


def read_img_metadata_tsv(path: Path) -> tuple[ImgLocation, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        return clean_img_metadata(csv.DictReader(handle, delimiter="\t"))


def _build_prepared(args: argparse.Namespace) -> None:
    model = build_deduplicated_model(read_prepared_jsonl(args.records_jsonl))
    locations = read_img_metadata_tsv(args.img_metadata) if args.img_metadata else ()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    validate_release(model, locations)
    write_marker_fasta(output / "markers.fasta", model.sequences)
    write_release_tables(output, model, locations)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    fetch = subparsers.add_parser("fetch", help="Fetch pinned sources atomically")
    fetch.add_argument("sources", nargs="+")
    fetch.add_argument("--config", type=Path, default=DEFAULT_SOURCES)
    fetch.add_argument("--output", type=Path, required=True)

    build = subparsers.add_parser(
        "build-prepared", help="Build from explicit prepared JSONL records"
    )
    build.add_argument("--records-jsonl", type=Path, required=True)
    build.add_argument("--img-metadata", type=Path)
    build.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "fetch":
        for name in args.sources:
            fetch_configured_source(name, args.output, config_path=args.config)
    else:
        _build_prepared(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
