#!/usr/bin/env python3
"""Deterministic, provenance-preserving SSU reference database builder.

Source acquisition and source-specific parsing live in :mod:`database_sources`.
Release validation and serialization live in :mod:`database_release_io`. This
module owns exact-sequence deduplication, taxonomy selection, and CLI orchestration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from database_contracts import (
    PR2_COMPARTMENTS,
    PROKARYOTIC_DOMAINS,
    DatabaseModel,
    PreferredTaxonomy,
    PreparedSourceRecord,
    ReleaseValidationError,
    SequenceRecord,
    SourceRecord,
    TaxonomyAssignment,
    TaxonomyError,
)
from database_sources import (
    DEFAULT_SOURCES,
    fetch_configured_source,
    normalize_sequence_for_hashing,
    sequence_identifier,
)
from taxonomy_utils import (
    common_value,
    lowest_common_ancestor as _shared_lca,
    taxonomy_path,
)


def _stable_id(prefix: str, *values: str) -> str:
    payload = "\0".join(values).encode("utf-8")
    return prefix + hashlib.sha256(payload).hexdigest()


def _taxonomy_tuple(value: Sequence[str] | str) -> tuple[str, ...]:
    try:
        return taxonomy_path(value)
    except ValueError as error:
        raise TaxonomyError(str(error)) from error


def _validate_assignment_taxonomy(
    value: Sequence[str] | str, method: str
) -> tuple[str, ...]:
    taxonomy = _taxonomy_tuple(value)
    if taxonomy[0] == "Unclassified" and method != "updated_reference_unclassified":
        raise TaxonomyError(
            "Unclassified taxonomy requires an unclassified assignment method: "
            f"{method!r}"
        )
    return taxonomy


def _reference_source(assignments: Iterable[TaxonomyAssignment]) -> str:
    return "+".join(
        sorted(
            {
                "IMG"
                if assignment.assignment_method.lower() != "native"
                else assignment.taxonomy_source
                for assignment in assignments
            }
        )
    )


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


def _centroid_summary(
    assignments: Iterable[TaxonomyAssignment],
) -> tuple[str, tuple[str, ...], str]:
    evidence = [assignment for assignment in assignments if assignment.centroid_name]
    if not evidence:
        return "", (), ""
    names = sorted({assignment.centroid_name for assignment in evidence})
    taxonomies = [
        assignment.centroid_taxonomy
        for assignment in evidence
        if assignment.centroid_taxonomy
    ]
    taxonomy = _shared_lca(taxonomies) if taxonomies else ()
    sources = sorted(
        {
            source
            for assignment in evidence
            for source in assignment.centroid_taxonomy_source.split("+")
            if source
        }
    )
    return (
        json.dumps(names, ensure_ascii=True, separators=(",", ":")),
        taxonomy,
        "+".join(sources) if taxonomy else "",
    )


def select_preferred_taxonomy(
    sequence_id: str, assignments: Iterable[TaxonomyAssignment]
) -> PreferredTaxonomy:
    candidates = list(assignments)
    if not candidates:
        raise TaxonomyError(f"No taxonomy assignments for {sequence_id}")
    for candidate in candidates:
        _validate_assignment_taxonomy(candidate.taxonomy, candidate.assignment_method)
    # Apply the source policy before testing for a domain conflict. Otherwise a
    # lower-priority, cluster-derived IMG assignment could turn an exact native
    # SILVA or PR2 match into a false cross-domain ambiguity.
    best_priority = min(taxonomy_priority(assignment) for assignment in candidates)
    best = [assignment for assignment in candidates if taxonomy_priority(assignment) == best_priority]
    centroid_names, centroid_taxonomy, centroid_taxonomy_source = _centroid_summary(best)
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
            reference_source=_reference_source(best),
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
            centroid_names=centroid_names,
            centroid_taxonomy=centroid_taxonomy,
            centroid_taxonomy_source=centroid_taxonomy_source,
        )
    taxonomy = (
        best[0].taxonomy
        if len(best) == 1
        else lowest_common_ancestor(assignment.taxonomy for assignment in best)
    )
    if not taxonomy:
        raise TaxonomyError(f"No common taxonomy domain for {sequence_id}")
    disagreement = len({assignment.taxonomy for assignment in best}) > 1
    taxonomy_sources = sorted({assignment.taxonomy_source for assignment in best})
    return PreferredTaxonomy(
        sequence_id=sequence_id,
        reference_source=_reference_source(best),
        taxonomy=taxonomy,
        taxonomy_source="+".join(taxonomy_sources),
        domain=taxonomy[0],
        compartment=common_value(
            (assignment.compartment for assignment in best), conflict="mixed"
        ),
        assignment_method="lowest_common_ancestor"
        if disagreement
        else best[0].assignment_method,
        cross_domain_conflict=cross_domain_conflict,
        centroid_names=centroid_names,
        centroid_taxonomy=centroid_taxonomy,
        centroid_taxonomy_source=centroid_taxonomy_source,
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
            taxonomy = _validate_assignment_taxonomy(
                record.taxonomy, record.assignment_method
            )
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

    by_identifier: dict[tuple[str, str], list[SourceRecord]] = defaultdict(list)
    for source_record in source_records:
        if source_record.reference_source.upper() == "IMG":
            by_identifier[
                (source_record.marker, source_record.source_identifier)
            ].append(source_record)
    assignments: list[TaxonomyAssignment] = []
    for row_number, row in enumerate(rows, 1):
        source_identifier = str(row.get("source_identifier", "")).strip()
        marker = str(row.get("marker", "")).strip()
        taxonomy_source = str(row.get("taxonomy_source", "")).strip()
        method = str(row.get("assignment_method", "")).strip()
        evidence = str(row.get("evidence", "")).strip()
        centroid_name = str(row.get("centroid_name", "")).strip()
        centroid_taxonomy_text = str(row.get("centroid_taxonomy", "")).strip()
        centroid_taxonomy_source = str(
            row.get("centroid_taxonomy_source", "")
        ).strip()
        if (marker, source_identifier) not in by_identifier:
            raise TaxonomyError(
                f"Unknown IMG marker/source_identifier on assignment row {row_number}: "
                f"{marker!r}/{source_identifier!r}"
            )
        if not taxonomy_source or method not in {
            "updated_reference_cluster",
            "updated_reference_derived",
            "updated_reference_unclassified",
        } or not evidence:
            raise TaxonomyError(
                "Derived IMG assignments require taxonomy_source, an updated-reference method, "
                "and evidence"
            )
        taxonomy = _validate_assignment_taxonomy(row.get("taxonomy", ""), method)
        centroid_taxonomy = (
            _taxonomy_tuple(centroid_taxonomy_text)
            if centroid_taxonomy_text
            else ()
        )
        if any(character in centroid_name for character in "\r\n\t|"):
            raise TaxonomyError(
                f"Derived IMG assignment has an invalid centroid name on row {row_number}"
            )
        if method == "updated_reference_cluster":
            if not centroid_name or not centroid_taxonomy or not centroid_taxonomy_source:
                raise TaxonomyError(
                    "Cluster-derived IMG assignments require centroid name, taxonomy, "
                    "and taxonomy source"
                )
            if centroid_taxonomy[: len(taxonomy)] != taxonomy:
                raise TaxonomyError(
                    "Cluster-derived IMG member taxonomy must be a prefix of the "
                    "centroid taxonomy"
                )
        elif method == "updated_reference_unclassified" and not centroid_name:
            raise TaxonomyError(
                "Unclassified IMG cluster assignments require a centroid name"
            )
        compartment = str(row.get("compartment", "")).strip()
        for source_record in by_identifier[(marker, source_identifier)]:
            assignment_id = _stable_id(
                "TAX_",
                source_record.source_record_id,
                ";".join(taxonomy),
                taxonomy_source,
                compartment,
                method,
                evidence,
                centroid_name,
                ";".join(centroid_taxonomy),
                centroid_taxonomy_source,
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
                    centroid_name,
                    centroid_taxonomy,
                    centroid_taxonomy_source,
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



def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    fetch = subparsers.add_parser("fetch", help="Fetch pinned sources atomically")
    fetch.add_argument("sources", nargs="+")
    fetch.add_argument("--config", type=Path, default=DEFAULT_SOURCES)
    fetch.add_argument("--output", type=Path, required=True)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "fetch":
        for name in args.sources:
            fetch_configured_source(name, args.output, config_path=args.config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
