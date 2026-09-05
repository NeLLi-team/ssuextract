#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import replace
from pathlib import Path

import duckdb

from hit_processing import HIT_FIELDS
from runtime_taxonomy import (
    BlastHit,
    RuntimeTaxonomyDecision,
    TaxonomyRecord,
    load_blast_hits,
    load_runtime_calibration,
    query_coverage,
    resolve_runtime_taxonomy,
)
from tree_schema import SUMMARY_TREE_FIELDS
from top_hit_reporting import (
    ReferenceRecord,
    load_query_sequences,
    load_reference_records,
    reference_record,
    write_top_hits,
)


SUMMARY_FIELDS = [
    "name",
    "sample",
    "model",
    "length",
    "coordinates",
    "strand",
    "sequence_type",
    "contig_name",
    "blast_sseqid",
    "blast_pident",
    "blast_length",
    "blast_bitscore",
    "is_assembled",
    "reference_source",
    "taxonomy",
    "taxonomy_source",
    "taxonomy_domain",
    "compartment",
    "taxonomy_assignment_method",
    "taxonomy_alternatives",
    "blast_tied_subjects",
    "blast_ties_truncated",
    "centroid_names",
    "centroid_taxonomy",
    "centroid_taxonomy_source",
    "reference_identifiers",
    "reference_versions",
    "query_sequence",
    *SUMMARY_TREE_FIELDS,
    "hit_length",
    "model_coverage",
    "fragment_count",
    "component_coordinates",
    "blast_candidate_taxonomy",
    "reference_taxonomy",
    "reference_taxonomy_source",
    "reference_taxonomy_assignment_method",
    "blast_query_coverage",
    "assignment_candidate_count",
    "assignment_unknown_count",
    "assignment_candidates_truncated",
]


def load_taxonomy_records(
    taxonomy_file: str | Path,
    subjects: set[str],
) -> dict[str, TaxonomyRecord]:
    if not subjects:
        return {}
    placeholders = ",".join("?" for _ in subjects)
    parameters = [str(Path(taxonomy_file)), *sorted(subjects)]
    connection = duckdb.connect(":memory:")
    try:
        connection.execute("SET threads = 1")
        columns = {
            str(row[0])
            for row in connection.execute(
                "DESCRIBE SELECT * FROM read_parquet(?)",
                [str(Path(taxonomy_file))],
            ).fetchall()
        }
        centroid_columns = {
            "centroid_names",
            "centroid_taxonomy",
            "centroid_taxonomy_source",
        }
        present_centroid_columns = centroid_columns.intersection(columns)
        if present_centroid_columns and present_centroid_columns != centroid_columns:
            raise ValueError("Preferred taxonomy has an incomplete centroid schema")
        centroid_fields = (
            ("centroid_names", "centroid_taxonomy", "centroid_taxonomy_source")
            if present_centroid_columns
            else (
                "'' AS centroid_names",
                "'' AS centroid_taxonomy",
                "'' AS centroid_taxonomy_source",
            )
        )
        query = f"""
            SELECT
                sequence_id,
                reference_source,
                taxonomy,
                taxonomy_source,
                domain,
                compartment,
                assignment_method,
                cross_domain_conflict,
                taxonomy_alternatives,
                {centroid_fields[0]},
                {centroid_fields[1]},
                {centroid_fields[2]}
            FROM read_parquet(?)
            WHERE sequence_id IN ({placeholders})
        """
        rows = connection.execute(query, parameters).fetchall()
    finally:
        connection.close()
    records: dict[str, TaxonomyRecord] = {}
    for row in rows:
        sequence_id = str(row[0])
        if sequence_id in records:
            raise ValueError(f"Duplicate preferred taxonomy row for {sequence_id}")
        records[sequence_id] = TaxonomyRecord(
            reference_source=str(row[1] or ""),
            taxonomy=str(row[2] or ""),
            taxonomy_source=str(row[3] or ""),
            domain=str(row[4] or ""),
            compartment=str(row[5] or ""),
            assignment_method=str(row[6] or ""),
            cross_domain_conflict=bool(row[7]),
            taxonomy_alternatives=str(row[8] or ""),
            centroid_names=str(row[9] or ""),
            centroid_taxonomy=str(row[10] or ""),
            centroid_taxonomy_source=str(row[11] or ""),
        )
        if records[sequence_id].centroid_names:
            try:
                names = json.loads(records[sequence_id].centroid_names)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"Invalid centroid_names JSON for {sequence_id}"
                ) from error
            if not isinstance(names, list) or not all(
                isinstance(name, str) and name and "|" not in name for name in names
            ):
                raise ValueError(f"Invalid centroid_names JSON for {sequence_id}")
            records[sequence_id] = replace(
                records[sequence_id], centroid_names="|".join(names)
            )
    missing = sorted(subjects - records.keys())
    if missing:
        preview = ", ".join(missing[:5])
        suffix = "" if len(missing) <= 5 else f" and {len(missing) - 5} more"
        raise ValueError(
            f"Missing taxonomy metadata for {len(missing)} BLAST subject(s): "
            f"{preview}{suffix}"
        )
    return records


def annotate_hits(
    hits_file: str | Path,
    m8_file: str | Path,
    output_file: str | Path,
    taxonomy_file: str | Path | None = None,
    max_targets: int = 500,
    *,
    query_fasta: str | Path | None = None,
    source_records_file: str | Path | None = None,
    top_hits_output: str | Path | None = None,
    top_hits: int = 5,
    marker: str = "16S",
    runtime_calibration: str | Path | None = None,
    reference_digest: str = "",
) -> None:
    if max_targets < 1:
        raise ValueError("max_targets must be positive")
    if top_hits < 1:
        raise ValueError("top_hits must be positive")
    blast_hits = load_blast_hits(m8_file)
    taxonomy_records: dict[str, TaxonomyRecord] = {}
    reference_records: dict[str, ReferenceRecord] = {}
    subjects = {
        hit.subject
        for query_hits in blast_hits.values()
        for hit in query_hits
    }
    if taxonomy_file is not None:
        taxonomy_records = load_taxonomy_records(taxonomy_file, subjects)
    if source_records_file is not None:
        reference_records = load_reference_records(source_records_file, subjects)
    query_sequences = load_query_sequences(query_fasta) if query_fasta else {}
    calibration = load_runtime_calibration(runtime_calibration, reference_digest)

    with Path(hits_file).open(newline="") as hits_handle:
        reader = csv.DictReader(hits_handle, delimiter="\t")
        if reader.fieldnames != HIT_FIELDS:
            raise ValueError(
                f"Unexpected hit-table columns in {hits_file}: {reader.fieldnames}"
            )
        hit_rows = list(reader)

    missing_sequences = sorted(
        row["name"]
        for row in hit_rows
        if query_fasta and row["name"] not in query_sequences
    )
    if missing_sequences:
        raise ValueError(
            "Extracted FASTA is missing hit sequence(s): " + ", ".join(missing_sequences)
        )

    decisions: dict[str, RuntimeTaxonomyDecision] = {}
    for row in hit_rows:
        query = row["name"]
        query_length = (
            len(query_sequences[query])
            if query in query_sequences
            else int(row["length"])
        )
        decisions[query] = resolve_runtime_taxonomy(
            blast_hits.get(query, ()),
            taxonomy_records,
            query_length=query_length,
            marker=marker,
            calibration=calibration,
            max_targets=max_targets,
        )

    if top_hits_output is not None:
        write_top_hits(
            top_hits_output,
            hit_rows,
            blast_hits,
            taxonomy_records,
            reference_records,
            query_sequences,
            top_hits,
            {
                query: {hit.subject for hit in decision.candidate_hits}
                for query, decision in decisions.items()
            },
        )

    with Path(output_file).open("w", newline="") as output_handle:
        writer = csv.DictWriter(
            output_handle,
            fieldnames=SUMMARY_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in hit_rows:
            query_hits = blast_hits.get(row["name"], [])
            decision = decisions[row["name"]]
            blast_hit = (
                decision.candidate_hits[0]
                if decision.candidate_hits
                else (query_hits[0] if query_hits else None)
            )
            selected_taxonomy = (
                taxonomy_records.get(blast_hit.subject) if blast_hit else None
            )
            blast_reference = (
                reference_record(blast_hit.subject, reference_records)
                if blast_hit
                else ReferenceRecord("", "", ())
            )
            equal_best = [
                hit
                for hit in query_hits
                if query_hits and hit.bit_score == query_hits[0].bit_score
            ]
            writer.writerow(
                {
                    "name": row["name"],
                    "sample": row["sample"],
                    "model": row["model"],
                    "length": row["length"],
                    "coordinates": row["coordinates"],
                    "strand": row["strand"],
                    "sequence_type": row["sequence_type"],
                    "contig_name": row["contig_name"],
                    "blast_sseqid": blast_hit.subject if blast_hit else "",
                    "centroid_names": (
                        selected_taxonomy.centroid_names if selected_taxonomy else ""
                    ),
                    "centroid_taxonomy": (
                        selected_taxonomy.centroid_taxonomy if selected_taxonomy else ""
                    ),
                    "centroid_taxonomy_source": (
                        selected_taxonomy.centroid_taxonomy_source
                        if selected_taxonomy
                        else ""
                    ),
                    "reference_identifiers": blast_reference.identifiers,
                    "reference_versions": blast_reference.versions,
                    "query_sequence": query_sequences.get(row["name"], ""),
                    "taxonomy_mode": "blast",
                    "blast_taxonomy": decision.taxonomy,
                    "blast_taxonomy_source": decision.taxonomy_source,
                    "blast_taxonomy_domain": decision.domain,
                    "blast_compartment": decision.compartment,
                    "blast_taxonomy_assignment_method": decision.assignment_method,
                    "blast_taxonomy_alternatives": decision.taxonomy_alternatives,
                    "blast_pident": (
                        str(round(blast_hit.percent_identity, 2)) if blast_hit else ""
                    ),
                    "blast_length": blast_hit.alignment_length if blast_hit else "",
                    "blast_bitscore": (
                        format(blast_hit.bit_score, "g") if blast_hit else ""
                    ),
                    "is_assembled": row["is_assembled"],
                    "reference_source": (
                        selected_taxonomy.reference_source if selected_taxonomy else ""
                    ),
                    "taxonomy": decision.taxonomy,
                    "taxonomy_source": decision.taxonomy_source,
                    "taxonomy_domain": decision.domain,
                    "compartment": decision.compartment,
                    "taxonomy_assignment_method": decision.assignment_method,
                    "taxonomy_alternatives": decision.taxonomy_alternatives,
                    "blast_tied_subjects": len(equal_best) if equal_best else "",
                    "blast_ties_truncated": (
                        str(len(equal_best) > max_targets).lower()
                        if equal_best
                        else ""
                    ),
                    "hit_length": row.get("hit_length", row["length"]),
                    "model_coverage": row.get("model_coverage", ""),
                    "fragment_count": row.get("fragment_count", ""),
                    "component_coordinates": row.get(
                        "component_coordinates", row["coordinates"]
                    ),
                    "blast_candidate_taxonomy": decision.candidate_taxonomy,
                    "reference_taxonomy": (
                        selected_taxonomy.taxonomy
                        or selected_taxonomy.domain
                        or "Unclassified"
                        if selected_taxonomy
                        else ""
                    ),
                    "reference_taxonomy_source": (
                        selected_taxonomy.taxonomy_source if selected_taxonomy else ""
                    ),
                    "reference_taxonomy_assignment_method": (
                        selected_taxonomy.assignment_method if selected_taxonomy else ""
                    ),
                    "blast_query_coverage": (
                        format(
                            query_coverage(
                                blast_hit,
                                len(query_sequences[row["name"]])
                                if row["name"] in query_sequences
                                else int(row["length"]),
                            ),
                            ".6g",
                        )
                        if blast_hit
                        else ""
                    ),
                    "assignment_candidate_count": len(decision.candidate_hits),
                    "assignment_unknown_count": decision.unknown_candidate_count,
                    "assignment_candidates_truncated": str(
                        decision.candidates_truncated
                    ).lower(),
                }
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add the highest-bit-score BLAST hit to an SSUextract hit table."
    )
    parser.add_argument("--hits", required=True)
    parser.add_argument("--m8", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--query-fasta")
    parser.add_argument("--top-hits-output")
    parser.add_argument(
        "--top-hits",
        type=int,
        default=5,
        help=(
            "Number of overall BLAST hits to report per query; equal-best and "
            "best-reference-source evidence is retained."
        ),
    )
    parser.add_argument(
        "--taxonomy-db",
        help="Preferred-taxonomy Parquet for manifest-driven databases.",
    )
    parser.add_argument(
        "--source-records-db",
        help="Source-record Parquet for public reference identifiers and versions.",
    )
    parser.add_argument(
        "--max-targets",
        type=int,
        default=500,
        help="Policy limit; BLAST must request one additional overflow target.",
    )
    parser.add_argument("--marker", required=True, choices=("16S", "18S"))
    parser.add_argument("--runtime-calibration")
    parser.add_argument("--reference-digest", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    annotate_hits(
        args.hits,
        args.m8,
        args.output,
        args.taxonomy_db,
        args.max_targets,
        query_fasta=args.query_fasta,
        source_records_file=args.source_records_db,
        top_hits_output=args.top_hits_output,
        top_hits=args.top_hits,
        marker=args.marker,
        runtime_calibration=args.runtime_calibration,
        reference_digest=args.reference_digest,
    )


if __name__ == "__main__":
    main()
