#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from annotate_hits import BlastHit, TaxonomyRecord, load_blast_hits, load_taxonomy_records
from top_hit_reporting import (
    ReferenceRecord,
    load_query_sequences,
    load_reference_records,
    reference_record,
)
from tree_schema import REFERENCE_FIELDS, TREE_ASSIGNMENT_FIELDS


MARKERS = ("16S", "18S")


def _hit_key(hit: BlastHit) -> tuple[float, float, int, str]:
    return (-hit.bit_score, -hit.percent_identity, -hit.alignment_length, hit.subject)


def _parse_assignments(values: list[str], label: str) -> dict[str, str]:
    assignments: dict[str, str] = {}
    for value in values:
        key, separator, item = value.partition("=")
        if not separator or key not in MARKERS or not item:
            raise ValueError(f"Invalid {label} assignment: {value!r}")
        if key in assignments:
            raise ValueError(f"Duplicate {label} assignment for {key}")
        assignments[key] = item
    missing = sorted(set(MARKERS) - assignments.keys())
    if missing:
        raise ValueError(f"Missing {label} assignment(s): {', '.join(missing)}")
    return assignments


def _global_route_hits(
    hits_by_marker: dict[str, list[BlastHit]],
    detected_marker: str,
    route_hits: int,
) -> list[tuple[str, BlastHit]]:
    candidates = [
        (marker, hit)
        for marker in MARKERS
        for hit in hits_by_marker.get(marker, [])
    ]
    candidates.sort(
        key=lambda item: (
            *_hit_key(item[1]),
            item[0] != detected_marker,
            item[0],
        )
    )
    unique: dict[str, tuple[str, BlastHit]] = {}
    for marker, hit in candidates:
        unique.setdefault(hit.subject, (marker, hit))
    return list(unique.values())[:route_hits]


def choose_marker(
    hits_by_marker: dict[str, list[BlastHit]],
    detected_marker: str,
    route_hits: int,
) -> tuple[str, str, dict[str, int], dict[str, float]]:
    if detected_marker not in MARKERS:
        raise ValueError(f"Unsupported detected marker for tree routing: {detected_marker}")
    routed = _global_route_hits(hits_by_marker, detected_marker, route_hits)
    votes = {
        marker: sum(observed == marker for observed, _hit in routed)
        for marker in MARKERS
    }
    best_scores = {
        marker: (
            hits_by_marker[marker][0].bit_score
            if hits_by_marker.get(marker)
            else float("-inf")
        )
        for marker in MARKERS
    }
    if not routed:
        return detected_marker, "detected_model_no_blast_hits", votes, best_scores
    if votes["16S"] != votes["18S"]:
        selected = max(MARKERS, key=lambda marker: votes[marker])
        return selected, "majority_global_top_hits", votes, best_scores
    if best_scores["16S"] != best_scores["18S"]:
        selected = max(MARKERS, key=lambda marker: best_scores[marker])
        return selected, "best_bitscore_tiebreak", votes, best_scores
    return detected_marker, "detected_model_tiebreak", votes, best_scores


def _task_key(sample: str, model: str, query: str) -> str:
    digest = hashlib.sha256(f"{sample}\0{model}\0{query}".encode()).hexdigest()
    return f"q_{digest[:16]}"


def _taxonomy_values(record: TaxonomyRecord | None) -> dict[str, str]:
    if record is None:
        return {
            "reference_source": "",
            "taxonomy": "",
            "taxonomy_source": "",
            "taxonomy_domain": "",
            "compartment": "",
            "taxonomy_assignment_method": "",
            "centroid_names": "",
            "centroid_taxonomy": "",
            "centroid_taxonomy_source": "",
        }
    return {
        "reference_source": record.reference_source,
        "taxonomy": record.taxonomy,
        "taxonomy_source": record.taxonomy_source,
        "taxonomy_domain": record.domain,
        "compartment": record.compartment,
        "taxonomy_assignment_method": record.assignment_method,
        "centroid_names": record.centroid_names,
        "centroid_taxonomy": record.centroid_taxonomy,
        "centroid_taxonomy_source": record.centroid_taxonomy_source,
    }


def _reference_row(
    leaf_id: str,
    rank: int,
    hit: BlastHit,
    taxonomy: TaxonomyRecord | None,
    reference: ReferenceRecord,
) -> dict[str, object]:
    return {
        "leaf_id": leaf_id,
        "blast_sseqid": hit.subject,
        "hit_rank": rank,
        "reference_identifiers": reference.identifiers,
        "reference_versions": reference.versions,
        **_taxonomy_values(taxonomy),
        "blast_pident": format(hit.percent_identity, "g"),
        "blast_length": hit.alignment_length,
        "blast_evalue": format(hit.evalue, "g"),
        "blast_bitscore": format(hit.bit_score, "g"),
    }


def _format_score(value: float) -> str:
    return "" if value == float("-inf") else format(value, ".10g")


def _skipped_assignment(
    *,
    query: str,
    sample: str,
    detected_model: str,
    selected_marker: str,
    selected_model: str,
    decision: str,
    votes: dict[str, int],
    best_scores: dict[str, float],
) -> dict[str, str]:
    row = dict.fromkeys(TREE_ASSIGNMENT_FIELDS, "")
    row.update(
        {
            "name": query,
            "sample": sample,
            "model": detected_model,
            "tree_model": selected_model,
            "tree_marker": selected_marker,
            "tree_route_decision": decision,
            "tree_route_16s_votes": str(votes["16S"]),
            "tree_route_18s_votes": str(votes["18S"]),
            "tree_route_16s_best_bitscore": _format_score(best_scores["16S"]),
            "tree_route_18s_best_bitscore": _format_score(best_scores["18S"]),
            "tree_assignment_method": "tree_skipped_insufficient_references",
            "tree_basis_neighbors": "0",
        }
    )
    return row


def prepare_tree_tasks(
    *,
    query_fasta: str | Path,
    blast_files: dict[str, str | Path],
    taxonomy_file: str | Path,
    source_records_file: str | Path,
    sample: str,
    detected_model: str,
    detected_marker: str,
    marker_models: dict[str, str],
    output_directory: str | Path,
    skipped_assignments_file: str | Path,
    reference_count: int = 100,
    route_hits: int = 100,
) -> list[Path]:
    if reference_count < 3:
        raise ValueError("reference_count must be at least 3")
    if route_hits < 1:
        raise ValueError("route_hits must be positive")
    query_sequences = load_query_sequences(query_fasta)
    hits_by_marker = {
        marker: load_blast_hits(blast_files[marker]) for marker in MARKERS
    }
    observed_queries = {
        query
        for marker_hits in hits_by_marker.values()
        for query in marker_hits
    }
    unknown_queries = sorted(observed_queries - query_sequences.keys())
    if unknown_queries:
        raise ValueError(
            "Tree-routing BLAST output contains unknown query sequence(s): "
            + ", ".join(unknown_queries)
        )
    subjects = {
        hit.subject
        for marker_hits in hits_by_marker.values()
        for hits in marker_hits.values()
        for hit in hits
    }
    taxonomy_records = load_taxonomy_records(taxonomy_file, subjects)
    reference_records = load_reference_records(source_records_file, subjects)
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    task_directories: list[Path] = []
    skipped_assignments: list[dict[str, str]] = []

    for query, sequence in query_sequences.items():
        query_hits = {
            marker: hits_by_marker[marker].get(query, []) for marker in MARKERS
        }
        selected_marker, decision, votes, best_scores = choose_marker(
            query_hits, detected_marker, route_hits
        )
        selected_hits = query_hits[selected_marker][:reference_count]
        if len(selected_hits) < 3:
            skipped_assignments.append(
                _skipped_assignment(
                    query=query,
                    sample=sample,
                    detected_model=detected_model,
                    selected_marker=selected_marker,
                    selected_model=marker_models[selected_marker],
                    decision=decision,
                    votes=votes,
                    best_scores=best_scores,
                )
            )
            continue
        key = _task_key(sample, detected_model, query)
        task_directory = output / key
        task_directory.mkdir()
        query_record = SeqRecord(Seq(sequence), id=query, description="")
        with (task_directory / "query.fna").open("w") as handle:
            SeqIO.write([query_record], handle, "fasta")

        rows = [
            _reference_row(
                f"REF{rank:04d}",
                rank,
                hit,
                taxonomy_records.get(hit.subject),
                reference_record(hit.subject, reference_records),
            )
            for rank, hit in enumerate(selected_hits, start=1)
        ]
        with (task_directory / "references.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=REFERENCE_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        (task_directory / "reference_ids.txt").write_text(
            "".join(f"{row['blast_sseqid']}\n" for row in rows)
        )
        payload = {
            "schema_version": 1,
            "query_key": key,
            "name": query,
            "sample": sample,
            "detected_model": detected_model,
            "detected_marker": detected_marker,
            "tree_model": marker_models[selected_marker],
            "tree_marker": selected_marker,
            "tree_route_decision": decision,
            "tree_route_hits": route_hits,
            "tree_route_16s_votes": votes["16S"],
            "tree_route_18s_votes": votes["18S"],
            "tree_route_16s_best_bitscore": (
                None if best_scores["16S"] == float("-inf") else best_scores["16S"]
            ),
            "tree_route_18s_best_bitscore": (
                None if best_scores["18S"] == float("-inf") else best_scores["18S"]
            ),
            "tree_reference_count": len(rows),
        }
        (task_directory / "task.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )
        task_directories.append(task_directory)

    with Path(skipped_assignments_file).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TREE_ASSIGNMENT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(skipped_assignments)
    return task_directories


def build_alignment_input(
    task_directory: str | Path,
    reference_fasta: str | Path,
    output_file: str | Path,
) -> None:
    task = Path(task_directory)
    with (task / "references.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != REFERENCE_FIELDS:
            raise ValueError("Unexpected tree-reference table columns")
        rows = list(reader)
    fetched: dict[str, SeqRecord] = {}
    with Path(reference_fasta).open() as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in fetched:
                raise ValueError(f"Duplicate fetched reference: {record.id}")
            fetched[record.id] = record
    expected = {row["blast_sseqid"] for row in rows}
    if set(fetched) != expected:
        missing = sorted(expected - fetched.keys())
        extra = sorted(fetched.keys() - expected)
        raise ValueError(
            f"Fetched tree references differ from selection; missing={missing[:5]}, "
            f"extra={extra[:5]}"
        )
    with (task / "query.fna").open() as handle:
        queries = list(SeqIO.parse(handle, "fasta"))
    if len(queries) != 1:
        raise ValueError("Tree task must contain exactly one query sequence")
    records = [SeqRecord(queries[0].seq, id="QUERY", description="")]
    records.extend(
        SeqRecord(
            fetched[row["blast_sseqid"]].seq,
            id=row["leaf_id"],
            description="",
        )
        for row in rows
    )
    with Path(output_file).open("w") as handle:
        SeqIO.write(records, handle, "fasta")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select tree-mode BLAST references and prepare cmalign input."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--query-fasta", required=True)
    prepare.add_argument("--blast", action="append", required=True)
    prepare.add_argument("--taxonomy-db", required=True)
    prepare.add_argument("--source-records-db", required=True)
    prepare.add_argument("--sample", required=True)
    prepare.add_argument("--detected-model", required=True)
    prepare.add_argument("--detected-marker", required=True, choices=MARKERS)
    prepare.add_argument("--marker-model", action="append", required=True)
    prepare.add_argument("--reference-count", type=int, default=100)
    prepare.add_argument("--route-hits", type=int, default=100)
    prepare.add_argument("--output-directory", required=True)
    prepare.add_argument("--skipped-assignments-output", required=True)

    alignment = subparsers.add_parser("alignment-input")
    alignment.add_argument("--task-directory", required=True)
    alignment.add_argument("--reference-fasta", required=True)
    alignment.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "prepare":
        prepare_tree_tasks(
            query_fasta=args.query_fasta,
            blast_files=_parse_assignments(args.blast, "BLAST file"),
            taxonomy_file=args.taxonomy_db,
            source_records_file=args.source_records_db,
            sample=args.sample,
            detected_model=args.detected_model,
            detected_marker=args.detected_marker,
            marker_models=_parse_assignments(args.marker_model, "marker-model"),
            output_directory=args.output_directory,
            skipped_assignments_file=args.skipped_assignments_output,
            reference_count=args.reference_count,
            route_hits=args.route_hits,
        )
    else:
        build_alignment_input(
            args.task_directory,
            args.reference_fasta,
            args.output,
        )


if __name__ == "__main__":
    main()
