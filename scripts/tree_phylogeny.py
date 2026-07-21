#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from taxonomy_utils import common_value, lowest_common_ancestor, taxonomy_path
from tree_schema import REFERENCE_FIELDS, TREE_ASSIGNMENT_FIELDS, TREE_NEIGHBOR_FIELDS

GAP_CHARACTERS = frozenset("-.~_")
VALID_NUCLEOTIDES = frozenset("ACGTRYSWKMBDHVN")


def trim_alignment(
    input_file: str | Path,
    output_file: str | Path,
    qc_file: str | Path,
    *,
    maximum_gap_fraction: float = 0.9,
) -> dict[str, object]:
    if not 0 <= maximum_gap_fraction < 1:
        raise ValueError("maximum_gap_fraction must be in [0, 1)")
    with Path(input_file).open() as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if len(records) < 4:
        raise ValueError("Tree alignment requires one query and at least 3 references")
    if len({record.id for record in records}) != len(records):
        raise ValueError("Tree alignment contains duplicate leaf identifiers")
    if "QUERY" not in {record.id for record in records}:
        raise ValueError("Tree alignment is missing the QUERY leaf")
    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1:
        raise ValueError("cmalign output contains unequal sequence lengths")
    input_columns = lengths.pop()
    if input_columns == 0:
        raise ValueError("cmalign output is empty")

    raw_sequences = [str(record.seq) for record in records]
    insert_columns = {
        index
        for index in range(input_columns)
        if any(
            sequence[index] == "." or sequence[index].islower()
            for sequence in raw_sequences
        )
    }

    normalized: list[str] = []
    for record, raw_sequence in zip(records, raw_sequences, strict=True):
        normalized_characters: list[str] = []
        for character in raw_sequence:
            if character in GAP_CHARACTERS:
                normalized_characters.append("-")
                continue
            nucleotide = character.upper()
            normalized_characters.append("T" if nucleotide == "U" else nucleotide)
        sequence = "".join(normalized_characters)
        invalid = sorted(set(sequence) - VALID_NUCLEOTIDES - {"-"})
        if invalid:
            raise ValueError(
                f"Unexpected aligned nucleotide character(s) for {record.id}: "
                + ", ".join(invalid)
            )
        normalized.append(sequence)

    keep_columns = [
        index
        for index in range(input_columns)
        if index not in insert_columns
        and sum(sequence[index] == "-" for sequence in normalized) / len(normalized)
        <= maximum_gap_fraction
    ]
    if not keep_columns:
        raise ValueError("Gap trimming removed every alignment column")
    trimmed_records = [
        SeqRecord(
            Seq("".join(sequence[index] for index in keep_columns)),
            id=record.id,
            description="",
        )
        for record, sequence in zip(records, normalized, strict=True)
    ]
    query = next(record for record in trimmed_records if record.id == "QUERY")
    query_sites = sum(character != "-" for character in str(query.seq))
    if query_sites == 0:
        raise ValueError("The query has no residues after alignment trimming")
    with Path(output_file).open("w") as handle:
        SeqIO.write(trimmed_records, handle, "fasta")
    qc = {
        "schema_version": 1,
        "sequence_count": len(records),
        "input_columns": input_columns,
        "retained_columns": len(keep_columns),
        "removed_columns": input_columns - len(keep_columns),
        "removed_insert_columns": len(insert_columns),
        "removed_high_gap_columns": (
            input_columns - len(insert_columns) - len(keep_columns)
        ),
        "maximum_gap_fraction": maximum_gap_fraction,
        "query_residue_columns": query_sites,
    }
    Path(qc_file).write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")
    return qc


def _read_reference_rows(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != REFERENCE_FIELDS:
            raise ValueError(f"Unexpected tree-reference columns: {reader.fieldnames}")
        rows = list(reader)
    if len({row["leaf_id"] for row in rows}) != len(rows):
        raise ValueError("Tree-reference table contains duplicate leaf identifiers")
    return rows


def _lineage(row: dict[str, str]) -> tuple[tuple[str, ...], str, str]:
    if row["taxonomy_domain"] == "ambiguous":
        return (), "", "ambiguous"
    candidates = [
        (row["taxonomy"], row["taxonomy_source"], "preferred_taxonomy"),
        (
            row["centroid_taxonomy"],
            row["centroid_taxonomy_source"],
            "centroid_taxonomy",
        ),
    ]
    parsed: list[tuple[tuple[str, ...], str, str]] = []
    for value, source, basis in candidates:
        if not value or value == "Unclassified":
            continue
        try:
            parsed.append((taxonomy_path(value), source, basis))
        except ValueError:
            continue
    if parsed:
        return max(parsed, key=lambda item: (len(item[0]), item[2]))
    domain = row["taxonomy_domain"]
    if domain and domain not in {"ambiguous", "Unclassified"}:
        return (domain,), row["taxonomy_source"], "domain"
    return (), "", "unclassified"


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".10g")


def classify_tree(
    *,
    tree_file: str | Path,
    references_file: str | Path,
    task_file: str | Path,
    assignment_output: str | Path,
    neighbors_output: str | Path,
    assignment_neighbors: int = 5,
    inference_model: str = "GTR+F+R4",
) -> dict[str, str]:
    if assignment_neighbors < 1:
        raise ValueError("assignment_neighbors must be positive")
    from ete4 import Tree

    task = json.loads(Path(task_file).read_text())
    if task.get("schema_version") != 1:
        raise ValueError("Unsupported tree-task schema")
    reference_rows = _read_reference_rows(references_file)
    by_leaf = {row["leaf_id"]: row for row in reference_rows}
    tree = Tree(Path(tree_file).read_text().strip())
    leaf_names = set(tree.leaf_names())
    expected = {"QUERY", *by_leaf}
    if leaf_names != expected:
        raise ValueError(
            "IQ-TREE leaves differ from the prepared alignment: "
            f"missing={sorted(expected - leaf_names)[:5]}, "
            f"extra={sorted(leaf_names - expected)[:5]}"
        )
    query = tree["QUERY"]
    neighbors: list[tuple[float, int, str, dict[str, str], tuple[str, ...], str, str]] = []
    for leaf_id, row in by_leaf.items():
        lineage, lineage_source, lineage_basis = _lineage(row)
        neighbors.append(
            (
                float(tree.get_distance(query, tree[leaf_id])),
                int(row["hit_rank"]),
                leaf_id,
                row,
                lineage,
                lineage_source,
                lineage_basis,
            )
        )
    neighbors.sort(key=lambda item: (item[0], item[1], item[2]))
    informative = [neighbor for neighbor in neighbors if neighbor[4]]
    basis = informative[:assignment_neighbors]
    if basis:
        boundary_distance = basis[-1][0]
        basis = [
            neighbor
            for neighbor in informative
            if neighbor[0] < boundary_distance
            or math.isclose(
                neighbor[0],
                boundary_distance,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        ]
    lineage_lca = lowest_common_ancestor(neighbor[4] for neighbor in basis)
    taxonomy = ";".join(lineage_lca) if lineage_lca else "Unclassified"
    domain = lineage_lca[0] if lineage_lca else "Unclassified"
    taxonomy_sources = sorted(
        {neighbor[5] for neighbor in basis if neighbor[5]}
    )
    taxonomy_source = "|".join(taxonomy_sources)
    compartment = common_value(
        (neighbor[3]["compartment"] for neighbor in basis), conflict="mixed"
    )
    nearest = neighbors[0]
    query_parent = query.parent
    parent_support = None if query_parent is None else query_parent.support
    edge_support = (
        None
        if query_parent is None or query_parent.is_root or parent_support is None
        else float(parent_support)
    )
    assignment = {
        "name": str(task["name"]),
        "sample": str(task["sample"]),
        "model": str(task["detected_model"]),
        "tree_model": str(task["tree_model"]),
        "tree_marker": str(task["tree_marker"]),
        "tree_route_decision": str(task["tree_route_decision"]),
        "tree_route_16s_votes": str(task["tree_route_16s_votes"]),
        "tree_route_18s_votes": str(task["tree_route_18s_votes"]),
        "tree_route_16s_best_bitscore": _format_number(
            task.get("tree_route_16s_best_bitscore")
        ),
        "tree_route_18s_best_bitscore": _format_number(
            task.get("tree_route_18s_best_bitscore")
        ),
        "tree_taxonomy": taxonomy,
        "tree_taxonomy_source": taxonomy_source,
        "tree_taxonomy_domain": domain,
        "tree_compartment": compartment,
        "tree_assignment_method": (
            "tree_nearest_named_lca"
            if basis
            else "tree_unclassified_no_named_neighbors"
        ),
        "tree_basis_neighbors": str(len(basis)),
        "tree_nearest_sseqid": nearest[3]["blast_sseqid"],
        "tree_nearest_reference_identifiers": nearest[3]["reference_identifiers"],
        "tree_nearest_distance": _format_number(nearest[0]),
        "tree_query_edge_support": _format_number(edge_support),
        "tree_inference_model": inference_model,
    }
    with Path(assignment_output).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TREE_ASSIGNMENT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(assignment)

    basis_leaf_ids = {neighbor[2] for neighbor in basis}
    with Path(neighbors_output).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TREE_NEIGHBOR_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for rank, neighbor in enumerate(neighbors, start=1):
            distance, _hit_rank, leaf_id, row, lineage, source, lineage_basis = neighbor
            writer.writerow(
                {
                    "name": task["name"],
                    "sample": task["sample"],
                    "model": task["detected_model"],
                    "tree_model": task["tree_model"],
                    "tree_marker": task["tree_marker"],
                    "tree_neighbor_rank": rank,
                    "tree_distance": _format_number(distance),
                    "used_for_assignment": str(leaf_id in basis_leaf_ids).lower(),
                    "tree_lineage": ";".join(lineage),
                    "tree_lineage_source": source,
                    "tree_lineage_basis": lineage_basis,
                    **row,
                }
            )
    return assignment


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Trim SSU covariance-model alignments and classify tree neighbors."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    trim = subparsers.add_parser("trim")
    trim.add_argument("--input", required=True)
    trim.add_argument("--output", required=True)
    trim.add_argument("--qc", required=True)
    trim.add_argument("--maximum-gap-fraction", type=float, default=0.9)

    classify = subparsers.add_parser("classify")
    classify.add_argument("--tree", required=True)
    classify.add_argument("--references", required=True)
    classify.add_argument("--task", required=True)
    classify.add_argument("--assignment-output", required=True)
    classify.add_argument("--neighbors-output", required=True)
    classify.add_argument("--assignment-neighbors", type=int, default=5)
    classify.add_argument("--inference-model", default="GTR+F+R4")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "trim":
        trim_alignment(
            args.input,
            args.output,
            args.qc,
            maximum_gap_fraction=args.maximum_gap_fraction,
        )
    else:
        classify_tree(
            tree_file=args.tree,
            references_file=args.references,
            task_file=args.task,
            assignment_output=args.assignment_output,
            neighbors_output=args.neighbors_output,
            assignment_neighbors=args.assignment_neighbors,
            inference_model=args.inference_model,
        )


if __name__ == "__main__":
    main()
