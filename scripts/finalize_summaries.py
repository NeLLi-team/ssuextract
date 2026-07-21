#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import glob
import re
from collections import Counter
from pathlib import Path

from annotate_hits import SUMMARY_FIELDS
from hit_processing import META_FIELDS
from top_hit_reporting import TOP_HIT_FIELDS
from tree_schema import TREE_ASSIGNMENT_FIELDS, TREE_NEIGHBOR_FIELDS


CATEGORY_MAPPING = {
    "Bacteria": "BacteriaSSU",
    "Mitochondria": "MitochondriaSSU",
    "Chloroplast": "PlastidSSU",
    "Holosporales": "HolosporalesSSU",
    "Patescibacteria": "PatescibacteriaSSU",
    "Legionellales": "LegionellalesSSU",
    "Rickettsiales": "RickettsialesSSU",
    "Dependentiae": "DependentiaeSSU",
    "Cyanobacteria": "CyanobacteriaSSU",
    "Archaea": "ArchaeaSSU",
    "Eukaryota": "EukaryotaSSU",
}

COMPARTMENT_MAPPING = {
    "mitochondrion": "MitochondriaSSU",
    "plastid": "PlastidSSU",
    "apicoplast": "PlastidSSU",
    "chromatophore": "PlastidSSU",
}


def load_metadata(pattern: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for filename in sorted(glob.glob(pattern)):
        with Path(filename).open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != META_FIELDS:
                raise ValueError(
                    f"Unexpected metadata columns in {filename}: {reader.fieldnames}"
                )
            rows.extend(reader)
    return rows


def load_summary_rows(pattern: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for filename in sorted(glob.glob(pattern)):
        with Path(filename).open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != SUMMARY_FIELDS:
                raise ValueError(
                    f"Unexpected summary columns in {filename}: {reader.fieldnames}"
                )
            rows.extend(reader)
    return sorted(
        rows,
        key=lambda row: (
            row["sample"],
            row["model"],
            row["contig_name"],
            int(row["coordinates"].split("-", maxsplit=1)[0]),
            int(row["coordinates"].split("-", maxsplit=1)[1]),
            row["strand"],
        ),
    )


def write_detailed_summary(
    rows: list[dict[str, str]], output_file: str | Path
) -> None:
    with Path(output_file).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=SUMMARY_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_top_hit_rows(pattern: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for filename in sorted(glob.glob(pattern)):
        with Path(filename).open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != TOP_HIT_FIELDS:
                raise ValueError(
                    f"Unexpected top-hit columns in {filename}: {reader.fieldnames}"
                )
            rows.extend(reader)
    return sorted(
        rows,
        key=lambda row: (
            row["sample"],
            row["model"],
            row["name"],
            int(row["hit_rank"]),
            row["blast_sseqid"],
        ),
    )


def write_top_hit_summary(
    rows: list[dict[str, str]], output_file: str | Path
) -> None:
    with Path(output_file).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TOP_HIT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _load_rows(
    pattern: str,
    expected_fields: list[str],
    label: str,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for filename in sorted(glob.glob(pattern)):
        with Path(filename).open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != expected_fields:
                raise ValueError(
                    f"Unexpected {label} columns in {filename}: {reader.fieldnames}"
                )
            rows.extend(reader)
    return rows


def load_tree_assignment_rows(pattern: str) -> list[dict[str, str]]:
    rows = _load_rows(pattern, TREE_ASSIGNMENT_FIELDS, "tree-assignment")
    keys = [(row["sample"], row["model"], row["name"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate tree assignment for one extracted query")
    return rows


def apply_tree_assignments(
    rows: list[dict[str, str]],
    assignments: list[dict[str, str]],
    taxonomy_mode: str,
) -> list[dict[str, str]]:
    if taxonomy_mode not in {"blast", "tree"}:
        raise ValueError(f"Unsupported taxonomy mode: {taxonomy_mode}")
    if taxonomy_mode == "blast":
        if assignments:
            raise ValueError("Tree assignments were supplied in BLAST taxonomy mode")
        return rows
    by_key = {
        (assignment["sample"], assignment["model"], assignment["name"]): assignment
        for assignment in assignments
    }
    expected = {(row["sample"], row["model"], row["name"]) for row in rows}
    if set(by_key) != expected:
        missing = sorted(expected - by_key.keys())
        extra = sorted(by_key.keys() - expected)
        raise ValueError(
            "Tree assignments do not match extracted queries; "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )
    merged: list[dict[str, str]] = []
    for row in rows:
        assignment = by_key[(row["sample"], row["model"], row["name"])]
        updated = dict(row)
        updated.update(
            {
                field: assignment[field]
                for field in TREE_ASSIGNMENT_FIELDS[3:]
            }
        )
        if assignment["tree_assignment_method"].startswith("tree_skipped_"):
            updated["taxonomy_mode"] = "blast"
        else:
            updated.update(
                {
                    "taxonomy_mode": "tree",
                    "taxonomy": assignment["tree_taxonomy"],
                    "taxonomy_source": assignment["tree_taxonomy_source"],
                    "taxonomy_domain": assignment["tree_taxonomy_domain"],
                    "compartment": assignment["tree_compartment"],
                    "taxonomy_assignment_method": assignment[
                        "tree_assignment_method"
                    ],
                    "taxonomy_alternatives": "",
                }
            )
        merged.append(updated)
    return merged


def load_tree_neighbor_rows(pattern: str) -> list[dict[str, str]]:
    rows = _load_rows(pattern, TREE_NEIGHBOR_FIELDS, "tree-neighbor")
    return sorted(
        rows,
        key=lambda row: (
            row["sample"],
            row["model"],
            row["name"],
            int(row["tree_neighbor_rank"]),
            row["leaf_id"],
        ),
    )


def write_tree_neighbors(
    rows: list[dict[str, str]], output_file: str | Path
) -> None:
    with Path(output_file).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TREE_NEIGHBOR_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_category_summary(
    rows: list[dict[str, str]],
    metadata: list[dict[str, str]],
    output_file: str | Path,
) -> None:
    samples = sorted({row["sample"] for row in metadata})
    counts = {sample: Counter() for sample in samples}
    categories_by_contig: dict[tuple[str, str], set[str]] = {}

    for row in rows:
        if not row["blast_sseqid"]:
            continue
        key = (row["sample"], row["contig_name"])
        taxonomy = row.get("taxonomy", "")
        taxonomy_domain = row.get("taxonomy_domain", "")
        compartment = row.get("compartment", "")
        tokens = (
            re.split(r"[-_;]", taxonomy)
            if taxonomy
            else re.split(r"[-_;|]", row["blast_sseqid"])
        )
        if taxonomy_domain:
            tokens.append(taxonomy_domain)
        row_categories: set[str] = set()
        for token in tokens:
            category = CATEGORY_MAPPING.get(token)
            if category:
                row_categories.add(category)
        compartment_category = COMPARTMENT_MAPPING.get(compartment)
        if compartment_category:
            row_categories.add(compartment_category)
        categories_by_contig.setdefault(key, set()).update(row_categories)

    for (sample, _contig), categories in categories_by_contig.items():
        for category in categories:
            counts[sample][category] += 1

    categories = [
        category
        for category in CATEGORY_MAPPING.values()
        if any(sample_counts[category] for sample_counts in counts.values())
    ]
    with Path(output_file).open("w", newline="") as handle:
        if not categories:
            handle.write("\n")
            for sample in samples:
                handle.write(f"{sample}\n")
            return

        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([""] + categories)
        for sample in samples:
            writer.writerow([sample] + [counts[sample][category] for category in categories])


def merge_m8_files(pattern: str, output_file: str | Path) -> None:
    output_path = Path(output_file).resolve()
    input_files = [
        Path(filename)
        for filename in sorted(glob.glob(pattern))
        if Path(filename).resolve() != output_path
    ]
    with output_path.open("w") as output_handle:
        for filename in input_files:
            with filename.open() as input_handle:
                for line in input_handle:
                    output_handle.write(line)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create deterministic SSUextract detailed and category summaries."
    )
    parser.add_argument("--summary-glob", default="*.summary.tsv")
    parser.add_argument("--metadata-glob", default="*.meta.tsv")
    parser.add_argument("--m8-glob", default="*.m8")
    parser.add_argument("--top-hits-glob", default="*.top_hits.tsv")
    parser.add_argument("--tree-assignment-glob", default="*.tree_assignment.tsv")
    parser.add_argument("--tree-neighbor-glob", default="*.tree_neighbors.tsv")
    parser.add_argument(
        "--taxonomy-mode", choices=("blast", "tree"), default="blast"
    )
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--category-output", required=True)
    parser.add_argument("--merged-m8-output", required=True)
    parser.add_argument("--top-hits-output", required=True)
    parser.add_argument("--tree-neighbor-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metadata = load_metadata(args.metadata_glob)
    rows = load_summary_rows(args.summary_glob)
    top_hit_rows = load_top_hit_rows(args.top_hits_glob)
    tree_assignments = load_tree_assignment_rows(args.tree_assignment_glob)
    tree_neighbors = load_tree_neighbor_rows(args.tree_neighbor_glob)
    rows = apply_tree_assignments(rows, tree_assignments, args.taxonomy_mode)
    write_detailed_summary(rows, args.summary_output)
    write_top_hit_summary(top_hit_rows, args.top_hits_output)
    write_tree_neighbors(tree_neighbors, args.tree_neighbor_output)
    write_category_summary(rows, metadata, args.category_output)
    merge_m8_files(args.m8_glob, args.merged_m8_output)


if __name__ == "__main__":
    main()
