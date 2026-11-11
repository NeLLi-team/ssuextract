#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import glob
import json
import re
import sqlite3
import sys
from collections import Counter
from collections.abc import Callable, Iterable, Iterator
from contextlib import closing
from itertools import chain
from pathlib import Path
from tempfile import TemporaryDirectory

from annotate_hits import SUMMARY_FIELDS
from hit_processing import META_FIELDS
from top_hit_reporting import TOP_HIT_FIELDS
from tree_schema import TREE_ASSIGNMENT_FIELDS, TREE_NEIGHBOR_FIELDS

csv.field_size_limit(sys.maxsize)


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

UNCLASSIFIED_CATEGORY = "Unclassified"

COMPARTMENT_MAPPING = {
    "mitochondrion": "MitochondriaSSU",
    "plastid": "PlastidSSU",
    "apicoplast": "PlastidSSU",
    "chromatophore": "PlastidSSU",
}


def _summary_sort_key(row: dict[str, str]) -> tuple[str | int, ...]:
    start, end = row["coordinates"].split("-", maxsplit=1)
    return (
        row["sample"],
        row["model"],
        row["contig_name"],
        int(start),
        int(end),
        row["strand"],
    )


def _top_hit_sort_key(row: dict[str, str]) -> tuple[str | int, ...]:
    return (
        row["sample"],
        row["model"],
        row["name"],
        int(row["hit_rank"]),
        row["blast_sseqid"],
    )


def _sort_rows(
    rows: Iterable[dict[str, str]],
    key: Callable[[dict[str, str]], tuple[str | int, ...]],
) -> Iterator[dict[str, str]]:
    """Sort on disk, preserving input order for equal keys."""
    iterator = iter(rows)
    first = next(iterator, None)
    if first is None:
        return
    columns = ", ".join(f"sort{index}" for index in range(len(key(first))))
    placeholders = ", ".join("?" for _ in range(len(key(first)) + 1))
    with (
        TemporaryDirectory(prefix="ssuextract-sort-", dir=Path.cwd()) as scratch,
        closing(sqlite3.connect(Path(scratch) / "rows.sqlite")) as connection,
    ):
        connection.execute("PRAGMA temp_store = FILE")
        connection.execute("PRAGMA cache_size = -2048")
        connection.execute(f"CREATE TABLE rows ({columns}, payload TEXT)")
        connection.executemany(
            f"INSERT INTO rows VALUES ({placeholders})",  # noqa: S608  # Only placeholder count varies; values are bound.
            ((*key(row), json.dumps(row)) for row in chain((first,), iterator)),
        )
        connection.commit()
        for (payload,) in connection.execute(
            f"SELECT payload FROM rows ORDER BY {columns}, rowid"  # noqa: S608  # Column names are generated numeric indices.
        ):
            yield json.loads(payload)


def write_detailed_summary(
    rows: Iterable[dict[str, str]], output_file: str | Path
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


def write_top_hit_summary(
    rows: Iterable[dict[str, str]], output_file: str | Path
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
) -> Iterator[dict[str, str]]:
    for filename in sorted(glob.glob(pattern)):
        with Path(filename).open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != expected_fields:
                raise ValueError(
                    f"Unexpected {label} columns in {filename}: {reader.fieldnames}"
                )
            yield from reader


def load_tree_assignment_rows(pattern: str) -> list[dict[str, str]]:
    rows = list(_load_rows(pattern, TREE_ASSIGNMENT_FIELDS, "tree-assignment"))
    keys = [(row["sample"], row["model"], row["name"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate tree assignment for one extracted query")
    return rows


def apply_tree_assignments(
    rows: Iterable[dict[str, str]],
    assignments: list[dict[str, str]],
    taxonomy_mode: str,
) -> Iterator[dict[str, str]]:
    """Apply tree calls as rows arrive and check the complete query key set."""
    if taxonomy_mode not in {"blast", "tree"}:
        raise ValueError(f"Unsupported taxonomy mode: {taxonomy_mode}")
    if taxonomy_mode == "blast":
        if assignments:
            raise ValueError("Tree assignments were supplied in BLAST taxonomy mode")
        yield from rows
        return
    by_key = {
        (assignment["sample"], assignment["model"], assignment["name"]): assignment
        for assignment in assignments
    }
    expected: set[tuple[str, str, str]] = set()
    for row in rows:
        query_key = (row["sample"], row["model"], row["name"])
        expected.add(query_key)
        if query_key not in by_key:
            continue
        assignment = by_key[query_key]
        updated = dict(row)
        updated.update(
            {
                field: assignment[field]
                for field in TREE_ASSIGNMENT_FIELDS[3:]
            }
        )
        if assignment["tree_assignment_method"].startswith("tree_skipped_"):
            updated.update(
                {
                    "taxonomy_mode": "blast",
                    "taxonomy": assignment["tree_route_blast_taxonomy"],
                    "taxonomy_source": assignment[
                        "tree_route_blast_taxonomy_source"
                    ],
                    "taxonomy_domain": assignment[
                        "tree_route_blast_taxonomy_domain"
                    ],
                    "compartment": assignment["tree_route_blast_compartment"],
                    "taxonomy_assignment_method": assignment[
                        "tree_route_blast_assignment_method"
                    ],
                    "taxonomy_alternatives": assignment[
                        "tree_taxonomy_alternatives"
                    ],
                }
            )
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
                    "taxonomy_alternatives": assignment[
                        "tree_taxonomy_alternatives"
                    ],
                }
            )
        yield updated
    if set(by_key) != expected:
        missing = sorted(expected - by_key.keys())
        extra = sorted(by_key.keys() - expected)
        raise ValueError(
            "Tree assignments do not match extracted queries; "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )


def load_tree_neighbor_rows(pattern: str) -> Iterator[dict[str, str]]:
    """Stream neighbors in query and numeric neighbor-rank order."""
    rows = _load_rows(pattern, TREE_NEIGHBOR_FIELDS, "tree-neighbor")
    return _sort_rows(
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
    rows: Iterable[dict[str, str]], output_file: str | Path
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
    rows: Iterable[dict[str, str]],
    metadata: Iterable[dict[str, str]],
    output_file: str | Path,
) -> None:
    samples = sorted({row["sample"] for row in metadata})
    counts = {sample: Counter() for sample in samples}
    categories_by_contig: dict[tuple[str, str], set[str]] = {}

    for row in rows:
        key = (row["sample"], row["contig_name"])
        taxonomy = row.get("taxonomy", "")
        taxonomy_domain = row.get("taxonomy_domain", "")
        compartment = row.get("compartment", "")
        tokens = re.split(r"[-_;]", taxonomy)
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
        if taxonomy == UNCLASSIFIED_CATEGORY:
            row_categories.add(UNCLASSIFIED_CATEGORY)
        categories_by_contig.setdefault(key, set()).update(row_categories)

    for (sample, _contig), categories in categories_by_contig.items():
        for category in categories:
            counts[sample][category] += 1

    categories = [
        category
        for category in [*CATEGORY_MAPPING.values(), UNCLASSIFIED_CATEGORY]
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
    metadata = _load_rows("*.meta.tsv", META_FIELDS, "metadata")
    rows = _sort_rows(
        _load_rows("*.summary.tsv", SUMMARY_FIELDS, "summary"), _summary_sort_key
    )
    top_hit_rows = _sort_rows(
        _load_rows("*.top_hits.tsv", TOP_HIT_FIELDS, "top-hit"), _top_hit_sort_key
    )
    tree_assignments = load_tree_assignment_rows("*.tree_assignment.tsv")
    tree_neighbors = load_tree_neighbor_rows("*.tree_neighbors.tsv")
    rows = apply_tree_assignments(rows, tree_assignments, args.taxonomy_mode)
    write_detailed_summary(rows, args.summary_output)
    write_top_hit_summary(top_hit_rows, args.top_hits_output)
    write_tree_neighbors(tree_neighbors, args.tree_neighbor_output)
    with Path(args.summary_output).open(newline="", encoding="utf-8") as handle:
        write_category_summary(
            csv.DictReader(handle, delimiter="\t"), metadata, args.category_output
        )
    merge_m8_files("*.m8", args.merged_m8_output)


if __name__ == "__main__":
    main()
