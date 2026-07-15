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
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--category-output", required=True)
    parser.add_argument("--merged-m8-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metadata = load_metadata(args.metadata_glob)
    rows = load_summary_rows(args.summary_glob)
    write_detailed_summary(rows, args.summary_output)
    write_category_summary(rows, metadata, args.category_output)
    merge_m8_files(args.m8_glob, args.merged_m8_output)


if __name__ == "__main__":
    main()
