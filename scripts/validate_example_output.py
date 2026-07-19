#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from pathlib import Path


KEY_FIELDS = ("sample", "contig_name", "model")
ANNOTATION_FIELDS = (
    "coordinates",
    "strand",
    "blast_sseqid",
    "centroid_names",
    "centroid_taxonomy",
    "centroid_taxonomy_source",
    "reference_source",
    "taxonomy",
    "taxonomy_source",
    "taxonomy_domain",
    "taxonomy_assignment_method",
)
CENTROID_FIELDS = (
    "centroid_names",
    "centroid_taxonomy",
    "centroid_taxonomy_source",
)
EXPECTED_MODEL_COUNTS = Counter({"RF00177": 9, "RF01960": 1})
CENTROID_CONTRACT_VERSION = (1, 0, 2)


def read_tsv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"{path} has no header")
            rows = list(reader)
    except OSError as error:
        raise ValueError(f"could not read {path}: {error}") from error
    return rows


def keyed_rows(
    rows: list[dict[str, str]], required_fields: tuple[str, ...], source: Path
) -> dict[tuple[str, ...], dict[str, str]]:
    if not rows:
        raise ValueError(f"{source} contains no annotation rows")
    missing_fields = [field for field in required_fields if field not in rows[0]]
    if missing_fields:
        raise ValueError(f"{source} is missing fields: {', '.join(missing_fields)}")

    keyed = {}
    for line_number, row in enumerate(rows, start=2):
        key = tuple(row[field] for field in KEY_FIELDS)
        if any(not value for value in key):
            raise ValueError(f"{source}:{line_number} has an empty example key")
        if key in keyed:
            raise ValueError(f"{source}:{line_number} repeats example key {key}")
        keyed[key] = row
    return keyed


def semantic_version(value: str) -> tuple[int, int, int]:
    fields = value.split(".")
    if len(fields) != 3 or any(not field.isdigit() for field in fields):
        raise ValueError(f"database version must use semantic x.y.z form: {value!r}")
    return tuple(int(field) for field in fields)


def validate_centroid_contract(
    rows: list[dict[str, str]],
    profile: str,
    database_version: str,
    source: Path,
) -> None:
    if semantic_version(database_version) < CENTROID_CONTRACT_VERSION:
        return

    derived = [
        row
        for row in rows
        if row["taxonomy_assignment_method"] == "updated_reference_cluster"
    ]
    native = [row for row in rows if row["taxonomy_assignment_method"] == "native"]
    if len(derived) + len(native) != len(rows):
        raise ValueError(f"{source} contains an unsupported assignment method")

    if profile == "curated":
        if derived or len(native) != 10:
            raise ValueError(
                f"{source} curated v{database_version} must contain 10 native rows"
            )
        if any(row[field] for row in rows for field in CENTROID_FIELDS):
            raise ValueError(
                f"{source} curated v{database_version} contains centroid evidence"
            )
        return

    if profile != "img":
        return
    if len(derived) != 6 or len(native) != 4:
        raise ValueError(
            f"{source} IMG v{database_version} must contain 6 cluster-derived "
            f"and 4 native rows; found {len(derived)} and {len(native)}"
        )
    if any(not row[field] for row in derived for field in CENTROID_FIELDS):
        raise ValueError(
            f"{source} IMG v{database_version} cluster-derived rows require "
            "nonempty centroid names, taxonomy, and taxonomy source"
        )
    if any(row[field] for row in native for field in CENTROID_FIELDS):
        raise ValueError(
            f"{source} IMG v{database_version} native rows contain centroid evidence"
        )
    if any(row["taxonomy"] != row["taxonomy_domain"] for row in derived):
        raise ValueError(
            f"{source} IMG v{database_version} cluster-derived member taxonomy "
            "must remain domain-capped"
        )
    deep_prokaryotic = [
        row
        for row in derived
        if row["taxonomy_domain"] in {"Bacteria", "Archaea"}
        and ";" in row["centroid_taxonomy"]
    ]
    if len(deep_prokaryotic) != 5:
        raise ValueError(
            f"{source} IMG v{database_version} must expose deeper centroid taxonomy "
            "for all five bacterial or archaeal cluster-derived rows"
        )


def validate_example(
    summary: Path,
    expectations: Path,
    profile: str,
    database_version: str,
) -> Counter[str]:
    all_expectations = read_tsv(expectations)
    expected_rows = [
        row
        for row in all_expectations
        if row.get("profile") == profile
        and row.get("database_version") == database_version
    ]
    if not expected_rows:
        supported_versions = sorted(
            {
                row.get("database_version", "")
                for row in all_expectations
                if row.get("profile") == profile
                and row.get("database_version", "")
            }
        )
        supported = ", ".join(f"v{version}" for version in supported_versions)
        if not supported:
            supported = "none"
        raise ValueError(
            f"this SSUextract checkout has no annotation contract for database "
            f"profile {profile!r} v{database_version}; supported versions: "
            f"{supported}. Update SSUextract before running this example with "
            "a different database release"
        )
    for row in expected_rows:
        for field in CENTROID_FIELDS:
            row.setdefault(field, "")

    expected = keyed_rows(
        expected_rows,
        (*KEY_FIELDS, *ANNOTATION_FIELDS),
        expectations,
    )
    validate_centroid_contract(
        list(expected.values()), profile, database_version, expectations
    )
    expected_model_counts = Counter(row["model"] for row in expected.values())
    if expected_model_counts != EXPECTED_MODEL_COUNTS:
        raise ValueError(
            "example annotation contract must contain exactly 9 RF00177 "
            "16S rRNA gene annotations and 1 RF01960 18S rRNA gene annotation; "
            f"found {dict(expected_model_counts)}"
        )
    observed = keyed_rows(
        read_tsv(summary),
        (*KEY_FIELDS, *ANNOTATION_FIELDS),
        summary,
    )
    validate_centroid_contract(
        list(observed.values()), profile, database_version, summary
    )

    missing = sorted(set(expected) - set(observed))
    extra = sorted(set(observed) - set(expected))
    changed = {
        key: {
            field: {"expected": expected[key][field], "observed": observed[key][field]}
            for field in ANNOTATION_FIELDS
            if expected[key][field] != observed[key][field]
        }
        for key in sorted(set(expected) & set(observed))
    }
    changed = {key: fields for key, fields in changed.items() if fields}
    if missing or extra or changed:
        raise ValueError(
            "annotation contract mismatch: "
            f"missing={missing}, extra={extra}, changed={changed}"
        )

    model_counts = Counter(row["model"] for row in observed.values())
    return model_counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate the bundled example against its biological contract."
    )
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--expectations", type=Path, required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--database-version", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        model_counts = validate_example(
            args.summary,
            args.expectations,
            args.profile,
            args.database_version,
        )
    except ValueError as error:
        raise SystemExit(f"Example validation failed: {error}") from error
    total = sum(model_counts.values())
    marker_16s = "16S rRNA gene" if model_counts["RF00177"] == 1 else "16S rRNA genes"
    marker_18s = "18S rRNA gene" if model_counts["RF01960"] == 1 else "18S rRNA genes"
    print(
        f"Example validation passed: {total} annotations "
        f"({model_counts['RF00177']} {marker_16s}, "
        f"{model_counts['RF01960']} {marker_18s}) with database profile "
        f"{args.profile} v{args.database_version}."
    )


if __name__ == "__main__":
    main()
