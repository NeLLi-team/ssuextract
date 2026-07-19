#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass, replace
from pathlib import Path

import duckdb

from hit_processing import HIT_FIELDS
from taxonomy_utils import common_value as _shared_common_value
from taxonomy_utils import lowest_common_ancestor, taxonomy_path


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
]


@dataclass(frozen=True)
class BlastHit:
    subject: str
    percent_identity: float
    alignment_length: int
    bit_score: float


@dataclass(frozen=True)
class TaxonomyRecord:
    reference_source: str
    taxonomy: str
    taxonomy_source: str
    domain: str
    compartment: str
    assignment_method: str
    cross_domain_conflict: bool
    taxonomy_alternatives: str
    centroid_names: str
    centroid_taxonomy: str
    centroid_taxonomy_source: str


def load_best_blast_hits(m8_file: str | Path) -> dict[str, list[BlastHit]]:
    best_hits: dict[str, list[BlastHit]] = {}
    with Path(m8_file).open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(
                    f"Malformed BLAST m8 row at {m8_file}:{line_number}: "
                    f"expected at least 12 fields, found {len(fields)}"
                )

            query = fields[0]
            hit = BlastHit(
                subject=fields[1],
                percent_identity=float(fields[2]),
                alignment_length=int(fields[3]),
                bit_score=float(fields[11]),
            )
            current = best_hits.get(query, [])
            if not current or hit.bit_score > current[0].bit_score:
                best_hits[query] = [hit]
            elif hit.bit_score == current[0].bit_score:
                current.append(hit)
    for query, hits in best_hits.items():
        hits.sort(
            key=lambda hit: (
                -hit.percent_identity,
                -hit.alignment_length,
                hit.subject,
            )
        )
        unique_subjects: dict[str, BlastHit] = {}
        for hit in hits:
            unique_subjects.setdefault(hit.subject, hit)
        best_hits[query] = list(unique_subjects.values())
    return best_hits


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


def lowest_common_taxonomy(taxonomies: list[str]) -> str:
    paths = []
    for taxonomy in taxonomies:
        if not taxonomy:
            continue
        try:
            paths.append(taxonomy_path(taxonomy))
        except ValueError:
            continue
    if not paths:
        return ""
    return ";".join(lowest_common_ancestor(paths))


def common_value(values: list[str]) -> str:
    return _shared_common_value(values, conflict="mixed")


def merged_taxonomy_alternatives(records: list[TaxonomyRecord]) -> str:
    alternatives: dict[str, dict[str, str]] = {}
    for record in records:
        if record.taxonomy_alternatives:
            try:
                parsed = json.loads(record.taxonomy_alternatives)
            except json.JSONDecodeError as error:
                raise ValueError("Invalid taxonomy_alternatives JSON") from error
            if not isinstance(parsed, list) or not all(
                isinstance(alternative, dict) for alternative in parsed
            ):
                raise ValueError("taxonomy_alternatives must be a JSON array of objects")
            candidates = parsed
        elif record.taxonomy or record.domain:
            candidates = [
                {
                    "taxonomy_source": record.taxonomy_source,
                    "taxonomy": record.taxonomy,
                    "domain": record.domain,
                    "compartment": record.compartment,
                    "assignment_method": record.assignment_method,
                }
            ]
        else:
            candidates = []
        for candidate in candidates:
            normalized = {str(key): str(value) for key, value in candidate.items()}
            key = json.dumps(normalized, ensure_ascii=True, sort_keys=True)
            alternatives[key] = normalized
    return json.dumps(
        [alternatives[key] for key in sorted(alternatives)],
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def annotate_hits(
    hits_file: str | Path,
    m8_file: str | Path,
    output_file: str | Path,
    taxonomy_file: str | Path | None = None,
    max_targets: int = 500,
) -> None:
    if max_targets < 1:
        raise ValueError("max_targets must be positive")
    best_hits = load_best_blast_hits(m8_file)
    taxonomy_records: dict[str, TaxonomyRecord] = {}
    if taxonomy_file is not None:
        subjects = {
            hit.subject
            for query_hits in best_hits.values()
            for hit in query_hits
        }
        taxonomy_records = load_taxonomy_records(taxonomy_file, subjects)
    with Path(hits_file).open(newline="") as hits_handle, Path(output_file).open(
        "w", newline=""
    ) as output_handle:
        reader = csv.DictReader(hits_handle, delimiter="\t")
        if reader.fieldnames != HIT_FIELDS:
            raise ValueError(
                f"Unexpected hit-table columns in {hits_file}: {reader.fieldnames}"
            )

        writer = csv.DictWriter(
            output_handle,
            fieldnames=SUMMARY_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            tied_hits = best_hits.get(row["name"], [])
            blast_hit = tied_hits[0] if tied_hits else None
            tied_taxonomies = [
                taxonomy_records[hit.subject]
                for hit in tied_hits
                if hit.subject in taxonomy_records
            ]
            blast_taxonomy = (
                taxonomy_records.get(blast_hit.subject) if blast_hit else None
            )
            ties_truncated = len(tied_hits) > max_targets
            concrete_domains = {
                record.domain
                for record in tied_taxonomies
                if record.domain not in {"", "ambiguous", "Unclassified"}
            }
            cross_domain_ambiguity = (
                len(concrete_domains) > 1
                or any(
                    record.cross_domain_conflict or record.domain == "ambiguous"
                    for record in tied_taxonomies
                )
            )
            if cross_domain_ambiguity:
                taxonomy_domain = "ambiguous"
                taxonomy = ""
                compartment = "mixed"
                assignment_method = (
                    tied_taxonomies[0].assignment_method
                    if len(tied_taxonomies) == 1
                    else "cross_domain_ambiguous_equal_best"
                )
                taxonomy_alternatives = merged_taxonomy_alternatives(tied_taxonomies)
            else:
                taxonomy_domain = common_value(
                    [record.domain for record in tied_taxonomies]
                )
                taxonomy = (
                    taxonomy_domain
                    if ties_truncated
                    else lowest_common_taxonomy(
                        [record.taxonomy for record in tied_taxonomies]
                    )
                )
                compartment = "" if ties_truncated else common_value(
                    [record.compartment for record in tied_taxonomies]
                )
                assignment_method = (
                    "truncated_equal_best_lca"
                    if ties_truncated and tied_taxonomies
                    else common_value(
                        [record.assignment_method for record in tied_taxonomies]
                    )
                )
                taxonomy_alternatives = common_value(
                    [record.taxonomy_alternatives for record in tied_taxonomies]
                )
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
                        blast_taxonomy.centroid_names if blast_taxonomy else ""
                    ),
                    "centroid_taxonomy": (
                        blast_taxonomy.centroid_taxonomy if blast_taxonomy else ""
                    ),
                    "centroid_taxonomy_source": (
                        blast_taxonomy.centroid_taxonomy_source
                        if blast_taxonomy
                        else ""
                    ),
                    "blast_pident": (
                        str(round(blast_hit.percent_identity, 2)) if blast_hit else ""
                    ),
                    "blast_length": blast_hit.alignment_length if blast_hit else "",
                    "blast_bitscore": (
                        format(blast_hit.bit_score, "g") if blast_hit else ""
                    ),
                    "is_assembled": row["is_assembled"],
                    "reference_source": common_value(
                        [record.reference_source for record in tied_taxonomies]
                    ),
                    "taxonomy": taxonomy,
                    "taxonomy_source": common_value(
                        [record.taxonomy_source for record in tied_taxonomies]
                    ),
                    "taxonomy_domain": taxonomy_domain,
                    "compartment": compartment,
                    "taxonomy_assignment_method": assignment_method,
                    "taxonomy_alternatives": taxonomy_alternatives,
                    "blast_tied_subjects": len(tied_hits) if tied_hits else "",
                    "blast_ties_truncated": (
                        str(ties_truncated).lower() if tied_hits else ""
                    ),
                }
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add the highest-bit-score BLAST hit to an SSUextract hit table."
    )
    parser.add_argument("--hits", required=True)
    parser.add_argument("--m8", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--taxonomy-db",
        help="Preferred-taxonomy Parquet for manifest-driven databases.",
    )
    parser.add_argument(
        "--max-targets",
        type=int,
        default=500,
        help="Policy limit; BLAST must request one additional overflow target.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    annotate_hits(args.hits, args.m8, args.output, args.taxonomy_db, args.max_targets)


if __name__ == "__main__":
    main()
