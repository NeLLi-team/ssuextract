"""Immutable models and input loading for IMG cluster classification."""

from __future__ import annotations

import ast
import csv
import json
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable, Mapping

from taxonomy_utils import taxonomy_path


BLAST_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qcovs",
    "bitscore",
)


@dataclass(frozen=True)
class BlastHit:
    query: str
    subject: str
    percent_identity: Decimal
    alignment_length: int
    query_length: int
    subject_length: int
    query_coverage: Decimal
    bit_score: Decimal


@dataclass(frozen=True)
class TaxonomyRecord:
    taxonomy: tuple[str, ...]
    taxonomy_source: str
    domain: str
    compartment: str
    cross_domain_conflict: bool = False
    taxonomy_alternatives: str = ""


@dataclass(frozen=True)
class Cluster:
    cluster_id: str
    centroid: str
    img_members: tuple[str, ...]


@dataclass(frozen=True)
class CalibrationStratum:
    status: str
    rank_cap: int | None
    reason: str


@dataclass(frozen=True)
class CalibrationData:
    rank_caps: dict[str, int]
    strata: dict[str, CalibrationStratum]


def _decimal(value: str, field: str, line_number: int) -> Decimal:
    try:
        number = Decimal(value)
    except InvalidOperation as error:
        raise ValueError(
            f"Invalid BLAST {field} on line {line_number}: {value!r}"
        ) from error
    if not number.is_finite():
        raise ValueError(f"Non-finite BLAST {field} on line {line_number}")
    return number


def parse_blast_hits(lines: Iterable[str]) -> dict[str, tuple[BlastHit, ...]]:
    """Parse eight-column BLAST outfmt 6 and reject multiple HSPs per subject."""

    by_query: dict[str, list[BlastHit]] = defaultdict(list)
    observed_pairs: set[tuple[str, str]] = set()
    for line_number, raw_line in enumerate(lines, 1):
        line = raw_line.rstrip("\r\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != len(BLAST_FIELDS):
            raise ValueError(
                f"Malformed BLAST row {line_number}: expected {len(BLAST_FIELDS)} "
                f"columns ({' '.join(BLAST_FIELDS)}), found {len(fields)}"
            )
        query, subject = fields[:2]
        if not query or not subject:
            raise ValueError(f"Empty BLAST query or subject on line {line_number}")
        pair = (query, subject)
        if pair in observed_pairs:
            raise ValueError(
                f"Multiple HSPs for BLAST query/subject {query!r}/{subject!r}; "
                "run BLAST with one HSP per subject"
            )
        observed_pairs.add(pair)
        percent_identity = _decimal(fields[2], "pident", line_number)
        try:
            alignment_length, query_length, subject_length = map(int, fields[3:6])
        except ValueError as error:
            raise ValueError(f"Invalid BLAST sequence length on line {line_number}") from error
        if min(alignment_length, query_length, subject_length) < 1:
            raise ValueError(f"Non-positive BLAST sequence length on line {line_number}")
        query_coverage = _decimal(fields[6], "qcovs", line_number)
        bit_score = _decimal(fields[7], "bitscore", line_number)
        if not Decimal("0") <= percent_identity <= Decimal("100"):
            raise ValueError(f"BLAST pident outside [0, 100] on line {line_number}")
        if not Decimal("0") <= query_coverage <= Decimal("100"):
            raise ValueError(f"BLAST qcovs outside [0, 100] on line {line_number}")
        if bit_score < 0:
            raise ValueError(f"Negative BLAST bitscore on line {line_number}")
        by_query[query].append(
            BlastHit(
                query,
                subject,
                percent_identity,
                alignment_length,
                query_length,
                subject_length,
                query_coverage,
                bit_score,
            )
        )

    result: dict[str, tuple[BlastHit, ...]] = {}
    for query, hits in by_query.items():
        result[query] = tuple(
            sorted(
                hits,
                key=lambda hit: (
                    -hit.bit_score,
                    -hit.percent_identity,
                    -hit.query_coverage,
                    hit.subject,
                ),
            )
        )
    return result


def _taxonomy_path(value: object) -> tuple[str, ...]:
    try:
        return taxonomy_path(str(value or ""))
    except ValueError as error:
        raise ValueError("Reference taxonomy must contain a non-empty domain") from error


def _boolean(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value or "").strip().lower()
    if text in {"", "false", "0"}:
        return False
    if text in {"true", "1"}:
        return True
    raise ValueError(f"Invalid boolean value: {value!r}")


def _taxonomy_record(row: Mapping[str, object], subject: str) -> TaxonomyRecord:
    source = str(row.get("taxonomy_source", "")).strip()
    domain = str(row.get("domain", "")).strip()
    cross_domain_conflict = _boolean(row.get("cross_domain_conflict", False))
    alternatives = str(row.get("taxonomy_alternatives", "") or "").strip()
    if not source or not domain:
        raise ValueError(f"Reference taxonomy row {subject!r} lacks source or domain")
    if cross_domain_conflict:
        if domain != "ambiguous" or not alternatives:
            raise ValueError(
                f"Ambiguous reference taxonomy row {subject!r} lacks explicit alternatives"
            )
        try:
            parsed_alternatives = json.loads(alternatives)
        except json.JSONDecodeError as error:
            raise ValueError(
                f"Ambiguous reference taxonomy row {subject!r} has invalid alternatives JSON"
            ) from error
        if not isinstance(parsed_alternatives, list) or not parsed_alternatives:
            raise ValueError(
                f"Ambiguous reference taxonomy row {subject!r} has no alternatives"
            )
        return TaxonomyRecord(
            taxonomy=(),
            taxonomy_source=source,
            domain=domain,
            compartment=str(row.get("compartment", "") or "").strip(),
            cross_domain_conflict=True,
            taxonomy_alternatives=_canonical_json(parsed_alternatives),
        )
    taxonomy = _taxonomy_path(row.get("taxonomy", ""))
    if taxonomy[0] != domain:
        raise ValueError(
            f"Reference taxonomy row {subject!r} has domain {domain!r} "
            f"but taxonomy begins with {taxonomy[0]!r}"
        )
    return TaxonomyRecord(
        taxonomy=taxonomy,
        taxonomy_source=source,
        domain=domain,
        compartment=str(row.get("compartment", "") or "").strip(),
        cross_domain_conflict=False,
        taxonomy_alternatives="",
    )


def load_taxonomy_tsv(
    path: str | Path, subjects: set[str]
) -> dict[str, TaxonomyRecord]:
    records: dict[str, TaxonomyRecord] = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sequence_id", "taxonomy", "taxonomy_source", "domain"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Taxonomy TSV lacks columns: {sorted(missing)}")
        for row in reader:
            subject = row["sequence_id"].strip()
            if subject not in subjects:
                continue
            if subject in records:
                raise ValueError(f"Duplicate reference taxonomy row for {subject}")
            records[subject] = _taxonomy_record(row, subject)
    return records


def load_taxonomy_parquet(
    path: str | Path, subjects: set[str]
) -> dict[str, TaxonomyRecord]:
    try:
        import duckdb
    except ImportError as error:  # pragma: no cover - project environment provides DuckDB
        raise RuntimeError("DuckDB is required to read taxonomy Parquet") from error

    if not subjects:
        return {}
    connection = duckdb.connect(":memory:")
    try:
        raw_thread_limit = (
            os.environ.get("SLURM_CPUS_PER_TASK")
            or os.environ.get("OMP_NUM_THREADS")
            or "1"
        )
        thread_limit = int(raw_thread_limit)
        if thread_limit < 1:
            raise ValueError("DuckDB thread limit must be positive")
        connection.execute("SET threads = ?", [thread_limit])
        connection.execute("CREATE TABLE wanted(sequence_id VARCHAR PRIMARY KEY)")
        connection.executemany(
            "INSERT INTO wanted VALUES (?)", [(subject,) for subject in sorted(subjects)]
        )
        rows = connection.execute(
            """
            SELECT p.sequence_id, p.taxonomy, p.taxonomy_source, p.domain, p.compartment
                 , p.cross_domain_conflict, p.taxonomy_alternatives
            FROM read_parquet(?) AS p
            INNER JOIN wanted AS w USING (sequence_id)
            ORDER BY p.sequence_id
            """,
            [str(Path(path))],
        ).fetchall()
    finally:
        connection.close()

    records: dict[str, TaxonomyRecord] = {}
    for (
        sequence_id,
        taxonomy,
        source,
        domain,
        compartment,
        cross_domain_conflict,
        taxonomy_alternatives,
    ) in rows:
        subject = str(sequence_id)
        if subject in records:
            raise ValueError(f"Duplicate reference taxonomy row for {subject}")
        records[subject] = _taxonomy_record(
            {
                "taxonomy": taxonomy,
                "taxonomy_source": source,
                "domain": domain,
                "compartment": compartment,
                "cross_domain_conflict": cross_domain_conflict,
                "taxonomy_alternatives": taxonomy_alternatives,
            },
            subject,
        )
    return records


def load_taxonomy(
    path: str | Path, subjects: set[str]
) -> dict[str, TaxonomyRecord]:
    source = Path(path)
    if source.suffix.lower() == ".parquet":
        return load_taxonomy_parquet(source, subjects)
    return load_taxonomy_tsv(source, subjects)


def load_calibration(path: str | Path) -> CalibrationData:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Could not read taxonomy calibration {path}: {error}") from error
    if not isinstance(data, dict) or data.get("schema_version") != 2:
        raise ValueError("taxonomy calibration schema_version must be 2")
    raw_caps = data.get("rank_caps")
    raw_strata = data.get("strata")
    if not isinstance(raw_caps, dict):
        raise ValueError("taxonomy calibration must contain a rank_caps object")
    if not isinstance(raw_strata, dict) or not raw_strata:
        raise ValueError("taxonomy calibration must contain a non-empty strata object")
    caps: dict[str, int] = {}
    for key, value in raw_caps.items():
        if not isinstance(key, str) or not key or type(value) is not int or value < 0:
            raise ValueError("taxonomy calibration rank caps must be non-negative integers")
        caps[key] = value
    strata: dict[str, CalibrationStratum] = {}
    calibrated_keys: set[str] = set()
    for key, raw_stratum in raw_strata.items():
        if not isinstance(key, str) or not key or not isinstance(raw_stratum, dict):
            raise ValueError("taxonomy calibration strata must be keyed objects")
        if not {"status", "rank_cap", "reason"}.issubset(raw_stratum):
            raise ValueError(
                "taxonomy calibration strata require status, rank_cap, and reason"
            )
        status = raw_stratum.get("status")
        rank_cap = raw_stratum.get("rank_cap")
        reason = raw_stratum.get("reason")
        if status == "calibrated":
            if type(rank_cap) is not int or rank_cap < 0 or reason != "":
                raise ValueError(
                    "calibrated taxonomy strata require a non-negative rank_cap "
                    "and an empty reason"
                )
            calibrated_keys.add(key)
        elif status == "failed":
            if rank_cap is not None or not isinstance(reason, str) or not reason:
                raise ValueError(
                    "failed taxonomy strata require a null rank_cap and a reason"
                )
        else:
            raise ValueError(
                "taxonomy calibration stratum status must be 'calibrated' or 'failed'"
            )
        strata[key] = CalibrationStratum(status, rank_cap, reason)
    if set(caps) != calibrated_keys or any(
        caps[key] != strata[key].rank_cap for key in calibrated_keys
    ):
        raise ValueError(
            "taxonomy calibration rank_caps must exactly match calibrated strata"
        )
    return CalibrationData(caps, strata)


def parse_clusters(lines: Iterable[str]) -> tuple[Cluster, ...]:
    """Read only cluster identity, centroid, and IMG membership columns."""

    csv.field_size_limit(sys.maxsize)
    reader = csv.DictReader(lines, delimiter="\t")
    required = {"cluster_id", "centroid", "sequences"}
    missing = required - set(reader.fieldnames or ())
    if missing:
        raise ValueError(f"Cluster table lacks columns: {sorted(missing)}")

    clusters: list[Cluster] = []
    cluster_ids: set[str] = set()
    member_clusters: dict[str, str] = {}
    for row_number, row in enumerate(reader, 2):
        cluster_id = (row.get("cluster_id") or "").strip()
        centroid = (row.get("centroid") or "").strip()
        if not cluster_id or not centroid:
            raise ValueError(f"Cluster row {row_number} lacks cluster_id or centroid")
        if cluster_id in cluster_ids:
            raise ValueError(f"Duplicate cluster_id {cluster_id!r}")
        cluster_ids.add(cluster_id)
        try:
            raw_members = ast.literal_eval(row.get("sequences") or "")
        except (SyntaxError, ValueError) as error:
            raise ValueError(
                f"Cluster {cluster_id!r} has an invalid sequences list"
            ) from error
        if not isinstance(raw_members, list) or not all(
            isinstance(member, str) for member in raw_members
        ):
            raise ValueError(f"Cluster {cluster_id!r} sequences must be a list of strings")
        img_members = tuple(sorted({member for member in raw_members if member.startswith("IMG_")}))
        for member in img_members:
            previous = member_clusters.get(member)
            if previous is not None and previous != cluster_id:
                raise ValueError(
                    f"IMG member {member!r} occurs in clusters {previous!r} and {cluster_id!r}"
                )
            member_clusters[member] = cluster_id
        clusters.append(Cluster(cluster_id, centroid, img_members))
    return tuple(sorted(clusters, key=lambda cluster: (cluster.cluster_id, cluster.centroid)))


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
