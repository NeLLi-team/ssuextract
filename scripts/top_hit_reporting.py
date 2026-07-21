from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import duckdb
from Bio import SeqIO

if TYPE_CHECKING:
    from annotate_hits import BlastHit, TaxonomyRecord


TOP_HIT_FIELDS = [
    "name",
    "sample",
    "model",
    "hit_rank",
    "selection_reason",
    "blast_sseqid",
    "reference_identifiers",
    "reference_versions",
    "reference_source",
    "taxonomy",
    "taxonomy_source",
    "taxonomy_domain",
    "compartment",
    "taxonomy_assignment_method",
    "taxonomy_alternatives",
    "centroid_names",
    "centroid_taxonomy",
    "centroid_taxonomy_source",
    "blast_pident",
    "blast_length",
    "blast_mismatch",
    "blast_gapopen",
    "blast_qstart",
    "blast_qend",
    "blast_sstart",
    "blast_send",
    "blast_evalue",
    "blast_bitscore",
    "query_length",
    "query_sequence",
]


@dataclass(frozen=True)
class ReferenceRecord:
    identifiers: str
    versions: str
    sources: tuple[str, ...]


def load_reference_records(
    source_records_file: str | Path,
    subjects: set[str],
) -> dict[str, ReferenceRecord]:
    if not subjects:
        return {}
    placeholders = ",".join("?" for _ in subjects)
    parameters = [str(Path(source_records_file)), *sorted(subjects)]
    connection = duckdb.connect(":memory:")
    try:
        connection.execute("SET threads = 1")
        rows = connection.execute(
            f"""
            SELECT sequence_id, reference_source, source_version, source_identifier
            FROM read_parquet(?)
            WHERE sequence_id IN ({placeholders})
            ORDER BY sequence_id, reference_source, source_version, source_identifier
            """,
            parameters,
        ).fetchall()
    finally:
        connection.close()

    identifiers: dict[str, set[str]] = {}
    versions: dict[str, set[str]] = {}
    for sequence_id, source, version, identifier in rows:
        sequence_id = str(sequence_id)
        source = str(source or "")
        version = str(version or "")
        identifier = str(identifier or "")
        for label, value in (
            ("reference source", source),
            ("source version", version),
            ("source identifier", identifier),
        ):
            if not value or any(character in value for character in "\r\n\t|"):
                raise ValueError(f"Invalid {label} for {sequence_id}")
        identifiers.setdefault(sequence_id, set()).add(f"{source}:{identifier}")
        versions.setdefault(sequence_id, set()).add(f"{source}:{version}")

    missing = sorted(subjects - identifiers.keys())
    if missing:
        preview = ", ".join(missing[:5])
        suffix = "" if len(missing) <= 5 else f" and {len(missing) - 5} more"
        raise ValueError(
            f"Missing source metadata for {len(missing)} BLAST subject(s): "
            f"{preview}{suffix}"
        )
    return {
        sequence_id: ReferenceRecord(
            identifiers="|".join(sorted(values)),
            versions="|".join(sorted(versions[sequence_id])),
            sources=tuple(
                sorted({value.split(":", maxsplit=1)[0] for value in values})
            ),
        )
        for sequence_id, values in identifiers.items()
    }


def load_query_sequences(fasta_file: str | Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    with Path(fasta_file).open() as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in sequences:
                raise ValueError(f"Duplicate extracted FASTA identifier: {record.id}")
            sequence = str(record.seq)
            if not sequence:
                raise ValueError(f"Empty extracted FASTA sequence: {record.id}")
            sequences[record.id] = sequence
    return sequences


def reference_record(
    subject: str,
    references: dict[str, ReferenceRecord],
) -> ReferenceRecord:
    return references.get(subject, ReferenceRecord(subject, "", ()))


def _top_hit_taxonomy(record: TaxonomyRecord | None) -> str:
    if record is None:
        return ""
    return record.taxonomy or record.domain or "Unclassified"


def select_reported_hits(
    hits: list[BlastHit],
    references: dict[str, ReferenceRecord],
    top_hits: int,
) -> list[tuple[int, BlastHit, str]]:
    reason_order = (
        "overall_top_n",
        "equal_best_assignment",
        "best_IMG",
        "best_PR2",
        "best_SILVA",
    )
    reasons: dict[str, set[str]] = {}
    for hit in hits[:top_hits]:
        reasons.setdefault(hit.subject, set()).add("overall_top_n")
    if hits:
        best_score = hits[0].bit_score
        for hit in hits:
            if hit.bit_score != best_score:
                break
            reasons.setdefault(hit.subject, set()).add("equal_best_assignment")
    for source in ("IMG", "PR2", "SILVA"):
        best = next(
            (
                hit
                for hit in hits
                if source in reference_record(hit.subject, references).sources
            ),
            None,
        )
        if best is not None:
            reasons.setdefault(best.subject, set()).add(f"best_{source}")
    return [
        (
            rank,
            hit,
            "|".join(
                reason for reason in reason_order if reason in reasons[hit.subject]
            ),
        )
        for rank, hit in enumerate(hits, start=1)
        if hit.subject in reasons
    ]


def write_top_hits(
    output_file: str | Path,
    hit_rows: list[dict[str, str]],
    blast_hits: dict[str, list[BlastHit]],
    taxonomy_records: dict[str, TaxonomyRecord],
    reference_records: dict[str, ReferenceRecord],
    query_sequences: dict[str, str],
    top_hits: int,
) -> None:
    with Path(output_file).open("w", newline="") as output_handle:
        writer = csv.DictWriter(
            output_handle,
            fieldnames=TOP_HIT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in hit_rows:
            query = row["name"]
            query_sequence = query_sequences.get(query, "")
            reported_hits = select_reported_hits(
                blast_hits.get(query, []), reference_records, top_hits
            )
            for rank, hit, selection_reason in reported_hits:
                taxonomy = taxonomy_records.get(hit.subject)
                reference = reference_record(hit.subject, reference_records)
                writer.writerow(
                    {
                        "name": query,
                        "sample": row["sample"],
                        "model": row["model"],
                        "hit_rank": rank,
                        "selection_reason": selection_reason,
                        "blast_sseqid": hit.subject,
                        "reference_identifiers": reference.identifiers,
                        "reference_versions": reference.versions,
                        "reference_source": (
                            taxonomy.reference_source if taxonomy else ""
                        ),
                        "taxonomy": _top_hit_taxonomy(taxonomy),
                        "taxonomy_source": (
                            taxonomy.taxonomy_source if taxonomy else ""
                        ),
                        "taxonomy_domain": taxonomy.domain if taxonomy else "",
                        "compartment": taxonomy.compartment if taxonomy else "",
                        "taxonomy_assignment_method": (
                            taxonomy.assignment_method if taxonomy else ""
                        ),
                        "taxonomy_alternatives": (
                            taxonomy.taxonomy_alternatives if taxonomy else ""
                        ),
                        "centroid_names": taxonomy.centroid_names if taxonomy else "",
                        "centroid_taxonomy": (
                            taxonomy.centroid_taxonomy if taxonomy else ""
                        ),
                        "centroid_taxonomy_source": (
                            taxonomy.centroid_taxonomy_source if taxonomy else ""
                        ),
                        "blast_pident": format(hit.percent_identity, "g"),
                        "blast_length": hit.alignment_length,
                        "blast_mismatch": hit.mismatches,
                        "blast_gapopen": hit.gap_opens,
                        "blast_qstart": hit.query_start,
                        "blast_qend": hit.query_end,
                        "blast_sstart": hit.subject_start,
                        "blast_send": hit.subject_end,
                        "blast_evalue": format(hit.evalue, "g"),
                        "blast_bitscore": format(hit.bit_score, "g"),
                        "query_length": str(
                            len(query_sequence) if query_sequence else row["length"]
                        ),
                        "query_sequence": query_sequence,
                    }
                )
