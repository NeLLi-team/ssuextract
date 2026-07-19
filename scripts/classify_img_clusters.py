#!/usr/bin/env python3
"""Classify IMG clusters from centroid BLAST hits against updated references.

The cluster table is used only for cluster identifiers, centroids, and member
propagation.  Its legacy taxonomy columns are deliberately ignored.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import tempfile
from collections import Counter
from decimal import Decimal
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence, TextIO

from atomic_io import replace_and_fsync
from img_classification_data import (
    BLAST_FIELDS,
    BlastHit,
    CalibrationData,
    CalibrationStratum,
    Cluster,
    TaxonomyRecord,
    load_calibration,
    load_taxonomy,
    load_taxonomy_parquet,
    load_taxonomy_tsv,
    parse_blast_hits,
    parse_clusters,
    _taxonomy_record,
)
import img_search_provenance
from taxonomy_utils import common_value as _shared_common_value
from taxonomy_utils import lowest_common_ancestor as _shared_lca


ASSIGNMENT_FIELDS = (
    "source_identifier",
    "taxonomy",
    "taxonomy_source",
    "assignment_method",
    "evidence",
    "compartment",
    "cluster_id",
    "centroid",
    "centroid_name",
    "centroid_taxonomy",
    "centroid_taxonomy_source",
    "evidence_id",
)
MIN_QUERY_COVERAGE = Decimal("80")
CANDIDATE_BITSCORE_FRACTION = Decimal("0.98")
BLAST_MAX_TARGETS = 500
BLAST_FETCH_TARGETS = BLAST_MAX_TARGETS + 1


def _decimal_text(value: Decimal) -> str:
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text or "0"


def classification_policy(
    *,
    max_targets: int = BLAST_MAX_TARGETS,
    blast_fetch_targets: int = BLAST_FETCH_TARGETS,
) -> dict[str, object]:
    return {
        "blast_fields": list(BLAST_FIELDS),
        "min_query_coverage_percent": _decimal_text(MIN_QUERY_COVERAGE),
        "candidate_bitscore_fraction": _decimal_text(CANDIDATE_BITSCORE_FRACTION),
        "max_targets": max_targets,
        "blast_fetch_targets": blast_fetch_targets,
        "species_requires_all_candidates_exact_and_agreeing": True,
    }


def _file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_portable_search_provenance(path: str | Path) -> dict[str, object]:
    source = Path(path)
    try:
        source_sha256 = _file_sha256(source)
        payload = json.loads(source.read_text(encoding="utf-8"))
        portable = img_search_provenance.portable_provenance(
            payload, source_sha256
        )
    except (OSError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"Invalid IMG search provenance {source}: {error}") from error
    if _file_sha256(source) != source_sha256:
        raise ValueError("IMG search provenance changed while it was being loaded")
    return portable


def _classification_output_hashes(
    assignments: Sequence[Mapping[str, str]],
    outcomes: Sequence[Mapping[str, object]],
) -> dict[str, dict[str, str]]:
    assignment_digest = hashlib.sha256()
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=ASSIGNMENT_FIELDS, delimiter="\t", lineterminator="\n"
    )

    def update_assignment_digest() -> None:
        assignment_digest.update(buffer.getvalue().encode("utf-8"))
        buffer.seek(0)
        buffer.truncate(0)

    writer.writeheader()
    update_assignment_digest()
    for row in assignments:
        writer.writerow(row)
        update_assignment_digest()

    outcome_digest = hashlib.sha256()
    for row in outcomes:
        outcome_digest.update((_canonical_json(dict(row)) + "\n").encode("utf-8"))
    return {
        "assignments_tsv": {"sha256": assignment_digest.hexdigest()},
        "outcomes_jsonl": {"sha256": outcome_digest.hexdigest()},
    }


def lowest_common_ancestor(taxonomies: Iterable[Sequence[str]]) -> tuple[str, ...]:
    return _shared_lca(taxonomies)


def _is_pr2_species_path(
    taxonomy: Sequence[str], records: Sequence[TaxonomyRecord]
) -> bool:
    return (
        len(taxonomy) == 9
        and taxonomy[0] == "Eukaryota"
        and {record.taxonomy_source for record in records} == {"PR2"}
    )


def _common_value(values: Iterable[str]) -> str:
    return _shared_common_value(values)


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def centroid_name(value: str) -> str:
    """Return a display name without SILVA lineage text embedded in the header."""

    name = value.strip()
    if name.startswith("REF_SILVA_"):
        name = name.split(";", 1)[0]
    if not name or any(character in name for character in "\r\n\t|"):
        raise ValueError("cluster centroid has an invalid display name")
    return name


def _classify_cluster(
    cluster: Cluster,
    hits: Sequence[BlastHit],
    taxonomy_records: Mapping[str, TaxonomyRecord],
    *,
    max_targets: int,
    marker: str | None = None,
    calibration_strata: Mapping[str, CalibrationStratum] | None = None,
    propagation_rank_cap: int | None = None,
) -> dict[str, object]:
    eligible = [hit for hit in hits if hit.query_coverage >= MIN_QUERY_COVERAGE]
    base: dict[str, object] = {
        "cluster_id": cluster.cluster_id,
        "centroid": cluster.centroid,
        "img_member_count": len(cluster.img_members),
        "returned_hit_count": len(hits),
        "eligible_hit_count": len(eligible),
        "coverage_filtered_hit_count": len(hits) - len(eligible),
        "policy": classification_policy(max_targets=max_targets),
    }
    if not cluster.img_members:
        return {**base, "classification_status": "unclassified", "reason": "no_img_members"}
    if not hits:
        return {**base, "classification_status": "unclassified", "reason": "no_hits"}
    if not eligible:
        return {
            **base,
            "classification_status": "unclassified",
            "reason": "no_hits_at_minimum_query_coverage",
        }

    best_bit_score = eligible[0].bit_score
    threshold = best_bit_score * CANDIDATE_BITSCORE_FRACTION
    candidates = [hit for hit in eligible if hit.bit_score >= threshold]
    missing_taxonomy = sorted(
        hit.subject for hit in candidates if hit.subject not in taxonomy_records
    )
    candidate_evidence = [
        {
            "subject": hit.subject,
            "percent_identity": _decimal_text(hit.percent_identity),
            "alignment_length": hit.alignment_length,
            "query_length": hit.query_length,
            "subject_length": hit.subject_length,
            "query_coverage": _decimal_text(hit.query_coverage),
            "bit_score": _decimal_text(hit.bit_score),
            "taxonomy": ";".join(taxonomy_records[hit.subject].taxonomy)
            if hit.subject in taxonomy_records
            else "",
            "taxonomy_source": taxonomy_records[hit.subject].taxonomy_source
            if hit.subject in taxonomy_records
            else "",
            "cross_domain_conflict": taxonomy_records[hit.subject].cross_domain_conflict
            if hit.subject in taxonomy_records
            else False,
            "taxonomy_alternatives": json.loads(
                taxonomy_records[hit.subject].taxonomy_alternatives
            )
            if hit.subject in taxonomy_records
            and taxonomy_records[hit.subject].taxonomy_alternatives
            else [],
        }
        for hit in candidates
    ]
    base.update(
        {
            "best_bit_score": _decimal_text(best_bit_score),
            "candidate_bit_score_threshold": _decimal_text(threshold),
            "candidate_count": len(candidates),
            "candidates": candidate_evidence,
        }
    )
    if missing_taxonomy:
        return {
            **base,
            "classification_status": "unclassified",
            "reason": "missing_reference_taxonomy",
            "missing_taxonomy_subjects": missing_taxonomy,
        }

    records = [taxonomy_records[hit.subject] for hit in candidates]
    ambiguous_subjects = sorted(
        hit.subject
        for hit in candidates
        if taxonomy_records[hit.subject].cross_domain_conflict
    )
    if ambiguous_subjects:
        return {
            **base,
            "classification_status": "unclassified",
            "reason": "ambiguous_reference_taxonomy",
            "ambiguous_taxonomy_subjects": ambiguous_subjects,
        }
    # BLAST applies max_target_seqs before our query-coverage filter.  Therefore
    # the overflow sentinel must be read from the raw returned-hit boundary: a
    # low-coverage sentinel can still hide later, equally scoring eligible hits.
    truncated = len(hits) > max_targets and hits[max_targets].bit_score >= threshold
    if truncated:
        domains = sorted({record.domain for record in records if record.domain})
        taxonomy = (domains[0],) if len(domains) == 1 else ()
    else:
        taxonomy = lowest_common_ancestor(record.taxonomy for record in records)

    if not taxonomy:
        return {
            **base,
            "classification_status": "unclassified",
            "reason": "no_common_domain",
            "truncated": truncated,
        }

    calibration_rank_cap: int | None = None
    required_strata: list[str] = []
    if calibration_strata is not None:
        if not marker:
            raise ValueError("marker is required when calibration strata are supplied")
        required_strata = [
            f"{marker}|{source}|{taxonomy[0]}"
            for source in sorted({record.taxonomy_source for record in records})
        ]
        missing_strata = [key for key in required_strata if key not in calibration_strata]
        if missing_strata:
            raise ValueError(
                "taxonomy calibration lacks stratum "
                + ", ".join(repr(key) for key in missing_strata)
            )
        failed_strata = [
            {
                "key": key,
                "reason": calibration_strata[key].reason,
            }
            for key in required_strata
            if calibration_strata[key].status == "failed"
        ]
        if failed_strata:
            return {
                **base,
                "classification_status": "unclassified",
                "reason": "calibration_stratum_failed",
                "truncated": truncated,
                "failed_calibration_strata": failed_strata,
            }
        cap_values = [calibration_strata[key].rank_cap for key in required_strata]
        if any(value is None for value in cap_values):
            raise ValueError("calibrated taxonomy stratum lacks rank_cap")
        calibration_rank_cap = min(int(value) for value in cap_values)

    candidates_exact_and_agreeing = all(
        hit.percent_identity == Decimal("100")
        and hit.query_coverage == Decimal("100")
        and hit.alignment_length == hit.query_length == hit.subject_length
        for hit in candidates
    ) and len({record.taxonomy for record in records}) == 1
    species_guard_applied = False
    if _is_pr2_species_path(taxonomy, records) and not candidates_exact_and_agreeing:
        taxonomy = taxonomy[:-1]
        species_guard_applied = True

    taxonomy_sources = "+".join(sorted({record.taxonomy_source for record in records}))
    exact_pr2_species = _is_pr2_species_path(taxonomy, records) and candidates_exact_and_agreeing
    if calibration_rank_cap is not None and not exact_pr2_species:
        taxonomy = taxonomy[: calibration_rank_cap + 1]
    centroid_taxonomy = taxonomy
    if propagation_rank_cap is not None:
        if propagation_rank_cap < 0:
            raise ValueError("propagation_rank_cap must be non-negative")
        taxonomy = taxonomy[: propagation_rank_cap + 1]
    compartment = "" if truncated else _common_value(record.compartment for record in records)
    result = {
        **base,
        "classification_status": "classified",
        "reason": "",
        "taxonomy": ";".join(taxonomy),
        "taxonomy_source": taxonomy_sources,
        "centroid_taxonomy": ";".join(centroid_taxonomy),
        "centroid_taxonomy_source": taxonomy_sources,
        "domain": taxonomy[0],
        "compartment": compartment,
        "assignment_method": "updated_reference_cluster",
        "truncated": truncated,
        "species_guard_applied": species_guard_applied,
        "species_called": _is_pr2_species_path(taxonomy, records),
        "calibration_rank_cap": calibration_rank_cap,
        "propagation_rank_cap": propagation_rank_cap,
    }
    evidence_payload = dict(result)
    evidence_payload.pop("evidence_id", None)
    result["evidence_id"] = "IMGEV_" + hashlib.sha256(
        _canonical_json(evidence_payload).encode("utf-8")
    ).hexdigest()
    return result


def classify_clusters(
    clusters: Sequence[Cluster],
    hits_by_query: Mapping[str, Sequence[BlastHit]],
    taxonomy_records: Mapping[str, TaxonomyRecord],
    *,
    max_targets: int = BLAST_MAX_TARGETS,
    marker: str | None = None,
    calibration_strata: Mapping[str, CalibrationStratum] | None = None,
    propagation_rank_cap: int | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, object]], dict[str, object]]:
    if max_targets < 1:
        raise ValueError("max_targets must be positive")
    aliases: dict[str, Cluster] = {}
    for cluster in clusters:
        for alias in {cluster.cluster_id, cluster.centroid}:
            previous = aliases.get(alias)
            if previous is not None and previous != cluster:
                raise ValueError(f"Cluster query alias {alias!r} is ambiguous")
            aliases[alias] = cluster
    unmatched_queries = sorted(set(hits_by_query) - set(aliases))
    if unmatched_queries:
        preview = ", ".join(unmatched_queries[:5])
        raise ValueError(f"BLAST queries do not match a cluster centroid or ID: {preview}")

    assignments: list[dict[str, str]] = []
    outcomes: list[dict[str, object]] = []
    qc: Counter[str] = Counter()
    reasons: Counter[str] = Counter()
    domains: Counter[str] = Counter()
    for cluster in sorted(clusters, key=lambda value: (value.cluster_id, value.centroid)):
        matched_queries = [
            alias for alias in {cluster.cluster_id, cluster.centroid} if alias in hits_by_query
        ]
        if len(matched_queries) > 1:
            raise ValueError(
                f"Cluster {cluster.cluster_id!r} has BLAST hits under both cluster ID and centroid"
            )
        hits = tuple(hits_by_query[matched_queries[0]]) if matched_queries else ()
        outcome = _classify_cluster(
            cluster,
            hits,
            taxonomy_records,
            max_targets=max_targets,
            marker=marker,
            calibration_strata=calibration_strata,
            propagation_rank_cap=propagation_rank_cap,
        )
        outcomes.append(outcome)
        qc["clusters_total"] += 1
        qc["img_members_total"] += len(cluster.img_members)
        qc["hits_total"] += len(hits)
        qc["coverage_filtered_hits"] += int(outcome["coverage_filtered_hit_count"])
        qc["eligible_hits"] += int(outcome["eligible_hit_count"])
        qc["candidate_hits"] += int(outcome.get("candidate_count", 0))
        status = str(outcome["classification_status"])
        qc[f"clusters_{status}"] += 1
        if status == "unclassified":
            reason = str(outcome["reason"])
            reasons[reason] += 1
            if cluster.img_members:
                evidence_id = "IMGEV_" + hashlib.sha256(
                    _canonical_json(outcome).encode("utf-8")
                ).hexdigest()
                evidence = (
                    "centroid_blast_unclassified "
                    f"evidence_id={evidence_id} reason={reason}"
                )
                for member in cluster.img_members:
                    assignments.append(
                        {
                            "source_identifier": member,
                            "taxonomy": "Unclassified",
                            "taxonomy_source": "SILVA+PR2",
                            "assignment_method": "updated_reference_unclassified",
                            "evidence": evidence,
                            "compartment": "",
                            "cluster_id": cluster.cluster_id,
                            "centroid": cluster.centroid,
                            "centroid_name": centroid_name(cluster.centroid),
                            "centroid_taxonomy": "",
                            "centroid_taxonomy_source": "",
                            "evidence_id": evidence_id,
                        }
                    )
                    qc["assignments_written"] += 1
            continue
        if bool(outcome["truncated"]):
            qc["clusters_truncated"] += 1
        if bool(outcome["species_called"]):
            qc["clusters_species_called"] += 1
        domains[str(outcome["domain"])] += 1
        evidence = (
            "centroid_blast_lca "
            f"evidence_id={outcome['evidence_id']} candidates={outcome['candidate_count']} "
            f"best_bitscore={outcome['best_bit_score']} truncated="
            f"{str(outcome['truncated']).lower()}"
        )
        for member in cluster.img_members:
            assignments.append(
                {
                    "source_identifier": member,
                    "taxonomy": str(outcome["taxonomy"]),
                    "taxonomy_source": str(outcome["taxonomy_source"]),
                    "assignment_method": "updated_reference_cluster",
                    "evidence": evidence,
                    "compartment": str(outcome["compartment"]),
                    "cluster_id": cluster.cluster_id,
                    "centroid": cluster.centroid,
                    "centroid_name": centroid_name(cluster.centroid),
                    "centroid_taxonomy": str(outcome["centroid_taxonomy"]),
                    "centroid_taxonomy_source": str(
                        outcome["centroid_taxonomy_source"]
                    ),
                    "evidence_id": str(outcome["evidence_id"]),
                }
            )
            qc["assignments_written"] += 1

    qc_output: dict[str, object] = dict(sorted(qc.items()))
    qc_output["classified_clusters_by_domain"] = dict(sorted(domains.items()))
    qc_output["unclassified_clusters_by_reason"] = dict(sorted(reasons.items()))
    return assignments, outcomes, qc_output


def _atomic_text_writer(path: Path) -> Iterator[TextIO]:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", text=True
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        replace_and_fsync(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _write_outputs(
    assignments_tsv: Path,
    assignments_jsonl: Path,
    outcomes_jsonl: Path,
    qc_json: Path,
    assignments: Sequence[Mapping[str, str]],
    outcomes: Sequence[Mapping[str, object]],
    qc: Mapping[str, object],
) -> None:
    from contextlib import contextmanager

    writer_context = contextmanager(_atomic_text_writer)
    with writer_context(assignments_tsv) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=ASSIGNMENT_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(assignments)
    with writer_context(assignments_jsonl) as handle:
        for row in assignments:
            handle.write(_canonical_json(dict(row)) + "\n")
    with writer_context(outcomes_jsonl) as handle:
        for row in outcomes:
            handle.write(_canonical_json(dict(row)) + "\n")
    with writer_context(qc_json) as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--taxonomy", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--assignments-tsv", type=Path, required=True)
    parser.add_argument("--assignments-jsonl", type=Path, required=True)
    parser.add_argument("--outcomes-jsonl", type=Path, required=True)
    parser.add_argument("--qc-json", type=Path, required=True)
    parser.add_argument(
        "--blast-fetch-targets",
        type=int,
        required=True,
        help="Upstream BLAST max_target_seqs; must include one overflow sentinel.",
    )
    parser.add_argument("--marker", choices=("16S", "18S"))
    parser.add_argument("--calibration", type=Path)
    parser.add_argument("--search-provenance", type=Path)
    parser.add_argument(
        "--propagation-rank-cap",
        type=int,
        help="Maximum zero-based taxonomy rank propagated from a cluster centroid.",
    )
    return parser


def _validate_search_input_bindings(
    portable_search: Mapping[str, object], bindings: Mapping[str, Path]
) -> None:
    files = portable_search.get("files")
    if not isinstance(files, dict):
        raise ValueError("IMG search provenance lacks file bindings")
    for role, path in bindings.items():
        record = files.get(role)
        expected = record.get("sha256") if isinstance(record, dict) else None
        if not isinstance(expected, str) or _file_sha256(path) != expected:
            raise ValueError(
                f"IMG classifier input does not match search provenance: {role}"
            )


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.blast_fetch_targets != BLAST_FETCH_TARGETS:
        raise ValueError(
            "--blast-fetch-targets must be "
            f"{BLAST_FETCH_TARGETS} for the {BLAST_MAX_TARGETS}-hit policy"
        )
    if bool(args.marker) != bool(args.calibration):
        raise ValueError("--marker and --calibration must be supplied together")
    portable_search = (
        load_portable_search_provenance(args.search_provenance)
        if args.search_provenance
        else None
    )
    if portable_search is not None and portable_search["marker"] != args.marker:
        raise ValueError("IMG search provenance marker does not match --marker")
    search_input_bindings = {
        "blast_m8": args.blast,
        "clusters": args.clusters,
        "preferred_taxonomy": args.taxonomy,
    }
    if portable_search is not None:
        _validate_search_input_bindings(portable_search, search_input_bindings)
    calibration_sha256: str | None = None
    if args.calibration:
        calibration_sha256 = _file_sha256(args.calibration)
        calibration = load_calibration(args.calibration)
        if _file_sha256(args.calibration) != calibration_sha256:
            raise ValueError("taxonomy calibration changed while it was being loaded")
    else:
        calibration = None
    with args.blast.open(encoding="utf-8") as handle:
        hits_by_query = parse_blast_hits(handle)
    with args.clusters.open(newline="", encoding="utf-8") as handle:
        clusters = parse_clusters(handle)
    subjects = {
        hit.subject for hits in hits_by_query.values() for hit in hits
    }
    taxonomy_records = load_taxonomy(args.taxonomy, subjects)
    assignments, outcomes, qc = classify_clusters(
        clusters,
        hits_by_query,
        taxonomy_records,
        max_targets=BLAST_MAX_TARGETS,
        marker=args.marker,
        calibration_strata=calibration.strata if calibration else None,
        propagation_rank_cap=args.propagation_rank_cap,
    )
    if portable_search is not None:
        _validate_search_input_bindings(portable_search, search_input_bindings)
    qc["classification_binding"] = {
        "schema_version": 1,
        "calibration": {
            "sha256": calibration_sha256,
            "schema_version": 2,
        }
        if calibration
        else None,
        "marker": args.marker,
        "propagation_rank_cap": args.propagation_rank_cap,
        "policy": classification_policy(
            max_targets=BLAST_MAX_TARGETS,
            blast_fetch_targets=args.blast_fetch_targets,
        ),
        "outputs": _classification_output_hashes(assignments, outcomes),
    }
    if portable_search is not None:
        qc["classification_binding"]["search_provenance"] = portable_search
    _write_outputs(
        args.assignments_tsv,
        args.assignments_jsonl,
        args.outcomes_jsonl,
        args.qc_json,
        assignments,
        outcomes,
        qc,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
