"""Resolve runtime BLAST taxonomy within explicit evidence limits."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from taxonomy_utils import common_value, lowest_common_ancestor, taxonomy_path


CLASSIFICATION_POLICY = "runtime_lca_v1"
STRATUM_BASIS = "predicted_candidate_source_set_and_lca_domain"
MINIMUM_QUERY_COVERAGE = 0.80
CANDIDATE_SCORE_FRACTION = 0.98
NATIVE_TAXONOMY_SOURCES = frozenset({"PR2", "SILVA"})
MARKERS = frozenset({"16S", "18S"})


@dataclass(frozen=True)
class BlastHit:
    subject: str
    percent_identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bit_score: float
    query_length: int | None = None
    subject_length: int | None = None


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


@dataclass(frozen=True)
class RankRule:
    rank_index: int
    min_candidate_identity: float


@dataclass(frozen=True)
class RuntimeCalibration:
    reference_search_contract_sha256: str
    max_targets: int
    rank_caps: dict[str, int]
    rank_rules: dict[str, tuple[RankRule, ...]]
    strata: dict[str, dict[str, object]]


@dataclass(frozen=True)
class RuntimeTaxonomyDecision:
    eligible_hits: tuple[BlastHit, ...]
    candidate_hits: tuple[BlastHit, ...]
    candidate_taxonomy: str
    taxonomy: str
    taxonomy_source: str
    domain: str
    compartment: str
    assignment_method: str
    taxonomy_alternatives: str
    candidates_truncated: bool
    unknown_candidate_count: int


def _hit_key(hit: BlastHit) -> tuple[float, float, int, str]:
    return (
        -hit.bit_score,
        -hit.percent_identity,
        -hit.alignment_length,
        hit.subject,
    )


def _finite_float(value: str, label: str, location: str) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"Invalid BLAST {label} at {location}: {value!r}") from error
    if not math.isfinite(number):
        raise ValueError(f"Invalid BLAST {label} at {location}: {value!r}")
    return number


def _integer(value: str, label: str, location: str) -> int:
    try:
        return int(value)
    except ValueError as error:
        raise ValueError(f"Invalid BLAST {label} at {location}: {value!r}") from error


def _validate_hit(hit: BlastHit, location: str) -> None:
    if not hit.subject:
        raise ValueError(f"Empty BLAST subject at {location}")
    if not all(
        math.isfinite(value)
        for value in (hit.percent_identity, hit.evalue, hit.bit_score)
    ):
        raise ValueError(f"BLAST numeric field is not finite at {location}")
    if not 0.0 <= hit.percent_identity <= 100.0:
        raise ValueError(f"BLAST percent identity is outside 0..100 at {location}")
    if hit.alignment_length < 1:
        raise ValueError(f"BLAST alignment length must be positive at {location}")
    if hit.mismatches < 0 or hit.gap_opens < 0:
        raise ValueError(f"BLAST mismatch and gap counts must be non-negative at {location}")
    if min(hit.query_start, hit.query_end, hit.subject_start, hit.subject_end) < 1:
        raise ValueError(f"BLAST coordinates must be positive at {location}")
    if hit.evalue < 0.0 or hit.bit_score < 0.0:
        raise ValueError(f"BLAST scores must be non-negative at {location}")
    if hit.query_length is not None:
        if hit.query_length < 1 or max(hit.query_start, hit.query_end) > hit.query_length:
            raise ValueError(f"BLAST query coordinates exceed qlen at {location}")
    if hit.subject_length is not None:
        if hit.subject_length < 1 or max(hit.subject_start, hit.subject_end) > hit.subject_length:
            raise ValueError(f"BLAST subject coordinates exceed slen at {location}")


def load_blast_hits(m8_file: str | Path) -> dict[str, list[BlastHit]]:
    """Read BLAST 12-column or ``std qlen slen`` output."""

    hits_by_query: dict[str, list[BlastHit]] = {}
    with Path(m8_file).open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            location = f"{m8_file}:{line_number}"
            if len(fields) not in {12, 14}:
                raise ValueError(
                    f"Malformed BLAST m8 row at {location}: expected 12 or 14 "
                    f"fields, found {len(fields)}"
                )
            query = fields[0]
            if not query:
                raise ValueError(f"Empty BLAST query at {location}")
            hit = BlastHit(
                subject=fields[1],
                percent_identity=_finite_float(fields[2], "percent identity", location),
                alignment_length=_integer(fields[3], "alignment length", location),
                mismatches=_integer(fields[4], "mismatch count", location),
                gap_opens=_integer(fields[5], "gap-open count", location),
                query_start=_integer(fields[6], "query start", location),
                query_end=_integer(fields[7], "query end", location),
                subject_start=_integer(fields[8], "subject start", location),
                subject_end=_integer(fields[9], "subject end", location),
                evalue=_finite_float(fields[10], "E-value", location),
                bit_score=_finite_float(fields[11], "bit score", location),
                query_length=(
                    _integer(fields[12], "query length", location)
                    if len(fields) == 14
                    else None
                ),
                subject_length=(
                    _integer(fields[13], "subject length", location)
                    if len(fields) == 14
                    else None
                ),
            )
            _validate_hit(hit, location)
            hits_by_query.setdefault(query, []).append(hit)
    for query, hits in hits_by_query.items():
        unique_subjects: dict[str, BlastHit] = {}
        for hit in sorted(hits, key=_hit_key):
            unique_subjects.setdefault(hit.subject, hit)
        hits_by_query[query] = list(unique_subjects.values())
    return hits_by_query


def query_coverage(hit: BlastHit, fallback_query_length: int) -> float:
    query_length = hit.query_length or fallback_query_length
    if query_length < 1:
        raise ValueError("Query length must be positive")
    covered = abs(hit.query_end - hit.query_start) + 1
    if covered > query_length:
        raise ValueError(
            f"BLAST query span exceeds query length for subject {hit.subject}"
        )
    coverage = covered / query_length
    if not math.isfinite(coverage) or not 0.0 <= coverage <= 1.0:
        raise ValueError(f"Invalid BLAST query coverage for subject {hit.subject}")
    return coverage


def _canonical_source_set(records: Sequence[TaxonomyRecord]) -> str:
    sources = {
        source.strip()
        for record in records
        for source in record.taxonomy_source.split("+")
        if source.strip()
    }
    return "+".join(sorted(sources))


def _record_lineage(record: TaxonomyRecord) -> tuple[str, ...]:
    value = record.taxonomy
    if value and value != "Unclassified":
        try:
            return taxonomy_path(value)
        except ValueError:
            return ()
    if record.domain and record.domain not in {"ambiguous", "Unclassified"}:
        return (record.domain,)
    return ()


def merged_taxonomy_alternatives(records: Sequence[TaxonomyRecord]) -> str:
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


def _stratum_parts(key: str) -> tuple[str, str, str]:
    parts = key.split("|")
    if len(parts) != 3 or any(not part for part in parts):
        raise ValueError(f"Invalid runtime calibration stratum: {key!r}")
    marker, source_set, domain = parts
    if marker not in MARKERS:
        raise ValueError(f"Invalid runtime calibration marker: {marker!r}")
    sources = source_set.split("+")
    if any(not source for source in sources) or sources != sorted(set(sources)):
        raise ValueError(f"Runtime calibration source set is not canonical: {source_set!r}")
    return marker, source_set, domain


def _maximum_nonexact_rank(source_set: str) -> int:
    if source_set == "PR2":
        return 7
    if source_set == "SILVA":
        return 5
    return 0


def load_runtime_calibration(
    path: str | Path | None,
    reference_digest: str,
) -> RuntimeCalibration | None:
    """Load a runtime-query calibration file, or return ``None`` when absent."""

    if path is None:
        return None
    payload = json.loads(Path(path).read_text())
    if not isinstance(payload, dict):
        raise ValueError("Runtime calibration must be a JSON object")
    if payload.get("schema_version") != 3:
        raise ValueError("Runtime calibration must use schema version 3")
    if payload.get("calibration_use_case") != "runtime_query":
        raise ValueError("Runtime calibration use case must be runtime_query")
    if payload.get("classification_policy") != CLASSIFICATION_POLICY:
        raise ValueError("Runtime calibration classification policy does not match")
    if payload.get("stratum_basis") != STRATUM_BASIS:
        raise ValueError("Runtime calibration stratum basis does not match")
    if payload.get("reference_search_contract_sha256") != reference_digest:
        raise ValueError("Runtime calibration reference search contract does not match")
    calibrated_max_targets = payload.get("max_targets")
    if (
        isinstance(calibrated_max_targets, bool)
        or not isinstance(calibrated_max_targets, int)
        or calibrated_max_targets < 1
    ):
        raise ValueError("Runtime calibration max_targets must be a positive integer")

    raw_caps = payload.get("rank_caps")
    raw_strata = payload.get("strata")
    raw_rules = payload.get("rank_rules", {})
    if not isinstance(raw_caps, dict) or not isinstance(raw_strata, dict):
        raise ValueError("Runtime calibration rank_caps and strata must be objects")
    if not isinstance(raw_rules, dict):
        raise ValueError("Runtime calibration rank_rules must be an object")
    expected_status = "calibrated" if raw_caps else "failed"
    if payload.get("status") != expected_status:
        raise ValueError(
            f"Runtime calibration root status must be {expected_status} "
            f"when rank_caps is {'nonempty' if raw_caps else 'empty'}"
        )

    rank_caps: dict[str, int] = {}
    strata: dict[str, dict[str, object]] = {}
    for key, raw in raw_strata.items():
        if not isinstance(key, str) or not isinstance(raw, dict):
            raise ValueError("Runtime calibration contains an invalid stratum")
        _marker, source_set, _domain = _stratum_parts(key)
        status = raw.get("status")
        reason = raw.get("reason")
        rank_cap = raw.get("rank_cap")
        if status not in {"calibrated", "failed"} or not isinstance(reason, str):
            raise ValueError(f"Runtime calibration stratum has invalid status: {key}")
        if status == "calibrated":
            if reason or not isinstance(rank_cap, int) or not 0 <= rank_cap <= 8:
                raise ValueError(f"Runtime calibration stratum is invalid: {key}")
            if rank_cap > _maximum_nonexact_rank(source_set):
                raise ValueError(
                    f"Runtime calibration rank cap can assign an unsupported "
                    f"nonexact rank: {key}"
                )
        elif not reason or rank_cap is not None:
            raise ValueError(f"Failed runtime calibration stratum is invalid: {key}")
        strata[key] = dict(raw)

    for key, cap in raw_caps.items():
        if not isinstance(key, str) or not isinstance(cap, int) or not 0 <= cap <= 8:
            raise ValueError("Runtime calibration contains an invalid rank cap")
        _marker, source_set, _domain = _stratum_parts(key)
        if cap > _maximum_nonexact_rank(source_set):
            raise ValueError(
                f"Runtime calibration rank cap can assign an unsupported "
                f"nonexact rank: {key}"
            )
        stratum = strata.get(key)
        if (
            stratum is None
            or stratum["status"] != "calibrated"
            or stratum["rank_cap"] != cap
        ):
            raise ValueError(f"Runtime calibration rank cap does not match stratum: {key}")
        rank_caps[key] = cap
    calibrated_keys = {
        key for key, stratum in strata.items() if stratum["status"] == "calibrated"
    }
    if calibrated_keys != set(rank_caps):
        raise ValueError("Runtime calibration calibrated strata and rank caps differ")

    rank_rules: dict[str, tuple[RankRule, ...]] = {}
    for key, raw_rule_list in raw_rules.items():
        if key not in rank_caps or not isinstance(raw_rule_list, list):
            raise ValueError(f"Runtime calibration rank rules have an invalid stratum: {key}")
        parsed: list[RankRule] = []
        for expected_index, raw_rule in enumerate(raw_rule_list):
            if not isinstance(raw_rule, dict):
                raise ValueError(f"Runtime calibration rank rule is invalid: {key}")
            rank_index = raw_rule.get("rank_index")
            threshold = raw_rule.get("min_candidate_identity")
            if rank_index != expected_index or not isinstance(threshold, (int, float)):
                raise ValueError(f"Runtime calibration rank rules are not contiguous: {key}")
            threshold = float(threshold)
            if not math.isfinite(threshold) or not 0.0 <= threshold <= 100.0:
                raise ValueError(f"Runtime calibration identity threshold is invalid: {key}")
            if parsed and threshold < parsed[-1].min_candidate_identity:
                raise ValueError(
                    f"Runtime calibration identity thresholds decrease with rank: {key}"
                )
            parsed.append(RankRule(rank_index, threshold))
        if len(parsed) != rank_caps[key] + 1:
            raise ValueError(f"Runtime calibration rank rules do not reach rank cap: {key}")
        rank_rules[key] = tuple(parsed)

    return RuntimeCalibration(
        reference_search_contract_sha256=reference_digest,
        max_targets=calibrated_max_targets,
        rank_caps=rank_caps,
        rank_rules=rank_rules,
        strata=strata,
    )


def _unclassified(
    *,
    eligible_hits: Sequence[BlastHit],
    candidate_hits: Sequence[BlastHit],
    candidate_taxonomy: str,
    assignment_method: str,
    taxonomy_alternatives: str = "[]",
    candidates_truncated: bool = False,
    domain: str = "Unclassified",
    unknown_candidate_count: int = 0,
) -> RuntimeTaxonomyDecision:
    return RuntimeTaxonomyDecision(
        eligible_hits=tuple(eligible_hits),
        candidate_hits=tuple(candidate_hits),
        candidate_taxonomy=candidate_taxonomy,
        taxonomy="" if domain == "ambiguous" else "Unclassified",
        taxonomy_source="",
        domain=domain,
        compartment="",
        assignment_method=assignment_method,
        taxonomy_alternatives=taxonomy_alternatives,
        candidates_truncated=candidates_truncated,
        unknown_candidate_count=unknown_candidate_count,
    )


def _is_exact_native(
    hit: BlastHit,
    record: TaxonomyRecord,
) -> bool:
    if (
        record.assignment_method != "native"
        or record.taxonomy_source not in NATIVE_TAXONOMY_SOURCES
        or not record.taxonomy
        or record.cross_domain_conflict
    ):
        return False
    if hit.query_length is None or hit.subject_length is None:
        return False
    return (
        hit.percent_identity == 100.0
        and hit.mismatches == 0
        and hit.gap_opens == 0
        and hit.alignment_length == hit.query_length == hit.subject_length
        and min(hit.query_start, hit.query_end) == 1
        and max(hit.query_start, hit.query_end) == hit.query_length
        and min(hit.subject_start, hit.subject_end) == 1
        and max(hit.subject_start, hit.subject_end) == hit.subject_length
    )


def resolve_runtime_taxonomy(
    hits: Sequence[BlastHit],
    taxonomy_records: Mapping[str, TaxonomyRecord],
    *,
    query_length: int,
    marker: str,
    calibration: RuntimeCalibration | None,
    max_targets: int,
) -> RuntimeTaxonomyDecision:
    """Resolve one query from eligible, near-best, unique BLAST subjects."""

    if marker not in MARKERS:
        raise ValueError(f"Unsupported runtime taxonomy marker: {marker}")
    if query_length < 1:
        raise ValueError("Query length must be positive")
    if max_targets < 1:
        raise ValueError("max_targets must be positive")
    if calibration is not None and max_targets != calibration.max_targets:
        raise ValueError(
            "Runtime taxonomy max_targets differs from the calibrated value: "
            f"caller={max_targets}, calibration={calibration.max_targets}"
        )

    unique: dict[str, BlastHit] = {}
    for hit in sorted(hits, key=_hit_key):
        _validate_hit(hit, f"subject {hit.subject}")
        if hit.query_length is not None and hit.query_length != query_length:
            raise ValueError(
                f"BLAST qlen differs from the query sequence for subject {hit.subject}"
            )
        query_coverage(hit, query_length)
        unique.setdefault(hit.subject, hit)
    ranked = list(unique.values())
    policy_hits = ranked[:max_targets]
    eligible_hits = [
        hit
        for hit in ranked
        if query_coverage(hit, query_length) >= MINIMUM_QUERY_COVERAGE
    ]
    policy_eligible_hits = [
        hit
        for hit in policy_hits
        if query_coverage(hit, query_length) >= MINIMUM_QUERY_COVERAGE
    ]
    if not policy_eligible_hits:
        return _unclassified(
            eligible_hits=eligible_hits,
            candidate_hits=(),
            candidate_taxonomy="Unclassified",
            assignment_method=(
                "candidate_set_truncated"
                if len(ranked) > max_targets
                else "no_eligible_blast_hits"
            ),
            candidates_truncated=len(ranked) > max_targets,
        )

    score_floor = policy_eligible_hits[0].bit_score * CANDIDATE_SCORE_FRACTION
    candidate_hits = [
        hit for hit in policy_eligible_hits if hit.bit_score >= score_floor
    ]
    candidates_truncated = (
        len(ranked) > max_targets and ranked[max_targets].bit_score >= score_floor
    )
    candidate_records = [
        taxonomy_records[hit.subject]
        for hit in candidate_hits
        if hit.subject in taxonomy_records
    ]
    alternatives = merged_taxonomy_alternatives(candidate_records)
    known = [
        (hit, taxonomy_records[hit.subject], _record_lineage(taxonomy_records[hit.subject]))
        for hit in candidate_hits
        if hit.subject in taxonomy_records
        and _record_lineage(taxonomy_records[hit.subject])
    ]
    known_lineages = [lineage for _hit, _record, lineage in known]
    unknown_candidate_count = sum(
        hit.subject not in taxonomy_records
        or (
            not _record_lineage(taxonomy_records[hit.subject])
            and not taxonomy_records[hit.subject].cross_domain_conflict
            and taxonomy_records[hit.subject].domain != "ambiguous"
        )
        for hit in candidate_hits
    )
    candidate_lca = lowest_common_ancestor(known_lineages)
    candidate_taxonomy = ";".join(candidate_lca) or "Unclassified"
    explicit_ambiguity = any(
        record.cross_domain_conflict or record.domain == "ambiguous"
        for record in candidate_records
    )
    known_domains = {lineage[0] for lineage in known_lineages}
    if explicit_ambiguity or len(known_domains) > 1:
        return _unclassified(
            eligible_hits=eligible_hits,
            candidate_hits=candidate_hits,
            candidate_taxonomy="",
            assignment_method="cross_domain_ambiguous_candidates",
            taxonomy_alternatives=alternatives,
            candidates_truncated=candidates_truncated,
            domain="ambiguous",
            unknown_candidate_count=unknown_candidate_count,
        )
    if candidates_truncated:
        return _unclassified(
            eligible_hits=eligible_hits,
            candidate_hits=candidate_hits,
            candidate_taxonomy=candidate_taxonomy,
            assignment_method="candidate_set_truncated",
            taxonomy_alternatives=alternatives,
            candidates_truncated=True,
            unknown_candidate_count=unknown_candidate_count,
        )
    if not candidate_lca:
        return _unclassified(
            eligible_hits=eligible_hits,
            candidate_hits=candidate_hits,
            candidate_taxonomy="Unclassified",
            assignment_method="no_known_candidate_taxonomy",
            taxonomy_alternatives=alternatives,
            unknown_candidate_count=unknown_candidate_count,
        )

    all_candidates_known = len(known) == len(candidate_hits)
    exact_lineages = {lineage for _hit, _record, lineage in known}
    if (
        all_candidates_known
        and len(exact_lineages) == 1
        and all(_is_exact_native(hit, record) for hit, record, _lineage in known)
    ):
        lineage = next(iter(exact_lineages))
        sources = _canonical_source_set([record for _hit, record, _lineage in known])
        compartment = common_value(
            (record.compartment for _hit, record, _lineage in known)
        )
        return RuntimeTaxonomyDecision(
            eligible_hits=tuple(eligible_hits),
            candidate_hits=tuple(candidate_hits),
            candidate_taxonomy=";".join(lineage),
            taxonomy=";".join(lineage),
            taxonomy_source=sources,
            domain=lineage[0],
            compartment=compartment,
            assignment_method="exact_native_match",
            taxonomy_alternatives=alternatives,
            candidates_truncated=False,
            unknown_candidate_count=0,
        )

    if calibration is None:
        return _unclassified(
            eligible_hits=eligible_hits,
            candidate_hits=candidate_hits,
            candidate_taxonomy=candidate_taxonomy,
            assignment_method="no_runtime_calibration",
            taxonomy_alternatives=alternatives,
            unknown_candidate_count=unknown_candidate_count,
        )
    source_set = _canonical_source_set([record for _hit, record, _lineage in known])
    stratum_key = f"{marker}|{source_set}|{candidate_lca[0]}"
    cap = calibration.rank_caps.get(stratum_key)
    if cap is None:
        stratum = calibration.strata.get(stratum_key)
        reason = "runtime_calibration_failed"
        if stratum is None:
            reason = "runtime_calibration_missing_stratum"
        return _unclassified(
            eligible_hits=eligible_hits,
            candidate_hits=candidate_hits,
            candidate_taxonomy=candidate_taxonomy,
            assignment_method=reason,
            taxonomy_alternatives=alternatives,
            unknown_candidate_count=unknown_candidate_count,
        )

    supported_cap = min(
        cap,
        len(candidate_lca) - 1,
        _maximum_nonexact_rank(source_set),
    )
    rules = calibration.rank_rules.get(stratum_key)
    if rules:
        minimum_identity = min(hit.percent_identity for hit in candidate_hits)
        supported_cap = -1
        for rule in rules:
            if rule.rank_index >= len(candidate_lca):
                break
            if minimum_identity < rule.min_candidate_identity:
                break
            supported_cap = rule.rank_index
        if supported_cap < 0:
            return _unclassified(
                eligible_hits=eligible_hits,
                candidate_hits=candidate_hits,
                candidate_taxonomy=candidate_taxonomy,
                assignment_method="runtime_calibration_identity_below_domain",
                taxonomy_alternatives=alternatives,
                unknown_candidate_count=unknown_candidate_count,
            )

    selected = candidate_lca[: supported_cap + 1]
    known_records = [record for _hit, record, _lineage in known]
    return RuntimeTaxonomyDecision(
        eligible_hits=tuple(eligible_hits),
        candidate_hits=tuple(candidate_hits),
        candidate_taxonomy=candidate_taxonomy,
        taxonomy=";".join(selected),
        taxonomy_source=source_set,
        domain=selected[0],
        compartment=common_value(record.compartment for record in known_records),
        assignment_method="runtime_calibrated_lca",
        taxonomy_alternatives=alternatives,
        candidates_truncated=False,
        unknown_candidate_count=unknown_candidate_count,
    )
