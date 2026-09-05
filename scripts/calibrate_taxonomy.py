#!/usr/bin/env python3
"""Calibrate IMG taxonomy rank caps by deterministic leave-one-reference-out BLAST."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import duckdb

import classify_img_clusters as classifier
import database_manager as manager
from atomic_io import replace_and_fsync


RANKS = {
    "SILVA": ("domain", "phylum", "class", "order", "family", "genus"),
    "PR2": (
        "domain",
        "supergroup",
        "division",
        "subdivision",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ),
}
MIN_CALIBRATION_CALLS = 100
HIGHER_RANK_MIN_PRECISION = 0.95
GENUS_MIN_PRECISION = 0.98
CALIBRATION_FETCH_TARGETS = classifier.BLAST_FETCH_TARGETS + 1
SEARCH_PROVENANCE_NAME = "search_provenance.json"
SEARCH_PROVENANCE_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class CalibrationPrediction:
    taxonomy: tuple[str, ...]
    candidate_sources: tuple[str, ...]
    candidate_count: int
    truncated: bool
    species_guard_applied: bool


def _run(command: list[str]) -> None:
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as error:
        raise RuntimeError(f"required executable not found: {command[0]}") from error
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"command failed with exit code {error.returncode}: {' '.join(command)}") from error


def select_calibration_rows(
    profile_directory: str | Path, samples_per_stratum: int
) -> list[dict[str, str]]:
    if samples_per_stratum < MIN_CALIBRATION_CALLS:
        raise ValueError(
            f"samples_per_stratum must be at least {MIN_CALIBRATION_CALLS}"
        )
    profile = Path(profile_directory)
    connection = duckdb.connect(":memory:")
    try:
        thread_text = (
            os.environ.get("SLURM_CPUS_PER_TASK")
            or os.environ.get("OMP_NUM_THREADS")
            or "1"
        )
        try:
            thread_limit = int(thread_text)
        except ValueError as error:
            raise ValueError(
                f"DuckDB thread limit must be an integer, found {thread_text!r}"
            ) from error
        if thread_limit < 1:
            raise ValueError("DuckDB thread limit must be positive")
        connection.execute("SET threads = ?", [thread_limit])
        rows = connection.execute(
            """
            WITH expanded AS (
                SELECT DISTINCT
                    p.sequence_id,
                    p.taxonomy,
                    p.taxonomy_source,
                    p.domain,
                    trim(marker_token.marker) AS marker
                FROM read_parquet(?) AS p
                INNER JOIN read_parquet(?) AS s USING (sequence_id)
                CROSS JOIN unnest(string_split(s.markers, ';')) AS marker_token(marker)
                WHERE p.taxonomy_source IN ('SILVA', 'PR2')
                  AND p.domain IN ('Bacteria', 'Archaea', 'Eukaryota')
                  AND trim(p.taxonomy) <> ''
                  AND p.cross_domain_conflict = FALSE
                  AND trim(marker_token.marker) IN ('16S', '18S')
            ), ranked AS (
                SELECT *, row_number() OVER (
                    PARTITION BY marker, taxonomy_source, domain
                    ORDER BY sequence_id
                ) AS sample_rank
                FROM expanded
            )
            SELECT sequence_id, taxonomy, taxonomy_source, domain, marker
            FROM ranked
            WHERE sample_rank <= ?
            ORDER BY marker, taxonomy_source, domain, sequence_id
            """,
            [
                str(profile / "tables" / "preferred_taxonomy.parquet"),
                str(profile / "tables" / "sequences.parquet"),
                samples_per_stratum,
            ],
        ).fetchall()
    finally:
        connection.close()
    return [
        {
            "sequence_id": str(row[0]),
            "taxonomy": str(row[1]),
            "taxonomy_source": str(row[2]),
            "domain": str(row[3]),
            "marker": str(row[4]),
        }
        for row in rows
    ]


def _query_ids_text(rows: Sequence[Mapping[str, str]]) -> str:
    sequence_ids = [row["sequence_id"] for row in rows]
    if len(sequence_ids) != len(set(sequence_ids)):
        raise RuntimeError("calibration query IDs are not unique")
    return "".join(f"{sequence_id}\n" for sequence_id in sequence_ids)


def _truth_tsv_text(rows: Sequence[Mapping[str, str]]) -> str:
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(
        handle,
        fieldnames=("sequence_id", "taxonomy", "taxonomy_source", "domain", "marker"),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _profile_manifest_sha256(profile: Path) -> str:
    manifest = profile / manager.MANIFEST_NAME
    if not manifest.is_file():
        raise RuntimeError(f"profile manifest is missing: {manifest}")
    return _sha256_file(manifest)


def _search_commands(
    profile: Path,
    marker: str,
    ids: Path,
    fasta: Path,
    blast_output: Path,
    threads: int,
) -> dict[str, list[str]]:
    database = str((profile / "blast" / marker).resolve())
    ids_path = str(ids.resolve())
    fasta_path = str(fasta.resolve())
    blast_output_path = str(blast_output.resolve())
    return {
        "blastdbcmd": [
            "blastdbcmd",
            "-db",
            database,
            "-entry_batch",
            ids_path,
            "-outfmt",
            "%f",
            "-out",
            fasta_path,
        ],
        "blastn": [
            "blastn",
            "-task",
            "blastn",
            "-query",
            fasta_path,
            "-db",
            database,
            "-evalue",
            "1e-5",
            "-max_target_seqs",
            str(CALIBRATION_FETCH_TARGETS),
            "-max_hsps",
            "1",
            "-num_threads",
            str(threads),
            "-outfmt",
            "6 qseqid sseqid pident length qlen slen qcovs bitscore",
            "-out",
            blast_output_path,
        ],
    }


def _search_provenance_payload(
    profile: Path,
    marker: str,
    ids: Path,
    truth: Path,
    fasta: Path,
    blast_output: Path,
    threads: int,
) -> dict[str, object]:
    files = (ids, truth, fasta, blast_output)
    missing = [path.name for path in files if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(
            f"BLAST for {marker} did not produce complete non-empty outputs: "
            f"{', '.join(missing)}"
        )
    return {
        "schema_version": SEARCH_PROVENANCE_SCHEMA_VERSION,
        "status": "complete",
        "profile_manifest_sha256": _profile_manifest_sha256(profile),
        "files": {path.name: _sha256_file(path) for path in files},
        "commands": _search_commands(
            profile, marker, ids, fasta, blast_output, threads
        ),
    }


def _write_search_provenance(
    profile: Path,
    marker: str,
    marker_directory: Path,
    ids: Path,
    truth: Path,
    fasta: Path,
    blast_output: Path,
    threads: int,
) -> None:
    payload = _search_provenance_payload(
        profile, marker, ids, truth, fasta, blast_output, threads
    )
    _write_json_durably(marker_directory / SEARCH_PROVENANCE_NAME, payload)


def _validate_complete_blast_results(
    blast_output: Path,
    marker_rows: Sequence[Mapping[str, str]],
    marker: str,
    context: str,
) -> None:
    if not blast_output.is_file() or blast_output.stat().st_size == 0:
        raise RuntimeError(f"{context} {marker}: BLAST output is missing or empty")
    with blast_output.open(encoding="utf-8") as handle:
        retained_hits = classifier.parse_blast_hits(handle)
    expected = {row["sequence_id"] for row in marker_rows}
    observed = set(retained_hits)
    missing_queries = sorted(expected - observed)
    extra_queries = sorted(observed - expected)
    if missing_queries or extra_queries:
        raise RuntimeError(
            f"{context} {marker}: BLAST query coverage differs "
            f"(missing={len(missing_queries)}, extra={len(extra_queries)})"
        )
    missing_self_hits = sorted(
        query
        for query in expected
        if not any(hit.subject == query for hit in retained_hits[query])
    )
    if missing_self_hits:
        raise RuntimeError(
            f"{context} {marker}: BLAST output lacks self hits for "
            f"{len(missing_self_hits)} queries"
        )


def _validate_reusable_blast(
    profile: Path,
    marker: str,
    marker_directory: Path,
    marker_rows: Sequence[Mapping[str, str]],
) -> Path:
    ids = marker_directory / "query_ids.txt"
    truth = marker_directory / "truth.tsv"
    fasta = marker_directory / "queries.fna"
    blast_output = marker_directory / "leave_one_out.m8"
    provenance_path = marker_directory / SEARCH_PROVENANCE_NAME
    required = (ids, truth, fasta, blast_output, provenance_path)
    missing = [path.name for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(
            f"cannot reuse BLAST for {marker}: missing or empty files: {', '.join(missing)}"
        )

    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as error:
        raise RuntimeError(
            f"cannot reuse BLAST for {marker}: invalid {SEARCH_PROVENANCE_NAME}"
        ) from error
    if not isinstance(provenance, dict) or set(provenance) != {
        "schema_version",
        "status",
        "profile_manifest_sha256",
        "files",
        "commands",
    }:
        raise RuntimeError(
            f"cannot reuse BLAST for {marker}: invalid provenance schema"
        )
    if (
        provenance["schema_version"] != SEARCH_PROVENANCE_SCHEMA_VERSION
        or provenance["status"] != "complete"
    ):
        raise RuntimeError(
            f"cannot reuse BLAST for {marker}: provenance is not a complete "
            f"schema-{SEARCH_PROVENANCE_SCHEMA_VERSION} search"
        )
    current_manifest_sha256 = _profile_manifest_sha256(profile)
    if provenance["profile_manifest_sha256"] != current_manifest_sha256:
        raise RuntimeError(
            f"cannot reuse BLAST for {marker}: profile manifest has changed"
        )
    expected_hashes = {
        path.name: _sha256_file(path) for path in (ids, truth, fasta, blast_output)
    }
    if provenance["files"] != expected_hashes:
        raise RuntimeError(
            f"cannot reuse BLAST for {marker}: retained search output hash mismatch"
        )
    commands = provenance.get("commands")
    blastn_command = commands.get("blastn") if isinstance(commands, dict) else None
    try:
        thread_index = blastn_command.index("-num_threads")
        recorded_threads = int(blastn_command[thread_index + 1])
    except (AttributeError, IndexError, TypeError, ValueError) as error:
        raise RuntimeError(
            f"cannot reuse BLAST for {marker}: search command lacks a valid thread count"
        ) from error
    if recorded_threads < 1:
        raise RuntimeError(
            f"cannot reuse BLAST for {marker}: search command has a non-positive thread count"
        )
    expected_commands = _search_commands(
        profile, marker, ids, fasta, blast_output, recorded_threads
    )
    if provenance["commands"] != expected_commands:
        raise RuntimeError(
            f"cannot reuse BLAST for {marker}: search command contract mismatch"
        )

    expected_ids = _query_ids_text(marker_rows)
    if ids.read_text(encoding="ascii") != expected_ids:
        raise RuntimeError(
            f"cannot reuse BLAST for {marker}: query_ids.txt does not match "
            "deterministic selection"
        )
    if truth.read_text(encoding="utf-8") != _truth_tsv_text(marker_rows):
        raise RuntimeError(
            f"cannot reuse BLAST for {marker}: truth.tsv does not match deterministic selection"
        )

    with tempfile.TemporaryDirectory(
        prefix="reuse-validation-", dir=marker_directory
    ) as temporary:
        regenerated = Path(temporary) / "queries.fna"
        _run(
            [
                "blastdbcmd",
                "-db",
                str(profile / "blast" / marker),
                "-entry_batch",
                str(ids),
                "-outfmt",
                "%f",
                "-out",
                str(regenerated),
            ]
        )
        if not regenerated.is_file() or regenerated.stat().st_size == 0:
            raise RuntimeError(
                f"cannot reuse BLAST for {marker}: regenerated query FASTA is empty"
            )
        if regenerated.read_bytes() != fasta.read_bytes():
            raise RuntimeError(
                f"cannot reuse BLAST for {marker}: regenerated query FASTA differs "
                "from retained FASTA"
            )

    _validate_complete_blast_results(
        blast_output,
        marker_rows,
        marker,
        "cannot reuse BLAST for",
    )
    return blast_output


def _write_queries(
    profile: Path,
    output: Path,
    rows: Sequence[Mapping[str, str]],
    threads: int,
    *,
    reuse_existing_blast: bool = False,
    reuse_blast_directory: str | Path | None = None,
) -> dict[str, Path]:
    if reuse_existing_blast and reuse_blast_directory is not None:
        raise ValueError(
            "reuse_existing_blast and reuse_blast_directory are mutually exclusive"
        )
    retained_root = (
        Path(reuse_blast_directory)
        if reuse_blast_directory is not None
        else output if reuse_existing_blast else None
    )
    by_marker: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        by_marker[row["marker"]].append(row)
    blast_outputs: dict[str, Path] = {}
    for marker, marker_rows in sorted(by_marker.items()):
        marker_directory = output / marker
        if retained_root is not None:
            blast_outputs[marker] = _validate_reusable_blast(
                profile, marker, retained_root / marker, marker_rows
            )
            continue
        marker_directory.mkdir(parents=True, exist_ok=True)
        ids = marker_directory / "query_ids.txt"
        truth = marker_directory / "truth.tsv"
        fasta = marker_directory / "queries.fna"
        blast_output = marker_directory / "leave_one_out.m8"
        provenance = marker_directory / SEARCH_PROVENANCE_NAME
        provenance.unlink(missing_ok=True)
        ids.write_text(_query_ids_text(marker_rows), encoding="ascii")
        truth.write_text(_truth_tsv_text(marker_rows), encoding="utf-8")
        commands = _search_commands(
            profile, marker, ids, fasta, blast_output, threads
        )
        _run(commands["blastdbcmd"])
        _run(commands["blastn"])
        _validate_complete_blast_results(
            blast_output,
            marker_rows,
            marker,
            "fresh BLAST search for",
        )
        _write_search_provenance(
            profile,
            marker,
            marker_directory,
            ids,
            truth,
            fasta,
            blast_output,
            threads,
        )
        blast_outputs[marker] = blast_output
    return blast_outputs


def wilson_lower_bound(correct: int, total: int, z: float = 1.959963984540054) -> float:
    if total == 0:
        return 0.0
    proportion = correct / total
    denominator = 1.0 + z * z / total
    center = proportion + z * z / (2.0 * total)
    spread = z * math.sqrt(
        proportion * (1.0 - proportion) / total + z * z / (4.0 * total * total)
    )
    return (center - spread) / denominator


def _prediction(
    query: str,
    hits: Sequence[classifier.BlastHit],
    taxonomy: Mapping[str, classifier.TaxonomyRecord],
) -> CalibrationPrediction:
    non_self = [hit for hit in hits if hit.subject != query]
    eligible, threshold, candidates = classifier.select_candidate_hits(non_self)
    if not eligible:
        return CalibrationPrediction((), (), 0, False, False)
    assert threshold is not None
    missing_subjects = sorted(
        {hit.subject for hit in candidates if hit.subject not in taxonomy}
    )
    if missing_subjects:
        preview = ", ".join(missing_subjects[:5])
        suffix = "" if len(missing_subjects) <= 5 else ", ..."
        raise RuntimeError(
            f"taxonomy is missing for {len(missing_subjects)} BLAST subject(s) "
            f"for query {query}: {preview}{suffix}"
        )
    records = [taxonomy[hit.subject] for hit in candidates]
    sources = tuple(sorted({record.taxonomy_source for record in records}))
    if any(record.cross_domain_conflict for record in records):
        return CalibrationPrediction((), sources, len(candidates), False, False)
    prediction, truncated, _exact, species_guard_applied = (
        classifier.resolve_candidate_taxonomy(
            non_self,
            candidates,
            records,
            threshold,
            max_targets=classifier.BLAST_MAX_TARGETS,
        )
    )
    return CalibrationPrediction(
        prediction,
        sources,
        len(candidates),
        truncated,
        species_guard_applied,
    )


def evaluate_calibration(
    profile_directory: str | Path,
    rows: Sequence[Mapping[str, str]],
    blast_outputs: Mapping[str, Path],
    *,
    samples_per_stratum_requested: int | None = None,
) -> dict[str, object]:
    if not rows:
        raise RuntimeError("calibration selected no rows")
    profile = Path(profile_directory)
    hits_by_marker = {}
    all_subjects: set[str] = set()
    for marker, path in blast_outputs.items():
        with path.open(encoding="utf-8") as handle:
            marker_hits = classifier.parse_blast_hits(handle)
        hits_by_marker[marker] = marker_hits
        all_subjects.update(
            hit.subject for hits in marker_hits.values() for hit in hits
        )
    taxonomy = classifier.load_taxonomy(
        profile / "tables" / "preferred_taxonomy.parquet", all_subjects
    )

    calls: Counter[tuple[str, str, str, int]] = Counter()
    truth_known_calls: Counter[tuple[str, str, str, int]] = Counter()
    truth_unknown_calls: Counter[tuple[str, str, str, int]] = Counter()
    correct: Counter[tuple[str, str, str, int]] = Counter()
    called_truth_mixtures: dict[
        tuple[str, str, str, int], Counter[str]
    ] = defaultdict(Counter)
    true_stratum_counts: Counter[tuple[str, str, str]] = Counter()
    marker_counts: Counter[str] = Counter()
    marker_truth_mixtures: dict[str, Counter[str]] = defaultdict(Counter)
    predicted_strata: set[tuple[str, str, str]] = set()
    for row in rows:
        marker = row["marker"]
        truth_source = row["taxonomy_source"]
        truth_domain = row["domain"]
        true_stratum = (marker, truth_source, truth_domain)
        true_class_key = "|".join(true_stratum)
        true_stratum_counts[true_stratum] += 1
        marker_counts[marker] += 1
        marker_truth_mixtures[marker][true_class_key] += 1
        truth = tuple(row["taxonomy"].split(";"))
        prediction = _prediction(
            row["sequence_id"],
            hits_by_marker.get(marker, {}).get(row["sequence_id"], ()),
            taxonomy,
        )
        if not prediction.taxonomy:
            continue
        predicted_domain = prediction.taxonomy[0]
        route_keys = classifier.calibration_stratum_keys(
            marker, prediction.candidate_sources, predicted_domain
        )
        for route_key in route_keys:
            route = tuple(route_key.split("|", 2))
            route_source = route[1]
            if route_source not in RANKS:
                raise RuntimeError(
                    f"unsupported predicted taxonomy source in calibration: {route_source!r}"
                )
            predicted_strata.add(route)
            rank_names = RANKS[route_source]
            calibrated_ranks = rank_names[:-1] if route_source == "PR2" else rank_names
            for index, _rank in enumerate(calibrated_ranks):
                if index >= len(prediction.taxonomy) or not prediction.taxonomy[index]:
                    continue
                key = (*route, index)
                calls[key] += 1
                called_truth_mixtures[key][true_class_key] += 1
                cross_domain = bool(truth) and truth[0] != predicted_domain
                comparable = (
                    cross_domain
                    or index == 0
                    or (
                        truth_source == route_source
                        and index < len(truth)
                        and all(truth[: index + 1])
                    )
                )
                if comparable:
                    truth_known_calls[key] += 1
                else:
                    truth_unknown_calls[key] += 1
                if comparable and prediction.taxonomy[: index + 1] == truth[: index + 1]:
                    correct[key] += 1

    metrics: list[dict[str, object]] = []
    rank_caps: dict[str, int] = {}
    stratum_results: dict[str, dict[str, object]] = {}
    all_strata = set(true_stratum_counts) | predicted_strata
    for stratum in sorted(all_strata):
        marker, source, domain = stratum
        if source not in RANKS:
            raise RuntimeError(
                f"unsupported taxonomy source in calibration stratum: {source!r}"
            )
        sampled = true_stratum_counts[stratum]
        evaluated = marker_counts[marker]
        accepted_cap = -1
        stratum_metrics: list[dict[str, object]] = []
        for index, rank in enumerate(RANKS[source][:-1] if source == "PR2" else RANKS[source]):
            key = (*stratum, index)
            called = calls[key]
            successes = correct[key]
            lower = wilson_lower_bound(successes, called)
            threshold = GENUS_MIN_PRECISION if rank == "genus" else HIGHER_RANK_MIN_PRECISION
            accepted = called >= MIN_CALIBRATION_CALLS and lower >= threshold
            if accepted and index == accepted_cap + 1:
                accepted_cap = index
            metric = {
                "marker": marker,
                "taxonomy_source": source,
                "domain": domain,
                "rank": rank,
                "rank_index": index,
                "sampled": sampled,
                "evaluated_queries": evaluated,
                "called": called,
                "truth_known_calls": truth_known_calls[key],
                "truth_unknown_calls": truth_unknown_calls[key],
                "correct": successes,
                "call_rate": called / evaluated if evaluated else 0.0,
                "precision": successes / called if called else 0.0,
                "wilson_95_lower": lower,
                "required_precision": threshold,
                "accepted": accepted,
                "called_true_class_mixture": dict(
                    sorted(called_truth_mixtures[key].items())
                ),
            }
            metrics.append(metric)
            stratum_metrics.append(metric)
        stratum_key = "|".join(stratum)
        if accepted_cap >= 0:
            status = "calibrated"
            rank_cap: int | None = accepted_cap
            reason = ""
            rank_caps[stratum_key] = accepted_cap
        else:
            status = "failed"
            rank_cap = None
            domain_key = (*stratum, 0)
            enough_calls = calls[domain_key] >= MIN_CALIBRATION_CALLS
            precise_enough = wilson_lower_bound(
                correct[domain_key], calls[domain_key]
            ) >= HIGHER_RANK_MIN_PRECISION
            if not enough_calls and not precise_enough:
                reason = "insufficient_domain_calls_and_precision_below_threshold"
            elif not enough_calls:
                reason = "insufficient_domain_calls"
            else:
                reason = "domain_precision_below_threshold"
        stratum_results[stratum_key] = {
            "marker": marker,
            "taxonomy_source": source,
            "domain": domain,
            "sampled": sampled,
            "evaluated_queries": evaluated,
            "sampled_true_class_mixture": dict(
                sorted(marker_truth_mixtures[marker].items())
            ),
            "status": status,
            "rank_cap": rank_cap,
            "reason": reason,
            "metrics": stratum_metrics,
        }
    requested = (
        samples_per_stratum_requested
        if samples_per_stratum_requested is not None
        else max(true_stratum_counts.values())
    )
    return {
        "schema_version": classifier.CALIBRATION_SCHEMA_VERSION,
        "status": "complete" if rank_caps else "failed",
        "calibration_use_case": classifier.CALIBRATION_USE_CASE,
        "stratum_basis": classifier.CALIBRATION_STRATUM_BASIS,
        "classification_policy": classifier.classification_policy(),
        "method": "deterministic_leave_one_reference_out",
        "samples_per_stratum_requested": requested,
        "minimum_calls": MIN_CALIBRATION_CALLS,
        "higher_rank_minimum_wilson_95_lower": HIGHER_RANK_MIN_PRECISION,
        "genus_minimum_wilson_95_lower": GENUS_MIN_PRECISION,
        "species_policy": "100 percent identity and coverage with unanimous exact candidates",
        "rank_caps": rank_caps,
        "strata": stratum_results,
        "metrics": metrics,
    }


def _write_json_durably(path: Path, value: Mapping[str, object]) -> None:
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(file_descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
        replace_and_fsync(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def run_calibration(
    profile_directory: str | Path,
    output_directory: str | Path,
    samples_per_stratum: int,
    threads: int,
    *,
    reuse_existing_blast: bool = False,
    reuse_blast_directory: str | Path | None = None,
) -> dict[str, object]:
    if threads < 1:
        raise ValueError("threads must be positive")
    profile = Path(profile_directory)
    manifest = manager.validate_profile_directory(profile, expected_profile="curated")
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    rows = select_calibration_rows(profile, samples_per_stratum)
    if not rows:
        raise RuntimeError("calibration selected no rows")
    blast_outputs = _write_queries(
        profile,
        output,
        rows,
        threads,
        reuse_existing_blast=reuse_existing_blast,
        reuse_blast_directory=reuse_blast_directory,
    )
    result = evaluate_calibration(
        profile,
        rows,
        blast_outputs,
        samples_per_stratum_requested=samples_per_stratum,
    )
    manifest_path = profile / manager.MANIFEST_NAME
    result["curated_profile"] = {
        "version": manifest["version"],
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
    }
    repository = Path(__file__).resolve().parents[1]
    result["implementation_sha256"] = {
        relative: _sha256_file(repository / relative)
        for relative in (
            "scripts/calibrate_taxonomy.py",
            "scripts/classify_img_clusters.py",
            "scripts/img_classification_data.py",
        )
    }
    result["search_bindings"] = {}
    for marker, blast_output in sorted(blast_outputs.items()):
        provenance_path = blast_output.parent / SEARCH_PROVENANCE_NAME
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        result["search_bindings"][marker] = {
            "search_provenance_sha256": _sha256_file(provenance_path),
            "schema_version": provenance["schema_version"],
            "status": provenance["status"],
            "profile_manifest_sha256": provenance["profile_manifest_sha256"],
            "files": provenance["files"],
        }
    _write_json_durably(output / "calibration.json", result)
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples-per-stratum", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--reuse-existing-blast",
        action="store_true",
        help=(
            "reuse retained leave_one_out.m8 only after validating deterministic "
            "IDs/truth, regenerated query FASTA, exact query coverage, and self hits"
        ),
    )
    parser.add_argument(
        "--reuse-blast-directory",
        type=Path,
        help=(
            "reuse validated retained searches from this directory while writing "
            "new calibration output"
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    result = run_calibration(
        args.profile,
        args.output,
        args.samples_per_stratum,
        args.threads,
        reuse_existing_blast=args.reuse_existing_blast,
        reuse_blast_directory=args.reuse_blast_directory,
    )
    print(json.dumps({"rank_caps": result["rank_caps"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
