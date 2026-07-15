#!/usr/bin/env python3
"""Build release-ready curated and IMG SSUextract database profiles."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import itertools
import json
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import assemble_database_profile as assembler
import build_database_release as builder
import classify_img_clusters as classifier
import database_sources
import img_search_provenance


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_CONFIG = REPO / "config" / "database_sources.json"
DEFAULT_NOTICES = REPO / "resources" / "database_notices"
CURATED_SOURCE_NAMES = (
    "silva_ssu_nr99_fasta",
    "silva_ssu_taxonomy",
    "pr2_ssu_fasta",
)
IMG_SOURCE_NAMES = (
    *CURATED_SOURCE_NAMES,
    "eukcensus_16s_fasta",
    "eukcensus_18s_fasta",
    "eukcensus_16s_clusters",
    "eukcensus_18s_clusters",
    "eukcensus_img_metadata",
)
SOURCE_TREE_FILES = (
    "config/database_sources.json",
    "config/img_search_execution.json",
    "pixi.lock",
    "pixi.toml",
    "resources/database_notices/EUKCENSUS_NOTICE.txt",
    "resources/database_notices/PR2_LICENSE.txt",
    "resources/database_notices/README.md",
    "resources/database_notices/SILVA_NOTICE.txt",
    "scripts/atomic_io.py",
    "scripts/assemble_database_profile.py",
    "scripts/build_database_profiles.py",
    "scripts/build_database_release.py",
    "scripts/calibrate_taxonomy.py",
    "scripts/classify_img_marker.sh",
    "scripts/classify_img_clusters.py",
    "scripts/img_classification_data.py",
    "scripts/database_contracts.py",
    "scripts/database_manager.py",
    "scripts/database_release_io.py",
    "scripts/database_sources.py",
    "scripts/extract_img_cluster_centroids.py",
    "scripts/img_chunked_search.py",
    "scripts/img_search_provenance.py",
    "scripts/search_img_marker.sh",
    "scripts/taxonomy_utils.py",
)
EVIDENCE_FIELDS = (
    "classification_status",
    "reason",
    "returned_hit_count",
    "eligible_hit_count",
    "coverage_filtered_hit_count",
    "policy",
    "best_bit_score",
    "candidate_bit_score_threshold",
    "candidate_count",
    "candidates",
    "missing_taxonomy_subjects",
    "ambiguous_taxonomy_subjects",
    "failed_calibration_strata",
    "truncated",
    "taxonomy",
    "taxonomy_source",
    "domain",
    "compartment",
    "assignment_method",
    "species_guard_applied",
    "species_called",
    "calibration_rank_cap",
    "propagation_rank_cap",
)
FORBIDDEN_EVIDENCE_KEYS = {
    "cluster_id",
    "centroid",
    "img_member_count",
    "img_members",
    "source_identifier",
    "taxon_oid",
    "project",
    "project_name",
    "email",
    "contact",
    "comments",
}


def load_source_catalog(path: str | Path) -> dict[str, object]:
    source = Path(path)
    try:
        catalog = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise builder.BuildError(f"Could not read source catalog {source}: {error}") from error
    if not isinstance(catalog, dict) or catalog.get("schema_version") != 1:
        raise builder.BuildError("source catalog schema_version must be 1")
    if not isinstance(catalog.get("sources"), dict):
        raise builder.BuildError("source catalog must contain a sources object")
    return catalog


def calibration_provenance(path: str | Path) -> dict[str, object]:
    source = Path(path)
    try:
        source_sha256 = _sha256(source)
        calibration = json.loads(source.read_text(encoding="utf-8"))
        classifier.load_calibration(source)
    except (OSError, json.JSONDecodeError, ValueError) as error:
        raise builder.BuildError(f"Invalid taxonomy calibration {source}: {error}") from error
    curated_profile = calibration.get("curated_profile")
    if not isinstance(curated_profile, dict):
        raise builder.BuildError("taxonomy calibration lacks curated_profile provenance")
    version = curated_profile.get("version")
    manifest_sha256 = curated_profile.get("manifest_sha256")
    if (
        not isinstance(version, str)
        or not version
        or not isinstance(manifest_sha256, str)
        or len(manifest_sha256) != 64
        or any(character not in "0123456789abcdef" for character in manifest_sha256)
    ):
        raise builder.BuildError("taxonomy calibration has invalid curated_profile provenance")
    if _sha256(source) != source_sha256:
        raise builder.BuildError("taxonomy calibration changed while it was being loaded")
    return {
        "sha256": source_sha256,
        "schema_version": calibration["schema_version"],
        "rank_caps": calibration["rank_caps"],
        "strata": calibration["strata"],
        "curated_profile": curated_profile,
    }


def validate_curated_manifest(
    path: str | Path, calibration: Mapping[str, object]
) -> dict[str, str]:
    source = Path(path)
    if not source.is_file():
        raise builder.BuildError(f"curated manifest is missing: {source}")
    actual_sha256 = _sha256(source)
    expected_sha256 = calibration["curated_profile"]["manifest_sha256"]
    if actual_sha256 != expected_sha256:
        raise builder.BuildError(
            "curated manifest SHA256 does not match taxonomy calibration: "
            f"{actual_sha256} != {expected_sha256}"
        )
    return {"sha256": actual_sha256}


def _load_json_object(path: Path, label: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise builder.BuildError(f"Invalid {label} {path}: {error}") from error
    if not isinstance(value, dict):
        raise builder.BuildError(f"{label} must be a JSON object: {path}")
    return value


def validate_classification_qc(
    assignment_path: str | Path,
    marker: str,
    calibration: Mapping[str, object],
    search_provenance: Mapping[str, object] | None = None,
) -> Path:
    assignment = Path(assignment_path).resolve()
    qc_path = assignment.parent / "qc.json"
    qc = _load_json_object(qc_path, "IMG classification QC")
    outcomes_path = assignment.parent / "outcomes.jsonl"
    if not assignment.is_file() or not outcomes_path.is_file():
        raise builder.BuildError(
            f"IMG classification outputs are missing for {marker}"
        )
    expected = {
        "schema_version": 1,
        "calibration": {
            "sha256": calibration["sha256"],
            "schema_version": calibration["schema_version"],
        },
        "marker": marker,
        "propagation_rank_cap": 0,
        "policy": classifier.classification_policy(),
        "outputs": {
            "assignments_tsv": {"sha256": _sha256(assignment)},
            "outcomes_jsonl": {"sha256": _sha256(outcomes_path)},
        },
    }
    if search_provenance is not None:
        expected["search_provenance"] = search_provenance
    binding = qc.get("classification_binding")
    if _canonical_json(binding) != _canonical_json(expected):
        raise builder.BuildError(
            f"IMG classification QC binding does not match build inputs for {marker}"
        )
    return qc_path


def validate_classification_search(
    assignment_path: str | Path,
    marker: str,
    curated_profile: str | Path,
    clusters: str | Path,
    source_fasta: str | Path,
    *,
    search_threads: int = 8,
) -> dict[str, object]:
    directory = Path(assignment_path).resolve().parent
    sidecar = directory / "search_provenance.json"
    try:
        payload = img_search_provenance.validate_search(
            marker,
            curated_profile,
            clusters,
            source_fasta,
            directory / "centroids.fna",
            directory / "centroids.m8",
            sidecar,
            search_threads,
        )
        return img_search_provenance.portable_provenance(
            payload, _sha256(sidecar)
        )
    except (OSError, RuntimeError, ValueError) as error:
        raise builder.BuildError(
            f"Invalid {marker} IMG search provenance: {error}"
        ) from error


def _candidate_strata(marker: str, outcome: Mapping[str, object]) -> list[str]:
    candidates = outcome.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise builder.BuildError("calibration-gated IMG outcome lacks candidates")
    domains: set[str] = set()
    sources: set[str] = set()
    for candidate in candidates:
        if not isinstance(candidate, dict):
            raise builder.BuildError("IMG outcome candidate must be an object")
        taxonomy = candidate.get("taxonomy")
        source = candidate.get("taxonomy_source")
        if (
            not isinstance(taxonomy, str)
            or not taxonomy
            or not isinstance(source, str)
            or not source
        ):
            raise builder.BuildError("IMG outcome candidate lacks taxonomy or source")
        domains.add(taxonomy.split(";", 1)[0])
        sources.add(source)
    if len(domains) != 1:
        raise builder.BuildError("calibration-gated IMG outcome lacks one candidate domain")
    domain = next(iter(domains))
    return [f"{marker}|{source}|{domain}" for source in sorted(sources)]


def validate_outcome_calibration(
    marker: str,
    outcome: Mapping[str, object],
    calibration: Mapping[str, object],
) -> None:
    status = outcome.get("classification_status")
    reason = outcome.get("reason")
    if status != "classified" and reason != "calibration_stratum_failed":
        return
    required_keys = _candidate_strata(marker, outcome)
    strata = calibration["strata"]
    missing = [key for key in required_keys if key not in strata]
    if missing:
        raise builder.BuildError(
            "IMG outcome requires missing calibration strata: " + ", ".join(missing)
        )
    failed = [
        {"key": key, "reason": strata[key]["reason"]}
        for key in required_keys
        if strata[key]["status"] == "failed"
    ]
    if status == "classified":
        if reason != "":
            raise builder.BuildError("classified IMG outcome has a non-empty reason")
        if failed:
            raise builder.BuildError("classified IMG outcome requires a failed stratum")
        expected_cap = min(strata[key]["rank_cap"] for key in required_keys)
        if type(outcome.get("calibration_rank_cap")) is not int or outcome.get(
            "calibration_rank_cap"
        ) != expected_cap:
            raise builder.BuildError("classified IMG outcome has wrong calibration rank cap")
        if (
            type(outcome.get("propagation_rank_cap")) is not int
            or outcome.get("propagation_rank_cap") != 0
        ):
            raise builder.BuildError("classified IMG outcome lacks the domain propagation cap")
        domain = required_keys[0].rsplit("|", 1)[1]
        if outcome.get("domain") != domain or outcome.get("taxonomy") != domain:
            raise builder.BuildError("classified IMG outcome exceeds the domain propagation cap")
        expected_sources = "+".join(key.split("|", 2)[1] for key in required_keys)
        if outcome.get("taxonomy_source") != expected_sources:
            raise builder.BuildError("classified IMG outcome taxonomy source is inconsistent")
        return
    evidence = outcome.get("failed_calibration_strata")
    if _canonical_json(evidence) != _canonical_json(failed) or not failed:
        raise builder.BuildError("failed IMG calibration outcome has inconsistent evidence")


def _source_paths(
    catalog: Mapping[str, object], source_directory: Path, source_names: Sequence[str]
) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for name in source_names:
        try:
            raw_entry = catalog["sources"][name]
        except KeyError as error:
            raise builder.BuildError(f"source catalog lacks {name!r}") from error
        if not isinstance(raw_entry, dict):
            raise builder.BuildError(f"source catalog entry {name!r} must be an object")
        filename = raw_entry.get("filename")
        if not isinstance(filename, str) or Path(filename).name != filename:
            raise builder.BuildError(f"source catalog entry {name!r} has an unsafe filename")
        path = source_directory / filename
        if not path.is_file():
            raise builder.BuildError(f"missing pinned source {name!r}: {path}")
        database_sources.validate_artifact(path, raw_entry)
        paths[name] = path
    return paths


def _source_version(catalog: Mapping[str, object], name: str) -> str:
    entry = catalog["sources"][name]
    version = entry.get("version") if isinstance(entry, dict) else None
    if not isinstance(version, str) or not version:
        raise builder.BuildError(f"source catalog entry {name!r} lacks a version")
    return version


def iter_curated_records(
    paths: Mapping[str, Path], catalog: Mapping[str, object]
) -> Iterator[builder.PreparedSourceRecord]:
    silva_version = _source_version(catalog, "silva_ssu_nr99_fasta")
    taxonomy_version = _source_version(catalog, "silva_ssu_taxonomy")
    if silva_version != taxonomy_version:
        raise builder.BuildError(
            "SILVA FASTA and taxonomy versions differ: "
            f"{silva_version!r} != {taxonomy_version!r}"
        )
    yield from builder.iter_curated_silva_records(
        paths["silva_ssu_nr99_fasta"],
        paths["silva_ssu_taxonomy"],
        source_version=silva_version,
    )
    yield from builder.iter_curated_pr2_records(
        paths["pr2_ssu_fasta"],
        source_version=_source_version(catalog, "pr2_ssu_fasta"),
    )


def iter_img_records(
    paths: Mapping[str, Path], catalog: Mapping[str, object]
) -> Iterator[builder.PreparedSourceRecord]:
    yield from builder.iter_img_records(
        paths["eukcensus_16s_fasta"],
        "16S",
        source_version=_source_version(catalog, "eukcensus_16s_fasta"),
    )
    yield from builder.iter_img_records(
        paths["eukcensus_18s_fasta"],
        "18S",
        source_version=_source_version(catalog, "eukcensus_18s_fasta"),
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True)


def _tool_version(command: Sequence[str]) -> str | None:
    try:
        result = subprocess.run(
            list(command), check=True, capture_output=True, text=True
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    lines = [
        line.strip()
        for line in (result.stdout + result.stderr).splitlines()
        if line.strip()
    ]
    return " | ".join(lines) or None


def _source_tree_provenance() -> dict[str, object]:
    files = []
    for relative_name in SOURCE_TREE_FILES:
        path = REPO / relative_name
        if not path.is_file():
            raise builder.BuildError(f"release source-tree file is missing: {relative_name}")
        files.append({"path": relative_name, "sha256": _sha256(path)})
    tree_sha256 = hashlib.sha256(_canonical_json(files).encode("utf-8")).hexdigest()
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain", "--", *SOURCE_TREE_FILES],
            cwd=REPO,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.splitlines()
    except (FileNotFoundError, subprocess.CalledProcessError):
        commit = None
        status = None
    try:
        duckdb_version = importlib.metadata.version("duckdb")
    except importlib.metadata.PackageNotFoundError:
        duckdb_version = None
    return {
        "repository": {
            "commit": commit,
            "relevant_worktree_dirty": None if status is None else bool(status),
            "source_tree_sha256": tree_sha256,
            "source_files": files,
        },
        "environment": {
            "python": platform.python_version(),
            "duckdb": duckdb_version,
            "pixi_lock_sha256": _sha256(REPO / "pixi.lock"),
            "tools": {
                "blastn": _tool_version(("blastn", "-version")),
                "blastdbcmd": _tool_version(("blastdbcmd", "-version")),
                "zstd": _tool_version(("zstd", "--version")),
            },
        },
        "command": {
            "entrypoint": "scripts/build_database_profiles.py",
            "profile_argument": "{profile}",
            "execution_environment": "pixi lockfile",
        },
    }


def _safe_evidence_payload(marker: str, outcome: Mapping[str, object]) -> dict[str, object]:
    status = outcome.get("classification_status")
    if status not in {"classified", "unclassified"}:
        raise builder.BuildError(f"IMG evidence outcome has invalid status: {status!r}")
    payload = {
        "schema_version": 1,
        "marker": marker,
        **{field: outcome[field] for field in EVIDENCE_FIELDS if field in outcome},
    }
    _assert_evidence_privacy(payload)
    return payload


def _assert_evidence_privacy(value: object) -> None:
    if isinstance(value, Mapping):
        forbidden = FORBIDDEN_EVIDENCE_KEYS.intersection(str(key) for key in value)
        if forbidden:
            raise builder.BuildError(
                f"IMG evidence contains forbidden identity/metadata keys: {sorted(forbidden)}"
            )
        for nested in value.values():
            _assert_evidence_privacy(nested)
    elif isinstance(value, (list, tuple)):
        for nested in value:
            _assert_evidence_privacy(nested)
    elif isinstance(value, str) and "IMG_" in value:
        raise builder.BuildError("IMG evidence contains a raw IMG identifier")


def _outcome_evidence_id(outcome: Mapping[str, object]) -> str:
    existing = outcome.get("evidence_id")
    if (
        isinstance(existing, str)
        and len(existing) == 70
        and existing.startswith("IMGEV_")
        and all(character in "0123456789abcdef" for character in existing[6:])
    ):
        return existing
    if existing is not None:
        raise builder.BuildError(f"invalid IMG outcome evidence_id: {existing!r}")
    return "IMGEV_" + hashlib.sha256(
        _canonical_json(dict(outcome)).encode("utf-8")
    ).hexdigest()


def create_img_evidence_catalog(
    assignment_paths: Sequence[Path],
    destination: Path,
    calibration: Mapping[str, object],
) -> tuple[dict[Path, dict[str, dict[str, str]]], dict[str, object]]:
    """Write a content-addressed evidence catalog with no IMG cluster identities."""

    mappings: dict[Path, dict[str, dict[str, str]]] = {}
    seen_evidence_ids: set[str] = set()
    input_files: list[dict[str, object]] = []
    destination.parent.mkdir(parents=True, exist_ok=True)
    body = destination.with_name(f".{destination.name}.body")
    try:
        with body.open("w", encoding="utf-8", newline="\n") as body_handle:
            for assignment_path in sorted(
                assignment_paths, key=lambda item: (item.parent.name, str(item))
            ):
                path = assignment_path.resolve()
                marker = path.parent.name
                if marker not in {"16S", "18S"}:
                    raise builder.BuildError(
                        "IMG assignment parent must identify marker 16S or 18S: "
                        f"{assignment_path}"
                    )
                outcomes_path = path.parent / "outcomes.jsonl"
                if not outcomes_path.is_file():
                    raise builder.BuildError(
                        f"IMG evidence outcomes are missing: {outcomes_path}"
                    )
                marker_mapping: dict[str, dict[str, str]] = {}
                with outcomes_path.open(encoding="utf-8") as handle:
                    for line_number, line in enumerate(handle, 1):
                        try:
                            outcome = json.loads(line)
                        except json.JSONDecodeError as error:
                            raise builder.BuildError(
                                "invalid IMG outcome JSON at "
                                f"{outcomes_path}:{line_number}"
                            ) from error
                        if not isinstance(outcome, dict):
                            raise builder.BuildError(
                                "IMG outcome must be an object at "
                                f"{outcomes_path}:{line_number}"
                            )
                        validate_outcome_calibration(marker, outcome, calibration)
                        old_id = _outcome_evidence_id(outcome)
                        payload = _safe_evidence_payload(marker, outcome)
                        content_sha256 = hashlib.sha256(
                            _canonical_json(payload).encode("utf-8")
                        ).hexdigest()
                        evidence_id = f"IMGEV_{content_sha256}"
                        if outcome["classification_status"] == "classified":
                            resolution = {
                                "evidence_id": evidence_id,
                                "taxonomy": str(outcome.get("taxonomy", "")),
                                "taxonomy_source": str(
                                    outcome.get("taxonomy_source", "")
                                ),
                                "assignment_method": str(
                                    outcome.get("assignment_method", "")
                                ),
                                "compartment": str(outcome.get("compartment", "")),
                            }
                        else:
                            resolution = {
                                "evidence_id": evidence_id,
                                "taxonomy": "Unclassified",
                                "taxonomy_source": "SILVA+PR2",
                                "assignment_method": "updated_reference_unclassified",
                                "compartment": "",
                            }
                        previous = marker_mapping.setdefault(old_id, resolution)
                        if previous != resolution:
                            raise builder.BuildError(
                                f"conflicting IMG evidence ID: {old_id}"
                            )
                        if evidence_id not in seen_evidence_ids:
                            seen_evidence_ids.add(evidence_id)
                            entry = {
                                "record_type": "evidence",
                                "evidence_id": evidence_id,
                                "content_sha256": content_sha256,
                                "decision": payload,
                            }
                            body_handle.write(_canonical_json(entry) + "\n")
                mappings[path] = marker_mapping
                input_files.extend(
                    [
                        {
                            "role": "assignment_table",
                            "marker": marker,
                            "bytes": path.stat().st_size,
                            "sha256": _sha256(path),
                        },
                        {
                            "role": "classification_outcomes",
                            "marker": marker,
                            "bytes": outcomes_path.stat().st_size,
                            "sha256": _sha256(outcomes_path),
                        },
                    ]
                )
        header = {
            "record_type": "catalog",
            "schema_version": 1,
            "evidence_records": len(seen_evidence_ids),
            "privacy": {
                "cluster_ids_packaged": False,
                "centroid_ids_packaged": False,
                "img_member_ids_packaged": False,
                "raw_cluster_tables_packaged": False,
            },
        }
        with destination.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(_canonical_json(header) + "\n")
            with body.open(encoding="utf-8") as body_handle:
                shutil.copyfileobj(body_handle, handle)
    finally:
        body.unlink(missing_ok=True)
    metadata = {
        "path": "EVIDENCE/img_taxonomy_evidence.jsonl",
        "records": len(seen_evidence_ids),
        "bytes": destination.stat().st_size,
        "sha256": _sha256(destination),
        "inputs": sorted(input_files, key=lambda row: (str(row["marker"]), str(row["role"]))),
    }
    return mappings, metadata


def read_assignment_rows(
    paths: Iterable[str | Path],
    evidence_maps: Mapping[Path, Mapping[str, Mapping[str, str]]] | None = None,
) -> Iterator[dict[str, str]]:
    for path_value in paths:
        path = Path(path_value).resolve()
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "source_identifier",
                "taxonomy",
                "taxonomy_source",
                "assignment_method",
                "evidence",
                "compartment",
            }
            missing = required - set(reader.fieldnames or ())
            if missing:
                raise builder.BuildError(
                    f"IMG assignment table {path} lacks columns: {sorted(missing)}"
                )
            for row in reader:
                if evidence_maps is not None:
                    old_id = str(row.get("evidence_id", "")).strip()
                    if not old_id:
                        raise builder.BuildError(
                            f"IMG assignment table {path} lacks evidence_id values"
                        )
                    try:
                        resolution = evidence_maps[path][old_id]
                    except KeyError as error:
                        raise builder.BuildError(
                            f"IMG assignment evidence_id {old_id!r} does not resolve in {path}"
                        ) from error
                    new_id = resolution["evidence_id"]
                    for field in (
                        "taxonomy",
                        "taxonomy_source",
                        "assignment_method",
                        "compartment",
                    ):
                        if str(row.get(field, "")).strip() != resolution[field]:
                            raise builder.BuildError(
                                f"IMG assignment {field} does not match evidence {old_id!r}"
                            )
                    evidence = str(row["evidence"])
                    if old_id not in evidence:
                        raise builder.BuildError(
                            f"IMG assignment evidence text does not contain {old_id!r}"
                        )
                    row["evidence"] = evidence.replace(old_id, new_id)
                    row["evidence_id"] = new_id
                yield row


def read_img_locations(
    path: str | Path, included_taxon_oids: set[str]
) -> tuple[tuple[builder.ImgLocation, ...], list[dict[str, object]]]:
    corrections: list[dict[str, object]] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        locations = builder.clean_img_metadata(
            csv.DictReader(handle, delimiter="\t"), corrections=corrections
        )
    filtered = tuple(
        location for location in locations if location.taxon_oid in included_taxon_oids
    )
    return filtered, corrections


def _release_files(profile: str, notices: Path) -> dict[str, Path]:
    names = ["SILVA_NOTICE.txt", "PR2_LICENSE.txt"]
    if profile == "img":
        names.append("EUKCENSUS_NOTICE.txt")
    files = {"README.md": notices / "README.md"}
    files.update({f"LICENSES/{name}": notices / name for name in names})
    return files


def _provenance_details(
    catalog: Mapping[str, object],
    source_names: Sequence[str],
    profile: str,
    corrections: Sequence[Mapping[str, object]] = (),
    extra: Mapping[str, object] | None = None,
    build_integrity: Mapping[str, object] | None = None,
) -> dict[str, object]:
    sources = catalog["sources"]
    selected = {name: sources[name] for name in source_names}
    integrity = json.loads(json.dumps(build_integrity or _source_tree_provenance()))
    integrity["command"]["profile_argument"] = profile
    silva_version = _source_version(catalog, "silva_ssu_nr99_fasta")
    pr2_version = _source_version(catalog, "pr2_ssu_fasta")
    details: dict[str, object] = {
        "taxonomy_policy_version": catalog.get("taxonomy_policy_version"),
        "source_inputs": selected,
        "sequence_normalization": "remove whitespace; uppercase; U to T",
        "exact_sequence_deduplication": True,
        "taxonomy_policy": {
            "Bacteria": f"SILVA {silva_version}",
            "Archaea": f"SILVA {silva_version}",
            "Eukaryota": f"PR2 {pr2_version}",
            "cross_domain_exact_sequence": "explicit ambiguity with native alternatives",
        },
        "raw_marker_fasta_packaged": False,
        "build_integrity": integrity,
    }
    if profile == "img":
        details.update(
            {
                "embedded_reference_records_retained": False,
                "img_taxonomy": "cluster-centroid BLAST against this release's curated profile",
                "cluster_propagation_rank_cap": "domain",
                "cluster_propagation_reason": (
                    "source cluster construction and within-cluster taxonomic coherence "
                    "are not documented"
                ),
                "cluster_tables_packaged": False,
                "img_metadata_columns": list(builder.IMG_LOCATION_COLUMNS),
                "metadata_corrections": list(corrections),
            }
        )
    if extra:
        details.update(extra)
    return details


def _assemble(
    model: builder.DatabaseModel,
    args: argparse.Namespace,
    profile: str,
    catalog: Mapping[str, object],
    source_names: Sequence[str],
    locations: Sequence[builder.ImgLocation] = (),
    corrections: Sequence[Mapping[str, object]] = (),
    extra_provenance: Mapping[str, object] | None = None,
    extra_release_files: Mapping[str, str | Path] | None = None,
) -> assembler.AssemblyResult:
    archive = args.archive_directory / f"ssuextract-db-{profile}-v{args.version}.tar.zst"
    return assembler.assemble_database_profile(
        model,
        args.output_root / profile,
        profile=profile,
        version=args.version,
        img_locations=locations,
        archive_path=archive,
        provenance_details=_provenance_details(
            catalog,
            source_names,
            profile,
            corrections,
            extra_provenance,
            getattr(args, "build_integrity", None),
        ),
        release_files={
            **_release_files(profile, args.notices),
            **(extra_release_files or {}),
        },
    )


def build_curated(args: argparse.Namespace) -> assembler.AssemblyResult:
    catalog = load_source_catalog(args.source_config)
    paths = _source_paths(catalog, args.source_directory, CURATED_SOURCE_NAMES)
    model = builder.build_deduplicated_model(iter_curated_records(paths, catalog))
    return _assemble(model, args, "curated", catalog, CURATED_SOURCE_NAMES)


def build_img(args: argparse.Namespace) -> assembler.AssemblyResult:
    assignment_paths = tuple(Path(path).resolve() for path in args.assignments)
    markers = [path.parent.name for path in assignment_paths]
    if sorted(markers) != ["16S", "18S"]:
        raise builder.BuildError(
            "IMG profile requires exactly one 16S and one 18S assignment table"
        )
    if getattr(args, "curated_manifest", None) is None:
        raise builder.BuildError("IMG profile requires a curated manifest")
    calibration = calibration_provenance(args.calibration)
    curated_manifest = validate_curated_manifest(args.curated_manifest, calibration)
    catalog = load_source_catalog(args.source_config)
    paths = _source_paths(catalog, args.source_directory, IMG_SOURCE_NAMES)
    curated_profile = Path(args.curated_manifest).resolve().parent
    classification_search = {
        marker: validate_classification_search(
            path,
            marker,
            curated_profile,
            paths[f"eukcensus_{marker.lower()}_clusters"],
            paths[f"eukcensus_{marker.lower()}_fasta"],
        )
        for path, marker in zip(assignment_paths, markers)
    }
    classification_qc = {
        marker: validate_classification_qc(
            path, marker, calibration, classification_search[marker]
        )
        for path, marker in zip(assignment_paths, markers)
    }
    with tempfile.TemporaryDirectory(prefix="ssuextract-img-evidence-") as temporary:
        evidence_path = Path(temporary) / "img_taxonomy_evidence.jsonl"
        evidence_maps, evidence_metadata = create_img_evidence_catalog(
            assignment_paths, evidence_path, calibration
        )
        prepared = itertools.chain(
            iter_curated_records(paths, catalog), iter_img_records(paths, catalog)
        )
        model = builder.build_deduplicated_model(prepared)
        derived = builder.ingest_derived_cluster_assignments(
            read_assignment_rows(assignment_paths, evidence_maps), model.source_records
        )
        model = builder.add_taxonomy_assignments(model, derived)
        taxon_oids = {
            record.taxon_oid
            for record in model.source_records
            if record.reference_source == "IMG" and record.taxon_oid is not None
        }
        locations, corrections = read_img_locations(
            paths["eukcensus_img_metadata"], taxon_oids
        )
        classification_inputs = {
            "calibration": calibration,
            "curated_manifest": curated_manifest,
            "search_provenance": classification_search,
            "evidence_catalog": evidence_metadata,
        }
        qc_files = {
            f"QC/{marker}_classification.json": path
            for marker, path in classification_qc.items()
        }
        qc_files["QC/taxonomy_calibration.json"] = args.calibration
        qc_files["EVIDENCE/img_taxonomy_evidence.jsonl"] = evidence_path
        return _assemble(
            model,
            args,
            "img",
            catalog,
            IMG_SOURCE_NAMES,
            locations,
            corrections,
            {"img_classification": classification_inputs},
            qc_files,
        )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", choices=("curated", "img"))
    parser.add_argument("--source-directory", type=Path, required=True)
    parser.add_argument("--source-config", type=Path, default=DEFAULT_SOURCE_CONFIG)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--archive-directory", type=Path, required=True)
    parser.add_argument("--version", default="1.0.0")
    parser.add_argument("--notices", type=Path, default=DEFAULT_NOTICES)
    parser.add_argument(
        "--assignments",
        type=Path,
        action="append",
        default=[],
        help="IMG assignment TSV; pass once per marker for the img profile.",
    )
    parser.add_argument("--calibration", type=Path)
    parser.add_argument("--curated-manifest", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.profile == "img" and (
        not args.assignments
        or args.calibration is None
        or args.curated_manifest is None
    ):
        raise SystemExit(
            "build-database-profiles: img requires --assignments, --calibration, "
            "and --curated-manifest"
        )
    args.output_root.mkdir(parents=True, exist_ok=True)
    args.archive_directory.mkdir(parents=True, exist_ok=True)
    # Snapshot the code/environment before long-running parsing or indexing so
    # provenance cannot drift if the worktree changes during an HPC build.
    args.build_integrity = _source_tree_provenance()
    result = build_curated(args) if args.profile == "curated" else build_img(args)
    output = {
        "profile": args.profile,
        "profile_directory": str(result.profile_directory),
        "archive": {
            "path": str(result.archive.path),
            "bytes": result.archive.bytes,
            "sha256": result.archive.sha256,
        },
        "counts": json.loads(
            (result.profile_directory / "provenance.json").read_text(encoding="utf-8")
        )["counts"],
    }
    json.dump(output, sys.stdout, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
