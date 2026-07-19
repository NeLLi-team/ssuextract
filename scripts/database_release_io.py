"""Release validation and atomic serialization for SSU database models."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Iterable, Sequence

from atomic_io import replace_and_fsync
from database_contracts import (
    BuildError,
    DatabaseModel,
    FastaFormatError,
    IMG_LOCATION_COLUMNS,
    ImgLocation,
    PREFERRED_COMPARTMENTS,
    PreferredTaxonomy,
    ReleaseValidationError,
    SequenceRecord,
    SourceRecord,
    TaxonomyAssignment,
)
from database_sources import (
    clean_img_metadata_row,
    sequence_identifier,
    validate_privacy_columns,
)


def _duplicates(values: Iterable[str]) -> set[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return duplicates


def validate_release(model: DatabaseModel, img_locations: Iterable[ImgLocation] = ()) -> None:
    """Apply release-blocking integrity, privacy, taxonomy, and sequence checks."""

    errors: list[str] = []
    sequence_ids = {record.sequence_id for record in model.sequences}
    duplicate_sequences = _duplicates(record.sequence_id for record in model.sequences)
    if duplicate_sequences:
        errors.append(f"duplicate subject IDs: {sorted(duplicate_sequences)}")
    for record in model.sequences:
        try:
            if sequence_identifier(record.sequence) != record.sequence_id:
                errors.append(f"invalid sequence hash for {record.sequence_id}")
        except FastaFormatError as error:
            errors.append(str(error))
        if not record.markers:
            errors.append(f"sequence has no marker membership: {record.sequence_id}")

    source_ids = {record.source_record_id for record in model.source_records}
    source_sequence_ids = {
        record.source_record_id: record.sequence_id for record in model.source_records
    }
    duplicate_sources = _duplicates(record.source_record_id for record in model.source_records)
    if duplicate_sources:
        errors.append(f"duplicate source record IDs: {sorted(duplicate_sources)}")
    for record in model.source_records:
        if record.sequence_id not in sequence_ids:
            errors.append(f"orphan source_records.sequence_id: {record.sequence_id}")
    expected_markers: dict[str, set[str]] = defaultdict(set)
    for record in model.source_records:
        expected_markers[record.sequence_id].add(record.marker)
    for record in model.sequences:
        if set(record.markers) != expected_markers.get(record.sequence_id, set()):
            errors.append(f"marker membership mismatch for {record.sequence_id}")

    assignment_ids = _duplicates(
        assignment.taxonomy_assignment_id for assignment in model.taxonomy_assignments
    )
    if assignment_ids:
        errors.append(f"duplicate taxonomy assignment IDs: {sorted(assignment_ids)}")
    pr2_assignment_suffixes = 0
    pr2_assignment_suffix_examples: list[str] = []
    for assignment in model.taxonomy_assignments:
        if assignment.source_record_id not in source_ids:
            errors.append(f"orphan taxonomy source_record_id: {assignment.source_record_id}")
        if assignment.sequence_id not in sequence_ids:
            errors.append(f"orphan taxonomy sequence_id: {assignment.sequence_id}")
        if source_sequence_ids.get(assignment.source_record_id) not in {
            None,
            assignment.sequence_id,
        }:
            errors.append(
                f"taxonomy/source sequence mismatch: {assignment.taxonomy_assignment_id}"
            )
        if assignment.taxonomy_source.upper() == "PR2" and any(
            ":" in taxon for taxon in assignment.taxonomy
        ):
            pr2_assignment_suffixes += 1
            if len(pr2_assignment_suffix_examples) < 5:
                pr2_assignment_suffix_examples.append(assignment.taxonomy_assignment_id)
        if assignment.centroid_name:
            if any(character in assignment.centroid_name for character in "\r\n\t|"):
                errors.append(
                    f"invalid centroid name: {assignment.taxonomy_assignment_id}"
                )
            if assignment.assignment_method == "updated_reference_cluster" and (
                not assignment.centroid_taxonomy
                or not assignment.centroid_taxonomy_source
                or assignment.centroid_taxonomy[: len(assignment.taxonomy)]
                != assignment.taxonomy
            ):
                errors.append(
                    "invalid cluster centroid evidence: "
                    f"{assignment.taxonomy_assignment_id}"
                )
        elif assignment.assignment_method in {
            "updated_reference_cluster",
            "updated_reference_unclassified",
        }:
            errors.append(
                f"missing cluster centroid name: {assignment.taxonomy_assignment_id}"
            )
    if pr2_assignment_suffixes:
        errors.append(
            f"PR2 compartment suffix remains in {pr2_assignment_suffixes} "
            "normalized taxonomy assignments; examples: "
            + ", ".join(pr2_assignment_suffix_examples)
        )

    preferred_ids = [record.sequence_id for record in model.preferred_taxonomy]
    if _duplicates(preferred_ids):
        errors.append("duplicate preferred taxonomy sequence IDs")
    missing_preferred = sequence_ids - set(preferred_ids)
    if missing_preferred:
        errors.append(f"missing preferred taxonomy for {len(missing_preferred)} sequences")
    if set(preferred_ids) - sequence_ids:
        errors.append("orphan preferred taxonomy sequence_id")
    preferred_domain_mismatches = 0
    preferred_domain_mismatch_examples: list[str] = []
    pr2_preferred_suffixes = 0
    pr2_preferred_suffix_examples: list[str] = []
    for record in model.preferred_taxonomy:
        if record.compartment not in PREFERRED_COMPARTMENTS:
            errors.append(
                f"invalid preferred taxonomy compartment for {record.sequence_id}: "
                f"{record.compartment!r}"
            )
        if record.taxonomy and record.domain != record.taxonomy[0]:
            preferred_domain_mismatches += 1
            if len(preferred_domain_mismatch_examples) < 5:
                preferred_domain_mismatch_examples.append(record.sequence_id)
        if "PR2" in record.taxonomy_source.upper().split("+") and any(
            ":" in taxon for taxon in record.taxonomy
        ):
            pr2_preferred_suffixes += 1
            if len(pr2_preferred_suffix_examples) < 5:
                pr2_preferred_suffix_examples.append(record.sequence_id)
        if record.cross_domain_conflict and not (
            record.domain == "ambiguous"
            and not record.taxonomy
            and record.compartment == "mixed"
            and record.assignment_method == "cross_domain_ambiguous_exact_sequence"
            and record.taxonomy_alternatives
        ):
            errors.append(f"unresolved cross-domain taxonomy conflict for {record.sequence_id}")
        if record.centroid_names:
            try:
                centroid_names = json.loads(record.centroid_names)
            except json.JSONDecodeError:
                centroid_names = None
            valid_centroid_names = (
                isinstance(centroid_names, list)
                and bool(centroid_names)
                and all(
                    isinstance(name, str)
                    and name
                    and not any(character in name for character in "\r\n\t|")
                    for name in centroid_names
                )
            )
            if valid_centroid_names:
                valid_centroid_names = centroid_names == sorted(set(centroid_names))
            if not valid_centroid_names:
                errors.append(f"invalid centroid names for {record.sequence_id}")
            if record.centroid_taxonomy and (
                not record.centroid_taxonomy_source
                or (
                    record.domain not in {"ambiguous", "Unclassified"}
                    and record.centroid_taxonomy[0] != record.domain
                )
            ):
                errors.append(f"invalid centroid taxonomy for {record.sequence_id}")
        elif record.centroid_taxonomy or record.centroid_taxonomy_source:
            errors.append(f"orphan centroid taxonomy for {record.sequence_id}")
    if preferred_domain_mismatches:
        errors.append(
            f"preferred taxonomy domain mismatch in {preferred_domain_mismatches} rows; "
            "examples: " + ", ".join(preferred_domain_mismatch_examples)
        )
    if pr2_preferred_suffixes:
        errors.append(
            f"PR2 compartment suffix remains in {pr2_preferred_suffixes} "
            "preferred taxonomy rows; examples: "
            + ", ".join(pr2_preferred_suffix_examples)
        )

    assignment_sources = {assignment.source_record_id for assignment in model.taxonomy_assignments}
    missing_img = [
        record.source_identifier
        for record in model.source_records
        if record.reference_source.upper() == "IMG" and record.source_record_id not in assignment_sources
    ]
    if missing_img:
        errors.append(f"IMG taxonomy is missing for {len(missing_img)} source records")

    try:
        validate_privacy_columns(IMG_LOCATION_COLUMNS)
        locations = tuple(img_locations)
        if _duplicates(location.taxon_oid for location in locations):
            errors.append("duplicate IMG location taxon_oid")
        for location in locations:
            clean_img_metadata_row(asdict(location))
    except ReleaseValidationError as error:
        errors.append(str(error))
    if errors:
        raise ReleaseValidationError("Release validation failed: " + "; ".join(errors))


def write_marker_fasta(
    path: str | Path, sequences: Iterable[SequenceRecord], *, line_width: int = 80
) -> Path:
    """Atomically write unique marker FASTA sorted by stable sequence ID."""

    if line_width < 1:
        raise ValueError("line_width must be positive")
    records = sorted(sequences, key=lambda record: record.sequence_id)
    duplicates = _duplicates(record.sequence_id for record in records)
    if duplicates:
        raise ReleaseValidationError(f"duplicate subject IDs: {sorted(duplicates)}")
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="ascii", newline="\n") as handle:
            for record in records:
                if sequence_identifier(record.sequence) != record.sequence_id:
                    raise ReleaseValidationError(f"invalid sequence/hash pair: {record.sequence_id}")
                handle.write(f">{record.sequence_id}\n")
                for offset in range(0, len(record.sequence), line_width):
                    handle.write(record.sequence[offset : offset + line_width] + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        replace_and_fsync(temporary, target)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return target


def _write_parquet(
    path: Path,
    create_sql: str,
    insert_sql: str,
    rows: Sequence[tuple[object, ...]],
) -> Path:
    try:
        import duckdb
    except ImportError as error:  # pragma: no cover - dependency is pinned by pixi
        raise BuildError("DuckDB is required to write release Parquet files") from error
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{path.name}.", dir=path.parent) as staging:
        temporary = Path(staging) / path.name
        connection = duckdb.connect(":memory:")
        try:
            # A single writer thread keeps Parquet bytes independent of the Slurm
            # allocation used for a release rebuild.
            connection.execute("SET threads = 1")
            connection.execute(create_sql)
            if rows:
                connection.executemany(insert_sql, rows)
            connection.execute(
                "COPY release_table TO ? (FORMAT PARQUET, COMPRESSION ZSTD)",
                [str(temporary)],
            )
            with temporary.open("rb") as handle:
                os.fsync(handle.fileno())
            replace_and_fsync(temporary, path)
        finally:
            connection.close()
    return path


def write_sequences_parquet(path: str | Path, records: Iterable[SequenceRecord]) -> Path:
    rows = [
        (
            record.sequence_id,
            len(record.sequence),
            hashlib.sha256(record.sequence.encode("ascii")).hexdigest(),
            ";".join(record.markers),
        )
        for record in sorted(records, key=lambda item: item.sequence_id)
    ]
    return _write_parquet(
        Path(path),
        "CREATE TABLE release_table(sequence_id VARCHAR, length BIGINT, sha256 VARCHAR, markers VARCHAR)",
        "INSERT INTO release_table VALUES (?, ?, ?, ?)",
        rows,
    )


def write_source_records_parquet(path: str | Path, records: Iterable[SourceRecord]) -> Path:
    rows = [
        (
            record.source_record_id,
            record.sequence_id,
            record.reference_source,
            record.source_version,
            record.source_identifier,
            record.original_header,
            record.marker,
            record.taxon_oid,
        )
        for record in sorted(records, key=lambda item: item.source_record_id)
    ]
    return _write_parquet(
        Path(path),
        """CREATE TABLE release_table(
            source_record_id VARCHAR, sequence_id VARCHAR, reference_source VARCHAR,
            source_version VARCHAR, source_identifier VARCHAR, original_header VARCHAR,
            marker VARCHAR, taxon_oid VARCHAR)""",
        "INSERT INTO release_table VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )


def write_taxonomy_assignments_parquet(
    path: str | Path, records: Iterable[TaxonomyAssignment]
) -> Path:
    rows = [
        (
            record.taxonomy_assignment_id,
            record.source_record_id,
            record.sequence_id,
            ";".join(record.taxonomy),
            record.taxonomy_source,
            record.domain,
            record.compartment,
            record.assignment_method,
            record.evidence,
            record.centroid_name,
            ";".join(record.centroid_taxonomy),
            record.centroid_taxonomy_source,
        )
        for record in sorted(records, key=lambda item: item.taxonomy_assignment_id)
    ]
    return _write_parquet(
        Path(path),
        """CREATE TABLE release_table(
            taxonomy_assignment_id VARCHAR, source_record_id VARCHAR, sequence_id VARCHAR,
            taxonomy VARCHAR, taxonomy_source VARCHAR, domain VARCHAR, compartment VARCHAR,
            assignment_method VARCHAR, evidence VARCHAR, centroid_name VARCHAR,
            centroid_taxonomy VARCHAR, centroid_taxonomy_source VARCHAR)""",
        "INSERT INTO release_table VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )


def write_preferred_taxonomy_parquet(
    path: str | Path, records: Iterable[PreferredTaxonomy]
) -> Path:
    rows = [
        (
            record.sequence_id,
            record.reference_source,
            ";".join(record.taxonomy),
            record.taxonomy_source,
            record.domain,
            record.compartment,
            record.assignment_method,
            record.cross_domain_conflict,
            record.taxonomy_alternatives,
            record.centroid_names,
            ";".join(record.centroid_taxonomy),
            record.centroid_taxonomy_source,
        )
        for record in sorted(records, key=lambda item: item.sequence_id)
    ]
    return _write_parquet(
        Path(path),
        """CREATE TABLE release_table(
            sequence_id VARCHAR, reference_source VARCHAR, taxonomy VARCHAR,
            taxonomy_source VARCHAR, domain VARCHAR, compartment VARCHAR,
            assignment_method VARCHAR, cross_domain_conflict BOOLEAN,
            taxonomy_alternatives VARCHAR, centroid_names VARCHAR,
            centroid_taxonomy VARCHAR, centroid_taxonomy_source VARCHAR)""",
        "INSERT INTO release_table VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )


def write_img_location_parquet(path: str | Path, records: Iterable[ImgLocation]) -> Path:
    validate_privacy_columns(IMG_LOCATION_COLUMNS)
    rows = [
        (record.taxon_oid, record.latitude, record.longitude)
        for record in sorted(records, key=lambda item: item.taxon_oid)
    ]
    return _write_parquet(
        Path(path),
        "CREATE TABLE release_table(taxon_oid VARCHAR, latitude DOUBLE, longitude DOUBLE)",
        "INSERT INTO release_table VALUES (?, ?, ?)",
        rows,
    )


def write_release_tables(
    output_directory: str | Path, model: DatabaseModel, img_locations: Iterable[ImgLocation] = ()
) -> None:
    output = Path(output_directory)
    locations = tuple(img_locations)
    validate_release(model, locations)
    write_sequences_parquet(output / "sequences.parquet", model.sequences)
    write_source_records_parquet(output / "source_records.parquet", model.source_records)
    write_taxonomy_assignments_parquet(
        output / "taxonomy_assignments.parquet", model.taxonomy_assignments
    )
    write_preferred_taxonomy_parquet(
        output / "preferred_taxonomy.parquet", model.preferred_taxonomy
    )
    write_img_location_parquet(output / "img_location.parquet", locations)
