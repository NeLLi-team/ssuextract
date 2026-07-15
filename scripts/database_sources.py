"""Pinned source acquisition and source-specific SSU record adapters."""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import math
import os
import re
import tempfile
import urllib.request
from pathlib import Path
from typing import Callable, Iterable, Iterator, Mapping, Sequence, TextIO

from atomic_io import replace_and_fsync
from database_contracts import (
    BuildError,
    FastaFormatError,
    FastaRecord,
    IMG_LOCATION_COLUMNS,
    ImgIdentifier,
    ImgLocation,
    PR2_COMPARTMENTS,
    PR2_RANKS,
    PROKARYOTIC_DOMAINS,
    PreparedSourceRecord,
    Pr2Header,
    ReleaseValidationError,
    SILVA_PROKARYOTIC_RANKS,
    SilvaHeader,
    SourceIntegrityError,
    TaxonomyError,
)
from taxonomy_utils import taxonomy_path


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCES = REPO / "config" / "database_sources.json"
FORBIDDEN_PII_TOKENS = frozenset(
    {"email", "contact", "name", "comment", "address", "phone", "institution"}
)
_IUPAC_NUCLEOTIDES = frozenset("ACGTURYSWKMBDHVN")
_IMG_ID = re.compile(r"^IMG_(?P<taxon_oid>[0-9]+)(?:[._|]|$)")


def load_source_config(path: str | Path = DEFAULT_SOURCES) -> dict[str, dict[str, object]]:
    """Read and minimally validate ``config/database_sources.json``."""

    source = Path(path)
    try:
        data = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise BuildError(f"Could not read source config {source}: {error}") from error
    if not isinstance(data, dict) or data.get("schema_version") != 1:
        raise BuildError("database source config schema_version must be 1")
    sources = data.get("sources")
    if not isinstance(sources, dict) or not sources:
        raise BuildError("database source config must contain a non-empty sources object")
    for name, entry in sources.items():
        if not isinstance(name, str) or not isinstance(entry, dict):
            raise BuildError("database source entries must be named JSON objects")
    return sources


def _new_md5():
    try:
        return hashlib.md5(usedforsecurity=False)
    except TypeError:  # pragma: no cover - Python builds without this keyword
        return hashlib.md5()


def validate_artifact(path: Path, source: Mapping[str, object]) -> None:
    expected_size = source.get("bytes")
    if not isinstance(expected_size, int) or expected_size < 0:
        raise SourceIntegrityError("source bytes must be a non-negative integer")
    actual_size = path.stat().st_size
    if actual_size != expected_size:
        raise SourceIntegrityError(
            f"Size mismatch for {path.name}: expected {expected_size}, found {actual_size}"
        )

    sha256 = hashlib.sha256()
    md5 = _new_md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            sha256.update(chunk)
            md5.update(chunk)
    checks = (("sha256", "SHA-256", sha256.hexdigest()), ("md5", "MD5", md5.hexdigest()))
    checked = False
    for key, label, actual in checks:
        expected = source.get(key)
        if expected is None:
            continue
        checked = True
        if not isinstance(expected, str) or actual.lower() != expected.lower():
            raise SourceIntegrityError(
                f"{label} mismatch for {path.name}: expected {expected}, found {actual}"
            )
    if not checked:
        raise SourceIntegrityError("source must pin at least one SHA-256 or MD5 checksum")


def fetch_artifact(
    source: Mapping[str, object],
    destination: str | Path,
    *,
    opener: Callable[..., object] = urllib.request.urlopen,
    timeout: int = 60,
) -> Path:
    """Fetch one pinned artifact and atomically install it after all checks pass."""

    target = Path(destination)
    url = source.get("url")
    if not isinstance(url, str) or not url:
        raise SourceIntegrityError("source URL must be a non-empty string")
    if not url.startswith("https://"):
        raise SourceIntegrityError("source URL must use HTTPS")
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        validate_artifact(target, source)
        return target

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".part", dir=target.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output, opener(url, timeout=timeout) as response:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                output.write(chunk)
            output.flush()
            os.fsync(output.fileno())
        validate_artifact(temporary, source)
        replace_and_fsync(temporary, target)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return target


def fetch_configured_source(
    name: str,
    destination_directory: str | Path,
    *,
    config_path: str | Path = DEFAULT_SOURCES,
    opener: Callable[..., object] = urllib.request.urlopen,
) -> Path:
    sources = load_source_config(config_path)
    if name not in sources:
        raise BuildError(f"Unknown configured source: {name}")
    entry = sources[name]
    filename = entry.get("filename")
    if not isinstance(filename, str) or Path(filename).name != filename:
        raise BuildError(f"Configured source {name!r} has an unsafe filename")
    return fetch_artifact(entry, Path(destination_directory) / filename, opener=opener)


def _open_fasta(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline=None)
    return path.open("rt", encoding="utf-8", newline=None)


def normalize_sequence_for_hashing(sequence: str) -> str:
    """Apply the only permitted hash normalization: whitespace, case, and U/T."""

    return "".join(sequence.split()).upper().replace("U", "T")


def validate_nucleotide_sequence(sequence: str, label: str = "sequence") -> str:
    normalized = normalize_sequence_for_hashing(sequence)
    if not normalized:
        raise FastaFormatError(f"Empty nucleotide {label}")
    invalid = sorted(set(normalized) - _IUPAC_NUCLEOTIDES)
    if invalid:
        raise FastaFormatError(f"Invalid nucleotide {label}: unexpected {''.join(invalid)!r}")
    return normalized


def sequence_identifier(sequence: str) -> str:
    normalized = validate_nucleotide_sequence(sequence)
    digest = hashlib.sha256(normalized.encode("ascii")).digest()
    # NCBI BLAST local identifiers are limited to 50 characters. Unpadded
    # base64url preserves all 256 digest bits in 43 safe characters.
    return "SSU_" + base64.urlsafe_b64encode(digest).rstrip(b"=").decode("ascii")


def iter_fasta(path: str | Path) -> Iterator[FastaRecord]:
    """Stream plain or gzip FASTA while validating every nucleotide record."""

    source = Path(path)
    with _open_fasta(source) as handle:
        yield from iter_fasta_lines(handle, label=str(source))


def iter_fasta_lines(lines: Iterable[str], label: str = "FASTA") -> Iterator[FastaRecord]:
    header: str | None = None
    chunks: list[str] = []
    record_number = 0
    for line_number, raw_line in enumerate(lines, 1):
        line = raw_line.rstrip("\r\n")
        if line.startswith(">"):
            if header is not None:
                sequence = "".join(chunks)
                validate_nucleotide_sequence(sequence, f"record {record_number} in {label}")
                yield FastaRecord(header, sequence)
            header = line[1:]
            record_number += 1
            chunks = []
            if not header:
                raise FastaFormatError(f"Empty FASTA header at {label}:{line_number}")
        else:
            if header is None:
                if line.strip():
                    raise FastaFormatError(f"Sequence before first header at {label}:{line_number}")
                continue
            chunks.append(line)
    if header is None:
        raise FastaFormatError(f"No FASTA records in {label}")
    sequence = "".join(chunks)
    validate_nucleotide_sequence(sequence, f"record {record_number} in {label}")
    yield FastaRecord(header, sequence)


def _parse_accession_coordinates(identifier: str) -> tuple[str, int | None, int | None]:
    parts = identifier.rsplit(".", 2)
    if len(parts) == 3 and parts[1].isdigit() and parts[2].isdigit():
        return parts[0], int(parts[1]), int(parts[2])
    return identifier, None, None


def parse_pr2_header(header: str) -> Pr2Header:
    """Parse the PR2 5.1.1 SSU ``taxo_long`` header contract."""

    text = header[1:] if header.startswith(">") else header
    fields = text.split("|")
    if len(fields) != 13:
        raise TaxonomyError(
            f"PR2 taxo_long header must have 13 pipe-delimited fields, found {len(fields)}"
        )
    source_identifier, gene, compartment = fields[:3]
    if not source_identifier or not gene:
        raise TaxonomyError("PR2 identifier and gene must be non-empty")
    coordinate_token, separator, orientation = source_identifier.rpartition("_")
    if not separator or not orientation:
        raise TaxonomyError(f"PR2 identifier lacks orientation suffix: {source_identifier}")
    accession, start, end = _parse_accession_coordinates(coordinate_token)
    if start is None or end is None or start < 0 or end < 0:
        raise TaxonomyError(f"PR2 identifier has invalid coordinates: {source_identifier}")
    taxonomy = tuple(field.strip() for field in fields[4:])
    if len(taxonomy) != len(PR2_RANKS) or not taxonomy[0]:
        raise TaxonomyError("PR2 header must provide a fixed nine-rank taxonomy with domain")
    return Pr2Header(
        source_identifier=source_identifier,
        accession=accession,
        start=start,
        end=end,
        orientation=orientation,
        gene=gene,
        compartment=compartment,
        taxonomy=taxonomy,
    )


def include_pr2_header(header: Pr2Header) -> bool:
    compartment = header.compartment.strip().lower()
    return header.taxonomy[0] == "Eukaryota" or compartment in PR2_COMPARTMENTS


def normalize_pr2_taxonomy(taxonomy: Sequence[str]) -> tuple[str, ...]:
    """Separate PR2 organellar domain suffixes from the host-domain taxonomy."""

    try:
        normalized = taxonomy_path(taxonomy)
    except ValueError as error:
        raise TaxonomyError(str(error)) from error
    if normalized[0].startswith("Eukaryota:"):
        normalized = ("Eukaryota", *normalized[1:])
    return normalized


def parse_silva_header(header: str) -> SilvaHeader:
    text = header[1:] if header.startswith(">") else header
    source_identifier, separator, lineage = text.partition(" ")
    taxonomy = tuple(part.strip() for part in lineage.split(";") if part.strip())
    if not source_identifier or not separator or not taxonomy:
        raise TaxonomyError("SILVA header must contain an identifier and semicolon taxonomy")
    accession, start, end = _parse_accession_coordinates(source_identifier)
    return SilvaHeader(source_identifier, accession, start, end, taxonomy)


def include_silva_header(header: SilvaHeader) -> bool:
    return bool(header.taxonomy and header.taxonomy[0] in PROKARYOTIC_DOMAINS)


def parse_silva_rank_table(lines: Iterable[str]) -> dict[str, str]:
    """Map canonical cumulative SILVA paths to their declared ranks."""

    rank_by_path: dict[str, str] = {}
    for line_number, raw_line in enumerate(lines, 1):
        line = raw_line.rstrip("\r\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            raise TaxonomyError(f"Malformed SILVA taxonomy row {line_number}")
        path = ";".join(part.strip() for part in fields[0].strip(";").split(";") if part.strip())
        rank = fields[2].strip().lower()
        if not path or not rank:
            raise TaxonomyError(f"Malformed SILVA taxonomy row {line_number}")
        previous = rank_by_path.get(path)
        if previous is not None and previous != rank:
            raise TaxonomyError(f"Conflicting SILVA rank for path {path!r}")
        rank_by_path[path] = rank
    return rank_by_path


def rank_silva_prokaryotic_path(
    taxonomy: Sequence[str] | str, rank_by_path: Mapping[str, str]
) -> tuple[str, ...]:
    """Project a SILVA lineage onto domain/phylum/class/order/family/genus."""

    lineage = (
        tuple(part.strip() for part in taxonomy.split(";") if part.strip())
        if isinstance(taxonomy, str)
        else tuple(str(part).strip() for part in taxonomy if str(part).strip())
    )
    if not lineage or lineage[0] not in PROKARYOTIC_DOMAINS:
        raise TaxonomyError("SILVA prokaryotic taxonomy must begin with Bacteria or Archaea")
    values = {rank: "" for rank in SILVA_PROKARYOTIC_RANKS}
    prefix: list[str] = []
    for taxon in lineage:
        prefix.append(taxon)
        rank = rank_by_path.get(";".join(prefix), "").lower()
        if rank in values and not values[rank]:
            values[rank] = taxon
    if not values["domain"]:
        raise TaxonomyError("SILVA taxonomy rank table did not identify the domain")
    return tuple(values[rank] for rank in SILVA_PROKARYOTIC_RANKS)


def load_silva_rank_table(path: str | Path) -> dict[str, str]:
    source = Path(path)
    if source.suffix.lower() == ".gz":
        with gzip.open(source, "rt", encoding="utf-8", newline=None) as handle:
            return parse_silva_rank_table(handle)
    with source.open("rt", encoding="utf-8", newline=None) as handle:
        return parse_silva_rank_table(handle)


def iter_curated_silva_records(
    fasta_path: str | Path,
    rank_table_path: str | Path,
    *,
    source_version: str,
) -> Iterator[PreparedSourceRecord]:
    """Yield rank-normalized bacterial and archaeal SILVA SSU records."""

    rank_by_path = load_silva_rank_table(rank_table_path)
    for record in iter_fasta(fasta_path):
        header = parse_silva_header(record.header)
        if not include_silva_header(header):
            continue
        yield PreparedSourceRecord(
            reference_source="SILVA",
            source_version=source_version,
            source_identifier=header.source_identifier,
            original_header=record.header,
            sequence=record.sequence,
            marker="16S",
            taxonomy=rank_silva_prokaryotic_path(header.taxonomy, rank_by_path),
            taxonomy_source="SILVA",
            assignment_method="native",
        )


def iter_curated_pr2_records(
    fasta_path: str | Path, *, source_version: str
) -> Iterator[PreparedSourceRecord]:
    """Yield PR2 eukaryotic and organellar SSU records with fixed-rank taxonomy."""

    marker_by_gene = {"16S_rRNA": "16S", "18S_rRNA": "18S"}
    for record in iter_fasta(fasta_path):
        header = parse_pr2_header(record.header)
        if not include_pr2_header(header):
            continue
        try:
            marker = marker_by_gene[header.gene]
        except KeyError as error:
            raise TaxonomyError(f"Unsupported PR2 SSU gene: {header.gene}") from error
        compartment = (
            header.compartment.strip().lower().replace("mitochondria", "mitochondrion")
            or "unknown"
        )
        yield PreparedSourceRecord(
            reference_source="PR2",
            source_version=source_version,
            source_identifier=header.source_identifier,
            original_header=record.header,
            sequence=record.sequence,
            marker=marker,
            taxonomy=normalize_pr2_taxonomy(header.taxonomy),
            taxonomy_source="PR2",
            compartment=compartment,
            assignment_method="native",
        )


def parse_img_identifier(header: str) -> ImgIdentifier:
    text = header[1:] if header.startswith(">") else header
    source_identifier = text.split(maxsplit=1)[0]
    match = _IMG_ID.match(source_identifier)
    if match is None:
        raise TaxonomyError(f"Not an IMG source identifier: {source_identifier!r}")
    return ImgIdentifier(source_identifier, match.group("taxon_oid"))


def iter_img_records(
    fasta_path: str | Path, marker: str, *, source_version: str
) -> Iterator[PreparedSourceRecord]:
    """Yield only IMG records, excluding every embedded legacy reference record."""

    if marker not in {"16S", "18S"}:
        raise TaxonomyError(f"Unsupported IMG marker: {marker}")
    for record in iter_fasta(fasta_path):
        if not record.header.startswith("IMG_"):
            continue
        parsed = parse_img_identifier(record.header)
        yield PreparedSourceRecord(
            reference_source="IMG",
            source_version=source_version,
            source_identifier=parsed.source_identifier,
            original_header=record.header,
            sequence=record.sequence,
            marker=marker,
            taxon_oid=parsed.taxon_oid,
        )


def _optional_coordinate(value: object, field: str, lower: float, upper: float) -> float | None:
    if value is None or (isinstance(value, str) and value.strip().lower() in {"", "na", "n/a", "null"}):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ReleaseValidationError(f"IMG {field} must be numeric or missing: {value!r}") from error
    if not math.isfinite(number) or not lower <= number <= upper:
        raise ReleaseValidationError(f"IMG {field} outside [{lower}, {upper}]: {number}")
    return number


def validate_privacy_columns(columns: Iterable[str]) -> None:
    normalized = {column.strip().lower() for column in columns}
    forbidden = sorted(
        column
        for column in normalized
        if any(token in re.split(r"[^a-z0-9]+", column) for token in FORBIDDEN_PII_TOKENS)
    )
    unexpected = sorted(normalized - set(IMG_LOCATION_COLUMNS))
    if forbidden:
        raise ReleaseValidationError(f"Forbidden PII columns: {', '.join(forbidden)}")
    if unexpected:
        raise ReleaseValidationError(
            f"IMG location output is restricted to {IMG_LOCATION_COLUMNS}; found {unexpected}"
        )


def clean_img_metadata_row(
    row: Mapping[str, object], *, corrections: list[dict[str, object]] | None = None
) -> ImgLocation:
    """Return only the privacy-reviewed IMG location allowlist."""

    normalized = {str(key).strip().lower(): value for key, value in row.items()}
    oid_value = normalized.get("taxon_oid", normalized.get("img genome id"))
    oid = "" if oid_value is None else str(oid_value).strip()
    if not oid.isdigit():
        raise ReleaseValidationError(f"IMG taxon_oid must contain only digits: {oid!r}")
    latitude_raw = normalized.get("latitude")
    longitude_raw = normalized.get("longitude")
    try:
        latitude = _optional_coordinate(latitude_raw, "latitude", -90.0, 90.0)
        longitude = _optional_coordinate(longitude_raw, "longitude", -180.0, 180.0)
    except ReleaseValidationError as original_error:
        if corrections is None:
            raise original_error
        # A plausible coordinate after swapping does not prove that the source
        # fields were transposed. Keep the record, but do not invent a location.
        latitude = None
        longitude = None
        corrections.append(
            {
                "taxon_oid": oid,
                "field": "latitude_longitude",
                "action": "dropped_invalid_pair",
                "reason": "invalid_or_out_of_range_coordinate",
            }
        )
    return ImgLocation(taxon_oid=oid, latitude=latitude, longitude=longitude)


def clean_img_metadata(
    rows: Iterable[Mapping[str, object]], *, corrections: list[dict[str, object]] | None = None
) -> tuple[ImgLocation, ...]:
    locations: dict[str, ImgLocation] = {}
    for row in rows:
        location = clean_img_metadata_row(row, corrections=corrections)
        previous = locations.get(location.taxon_oid)
        if previous is not None and previous != location:
            raise ReleaseValidationError(f"Conflicting IMG locations for taxon_oid {location.taxon_oid}")
        locations[location.taxon_oid] = location
    return tuple(locations[key] for key in sorted(locations))
