"""Typed contracts shared by SSU reference database build stages."""

from __future__ import annotations

from dataclasses import dataclass


PR2_RANKS = (
    "domain",
    "supergroup",
    "division",
    "subdivision",
    "class",
    "order",
    "family",
    "genus",
    "species",
)
SILVA_PROKARYOTIC_RANKS = ("domain", "phylum", "class", "order", "family", "genus")
PROKARYOTIC_DOMAINS = frozenset({"Bacteria", "Archaea"})
PR2_COMPARTMENTS = frozenset(
    {
        "nucleus",
        "nucleomorph",
        "plastid",
        "apicoplast",
        "chromatophore",
        "mitochondrion",
        "mitochondria",
    }
)
IMG_LOCATION_COLUMNS = ("taxon_oid", "latitude", "longitude")


class BuildError(RuntimeError):
    """Base class for deterministic database build failures."""


class SourceIntegrityError(BuildError):
    """A fetched source did not match its pinned size or digest."""


class FastaFormatError(BuildError):
    """A FASTA record is empty, malformed, or not nucleotide sequence."""


class TaxonomyError(BuildError):
    """A taxonomy header or explicit assignment is invalid."""


class ReleaseValidationError(BuildError):
    """The prepared release violates one or more release invariants."""


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str


@dataclass(frozen=True)
class Pr2Header:
    source_identifier: str
    accession: str
    start: int
    end: int
    orientation: str
    gene: str
    compartment: str
    taxonomy: tuple[str, ...]

    @property
    def coordinates(self) -> tuple[int, int]:
        return self.start, self.end


@dataclass(frozen=True)
class SilvaHeader:
    source_identifier: str
    accession: str
    start: int | None
    end: int | None
    taxonomy: tuple[str, ...]

    @property
    def coordinates(self) -> tuple[int, int] | None:
        if self.start is None or self.end is None:
            return None
        return self.start, self.end


@dataclass(frozen=True)
class ImgIdentifier:
    source_identifier: str
    taxon_oid: str


@dataclass(frozen=True)
class PreparedSourceRecord:
    reference_source: str
    source_version: str
    source_identifier: str
    original_header: str
    sequence: str
    marker: str
    taxonomy: tuple[str, ...] = ()
    taxonomy_source: str = ""
    compartment: str = ""
    assignment_method: str = "native"
    evidence: str = ""
    taxon_oid: str | None = None


@dataclass(frozen=True)
class SequenceRecord:
    sequence_id: str
    sequence: str
    markers: tuple[str, ...] = ()


@dataclass(frozen=True)
class SourceRecord:
    source_record_id: str
    sequence_id: str
    reference_source: str
    source_version: str
    source_identifier: str
    original_header: str
    marker: str
    taxon_oid: str | None = None


@dataclass(frozen=True)
class TaxonomyAssignment:
    taxonomy_assignment_id: str
    source_record_id: str
    sequence_id: str
    taxonomy: tuple[str, ...]
    taxonomy_source: str
    domain: str
    compartment: str
    assignment_method: str
    evidence: str = ""


@dataclass(frozen=True)
class PreferredTaxonomy:
    sequence_id: str
    reference_source: str
    taxonomy: tuple[str, ...]
    taxonomy_source: str
    domain: str
    compartment: str
    assignment_method: str
    cross_domain_conflict: bool
    taxonomy_alternatives: str = ""


@dataclass(frozen=True)
class ImgLocation:
    taxon_oid: str
    latitude: float | None
    longitude: float | None


@dataclass(frozen=True)
class DatabaseModel:
    sequences: tuple[SequenceRecord, ...]
    source_records: tuple[SourceRecord, ...]
    taxonomy_assignments: tuple[TaxonomyAssignment, ...]
    preferred_taxonomy: tuple[PreferredTaxonomy, ...]
