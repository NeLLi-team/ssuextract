#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from Bio import SeqIO
from Bio.Seq import Seq


HIT_FIELDS = [
    "name",
    "sample",
    "model",
    "length",
    "coordinates",
    "strand",
    "sequence_type",
    "contig_name",
    "is_assembled",
    "hit_length",
    "model_coverage",
    "fragment_count",
    "component_coordinates",
]

META_FIELDS = ["sample", "model"]

ACCEPTED_HIT_FIELDS = [
    "subject",
    "model",
    "model_accession",
    "model_from",
    "model_to",
    "sequence_from",
    "sequence_to",
    "strand",
    "bit_score",
    "e_value",
]

# Explicit allowlist for the two bundled Rfam models whose overlapping hits are
# resolved. This is not a complete Rfam clan registry.
RESOLVED_SSU_MODELS = {
    "RF00177": "SSU_rRNA_bacteria",
    "RF01960": "SSU_rRNA_eukarya",
}
RESOLVED_SSU_ACCESSIONS_BY_NAME = {
    name: accession for accession, name in RESOLVED_SSU_MODELS.items()
}


@dataclass(frozen=True)
class CmHit:
    subject: str
    model: str
    model_accession: str
    model_from: int
    model_to: int
    sequence_from: int
    sequence_to: int
    strand: str
    bit_score: float
    e_value: float
    included: bool

    @property
    def model_start(self) -> int:
        return min(self.model_from, self.model_to)

    @property
    def model_end(self) -> int:
        return max(self.model_from, self.model_to)

    @property
    def sequence_start(self) -> int:
        return min(self.sequence_from, self.sequence_to)

    @property
    def sequence_end(self) -> int:
        return max(self.sequence_from, self.sequence_to)


@dataclass(frozen=True)
class CmModel:
    name: str
    accession: str
    length: int


@dataclass(frozen=True)
class ExtractionRegion:
    subject: str
    start: int
    end: int
    strand: str
    sequence_type: str
    is_assembled: bool
    hit_length: int | None = None
    model_coverage: float = 0.0
    component_coordinates: tuple[tuple[int, int], ...] = ()

    def __post_init__(self) -> None:
        if self.hit_length is None:
            object.__setattr__(self, "hit_length", self.end - self.start + 1)
        if not self.component_coordinates:
            object.__setattr__(
                self, "component_coordinates", ((self.start, self.end),)
            )

    @property
    def fragment_count(self) -> int:
        return len(self.component_coordinates)


@dataclass(frozen=True)
class ExtractedRecord:
    region: ExtractionRegion
    sequence: str

    @property
    def name(self) -> str:
        region = self.region
        return (
            f"{region.subject}|{region.start}-{region.end}|"
            f"strand_{region.strand}|{region.sequence_type}"
        )


def _validated_model_accession(
    name: str, accession: str | None, location: str
) -> str:
    expected_accession = RESOLVED_SSU_ACCESSIONS_BY_NAME.get(name)
    if expected_accession is not None and accession != expected_accession:
        observed = accession if accession is not None else "missing"
        raise ValueError(
            f"Bundled SSU model {name!r} has accession {observed!r} at "
            f"{location}; expected {expected_accession!r}"
        )
    if accession in RESOLVED_SSU_MODELS:
        expected_name = RESOLVED_SSU_MODELS[accession]
        if name != expected_name:
            raise ValueError(
                f"Bundled SSU accession {accession!r} has model name {name!r} "
                f"at {location}; expected {expected_name!r}"
            )
    return accession if accession is not None else name


def read_covariance_model(model_file: str | Path) -> CmModel:
    values: dict[str, str] = {}
    model_records = 0
    first_record_complete = False
    with Path(model_file).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            if fields[0].startswith("INFERNAL"):
                model_records += 1
            if fields[0] == "//":
                first_record_complete = True
                continue
            # The embedded HMMER3 filter follows the first // and repeats NAME/ACC.
            if first_record_complete:
                continue
            if fields[0] not in {"NAME", "ACC", "CLEN"}:
                continue
            if len(fields) != 2:
                raise ValueError(
                    f"Malformed {fields[0]} field at {model_file}:{line_number}"
                )
            if fields[0] in values:
                raise ValueError(
                    f"Duplicate {fields[0]} field at {model_file}:{line_number}"
                )
            values[fields[0]] = fields[1]

    if model_records != 1:
        raise ValueError(
            f"Covariance-model file {model_file} contains {model_records} models; "
            "exactly one is required"
        )
    if not first_record_complete:
        raise ValueError(f"Covariance model {model_file} lacks its // terminator")
    missing = sorted({"NAME", "CLEN"} - values.keys())
    if missing:
        raise ValueError(
            f"Covariance model {model_file} lacks required fields: "
            + ", ".join(missing)
        )
    try:
        length = int(values["CLEN"])
    except ValueError as error:
        raise ValueError(
            f"Invalid CLEN field in covariance model {model_file}: "
            f"{values['CLEN']!r}"
        ) from error
    if length <= 0:
        raise ValueError(
            f"Invalid CLEN field in covariance model {model_file}: {length}"
        )
    accession = _validated_model_accession(
        values["NAME"], values.get("ACC"), str(model_file)
    )
    if accession in RESOLVED_SSU_MODELS and Path(model_file).stem != accession:
        raise ValueError(
            f"Bundled SSU accession {accession!r} must use filename "
            f"{accession}.cm, not {Path(model_file).name!r}"
        )
    return CmModel(values["NAME"], accession, length)


def _build_cm_hit(
    *,
    subject: str,
    model: str,
    model_accession: str | None,
    model_from: str,
    model_to: str,
    sequence_from: str,
    sequence_to: str,
    strand: str,
    bit_score: str,
    e_value: str,
    included: bool,
    location: str,
    label: str,
) -> CmHit:
    if not subject or not model or model_accession == "":
        raise ValueError(f"{label} identity is empty at {location}")
    accession = _validated_model_accession(model, model_accession, location)
    try:
        parsed_model_from = int(model_from)
        parsed_model_to = int(model_to)
        parsed_sequence_from = int(sequence_from)
        parsed_sequence_to = int(sequence_to)
        parsed_bit_score = float(bit_score)
        parsed_e_value = float(e_value)
    except ValueError as error:
        raise ValueError(f"{label} numeric field is invalid at {location}") from error
    if min(
        parsed_model_from,
        parsed_model_to,
        parsed_sequence_from,
        parsed_sequence_to,
    ) <= 0:
        raise ValueError(f"{label} coordinate is not positive at {location}")
    if strand not in {"+", "-"}:
        raise ValueError(f"{label} strand is invalid at {location}")
    if not math.isfinite(parsed_bit_score):
        raise ValueError(f"{label} bit score is invalid at {location}")
    if not math.isfinite(parsed_e_value) or parsed_e_value < 0:
        raise ValueError(f"{label} E-value is invalid at {location}")
    return CmHit(
        subject=subject,
        model=model,
        model_accession=accession,
        model_from=parsed_model_from,
        model_to=parsed_model_to,
        sequence_from=parsed_sequence_from,
        sequence_to=parsed_sequence_to,
        strand=strand,
        bit_score=parsed_bit_score,
        e_value=parsed_e_value,
        included=included,
    )


def parse_cmsearch_tblout(tblout: str | Path) -> list[CmHit]:
    hits: list[CmHit] = []
    with Path(tblout).open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split(maxsplit=17)
            if len(fields) < 17:
                raise ValueError(
                    f"Malformed cmsearch tblout row at {tblout}:{line_number}: "
                    f"expected at least 17 fields, found {len(fields)}"
                )

            hits.append(
                _build_cm_hit(
                    subject=fields[0],
                    model=fields[2],
                    model_accession=None if fields[3] == "-" else fields[3],
                    model_from=fields[5],
                    model_to=fields[6],
                    sequence_from=fields[7],
                    sequence_to=fields[8],
                    strand=fields[9],
                    bit_score=fields[14],
                    e_value=fields[15],
                    included=fields[16] == "!",
                    location=f"{tblout}:{line_number}",
                    label="cmsearch hit",
                )
            )
    return hits


def _same_competing_clan(left: CmHit, right: CmHit) -> bool:
    return (
        left.model_accession != right.model_accession
        and left.model_accession in RESOLVED_SSU_MODELS
        and right.model_accession in RESOLVED_SSU_MODELS
    )


def _sequence_intervals_overlap(left: CmHit, right: CmHit) -> bool:
    return (
        left.subject == right.subject
        and left.strand == right.strand
        and left.sequence_start <= right.sequence_end
        and right.sequence_start <= left.sequence_end
    )


def _hit_output_order(hit: CmHit) -> tuple[object, ...]:
    return (
        hit.model_accession,
        hit.subject,
        hit.sequence_start,
        hit.sequence_end,
        hit.strand,
        hit.model_start,
        hit.model_end,
    )


def resolve_competing_model_hits(hits: Iterable[CmHit]) -> list[CmHit]:
    """Retain the best same-clan model explanation for each overlapping locus."""

    included = [hit for hit in hits if hit.included]
    if len(set(included)) != len(included):
        raise ValueError("Duplicate included cmsearch hit")

    accepted: list[CmHit] = []
    for candidate in sorted(
        included,
        key=lambda hit: (hit.e_value, -hit.bit_score, *_hit_output_order(hit)),
    ):
        unresolved_overlaps = [
            winner
            for winner in accepted
            if candidate.model_accession != winner.model_accession
            and _sequence_intervals_overlap(candidate, winner)
            and not _same_competing_clan(candidate, winner)
        ]
        if unresolved_overlaps:
            models = sorted(
                {
                    candidate.model_accession,
                    *(hit.model_accession for hit in unresolved_overlaps),
                }
            )
            raise ValueError(
                "Cannot resolve an overlap between models without an explicit "
                "competition rule: "
                + ", ".join(models)
            )
        competitors = [
            winner
            for winner in accepted
            if _same_competing_clan(candidate, winner)
            and _sequence_intervals_overlap(candidate, winner)
        ]
        if not competitors:
            accepted.append(candidate)
            continue
        if any(
            candidate.e_value == winner.e_value
            and candidate.bit_score == winner.bit_score
            for winner in competitors
        ):
            models = sorted(
                {candidate.model_accession, *(winner.model_accession for winner in competitors)}
            )
            raise ValueError(
                "Indistinguishable competing SSU model hits for "
                f"{candidate.subject}:{candidate.sequence_start}-{candidate.sequence_end} "
                f"on strand {candidate.strand}: {', '.join(models)}"
            )

    return sorted(accepted, key=_hit_output_order)


def write_accepted_hits(hits: Iterable[CmHit], output: str | Path) -> None:
    with Path(output).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=ACCEPTED_HIT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for hit in hits:
            writer.writerow(
                {
                    "subject": hit.subject,
                    "model": hit.model,
                    "model_accession": hit.model_accession,
                    "model_from": hit.model_from,
                    "model_to": hit.model_to,
                    "sequence_from": hit.sequence_from,
                    "sequence_to": hit.sequence_to,
                    "strand": hit.strand,
                    "bit_score": format(hit.bit_score, ".17g"),
                    "e_value": format(hit.e_value, ".17g"),
                }
            )


def read_accepted_hits(path: str | Path) -> list[CmHit]:
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ACCEPTED_HIT_FIELDS:
            raise ValueError(
                f"Unexpected accepted-hit columns in {path}: {reader.fieldnames}"
            )
        hits: list[CmHit] = []
        for row in reader:
            location = f"{path}:{reader.line_num}"
            hits.append(
                _build_cm_hit(
                    subject=row["subject"],
                    model=row["model"],
                    model_accession=row["model_accession"],
                    model_from=row["model_from"],
                    model_to=row["model_to"],
                    sequence_from=row["sequence_from"],
                    sequence_to=row["sequence_to"],
                    strand=row["strand"],
                    bit_score=row["bit_score"],
                    e_value=row["e_value"],
                    included=True,
                    location=location,
                    label="Accepted hit",
                )
            )
    if len(set(hits)) != len(hits):
        raise ValueError(f"Accepted-hit table contains duplicate rows: {path}")
    return hits


def select_model_hits(hits: Iterable[CmHit], model: CmModel) -> list[CmHit]:
    selected: list[CmHit] = []
    for hit in hits:
        accession_matches = hit.model_accession == model.accession
        name_matches = hit.model == model.name
        if accession_matches != name_matches:
            raise ValueError(
                "Accepted-hit model identity disagrees with covariance model "
                f"{model.name!r}/{model.accession!r}: "
                f"{hit.model!r}/{hit.model_accession!r}"
            )
        if accession_matches:
            selected.append(hit)
    return selected


def _union_length(intervals: Iterable[tuple[int, int]]) -> int:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return sum(end - start + 1 for start, end in merged)


def _region_from_hits(hits: Sequence[CmHit], model_length: int) -> ExtractionRegion:
    sequence_intervals = tuple(
        (hit.sequence_start, hit.sequence_end) for hit in hits
    )
    model_support = _union_length(
        (hit.model_start, hit.model_end) for hit in hits
    )
    assembled = len(hits) > 1
    return ExtractionRegion(
        subject=hits[0].subject,
        start=min(start for start, _ in sequence_intervals),
        end=max(end for _, end in sequence_intervals),
        strand=hits[0].strand,
        sequence_type="assembled" if assembled else "simple",
        is_assembled=assembled,
        hit_length=_union_length(sequence_intervals),
        model_coverage=model_support / model_length,
        component_coordinates=sequence_intervals,
    )


def _is_full_model_hit(hit: CmHit, model_length: int) -> bool:
    return hit.model_start == 1 and hit.model_end == model_length


def _sequence_order(hit: CmHit) -> tuple[int, int]:
    if hit.strand == "+":
        return (hit.sequence_start, hit.sequence_end)
    return (-hit.sequence_end, -hit.sequence_start)


def _fragments_can_merge(left: CmHit, right: CmHit) -> bool:
    if left.model_end >= right.model_start:
        return False
    model_gap = right.model_start - left.model_end - 1

    if left.strand == "+":
        if left.sequence_end >= right.sequence_start:
            return False
        sequence_gap = right.sequence_start - left.sequence_end - 1
    else:
        if right.sequence_end >= left.sequence_start:
            return False
        sequence_gap = left.sequence_start - right.sequence_end - 1
    return sequence_gap <= model_gap


def _full_hit_blocks_merge(
    left: CmHit, right: CmHit, full_hits: Sequence[CmHit]
) -> bool:
    start = min(left.sequence_start, right.sequence_start)
    end = max(left.sequence_end, right.sequence_end)
    return any(
        hit.sequence_start <= end and hit.sequence_end >= start for hit in full_hits
    )


def resolve_extraction_regions(
    hits: Iterable[CmHit], model_length: int
) -> list[ExtractionRegion]:
    if model_length <= 0:
        raise ValueError("model_length must be positive")

    by_subject: dict[str, list[CmHit]] = defaultdict(list)
    for hit in hits:
        if hit.included:
            if hit.model_start < 1 or hit.model_end > model_length:
                raise ValueError(
                    f"Model interval for {hit.subject} is outside 1-{model_length}: "
                    f"{hit.model_start}-{hit.model_end}"
                )
            by_subject[hit.subject].append(hit)

    regions: list[ExtractionRegion] = []
    for subject, subject_hits in by_subject.items():
        full_hits = [
            hit for hit in subject_hits if _is_full_model_hit(hit, model_length)
        ]
        regions.extend(_region_from_hits([hit], model_length) for hit in full_hits)

        for strand in ("+", "-"):
            partial_hits = sorted(
                (
                    hit
                    for hit in subject_hits
                    if hit.strand == strand
                    and not _is_full_model_hit(hit, model_length)
                ),
                key=_sequence_order,
            )
            current: list[CmHit] = []
            for hit in partial_hits:
                if current and (
                    not _fragments_can_merge(current[-1], hit)
                    or _full_hit_blocks_merge(current[-1], hit, full_hits)
                ):
                    regions.append(_region_from_hits(current, model_length))
                    current = []
                current.append(hit)
            if current:
                regions.append(_region_from_hits(current, model_length))

    return sorted(
        regions,
        key=lambda region: (
            region.subject,
            region.start,
            region.end,
            region.strand,
            region.sequence_type,
        ),
    )


def extract_regions(
    fasta_file: str | Path,
    regions: Iterable[ExtractionRegion],
    minimum_length: int,
) -> list[ExtractedRecord]:
    if minimum_length < 0:
        raise ValueError("minimum_length must be non-negative")

    region_list = list(regions)
    target_ids = {region.subject for region in region_list}
    sequences: dict[str, str] = {}
    with Path(fasta_file).open() as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            if record.id not in target_ids:
                continue
            if record.id in sequences:
                raise ValueError(f"Duplicate FASTA identifier: {record.id}")
            sequences[record.id] = str(record.seq)

    missing = sorted(target_ids - sequences.keys())
    if missing:
        raise ValueError(
            "cmsearch subjects missing from FASTA: " + ", ".join(missing)
        )

    extracted: list[ExtractedRecord] = []
    for region in region_list:
        source = sequences[region.subject]
        if region.start < 1 or region.end > len(source) or region.start > region.end:
            raise ValueError(
                f"Invalid 1-based inclusive interval for {region.subject}: "
                f"{region.start}-{region.end} (sequence length {len(source)})"
            )

        sequence = source[region.start - 1 : region.end]
        if region.strand == "-":
            sequence = str(Seq(sequence).reverse_complement())
        if region.hit_length >= minimum_length:
            extracted.append(ExtractedRecord(region=region, sequence=sequence))
    return extracted


def write_extraction_outputs(
    records: Iterable[ExtractedRecord],
    sample: str,
    model: str,
    fasta_output: str | Path,
    hits_output: str | Path,
    metadata_output: str | Path,
) -> None:
    record_list = list(records)
    with Path(fasta_output).open("w") as fasta_handle:
        for record in record_list:
            output_name = f"{sample}|{model}|{record.name}"
            fasta_handle.write(f">{output_name}\n{record.sequence}\n")

    with Path(hits_output).open("w", newline="") as hits_handle:
        writer = csv.DictWriter(
            hits_handle,
            fieldnames=HIT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for record in record_list:
            region = record.region
            output_name = f"{sample}|{model}|{record.name}"
            writer.writerow(
                {
                    "name": output_name,
                    "sample": sample,
                    "model": model,
                    "length": len(record.sequence),
                    "coordinates": f"{region.start}-{region.end}",
                    "strand": region.strand,
                    "sequence_type": region.sequence_type,
                    "contig_name": region.subject,
                    "is_assembled": str(region.is_assembled),
                    "hit_length": region.hit_length,
                    "model_coverage": format(region.model_coverage, ".6g"),
                    "fragment_count": region.fragment_count,
                    "component_coordinates": ",".join(
                        f"{start}-{end}"
                        for start, end in region.component_coordinates
                    ),
                }
            )

    with Path(metadata_output).open("w", newline="") as metadata_handle:
        writer = csv.DictWriter(
            metadata_handle,
            fieldnames=META_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"sample": sample, "model": model})
