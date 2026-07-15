#!/usr/bin/env python3

from __future__ import annotations

import csv
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
]

META_FIELDS = ["sample", "model"]


@dataclass(frozen=True)
class CmHit:
    subject: str
    model: str
    model_from: int
    model_to: int
    sequence_from: int
    sequence_to: int
    strand: str
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
class ExtractionRegion:
    subject: str
    start: int
    end: int
    strand: str
    sequence_type: str
    is_assembled: bool


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


def read_model_length(model_file: str | Path) -> int:
    with Path(model_file).open() as handle:
        for line in handle:
            fields = line.split()
            if fields and fields[0] == "CLEN":
                return int(fields[1])
    raise ValueError(f"No CLEN field found in covariance model: {model_file}")


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

            strand = fields[9]
            if strand not in {"+", "-"}:
                raise ValueError(
                    f"Invalid strand at {tblout}:{line_number}: {strand!r}"
                )

            hits.append(
                CmHit(
                    subject=fields[0],
                    model=fields[2],
                    model_from=int(fields[5]),
                    model_to=int(fields[6]),
                    sequence_from=int(fields[7]),
                    sequence_to=int(fields[8]),
                    strand=strand,
                    included=fields[16] == "!",
                )
            )
    return hits


def _simple_region(hit: CmHit) -> ExtractionRegion:
    return ExtractionRegion(
        subject=hit.subject,
        start=hit.sequence_start,
        end=hit.sequence_end,
        strand=hit.strand,
        sequence_type="simple",
        is_assembled=False,
    )


def _is_full_model_hit(hit: CmHit, model_length: int) -> bool:
    return hit.model_start == 1 and hit.model_end == model_length


def _fragments_are_collinear(hits: Sequence[CmHit]) -> bool:
    if len({hit.strand for hit in hits}) != 1:
        return False

    model_endpoints: list[tuple[int, int]] = []
    sequence_endpoints: list[tuple[int, int]] = []
    for index, hit in enumerate(hits):
        start_marker = index * 2
        end_marker = start_marker + 1
        model_endpoints.extend(
            [
                (hit.model_from, start_marker),
                (hit.model_to, end_marker),
            ]
        )

        if hit.strand == "+":
            sequence_endpoints.extend(
                [
                    (hit.sequence_from, start_marker),
                    (hit.sequence_to, end_marker),
                ]
            )
        else:
            sequence_endpoints.extend(
                [
                    (-hit.sequence_from, start_marker),
                    (-hit.sequence_to, end_marker),
                ]
            )

    model_positions = [position for position, _ in model_endpoints]
    sequence_positions = [position for position, _ in sequence_endpoints]
    if len(set(model_positions)) != len(model_positions):
        return False
    if len(set(sequence_positions)) != len(sequence_positions):
        return False

    model_order = [marker for _, marker in sorted(model_endpoints)]
    sequence_order = [marker for _, marker in sorted(sequence_endpoints)]
    return model_order == sequence_order


def resolve_extraction_regions(
    hits: Iterable[CmHit], model_length: int
) -> list[ExtractionRegion]:
    by_subject: dict[str, list[CmHit]] = defaultdict(list)
    for hit in hits:
        if hit.included:
            by_subject[hit.subject].append(hit)

    regions: list[ExtractionRegion] = []
    for subject, subject_hits in by_subject.items():
        full_hits = [
            hit for hit in subject_hits if _is_full_model_hit(hit, model_length)
        ]
        partial_hits = [
            hit for hit in subject_hits if not _is_full_model_hit(hit, model_length)
        ]

        regions.extend(_simple_region(hit) for hit in full_hits)
        if len(partial_hits) == 1:
            regions.append(_simple_region(partial_hits[0]))
        elif len(partial_hits) > 1 and _fragments_are_collinear(partial_hits):
            regions.append(
                ExtractionRegion(
                    subject=subject,
                    start=min(hit.sequence_start for hit in partial_hits),
                    end=max(hit.sequence_end for hit in partial_hits),
                    strand=partial_hits[0].strand,
                    sequence_type="assembled",
                    is_assembled=True,
                )
            )
        else:
            regions.extend(_simple_region(hit) for hit in partial_hits)

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
        if len(sequence) >= minimum_length:
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
            fasta_handle.write(f">{record.name}\n{record.sequence}\n")

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
            writer.writerow(
                {
                    "name": record.name,
                    "sample": sample,
                    "model": model,
                    "length": len(record.sequence),
                    "coordinates": f"{region.start}-{region.end}",
                    "strand": region.strand,
                    "sequence_type": region.sequence_type,
                    "contig_name": region.subject,
                    "is_assembled": str(region.is_assembled),
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
