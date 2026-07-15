#!/usr/bin/env python3

import argparse
from pathlib import Path

from hit_processing import ExtractionRegion, extract_regions


def parse_seqmap(seqmap_file: str) -> list[ExtractionRegion]:
    regions: list[ExtractionRegion] = []
    for line_number, line in enumerate(Path(seqmap_file).read_text().splitlines(), start=1):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != 5:
            raise ValueError(
                f"Malformed seqmap row at {seqmap_file}:{line_number}: "
                f"expected 5 fields, found {len(fields)}"
            )
        subject = fields[0].split("_d_", maxsplit=1)[0]
        start, end = sorted((int(fields[1]), int(fields[2])))
        regions.append(
            ExtractionRegion(
                subject=subject,
                start=start,
                end=end,
                strand=fields[3],
                sequence_type=fields[4],
                is_assembled=fields[4] == "assembled",
            )
        )
    return regions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract sequences listed in a legacy SSUextract seqmap."
    )
    parser.add_argument("fasta")
    parser.add_argument("seqmap")
    parser.add_argument("output")
    parser.add_argument("minimum_length", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    records = extract_regions(
        args.fasta,
        parse_seqmap(args.seqmap),
        args.minimum_length,
    )
    with Path(args.output).open("w") as handle:
        for record in records:
            handle.write(f">{record.name}\n{record.sequence}\n")


if __name__ == "__main__":
    main()

