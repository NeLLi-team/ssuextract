#!/usr/bin/env python3

import argparse
from pathlib import Path

from hit_processing import (
    extract_regions,
    read_covariance_model,
    read_accepted_hits,
    resolve_extraction_regions,
    select_model_hits,
    write_extraction_outputs,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve cmsearch hits and extract 1-based inclusive FASTA intervals."
    )
    parser.add_argument("--model-file", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--accepted-hits", required=True)
    parser.add_argument("--minimum-length", type=int, required=True)
    parser.add_argument("--fasta-output", required=True)
    parser.add_argument("--hits-output", required=True)
    parser.add_argument("--metadata-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if Path(args.model_file).stem != args.model:
        raise ValueError(
            f"Model label {args.model!r} does not match covariance-model file "
            f"{Path(args.model_file).name!r}"
        )
    model = read_covariance_model(args.model_file)
    hits = select_model_hits(read_accepted_hits(args.accepted_hits), model)
    regions = resolve_extraction_regions(hits, model.length)
    records = extract_regions(args.fasta, regions, args.minimum_length)
    write_extraction_outputs(
        records=records,
        sample=args.sample,
        model=args.model,
        fasta_output=args.fasta_output,
        hits_output=args.hits_output,
        metadata_output=args.metadata_output,
    )


if __name__ == "__main__":
    main()
