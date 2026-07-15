#!/usr/bin/env python3

import argparse

from hit_processing import (
    extract_regions,
    parse_cmsearch_tblout,
    read_model_length,
    resolve_extraction_regions,
    write_extraction_outputs,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve cmsearch hits and extract 1-based inclusive FASTA intervals."
    )
    parser.add_argument("--cmsearch", required=True)
    parser.add_argument("--model-file", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--minimum-length", type=int, required=True)
    parser.add_argument("--fasta-output", required=True)
    parser.add_argument("--hits-output", required=True)
    parser.add_argument("--metadata-output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    model_length = read_model_length(args.model_file)
    hits = parse_cmsearch_tblout(args.cmsearch)
    regions = resolve_extraction_regions(hits, model_length)
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

