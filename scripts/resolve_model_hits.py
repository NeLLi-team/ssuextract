#!/usr/bin/env python3

from __future__ import annotations

import argparse

from hit_processing import (
    parse_cmsearch_tblout,
    resolve_competing_model_hits,
    write_accepted_hits,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve overlapping hits from homologous SSU covariance models."
    )
    parser.add_argument("--cmsearch", action="append", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    hits = [
        hit
        for path in args.cmsearch
        for hit in parse_cmsearch_tblout(path)
    ]
    write_accepted_hits(resolve_competing_model_hits(hits), args.output)


if __name__ == "__main__":
    main()
