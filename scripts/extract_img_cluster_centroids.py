#!/usr/bin/env python3
"""Extract cluster centroid sequences from a composite eukcensus FASTA."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Iterable, Sequence

import build_database_release as builder


def read_centroids(lines: Iterable[str]) -> dict[str, str]:
    csv.field_size_limit(sys.maxsize)
    reader = csv.DictReader(lines, delimiter="\t")
    required = {"cluster_id", "centroid"}
    missing = required - set(reader.fieldnames or ())
    if missing:
        raise ValueError(f"cluster table lacks columns: {sorted(missing)}")
    by_header: dict[str, str] = {}
    cluster_ids: set[str] = set()
    for row_number, row in enumerate(reader, 2):
        cluster_id = (row.get("cluster_id") or "").strip()
        centroid = (row.get("centroid") or "").strip()
        if not cluster_id or not centroid:
            raise ValueError(f"cluster row {row_number} lacks cluster_id or centroid")
        if cluster_id in cluster_ids:
            raise ValueError(f"duplicate cluster_id {cluster_id!r}")
        if centroid in by_header:
            raise ValueError(f"centroid {centroid!r} belongs to more than one cluster")
        cluster_ids.add(cluster_id)
        by_header[centroid] = cluster_id
    return by_header


def extract_centroids(
    cluster_path: str | Path, fasta_path: str | Path, output_path: str | Path
) -> dict[str, int]:
    with Path(cluster_path).open(newline="", encoding="utf-8") as handle:
        centroids = read_centroids(handle)
    sequences: dict[str, str] = {}
    for record in builder.iter_fasta(fasta_path):
        cluster_id = centroids.get(record.header)
        if cluster_id is None:
            continue
        if cluster_id in sequences:
            raise ValueError(f"multiple FASTA records match cluster {cluster_id!r}")
        sequences[cluster_id] = builder.normalize_sequence_for_hashing(record.sequence)
    missing = sorted(set(centroids.values()) - sequences.keys())
    if missing:
        preview = ", ".join(missing[:10])
        raise ValueError(f"missing FASTA sequence for {len(missing)} centroids: {preview}")
    records = [builder.FastaRecord(cluster_id, sequences[cluster_id]) for cluster_id in sequences]
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="ascii", newline="\n") as handle:
        for record in sorted(records, key=lambda item: item.header):
            handle.write(f">{record.header}\n{record.sequence}\n")
    return {"clusters": len(centroids), "centroids_written": len(records)}


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    report = extract_centroids(args.clusters, args.fasta, args.output)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
