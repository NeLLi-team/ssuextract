#!/usr/bin/env python3
"""Extract cluster centroid sequences from a composite eukcensus FASTA."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence

from database_sources import iter_fasta, normalize_sequence_for_hashing
from img_classification_data import parse_clusters


def read_centroids(lines: Iterable[str]) -> dict[str, str]:
    """Return centroid headers only for clusters with IMG members."""
    by_header: dict[str, str] = {}
    for cluster in parse_clusters(lines):
        if not cluster.img_members:
            continue
        if cluster.centroid in by_header:
            raise ValueError(
                f"centroid {cluster.centroid!r} belongs to more than one cluster"
            )
        by_header[cluster.centroid] = cluster.cluster_id
    return by_header


def extract_centroids(
    cluster_path: str | Path, fasta_path: str | Path, output_path: str | Path
) -> dict[str, int]:
    with Path(cluster_path).open(newline="", encoding="utf-8") as handle:
        centroids = read_centroids(handle)
    sequences: dict[str, str] = {}
    for record in iter_fasta(fasta_path):
        cluster_id = centroids.get(record.header)
        if cluster_id is None:
            continue
        if cluster_id in sequences:
            raise ValueError(f"multiple FASTA records match cluster {cluster_id!r}")
        sequences[cluster_id] = normalize_sequence_for_hashing(record.sequence)
    missing = sorted(set(centroids.values()) - sequences.keys())
    if missing:
        preview = ", ".join(missing[:10])
        raise ValueError(f"missing FASTA sequence for {len(missing)} centroids: {preview}")
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="ascii", newline="\n") as handle:
        for cluster_id, sequence in sorted(sequences.items()):
            handle.write(f">{cluster_id}\n{sequence}\n")
    return {"clusters": len(centroids), "centroids_written": len(sequences)}


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
