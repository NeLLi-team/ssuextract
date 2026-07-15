#!/usr/bin/env bash
set -euo pipefail

marker=${1:?usage: classify_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE CALIBRATION OUTPUT_DIRECTORY SEARCH_THREADS}
clusters=${2:?usage: classify_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE CALIBRATION OUTPUT_DIRECTORY SEARCH_THREADS}
source_fasta=${3:?usage: classify_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE CALIBRATION OUTPUT_DIRECTORY SEARCH_THREADS}
curated_profile=${4:?usage: classify_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE CALIBRATION OUTPUT_DIRECTORY SEARCH_THREADS}
calibration=${5:?usage: classify_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE CALIBRATION OUTPUT_DIRECTORY SEARCH_THREADS}
output_directory=${6:?usage: classify_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE CALIBRATION OUTPUT_DIRECTORY SEARCH_THREADS}
search_threads=${7:?usage: classify_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE CALIBRATION OUTPUT_DIRECTORY SEARCH_THREADS}

if [[ "${marker}" != 16S && "${marker}" != 18S ]]; then
    printf 'unsupported marker: %s\n' "${marker}" >&2
    exit 2
fi
if [[ ! "${search_threads}" =~ ^[1-8]$ ]]; then
    printf 'search threads must be an integer from 1 through 8\n' >&2
    exit 2
fi

mkdir -p "${output_directory}"
queries="${output_directory}/centroids.fna"
blast_output="${output_directory}/centroids.m8"
taxonomy="${curated_profile}/tables/preferred_taxonomy.parquet"
search_provenance="${output_directory}/search_provenance.json"

python scripts/img_search_provenance.py validate \
    --marker "${marker}" \
    --curated-profile "${curated_profile}" \
    --clusters "${clusters}" \
    --source-fasta "${source_fasta}" \
    --centroids-fasta "${queries}" \
    --blast-m8 "${blast_output}" \
    --sidecar "${search_provenance}" \
    --threads "${search_threads}"

python scripts/classify_img_clusters.py \
    --blast "${blast_output}" \
    --blast-fetch-targets 501 \
    --taxonomy "${taxonomy}" \
    --clusters "${clusters}" \
    --marker "${marker}" \
    --calibration "${calibration}" \
    --search-provenance "${search_provenance}" \
    --propagation-rank-cap 0 \
    --assignments-tsv "${output_directory}/assignments.tsv" \
    --assignments-jsonl "${output_directory}/assignments.jsonl" \
    --outcomes-jsonl "${output_directory}/outcomes.jsonl" \
    --qc-json "${output_directory}/qc.json"

test -s "${output_directory}/assignments.tsv"
test -s "${output_directory}/outcomes.jsonl"
test -s "${output_directory}/qc.json"
