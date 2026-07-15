#!/usr/bin/env bash
set -euo pipefail

marker=${1:?usage: search_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE OUTPUT_DIRECTORY THREADS}
clusters=${2:?usage: search_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE OUTPUT_DIRECTORY THREADS}
source_fasta=${3:?usage: search_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE OUTPUT_DIRECTORY THREADS}
curated_profile=${4:?usage: search_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE OUTPUT_DIRECTORY THREADS}
output_directory=${5:?usage: search_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE OUTPUT_DIRECTORY THREADS}
threads=${6:?usage: search_img_marker.sh MARKER CLUSTERS SOURCE_FASTA CURATED_PROFILE OUTPUT_DIRECTORY THREADS}

if [[ "${marker}" != 16S && "${marker}" != 18S ]]; then
    printf 'unsupported marker: %s\n' "${marker}" >&2
    exit 2
fi
if [[ ! "${threads}" =~ ^[1-8]$ ]]; then
    printf 'threads must be an integer from 1 through 8\n' >&2
    exit 2
fi

mkdir -p "${output_directory}"
queries="${output_directory}/centroids.fna"
blast_output="${output_directory}/centroids.m8"
search_provenance="${output_directory}/search_provenance.json"

# A fresh search invalidates both its completion marker and prior classifications.
rm -f \
    "${search_provenance}" \
    "${output_directory}/assignments.tsv" \
    "${output_directory}/assignments.jsonl" \
    "${output_directory}/outcomes.jsonl" \
    "${output_directory}/qc.json"

python scripts/extract_img_cluster_centroids.py \
    --clusters "${clusters}" \
    --fasta "${source_fasta}" \
    --output "${queries}" \
    --report "${output_directory}/centroid-extraction.json"

python scripts/img_search_provenance.py run \
    --marker "${marker}" \
    --curated-profile "${curated_profile}" \
    --clusters "${clusters}" \
    --source-fasta "${source_fasta}" \
    --centroids-fasta "${queries}" \
    --blast-m8 "${blast_output}" \
    --sidecar "${search_provenance}" \
    --threads "${threads}"

test -s "${queries}"
test -s "${blast_output}"
test -s "${search_provenance}"
