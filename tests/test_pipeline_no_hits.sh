#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(git rev-parse --show-toplevel)
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/ssuextract-no-hit.XXXXXX")
trap 'rm -rf "${test_dir}"' EXIT

mkdir -p "${test_dir}/input"
printf '>no_hit description retained\nACGTACGTACGT\n' > "${test_dir}/input/no_hit.fna"
mkdir -p "${test_dir}/database"
printf '>reference;Bacteria\nACGTACGTACGTACGT\n' > "${test_dir}/database/reference.fna"
makeblastdb \
    -in "${test_dir}/database/reference.fna" \
    -dbtype nucl \
    -out "${test_dir}/database/silva-138-1_pr2-4-12" \
    > "${test_dir}/makeblastdb.stdout" \
    2> "${test_dir}/makeblastdb.stderr"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

nextflow \
    -log "${test_dir}/nextflow.log" \
    run "${repo_dir}/main.nf" \
    -ansi-log false \
    -work-dir "${test_dir}/work" \
    --querydir "${test_dir}/input" \
    --modeldir "${repo_dir}/resources/models" \
    --database_path "${test_dir}/database" \
    --outdir "${test_dir}/out" \
    --min_extract_length 0 \
    --threads_per_job 1 \
    > "${test_dir}/stdout.txt" \
    2> "${test_dir}/stderr.txt"

test "$(wc -l < "${test_dir}/out/cmsearch_summary.tsv")" -eq 1
test "$(wc -l < "${test_dir}/out/cmsearch_summary.tab")" -eq 2
test "$(sed -n '2p' "${test_dir}/out/cmsearch_summary.tab")" = "no_hit"
test "$(find "${test_dir}/out/stats" -name '*.hits.tsv' -type f | wc -l)" -eq 2
test "$(find "${test_dir}/out/out" "${test_dir}/out/stats" -name '*.fna' -type f | wc -l)" -eq 0
