#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(git rev-parse --show-toplevel)
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/ssuextract-no-hit.XXXXXX")
cleanup() {
    status=$?
    if [[ ${status} -eq 0 ]]; then
        rm -rf "${test_dir}"
        return
    fi
    echo "No-hit integration failed; temporary evidence: ${test_dir}" >&2
    for log in \
        stdout.txt \
        stderr.txt \
        nextflow.log \
        unsupported.stdout.txt \
        unsupported.stderr.txt \
        unsupported.nextflow.log; do
        if [[ -s "${test_dir}/${log}" ]]; then
            echo "--- ${log} ---" >&2
            sed -n '1,400p' "${test_dir}/${log}" >&2
        fi
    done
}
trap cleanup EXIT

mkdir -p "${test_dir}/input"
query_fasta="${test_dir}/input/single_sample.fna"
printf '>no_hit description retained\nACGTACGTACGT\n' > "${query_fasta}"
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

(
    cd "${test_dir}"
    nextflow \
        -log "${test_dir}/nextflow.log" \
        run "${repo_dir}/main.nf" \
        -ansi-log false \
        -work-dir "${test_dir}/work" \
        --querydir "${query_fasta}" \
        --modeldir "${repo_dir}/resources/models" \
        --database_path "${test_dir}/database" \
        --min_extract_length 0 \
        --threads_per_job 1 \
        > "${test_dir}/stdout.txt" \
        2> "${test_dir}/stderr.txt"
)

outdir="${test_dir}/results/single_sample"

test -d "${outdir}"
test ! -e "${test_dir}/results/single_sample.fna"
test "$(wc -l < "${outdir}/cmsearch_summary.tsv")" -eq 1
test "$(wc -l < "${outdir}/blast_top_hits.tsv")" -eq 1
test "$(wc -l < "${outdir}/tree_nearest_neighbors.tsv")" -eq 1
test "$(wc -l < "${outdir}/cmsearch_summary.tab")" -eq 2
test "$(sed -n '2p' "${outdir}/cmsearch_summary.tab")" = "single_sample"
test "$(find "${outdir}/stats" -name '*.hits.tsv' -type f | wc -l)" -eq 2
test "$(find "${outdir}/out" "${outdir}/stats" -name '*.fna' -type f | wc -l)" -eq 0
test "$(find "${outdir}/m8" -name '*.top_hits.tsv' -type f | wc -l)" -eq 2

unsupported="${test_dir}/input/single_sample.txt"
cp "${query_fasta}" "${unsupported}"
if (
    cd "${test_dir}"
    nextflow \
        -log "${test_dir}/unsupported.nextflow.log" \
        run "${repo_dir}/main.nf" \
        -ansi-log false \
        -work-dir "${test_dir}/unsupported-work" \
        --querydir "${unsupported}" \
        --modeldir "${repo_dir}/resources/models" \
        --database_path "${test_dir}/database" \
        --outdir "${test_dir}/unsupported-results" \
        > "${test_dir}/unsupported.stdout.txt" \
        2> "${test_dir}/unsupported.stderr.txt"
); then
    echo "Unsupported single-file suffix unexpectedly started the pipeline" >&2
    exit 1
fi
grep -Fq -- \
    "--querydir file must end in .fna, .fa, or .fasta" \
    "${test_dir}/unsupported.stdout.txt" \
    "${test_dir}/unsupported.stderr.txt"
if grep -Fq "Submitted process" \
    "${test_dir}/unsupported.stdout.txt" \
    "${test_dir}/unsupported.stderr.txt"; then
    echo "Unsupported single-file suffix submitted a process" >&2
    exit 1
fi
