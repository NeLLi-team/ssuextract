#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
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
        unsupported.nextflow.log \
        duplicate.stdout.txt \
        duplicate.stderr.txt \
        duplicate.nextflow.log; do
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
profile_dir="${test_dir}/database/curated"
mkdir -p "${profile_dir}/blast" "${profile_dir}/metadata"
printf '>ref16\nACGTACGTACGTACGT\n' > "${test_dir}/database/16S.fna"
printf '>ref18\nACGTACGTACGTACGT\n' > "${test_dir}/database/18S.fna"
for marker in 16S 18S; do
    makeblastdb \
        -in "${test_dir}/database/${marker}.fna" \
        -dbtype nucl \
        -blastdb_version 5 \
        -out "${profile_dir}/blast/${marker}" \
        > "${test_dir}/makeblastdb.${marker}.stdout" \
        2> "${test_dir}/makeblastdb.${marker}.stderr"
done
python3 - "${profile_dir}" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

import duckdb


profile = Path(sys.argv[1])
taxonomy = profile / "metadata" / "preferred_taxonomy.parquet"
source_records = profile / "metadata" / "source_records.parquet"
connection = duckdb.connect(":memory:")
connection.execute(
    """
    CREATE TABLE taxonomy (
        sequence_id VARCHAR,
        reference_source VARCHAR,
        taxonomy VARCHAR,
        taxonomy_source VARCHAR,
        domain VARCHAR,
        compartment VARCHAR,
        assignment_method VARCHAR,
        cross_domain_conflict BOOLEAN,
        taxonomy_alternatives VARCHAR,
        centroid_names VARCHAR,
        centroid_taxonomy VARCHAR,
        centroid_taxonomy_source VARCHAR
    )
    """
)
connection.executemany(
    "INSERT INTO taxonomy VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
    [
        ("ref16", "SILVA", "Bacteria", "SILVA", "Bacteria", "", "native", False, "", "", "", ""),
        ("ref18", "PR2", "Eukaryota", "PR2", "Eukaryota", "nucleus", "native", False, "", "", "", ""),
    ],
)
connection.execute("COPY taxonomy TO ? (FORMAT PARQUET)", [str(taxonomy)])
connection.execute(
    """
    CREATE TABLE source_records (
        source_record_id VARCHAR,
        sequence_id VARCHAR,
        reference_source VARCHAR,
        source_version VARCHAR,
        source_identifier VARCHAR,
        original_header VARCHAR,
        marker VARCHAR,
        taxon_oid VARCHAR
    )
    """
)
connection.executemany(
    "INSERT INTO source_records VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
    [
        ("SRC_ref16", "ref16", "SILVA", "138.2", "AB16.1.1537", "AB16.1.1537 Bacteria", "16S", None),
        ("SRC_ref18", "ref18", "PR2", "5.1.1", "AB18.1.1795_U", "AB18.1.1795_U|18S_rRNA", "18S", None),
    ],
)
connection.execute("COPY source_records TO ? (FORMAT PARQUET)", [str(source_records)])
connection.close()

artifacts = []
for path in (
    sorted((profile / "blast").glob("16S.*"))
    + sorted((profile / "blast").glob("18S.*"))
    + [taxonomy, source_records]
):
    artifacts.append(
        {
            "path": path.relative_to(profile).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    )
manifest = {
    "schema_version": 1,
    "profile": "curated",
    "version": "test-1",
    "artifacts": artifacts,
    "blast_databases": {
        "16S": {"prefix": "blast/16S"},
        "18S": {"prefix": "blast/18S"},
    },
    "taxonomy_database": {
        "preferred": "metadata/preferred_taxonomy.parquet",
        "source_records": "metadata/source_records.parquet",
    },
}
(profile / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
PY

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
        --query "${query_fasta}" \
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
test "$(find "${outdir}/pipeline_info" -name 'pipeline_dag_*.html' -type f | wc -l)" -eq 0
test "$(find "${outdir}/pipeline_info" -name 'execution_trace_*.txt' -type f | wc -l)" -eq 1
test "$(find "${outdir}/pipeline_info" -name 'execution_timeline_*.html' -type f | wc -l)" -eq 1
test "$(find "${outdir}/pipeline_info" -name 'execution_report_*.html' -type f | wc -l)" -eq 1

unsupported="${test_dir}/input/single_sample.txt"
cp "${query_fasta}" "${unsupported}"
if (
    cd "${test_dir}"
    nextflow \
        -log "${test_dir}/unsupported.nextflow.log" \
        run "${repo_dir}/main.nf" \
        -ansi-log false \
        -work-dir "${test_dir}/unsupported-work" \
        --query "${unsupported}" \
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
    "--query file must end in .fna, .fa, or .fasta" \
    "${test_dir}/unsupported.stdout.txt" \
    "${test_dir}/unsupported.stderr.txt"
if grep -Fq "Submitted process" \
    "${test_dir}/unsupported.stdout.txt" \
    "${test_dir}/unsupported.stderr.txt"; then
    echo "Unsupported single-file suffix submitted a process" >&2
    exit 1
fi

mkdir "${test_dir}/duplicate-input"
cp "${query_fasta}" "${test_dir}/duplicate-input/a.fna"
cp "${query_fasta}" "${test_dir}/duplicate-input/a.fa"
cp "${query_fasta}" "${test_dir}/duplicate-input/unique.fasta"
if (
    cd "${test_dir}"
    nextflow \
        -log "${test_dir}/duplicate.nextflow.log" \
        run "${repo_dir}/main.nf" \
        -ansi-log false \
        -work-dir "${test_dir}/duplicate-work" \
        --query "${test_dir}/duplicate-input" \
        --database_path "${test_dir}/database" \
        --outdir "${test_dir}/duplicate-results" \
        > "${test_dir}/duplicate.stdout.txt" \
        2> "${test_dir}/duplicate.stderr.txt"
); then
    echo "Duplicate sample identifiers unexpectedly started the pipeline" >&2
    exit 1
fi
for expected in "Duplicate sample identifiers" "a.fna" "a.fa"; do
    grep -Fq -- "${expected}" \
        "${test_dir}/duplicate.stdout.txt" \
        "${test_dir}/duplicate.stderr.txt"
done
if grep -Fq "Submitted process" \
    "${test_dir}/duplicate.stdout.txt" \
    "${test_dir}/duplicate.stderr.txt"; then
    echo "Duplicate sample identifiers submitted a process" >&2
    exit 1
fi
