#!/usr/bin/env bash
set -euo pipefail

repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/ssuextract-tree.XXXXXX")
cleanup() {
    status=$?
    if [[ ${status} -eq 0 ]]; then
        rm -rf "${test_dir}"
        return
    fi
    echo "Tree integration failed; temporary evidence: ${test_dir}" >&2
    for log in \
        first.stdout first.stderr \
        resume.stdout resume.stderr \
        no-hit.stdout no-hit.stderr; do
        if [[ -s "${test_dir}/${log}" ]]; then
            echo "--- ${log} ---" >&2
            sed -n '1,400p' "${test_dir}/${log}" >&2
        fi
    done
}
trap cleanup EXIT

mkdir -p \
    "${test_dir}/query" \
    "${test_dir}/database/curated/blast" \
    "${test_dir}/database/curated/metadata"

python3 - \
    "${repo}/data/example/LKH565_P11_Ci.fna" \
    "${repo}/data/example/LKH462_P08_Rh.fna" \
    "${test_dir}/combined.fna" \
    "${test_dir}/query" <<'PY'
import sys
from pathlib import Path

from Bio import SeqIO


targets = {
    "LKH565_P11_Ci|NODE_19_length_16319_cov_4.794085": "bacterial_locus",
    "LKH462_P08_Rh|NODE_382_length_10603_cov_3.207614": "eukaryotic_locus",
}
records = [
    record
    for source in map(Path, sys.argv[1:3])
    for record in SeqIO.parse(source, "fasta")
    if record.id in targets
]
if {record.id for record in records} != set(targets):
    raise SystemExit("could not construct tree-classification loci")
SeqIO.write(records, sys.argv[3], "fasta")
for record in records:
    SeqIO.write([record], Path(sys.argv[4]) / f"{targets[record.id]}.fna", "fasta")
PY

for marker in 16S 18S; do
    if [[ "${marker}" == "16S" ]]; then
        model=RF00177
    else
        model=RF01960
    fi
    cmsearch \
        --anytrunc \
        --cpu 1 \
        -o /dev/null \
        --tblout "${test_dir}/${model}.out" \
        "${repo}/resources/models/${model}.cm" \
        "${test_dir}/combined.fna"
done

python3 "${repo}/scripts/resolve_model_hits.py" \
    --cmsearch "${test_dir}/RF00177.out" \
    --cmsearch "${test_dir}/RF01960.out" \
    --output "${test_dir}/accepted-hits.tsv"

for marker in 16S 18S; do
    if [[ "${marker}" == "16S" ]]; then
        model=RF00177
    else
        model=RF01960
    fi
    python3 "${repo}/scripts/extract_hits.py" \
        --model-file "${repo}/resources/models/${model}.cm" \
        --fasta "${test_dir}/combined.fna" \
        --sample tree-test \
        --model "${model}" \
        --accepted-hits "${test_dir}/accepted-hits.tsv" \
        --minimum-length 500 \
        --fasta-output "${test_dir}/${marker}.query.fna" \
        --hits-output "${test_dir}/${marker}.hits.tsv" \
        --metadata-output "${test_dir}/${marker}.meta.tsv"
done

python3 - \
    "${test_dir}/16S.query.fna" \
    "${test_dir}/18S.query.fna" \
    "${test_dir}/database/curated" <<'PY'
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import duckdb
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


profile = Path(sys.argv[3])
lineages = {
    "16S": "Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonas",
    "18S": "Eukaryota;Amoebozoa;Discosea;Echinamoebida;Echinamoeba",
}
sources = {"16S": ("SILVA", "138.2"), "18S": ("PR2", "5.1.1")}
taxonomy_rows = []
source_rows = []


def mutate(sequence: str, offset: int, count: int) -> str:
    values = list(sequence.upper().replace("U", "T"))
    replacements = {"A": "C", "C": "G", "G": "T", "T": "A"}
    candidates = [
        index
        for index in range(offset, len(values))
        if values[index] in replacements
    ]
    if len(candidates) < count:
        raise RuntimeError("reference sequence lacks enough mutable nucleotides")
    step = max(1, len(candidates) // count)
    selected = candidates[::step][:count]
    for index in selected:
        values[index] = replacements[values[index]]
    return "".join(values)


for marker, query_path in zip(("16S", "18S"), map(Path, sys.argv[1:3])):
    query_records = list(SeqIO.parse(query_path, "fasta"))
    if len(query_records) != 1:
        raise SystemExit(f"expected one extracted {marker} query")
    sequence = str(query_records[0].seq)
    records = []
    source, version = sources[marker]
    for rank, offset in enumerate((0, 100, 300), start=1):
        sequence_id = f"ref{marker}_{rank}"
        reference_sequence = (
            sequence if rank == 1 else mutate(sequence, offset, 12 * (rank - 1))
        )
        records.append(SeqRecord(Seq(reference_sequence), id=sequence_id, description=""))
        taxonomy_rows.append(
            (
                sequence_id,
                source,
                lineages[marker],
                source,
                lineages[marker].split(";", 1)[0],
                "nucleus" if marker == "18S" else "",
                "native",
                False,
                "",
                "",
                "",
                "",
            )
        )
        source_rows.append(
            (
                f"SRC_{sequence_id}",
                sequence_id,
                source,
                version,
                f"{source}_{marker}_{rank}",
                f"{source}_{marker}_{rank}",
                marker,
                None,
            )
        )
    fasta = profile / "blast" / f"{marker}.fna"
    SeqIO.write(records, fasta, "fasta")
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            str(fasta),
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-blastdb_version",
            "5",
            "-out",
            str(profile / "blast" / marker),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    fasta.unlink()

taxonomy = profile / "metadata" / "preferred_taxonomy.parquet"
sources_path = profile / "metadata" / "source_records.parquet"
connection = duckdb.connect(":memory:")
connection.execute(
    """
    CREATE TABLE taxonomy (
        sequence_id VARCHAR, reference_source VARCHAR, taxonomy VARCHAR,
        taxonomy_source VARCHAR, domain VARCHAR, compartment VARCHAR,
        assignment_method VARCHAR, cross_domain_conflict BOOLEAN,
        taxonomy_alternatives VARCHAR, centroid_names VARCHAR,
        centroid_taxonomy VARCHAR, centroid_taxonomy_source VARCHAR
    )
    """
)
connection.executemany(
    "INSERT INTO taxonomy VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
    taxonomy_rows,
)
connection.execute("COPY taxonomy TO ? (FORMAT PARQUET)", [str(taxonomy)])
connection.execute(
    """
    CREATE TABLE source_records (
        source_record_id VARCHAR, sequence_id VARCHAR, reference_source VARCHAR,
        source_version VARCHAR, source_identifier VARCHAR, original_header VARCHAR,
        marker VARCHAR, taxon_oid VARCHAR
    )
    """
)
connection.executemany(
    "INSERT INTO source_records VALUES (?, ?, ?, ?, ?, ?, ?, ?)", source_rows
)
connection.execute("COPY source_records TO ? (FORMAT PARQUET)", [str(sources_path)])
connection.close()

artifacts = []
for path in sorted((profile / "blast").iterdir()) + [taxonomy, sources_path]:
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
    "version": "tree-test",
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

run_tree_pipeline() {
    nextflow run "${repo}/main.nf" \
        --querydir "${test_dir}/query" \
        --modeldir "${repo}/resources/models" \
        --database_path "${test_dir}/database" \
        --database_profile curated \
        --outdir "${test_dir}/results" \
        --threads_per_job 1 \
        --tree_classification \
        --tree_reference_count 3 \
        --tree_assignment_neighbors 2 \
        -ansi-log false \
        -work-dir "${test_dir}/work" \
        "$@"
}

run_tree_pipeline >"${test_dir}/first.stdout" 2>"${test_dir}/first.stderr"
run_tree_pipeline -resume >"${test_dir}/resume.stdout" 2>"${test_dir}/resume.stderr"

python3 - \
    "${test_dir}/results/cmsearch_summary.tsv" \
    "${test_dir}/results/tree_nearest_neighbors.tsv" \
    "${test_dir}/results/phylogeny" \
    "${test_dir}/results/cmsearch_summary.tab" <<'PY'
import csv
import json
import sys
from collections import Counter
from pathlib import Path


with Path(sys.argv[1]).open(newline="") as handle:
    summaries = list(csv.DictReader(handle, delimiter="\t"))
with Path(sys.argv[2]).open(newline="") as handle:
    neighbors = list(csv.DictReader(handle, delimiter="\t"))
with Path(sys.argv[4]).open(newline="") as handle:
    categories = list(csv.DictReader(handle, delimiter="\t"))
if len(summaries) != 2:
    raise SystemExit(f"expected two tree-classified loci, found {len(summaries)}")
if {row["taxonomy_mode"] for row in summaries} != {"tree"}:
    raise SystemExit("tree mode did not become the selected taxonomy mode")
expected = {
    "RF00177": "Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonas",
    "RF01960": "Eukaryota;Amoebozoa;Discosea;Echinamoebida;Echinamoeba",
}
if {row["model"]: row["taxonomy"] for row in summaries} != expected:
    raise SystemExit("tree taxonomy does not match the synthetic named neighborhood")
if any(not row["blast_taxonomy"] or not row["query_sequence"] for row in summaries):
    raise SystemExit("tree summary did not retain BLAST taxonomy or query sequence")
if len(categories) != 2:
    raise SystemExit(f"expected two sample category rows, found {len(categories)}")
if sum(int(row.get("BacteriaSSU", 0)) for row in categories) != 1:
    raise SystemExit("tree taxonomy did not produce one bacterial category count")
if sum(int(row.get("EukaryotaSSU", 0)) for row in categories) != 1:
    raise SystemExit("tree taxonomy did not produce one eukaryotic category count")
if len(neighbors) != 6:
    raise SystemExit(f"expected six tree-neighbor rows, found {len(neighbors)}")
if set(Counter(row["name"] for row in neighbors).values()) != {3}:
    raise SystemExit("each query must retain all three synthetic tree neighbors")
if set(
    Counter(
        row["name"] for row in neighbors if row["used_for_assignment"] == "true"
    ).values()
) != {2}:
    raise SystemExit("each tree assignment must use two named neighbors")
for query_rows in (
    [row for row in neighbors if row["name"] == name]
    for name in {row["name"] for row in neighbors}
):
    distances = [float(row["tree_distance"]) for row in query_rows]
    if distances != sorted(distances):
        raise SystemExit("tree neighbors are not sorted by patristic distance")

tree_directories = [path for path in Path(sys.argv[3]).glob("*/*/*") if path.is_dir()]
if len(tree_directories) != 2:
    raise SystemExit(f"expected two per-query tree directories, found {len(tree_directories)}")
required = {
    "alignment_qc.json",
    "cmalign.trimmed.fna",
    "iqtree.treefile",
    "references.fna",
    "references.tsv",
    "task.json",
    "tool_versions.txt",
}
for directory in tree_directories:
    names = {path.name for path in directory.iterdir() if path.stat().st_size > 0}
    missing = required - names
    if missing:
        raise SystemExit(f"tree directory {directory} lacks {sorted(missing)}")
    qc = json.loads((directory / "alignment_qc.json").read_text())
    if qc["sequence_count"] != 4 or qc["retained_columns"] <= 0:
        raise SystemExit(f"invalid tree alignment QC in {directory}")
    if qc["removed_columns"] != (
        qc["removed_insert_columns"] + qc["removed_high_gap_columns"]
    ):
        raise SystemExit(f"inconsistent alignment masking QC in {directory}")
PY

mkdir "${test_dir}/no-hit-query"
printf '>no_tree_hit\nACGTACGTACGT\n' > "${test_dir}/no-hit-query/no_tree_hit.fna"
nextflow run "${repo}/main.nf" \
    --querydir "${test_dir}/no-hit-query" \
    --modeldir "${repo}/resources/models" \
    --database_path "${test_dir}/database" \
    --database_profile curated \
    --outdir "${test_dir}/no-hit-results" \
    --threads_per_job 1 \
    --tree_classification \
    --tree_reference_count 3 \
    --tree_assignment_neighbors 2 \
    -ansi-log false \
    -work-dir "${test_dir}/no-hit-work" \
    >"${test_dir}/no-hit.stdout" \
    2>"${test_dir}/no-hit.stderr"

test "$(wc -l < "${test_dir}/no-hit-results/cmsearch_summary.tsv")" -eq 1
test "$(wc -l < "${test_dir}/no-hit-results/blast_top_hits.tsv")" -eq 1
test "$(wc -l < "${test_dir}/no-hit-results/tree_nearest_neighbors.tsv")" -eq 1
test "$(wc -l < "${test_dir}/no-hit-results/cmsearch_summary.tab")" -eq 2
