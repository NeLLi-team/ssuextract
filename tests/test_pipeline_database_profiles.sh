#!/usr/bin/env bash
set -euo pipefail

repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/ssuextract-db-profiles.XXXXXX")
trap 'rm -rf "${test_dir}"' EXIT

mkdir -p \
    "${test_dir}/query" \
    "${test_dir}/database/curated/blast" \
    "${test_dir}/database/curated/metadata"

integration_fasta="${test_dir}/model-competition.fna"
python3 - \
    "${repo}/data/example/LKH565_P11_Ci.fna" \
    "${repo}/data/example/LKH462_P08_Rh.fna" \
    "${integration_fasta}" \
    "${test_dir}/query" <<'PY'
import sys
from pathlib import Path

from Bio import SeqIO


targets = {
    "LKH565_P11_Ci|NODE_19_length_16319_cov_4.794085",
    "LKH462_P08_Rh|NODE_382_length_10603_cov_3.207614",
}
records = [
    record
    for source in map(Path, sys.argv[1:3])
    for record in SeqIO.parse(source, "fasta")
    if record.id in targets
]
if {record.id for record in records} != targets:
    raise SystemExit("could not construct the two-locus model competition fixture")
SeqIO.write(records, sys.argv[3], "fasta")
query_dir = Path(sys.argv[4])
for record in records:
    sample = (
        "bacterial_locus"
        if "NODE_19_" in record.id
        else "eukaryotic_locus"
    )
    SeqIO.write([record], query_dir / f"{sample}.fna", "fasta")
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
        "${integration_fasta}"
done

python3 "${repo}/scripts/resolve_model_hits.py" \
    --cmsearch "${test_dir}/RF00177.out" \
    --cmsearch "${test_dir}/RF01960.out" \
    --output "${test_dir}/accepted-hits.tsv"

PYTHONPATH="${repo}/scripts" python3 - \
    "${test_dir}/RF00177.out" \
    "${test_dir}/RF01960.out" \
    "${test_dir}/accepted-hits.tsv" <<'PY'
import sys

from hit_processing import parse_cmsearch_tblout, read_accepted_hits


raw = [
    hit
    for path in sys.argv[1:3]
    for hit in parse_cmsearch_tblout(path)
    if hit.included
]
accessions = {hit.model_accession for hit in raw}
if accessions != {"RF00177", "RF01960"}:
    raise SystemExit(f"bundled models emitted wrong accessions: {sorted(accessions)}")

by_subject = {}
for hit in raw:
    by_subject.setdefault(hit.subject, set()).add(hit.model_accession)
if len(by_subject) != 2 or any(
    models != {"RF00177", "RF01960"}
    for models in by_subject.values()
):
    raise SystemExit(f"fixture did not produce cross-model competition: {by_subject}")

accepted = read_accepted_hits(sys.argv[3])
observed = {(hit.subject, hit.model_accession) for hit in accepted}
expected = {
    ("LKH565_P11_Ci|NODE_19_length_16319_cov_4.794085", "RF00177"),
    ("LKH462_P08_Rh|NODE_382_length_10603_cov_3.207614", "RF01960"),
}
if observed != expected:
    raise SystemExit(f"wrong combined model resolution: {sorted(observed)}")
PY

for marker in 16S 18S; do
    subject=ref${marker%S}
    if [[ "${marker}" == "16S" ]]; then
        model=RF00177
    else
        model=RF01960
    fi
    python3 "${repo}/scripts/extract_hits.py" \
        --model-file "${repo}/resources/models/${model}.cm" \
        --fasta "${integration_fasta}" \
        --sample database-test \
        --model "${model}" \
        --accepted-hits "${test_dir}/accepted-hits.tsv" \
        --minimum-length 500 \
        --fasta-output "${test_dir}/${model}.fna" \
        --hits-output "${test_dir}/${model}.hits.tsv" \
        --metadata-output "${test_dir}/${model}.meta.tsv"
    [[ "$(grep -c '^>' "${test_dir}/${model}.fna")" -eq 1 ]]
    awk -v subject="${subject}" '
        /^>/ { print ">" subject; next }
        { print }
    ' "${test_dir}/${model}.fna" \
        > "${test_dir}/database/curated/blast/${marker}.fna"
    makeblastdb \
        -in "${test_dir}/database/curated/blast/${marker}.fna" \
        -dbtype nucl \
        -blastdb_version 5 \
        -out "${test_dir}/database/curated/blast/${marker}" \
        >/dev/null
    rm "${test_dir}/database/curated/blast/${marker}.fna"
done

mkdir -p "${test_dir}/legacy-database"
{
    awk '/^>/ { print ">ref16"; next } { print }' "${test_dir}/RF00177.fna"
    awk '/^>/ { print ">ref18"; next } { print }' "${test_dir}/RF01960.fna"
} > "${test_dir}/legacy.fna"
makeblastdb \
    -in "${test_dir}/legacy.fna" \
    -dbtype nucl \
    -blastdb_version 4 \
    -out "${test_dir}/legacy-database/silva-138-1_pr2-4-12" \
    >/dev/null

if nextflow run "${repo}/main.nf" \
    --query "${test_dir}/query" \
    --modeldir "${repo}/resources/models" \
    --database_path "${test_dir}/legacy-database" \
    --database_profile img \
    --outdir "${test_dir}/legacy-fallback-results" \
    --threads_per_job 1 \
    -work-dir "${test_dir}/legacy-fallback-work" \
    >"${test_dir}/legacy-fallback.log" 2>&1; then
    echo "IMG profile unexpectedly accepted the deprecated curated-only database" >&2
    exit 1
fi
grep -F "Database profile 'img' is not installed" "${test_dir}/legacy-fallback.log"

if nextflow run "${repo}/main.nf" \
    --query "${test_dir}/query" \
    --modeldir "${repo}/resources/models" \
    --database_path "${test_dir}/legacy-database" \
    --database_profile curated \
    --outdir "${test_dir}/legacy-tree-results" \
    --threads_per_job 1 \
    --tree_classification \
    -work-dir "${test_dir}/legacy-tree-work" \
    >"${test_dir}/legacy-tree.log" 2>&1; then
    echo "Tree mode unexpectedly accepted the deprecated database" >&2
    exit 1
fi
grep -F -- \
    "--tree_classification requires a managed curated or img database profile" \
    "${test_dir}/legacy-tree.log"

python3 - "${test_dir}/database/curated" <<'PY'
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
        taxonomy_alternatives VARCHAR
    )
    """
)
connection.executemany(
    "INSERT INTO taxonomy VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
    [
        ("ref16", "SILVA", "Bacteria", "SILVA", "Bacteria", "", "native", False, ""),
        ("ref18", "PR2", "Eukaryota", "PR2", "Eukaryota", "nucleus", "native", False, ""),
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
for path in sorted((profile / "blast").glob("16S.*")) + sorted(
    (profile / "blast").glob("18S.*")
) + [taxonomy, source_records]:
    relative = path.relative_to(profile).as_posix()
    artifacts.append(
        {
            "path": relative,
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

cp -a "${test_dir}/database" "${test_dir}/shell-safe-database"
python3 - "${test_dir}/shell-safe-database/curated" <<'PY'
import json
import subprocess
import sys
from pathlib import Path


profile = Path(sys.argv[1])
blast = profile / "blast"
hostile_prefix = "16S'$(:>PWNED_DB)"
source_fasta = blast / "16S-hostile-source.fna"
subprocess.run(
    [
        "blastdbcmd",
        "-db",
        str(blast / "16S"),
        "-entry",
        "all",
        "-outfmt",
        "%f",
        "-out",
        str(source_fasta),
    ],
    check=True,
)
for path in sorted(blast.glob("16S.*")):
    if path != source_fasta:
        path.unlink()
subprocess.run(
    [
        "makeblastdb",
        "-in",
        str(source_fasta),
        "-dbtype",
        "nucl",
        "-blastdb_version",
        "5",
        "-out",
        str(blast / hostile_prefix),
    ],
    check=True,
    stdout=subprocess.DEVNULL,
)
source_fasta.unlink()

taxonomy = profile / "metadata" / "preferred_taxonomy.parquet"
hostile_taxonomy = taxonomy.with_name("preferred`touch PWNED_TAX`.parquet")
taxonomy.rename(hostile_taxonomy)
source_records = profile / "metadata" / "source_records.parquet"
hostile_source_records = source_records.with_name(
    "source`touch PWNED_SOURCE`.parquet"
)
source_records.rename(hostile_source_records)

manifest_path = profile / "manifest.json"
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
manifest["blast_databases"]["16S"]["prefix"] = f"blast/{hostile_prefix}"
manifest["taxonomy_database"]["preferred"] = (
    f"metadata/{hostile_taxonomy.name}"
)
manifest["taxonomy_database"]["source_records"] = (
    f"metadata/{hostile_source_records.name}"
)
manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
PY

if ! nextflow run "${repo}/main.nf" \
    --query "${test_dir}/query" \
    --modeldir "${repo}/resources/models" \
    --database_path "${test_dir}/shell-safe-database" \
    --database_profile curated \
    --outdir "${test_dir}/shell-safe-results" \
    --threads_per_job 1 \
    -work-dir "${test_dir}/shell-safe-work" \
    >"${test_dir}/shell-safe.log" 2>&1; then
    cat "${test_dir}/shell-safe.log" >&2
    exit 1
fi

if find "${test_dir}/shell-safe-work" -type f \
    \( -name PWNED_DB -o -name PWNED_TAX -o -name PWNED_SOURCE \) \
    -print -quit | grep -q .; then
    echo "Manifest-derived path executed shell syntax" >&2
    exit 1
fi

cp -a "${test_dir}/database" "${test_dir}/traversal-database"
cp \
    "${test_dir}/database/curated/metadata/preferred_taxonomy.parquet" \
    "${test_dir}/traversal-database/outside.parquet"
python3 - "${test_dir}/traversal-database/curated/manifest.json" <<'PY'
import json
import sys
from pathlib import Path


path = Path(sys.argv[1])
manifest = json.loads(path.read_text(encoding="utf-8"))
manifest["taxonomy_database"]["preferred"] = "../outside.parquet"
path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
PY

if nextflow run "${repo}/main.nf" \
    --query "${test_dir}/query" \
    --modeldir "${repo}/resources/models" \
    --database_path "${test_dir}/traversal-database" \
    --database_profile curated \
    --outdir "${test_dir}/traversal-results" \
    --threads_per_job 1 \
    -work-dir "${test_dir}/traversal-work" \
    >"${test_dir}/traversal.log" 2>&1; then
    echo "Escaping taxonomy path unexpectedly passed profile validation" >&2
    exit 1
fi
grep -F "escapes database profile" "${test_dir}/traversal.log"

cp -a "${test_dir}/database" "${test_dir}/missing-source-database"
python3 - "${test_dir}/missing-source-database/curated/manifest.json" <<'PY'
import json
import sys
from pathlib import Path


path = Path(sys.argv[1])
manifest = json.loads(path.read_text(encoding="utf-8"))
del manifest["taxonomy_database"]["source_records"]
path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
PY

if nextflow run "${repo}/main.nf" \
    --query "${test_dir}/query" \
    --modeldir "${repo}/resources/models" \
    --database_path "${test_dir}/missing-source-database" \
    --database_profile curated \
    --outdir "${test_dir}/missing-source-results" \
    --threads_per_job 1 \
    -work-dir "${test_dir}/missing-source-work" \
    >"${test_dir}/missing-source.log" 2>&1; then
    echo "Manifest without source_records unexpectedly passed validation" >&2
    exit 1
fi
grep -F "missing taxonomy_database.source_records" "${test_dir}/missing-source.log"

nextflow run "${repo}/main.nf" \
    --query "${test_dir}/query" \
    --modeldir "${repo}/resources/models" \
    --database_path "${test_dir}/database" \
    --database_profile curated \
    --outdir "${test_dir}/results" \
    --threads_per_job 1 \
    -work-dir "${test_dir}/work" \
    >/dev/null

python3 - \
    "${test_dir}/results/cmsearch_summary.tsv" \
    "${test_dir}/results/blast_top_hits.tsv" <<'PY'
import csv
import sys
from pathlib import Path


with Path(sys.argv[1]).open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
with Path(sys.argv[2]).open(newline="") as handle:
    top_hits = list(csv.DictReader(handle, delimiter="\t"))
observed = {
    (
        row["sample"],
        row["model"],
        row["blast_sseqid"],
        row["taxonomy_domain"],
        row["taxonomy_source"],
    )
    for row in rows
}
expected = {
    ("bacterial_locus", "RF00177", "ref16", "Bacteria", "SILVA"),
    ("eukaryotic_locus", "RF01960", "ref18", "Eukaryota", "PR2"),
}
if len(rows) != 2 or observed != expected:
    raise SystemExit(
        f"wrong resolved model/database routes: expected {expected}, found {observed}"
    )
if any(not row["query_sequence"] for row in rows):
    raise SystemExit("summary is missing extracted query sequences")
if {row["reference_identifiers"] for row in rows} != {
    "SILVA:AB16.1.1537",
    "PR2:AB18.1.1795_U",
}:
    raise SystemExit("summary is missing public reference identifiers")
if len(top_hits) != 2 or any(not row["taxonomy"] for row in top_hits):
    raise SystemExit("metadata-joined top-hit table is incomplete")
if any(not row["query_sequence"] for row in top_hits):
    raise SystemExit("top-hit table is missing extracted query sequences")
PY
