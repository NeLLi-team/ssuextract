# SSUextract

SSUextract is a Nextflow pipeline for extracting SSU rRNA sequences from contigs with covariance models and annotating the hits with BLAST.

## What works today

- `pixi` manages the runtime environment.
- `nextflow` runs the pipeline itself.
- The default annotation database is SILVA/PR2.
- `pixi run ssuextract` now checks for the database and downloads it on first use.

## Quick Start

```bash
git clone https://github.com/NeLLi-team/ssuextract.git
cd ssuextract

# First-time setup. Prompts for the database location if no local config exists.
pixi run setup

# One-file smoke test on the smaller bundled example.
pixi run example

# Full bundled example.
pixi run ssuextract
```

## First-Run Database Setup

The runtime expects a BLAST database with the `silva-138-1_pr2-4-12` prefix.

On the first `pixi run setup`, `pixi run ssuextract`, or `pixi run example`:

1. SSUextract checks `config/local.config` for `database_path`.
2. If no local database path is configured, it prompts for one.
3. Pressing Enter accepts the default: `./resources/database`.
4. Missing database files are downloaded automatically.

The selected path is saved to `config/local.config`, which is gitignored.

If you want to bypass the saved path for a single run, pass `--database_path` directly:

```bash
pixi run ssuextract -- --database_path /path/to/database
```

## Commands

```bash
# Resolve database path and download the default database if needed
pixi run setup
pixi run download-db

# Smoke test on the smaller bundled file
pixi run example

# Run the full bundled example in data/example
pixi run ssuextract

# Preview the Nextflow graph without executing jobs
pixi run dryrun

# Run with custom arguments
pixi run ssuextract -- --querydir data/my_dataset --threads_per_job 4

# Development helpers
pixi run -e dev run-verbose
pixi run -e dev report
pixi run -e dev dag
```

## One-File Smoke Test

The pipeline consumes a directory of FASTA files, not a single FASTA path.

`pixi run example` stages the smaller bundled file `data/example/LKH565_P11_Ci.fna` into a temporary directory and runs:

```bash
nextflow run main.nf \
    --querydir <temporary_dir> \
    --outdir results/smoke \
    --threads_per_job 1
```

If you want to do the same manually:

```bash
tmpdir=$(mktemp -d /tmp/ssuextract-smoke.XXXXXX)
cp data/example/LKH565_P11_Ci.fna "$tmpdir/"
pixi run nextflow run main.nf --querydir "$tmpdir" --outdir results/smoke_test --threads_per_job 1
```

## Custom Runs

Run on your own assemblies by putting `.fna`, `.fa`, or `.fasta` files in a directory and pointing `--querydir` at it.

```bash
mkdir -p data/my_dataset
cp /path/to/*.fna data/my_dataset/
pixi run ssuextract -- --querydir data/my_dataset --threads_per_job 4
```

You can also call Nextflow directly if you prefer:

```bash
pixi run nextflow run main.nf --querydir data/my_dataset --threads_per_job 4
```

If you configured a non-default database location, either keep using `pixi run ssuextract` or pass `--database_path` explicitly when calling Nextflow directly.

## Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--querydir` | `data/example` | Directory containing `.fna`, `.fa`, or `.fasta` files |
| `--modeldir` | `resources/models` | Directory containing covariance models |
| `--outdir` | `results/<dataset>` | Output directory |
| `--database_path` | `resources/database` | Directory containing the BLAST database files |
| `--min_extract_length` | `500` | Minimum extracted sequence length |
| `--threads_per_job` | `2` | Threads per process |

## Project Layout

```text
ssuextract/
├── main.nf
├── nextflow.config
├── pixi.toml
├── config/
│   ├── base.config
│   ├── environment.yml
│   └── local.config      # created on first setup, gitignored
├── data/
│   └── example/
├── resources/
│   ├── models/
│   └── database/         # created on setup, gitignored
├── results/              # generated outputs, gitignored
└── scripts/
    ├── pipeline_cli.sh   # pixi-facing setup/run wrapper
    ├── get_cmstats.py
    ├── get_cmsequences.py
    ├── cmprocessing.py
    └── get_table.py
```

## Outputs

Main outputs are written under the selected `outdir`.

Important files:

- `cmsearch_summary.tsv`: per-hit summary table with coordinates, lengths, and top BLAST hit
- `cmsearch_summary.tab`: per-sample category counts
- `pipeline_info/`: Nextflow reports, trace, and DAG

## Notes

- The current documented and supported database path is SILVA/PR2.
- `pixi run ssuextract` and `pixi run example` both auto-bootstrap the database if it is missing.
- For non-interactive runs with no saved config, SSUextract falls back to `./resources/database`.
