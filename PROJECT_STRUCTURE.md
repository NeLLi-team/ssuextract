# SSUextract Project Structure

This repository is intentionally small. The main workflow lives in Nextflow, and `pixi` provides the entrypoints users actually run.

## Layout

```text
ssuextract/
├── main.nf                 # Nextflow workflow
├── nextflow.config         # Runtime defaults and profiles
├── pixi.toml               # Environment and user-facing tasks
├── config/
│   ├── base.config         # Process resources
│   ├── environment.yml     # Conda profile definition
│   └── local.config        # Optional local database path override, gitignored
├── data/
│   └── example/            # Bundled example assemblies
├── resources/
│   ├── models/             # RF00177.cm and RF01960.cm
│   └── database/           # Downloaded BLAST database, gitignored
├── results/                # Pipeline outputs, gitignored
└── scripts/
    ├── pipeline_cli.sh     # Setup/run/smoke wrapper used by pixi
    ├── get_cmstats.py      # Parse cmsearch hits
    ├── get_cmsequences.py  # Extract hit sequences
    ├── cmprocessing.py     # Summarize BLAST categories
    └── get_table.py        # Build final TSV output
```

## Rationale

- Nextflow stays focused on batch execution.
- Interactive setup happens in `scripts/pipeline_cli.sh`, not inside `main.nf`.
- The database location is persisted in `config/local.config` so users are not prompted every run.
- The bundled smoke test uses the smaller example file without forcing a permanent extra dataset into the repo.

## Expected Run Paths

```bash
pixi run setup
pixi run example
pixi run ssuextract
pixi run ssuextract -- --querydir data/my_dataset
```

## Current Constraints

- The supported annotation database is the SILVA/PR2 BLAST database with prefix `silva-138-1_pr2-4-12`.
- The pipeline reads directories of FASTA files rather than single FASTA file paths.
- Output and downloaded resources are intentionally gitignored.
