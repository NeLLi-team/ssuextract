# SSUextract project structure

```text
ssuextract/
├── main.nf                       # Four-stage Nextflow workflow
├── nextflow.config               # Parameters, profiles, reports, and manifest
├── pixi.toml                     # Environment constraints and commands
├── pixi.lock                     # Resolved cross-platform environment
├── config/
│   ├── base.config               # Per-process resource requests
│   ├── database_catalog.json     # Installable profile archives
│   ├── database_sources.json     # Pinned source versions, URLs, and hashes
│   ├── model_markers.json        # Covariance-model to 16S/18S mapping
│   ├── environment.yml           # Optional Nextflow Conda environment
│   └── local.config              # Local database path; generated and ignored
├── data/example/                 # Bundled example assemblies
├── resources/
│   ├── models/                   # RF00177 and RF01960 covariance models
│   └── database_notices/         # Source licenses and attribution
├── scripts/
│   ├── pipeline_cli.sh           # Database setup and Pixi-facing runner
│   ├── database_manager.py       # Profile download, validation, and resolution
│   ├── build_database_profiles.py # Reproducible source-to-release driver
│   ├── build_database_release.py # Source parsing and database construction
│   ├── database_sources.py       # Source-specific parsing
│   ├── database_contracts.py     # Taxonomy and record contracts
│   ├── database_release_io.py    # BLAST, Parquet, evidence, and manifest output
│   ├── assemble_database_profile.py # Validated profile/archive publication
│   ├── calibrate_taxonomy.py     # Leave-one-reference-out rank calibration
│   ├── classify_img_clusters.py  # Conservative IMG cluster classification
│   ├── classify_img_marker.sh    # Production centroid-classification contract
│   ├── hit_processing.py         # Typed hit parsing and sequence extraction
│   ├── extract_hits.py           # Extraction command-line interface
│   ├── annotate_hits.py          # BLAST annotation and taxonomy resolution
│   ├── finalize_summaries.py     # Deterministic final reports
│   ├── get_cmsequences.py        # Compatibility wrapper for legacy callers
│   └── check_version.py          # Release-version consistency gate
└── tests/                         # Unit and self-contained integration tests
```

`main.nf` carries sample and model identifiers as explicit channel metadata. Raw
search output is parsed once into hit and extraction records. Those records feed
the FASTA, detailed summary, category summary, and annotation outputs, avoiding
independent reconstruction of biological state in separate scripts.

Database construction is separate from runtime installation. The source catalog
pins every input; the build modules normalize sequences and taxonomy, write
marker-specific BLAST databases and compressed metadata, then assemble one
checksum-covered archive. Runtime setup accepts only a catalog entry whose
archive and internal manifest both validate.

Run all automated checks with:

```bash
pixi run test
```

The integration tests create their own FASTA inputs and BLAST databases, use
absolute paths, exercise no-hit and marker-routing behavior, and remove their
temporary files afterward.
