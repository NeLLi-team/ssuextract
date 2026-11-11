# SSUextract project structure

```text
ssuextract/
├── main.nf                       # Nextflow workflow with seven processes
├── nextflow.config               # Parameters, reports, and manifest
├── pixi.toml                     # Environment constraints and commands
├── pixi.lock                     # Resolved cross-platform environment
├── mkdocs.yml                    # Documentation site configuration
├── README.md, CHANGELOG.md, LICENSE
├── .github/workflows/            # Test and documentation workflows
├── config/
│   ├── base.config               # Per-process resource requests
│   ├── database_catalog.json     # Installable profile archives
│   ├── database_sources.json     # Pinned source versions, URLs, and hashes
│   ├── model_markers.json        # Covariance-model to 16S/18S mapping
│   ├── img_search_execution.json # Chunk and thread layout of the IMG search
│   ├── empty.tree_assignment.tsv # Header-only tree table for BLAST-mode runs
│   ├── empty.tree_neighbors.tsv  # Header-only neighbor table for BLAST-mode runs
│   └── local.config              # Local database path and profile; generated and ignored
├── data/example/                 # Bundled example assemblies and expected annotations
├── resources/
│   ├── models/                   # RF00177 and RF01960 covariance models
│   └── database_notices/         # Source licenses and attribution
├── scripts/
│   ├── pipeline_cli.sh           # Database setup and Pixi-facing runner
│   ├── database_manager.py       # Profile validation and installation
│   ├── database_download.py      # Resumable, checksum-verified archive download
│   ├── database_updates.py       # Zenodo release discovery and verification
│   ├── atomic_io.py              # Durable file replacement
│   ├── build_database_profiles.py # Reproducible source-to-release driver
│   ├── build_database_release.py # Deduplication and preferred-taxonomy selection
│   ├── database_sources.py       # Source download and source-specific parsing
│   ├── database_contracts.py     # Taxonomy and record contracts
│   ├── database_release_io.py    # Release validation, marker FASTA, and Parquet output
│   ├── assemble_database_profile.py # BLAST indexes, manifest, and archive publication
│   ├── calibrate_taxonomy.py     # Leave-one-reference-out rank calibration
│   ├── extract_img_cluster_centroids.py # Centroid FASTA for the IMG search
│   ├── img_chunked_search.py     # Restartable, chunked centroid BLAST search
│   ├── img_search_provenance.py  # Search completion contract
│   ├── classify_img_clusters.py  # Conservative IMG cluster classification
│   ├── img_classification_data.py # Classification records and input loading
│   ├── classify_img_marker.sh    # Production centroid-classification contract
│   ├── hit_processing.py         # Typed hit parsing and sequence extraction
│   ├── resolve_model_hits.py     # Cross-model overlap resolution
│   ├── extract_hits.py           # Extraction command-line interface
│   ├── annotate_hits.py          # BLAST annotation and taxonomy resolution
│   ├── runtime_taxonomy.py       # Shared query evidence and calibrated assignment
│   ├── taxonomy_utils.py         # Ranked taxonomy paths and lowest common ancestor
│   ├── top_hit_reporting.py      # Ranked reference evidence
│   ├── tree_reference_selection.py # Marker routing and tree inputs
│   ├── tree_phylogeny.py         # Alignment trimming and tree evidence
│   ├── tree_schema.py            # Shared tree output fields
│   ├── finalize_summaries.py     # Deterministic final reports
│   ├── validate_example_output.py # Expected-annotation check for the example
│   └── check_version.py          # Release-version consistency gate
├── tests/                        # Unit and self-contained integration tests
├── docs/                         # MkDocs pages, figures, and their data tables
├── docs-env/                     # Pixi environment for the documentation build
└── notebooks/                    # Figure scripts and notebooks behind the docs
```

`main.nf` carries sample and model identifiers in channel metadata. Hit and
extraction records supply the FASTA, detailed summary, category summary, and
annotation outputs. The seven
processes are `CMSEARCH`, `RESOLVE_MODEL_HITS`, `EXTRACT_HITS`, `BLAST_ANNOTATE`,
`PREPARE_TREE_TASKS`, `TREE_CLASSIFY`, and `FINALIZE_SUMMARIES`. The two tree
processes run only with tree classification enabled.

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
