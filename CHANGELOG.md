# Changelog

## [Unreleased]

### Added

- Add an `Unclassified` column to `cmsearch_summary.tab`. It counts the contigs
  with at least one query whose selected taxonomy is `Unclassified`.
- Add `pixi run setup --help`.

### Fixed

- Record the Infernal version in each per-query `tool_versions.txt`. The file
  held the `cmalign` banner line and no version.
- List each taxonomy alternative once in `tree_taxonomy_alternatives`. Bacterial
  and archaeal alternatives appeared twice in tree mode.
- Show the validation error when a database profile is not ready. A missing
  profile is installed. A profile with a missing or changed artifact is
  replaced, but never under an explicit `--database_path`. Every other
  validation failure stops the command without a download.
- Keep a valid profile when `install` runs without `--force`, including a
  profile restored from an interrupted replacement.
- Reject unknown `setup` arguments, a database option without a value, and a
  profile name outside the catalog with exit status 2 before any change.
- Save `config/local.config` only after the database is ready and only when
  its content changes.
- Exit with a nonzero status when `setup --update` cannot reach the release
  record.
- Expand `~` at the database path prompt and save an absolute path.
- Reject a manifest artifact that is a symbolic link, including a link to a
  file inside the profile.
- Download an archive once when the server ignores the byte range and the full
  response fails its checksum.
- Report an unusable database root as a `database-manager` error instead of a
  Python traceback.
- Declare Python 3.14 as the minimum version. Database archives are `.tar.zst`
  files, which `tarfile` reads from Python 3.14.
- Validate the inputs of `img_chunked_search.py prepare --force` before any
  deletion, and delete only a chunk directory that holds an earlier plan.
- Record the `-title` argument in the `makeblastdb` command template of the
  database provenance.

### Changed

- Pass `nextflow lint` and run under the strict syntax parser
  (`NXF_SYNTAX_PARSER=v2`). `process.resourceLimits` applies `--max_cpus`,
  `--max_memory`, and `--max_time`.
- Write the per-run DAG report as `.html`; it needs no Graphviz.
- Require the `qlen` and `slen` columns in BLAST input and the taxonomy and
  source-record tables in `annotate_hits.py`.
- `database_manager.py validate` exits with status 3 for a profile that is not
  installed and 4 for a damaged artifact.
- `pixi run dryrun` uses the saved database selection.
- Compare competing model hits only within one contig and strand. The step
  took 8.6 s for 16,000 hits and now takes 0.03 s, with identical output.
- Load database release tables in bulk. In a comparison on five tables of
  20,000 rows, the write time fell from 122 s to 0.7 s and the Parquet files
  were byte-identical. Values of a type other than the column type are now
  rejected instead of cast.
- Read the reference taxonomy table for IMG classification in one pass.

### Removed

- Remove support for the unprofiled `silva-138-1_pr2-4-12` database layout,
  the `--querydir` migration message, the two old `config/local.config`
  formats, and 12-column BLAST input.
- Remove `database_manager.py resolve` with `--markers` and
  `--marker-override`.
- Remove the Nextflow profiles `debug`, `conda`, `local`, `slurm`, and `test`
  with `config/environment.yml`.
- Remove the Pixi tasks `download-db`, `dag`, `run-verbose`, and `report` and
  the `dev` environment.
- Remove `scripts/get_cmsequences.py`, `docs/requirements.txt`, and
  `docs/releases/1.1.0.md`. Nothing used them.

## [1.2.2] - 2026-09-04

### Fixed

- Limit fragment joins to nearby, ordered model spans. Apply the minimum-length
  filter to supported bases, and report model coverage and fragment coordinates.
- Bind resumed tasks to database contents and runtime code. Record these hashes
  with the search policy and target limit in run provenance.
- Require query coverage and candidate agreement for taxonomy calls. Report
  reference labels and supported query assignments in separate fields.
- Retain ambiguous and unknown nearest tree references, limit tree calls to
  supported route taxonomy, and count successful alternate marker routes.
- Use explicit BLASTN search settings and honor the CPU limit in tree tasks.
- Count centroid calibration errors by predicted source and domain. Reject
  calibration artifacts for a different policy, profile, or use case.
- Include sample and model identifiers in extracted sequence identifiers.

### Changed

- Require profile-bound runtime calibration for nonexact taxonomy calls. Without
  it, only unanimous, full-length exact matches to native references can assign
  taxonomy. Released database profiles v1.0.2 have no runtime calibration;
  other queries report `Unclassified` and retain reference evidence.
- Add schema-3 calibration support with optional identity thresholds per rank.
  Tests of fitted rules and replacement database builds remain in progress.
- Document that model coverage does not measure genome completeness or
  contamination. The pipeline has no dedicated chimera detector.

## [1.2.1] - 2026-07-22

### Added

- Accept one `.fna`, `.fa`, or `.fasta` file as input. When `--outdir` is
  omitted, a single file writes to `results/<file-stem>`.
- Accept `-q` as the short form of `--query` through the Pixi command wrapper.

### Changed

- Rename the FASTA file or directory input option from `--querydir` to
  `--query`. The removed name exits with a migration message.

## [1.2.0] - 2026-07-20

### Added

- Write `blast_top_hits.tsv` with ranked BLAST evidence, public source
  identifiers, source versions, preferred taxonomy, centroid evidence, and
  alignment statistics.
- Add `--top_hits` to set the number of overall subjects reported per extracted
  query. Equal-best assignment subjects and the best subject from each available
  reference source remain in the table outside that cutoff.
- Append public reference identifiers, reference versions, and the extracted
  query sequence to `cmsearch_summary.tsv`.
- Add optional `--tree_classification`. The mode searches both marker indexes,
  aligns the query with 100 selected references using the matching covariance
  model, infers an IQ-TREE 3 tree, and assigns the lowest common taxonomy of
  the nearest named references.
- Write distance-ranked tree evidence to `tree_nearest_neighbors.tsv` and keep
  the BLAST assignment in separate summary fields when tree taxonomy is used.

### Fixed

- Retain the BLAST assignment for an individual tree-mode query when fewer than
  three references are available, without stopping other queries.
- Run IQ-TREE with one inference thread to avoid thread-dependent ordering of
  near-equal branch distances.

## [1.1.1] - 2026-07-19

### Fixed

- Report the matched database sequence separately from the centroid used to classify
  an IMG cluster.
- Append centroid names, calibrated centroid taxonomy, and centroid taxonomy source
  to detailed output tables without changing existing column positions.
- Publish database version 1.0.2 with the centroid evidence required by these output
  fields. Sequence content and preferred member assignments are unchanged from
  database version 1.0.1.

## [1.1.0] - 2026-07-11

### Fixed

- Extract complete 1-based inclusive `cmsearch` intervals on both strands. Earlier
  versions omitted one nucleotide from every extracted interval.
- Retain FASTA records whose headers contain descriptions.
- Retain the final record in legacy sequence maps without a trailing newline.
- Produce valid empty summaries when no covariance-model hits pass filtering.
- Accept absolute query, model, database, and output paths.

### Changed

- Replace the multi-script temporary-directory pipeline with typed hit parsing and
  four explicit Nextflow processes.
- Parse each BLAST result once and select the highest-bit-score annotation.
- Generate detailed and category summaries from one canonical metadata stream.
- Sort final tables and merged BLAST output deterministically.
- Publish only owned artifacts; input assemblies are no longer copied into results.
- Validate identifiers, numeric parameters, covariance-model lengths, sequence
  references, and BLAST database files before or during processing.
- Track the resolved Pixi environment and test Linux in continuous integration.

### Removed

- Remove the legacy `cmprocessing.py`, `cmsearchout_extract_by_position_size.py`,
  `get_cmstats.py`, `get_table.py`, and `rename_fnaheaders.py` stages.

## [1.0.0] - 2025-06-06

- Migrate the pipeline from Snakemake to Nextflow.
- Add the RF00177 and RF01960 covariance models.
- Preserve the former Snakemake implementation on the `ssuextract-snk` branch.
