# Changelog

## [1.3.0] - 2026-11-10

### Added

- `Unclassified` column in `cmsearch_summary.tab`.
- `pixi run setup --help` and `SSUEXTRACT_SKIP_UPDATE_CHECK=1`.

### Changed

- `cmsearch_summary.tsv`: remove `blast_tied_subjects` and
  `blast_ties_truncated`; write booleans as `true`/`false`; use one
  `blast_pident` format in both tables.
- Resolve relative command-line paths against the launch directory.
- Stop writing a per-run DAG report.
- Tree mode reuses the annotation BLAST result; tree processes request 8 GB
  and 2 hours.
- Runs check database file sizes and ask Zenodo at most once a day;
  `pixi run setup` keeps the full check.
- Stream the final tables and resolve competing hits per contig, for lower
  memory and run time on large inputs.
- Support the strict Nextflow syntax (`NXF_SYNTAX_PARSER=v2`);
  `process.resourceLimits` applies the `--max_*` limits.
- Require `qlen` and `slen` in BLAST input and Python 3.14.
- Remeasure the performance page with database v1.0.2.

### Fixed

- Record the Infernal version in `tool_versions.txt`; list each entry once in
  `tree_taxonomy_alternatives`.
- Stop on duplicate sample names, malformed taxonomy, and duplicate extracted
  names.
- Database setup: show validation errors; install only a missing profile and
  replace only a damaged one, never under `--database_path`; reject invalid
  `setup` arguments; save the selection after the database is ready; fail
  `setup --update` when Zenodo is unreachable.
- Reject symbolic-link artifacts; download an archive once after an ignored
  byte range.
- Report 2,172 IMG taxa without complete coordinates (was 2,171).

### Database build

- Write tables in bulk, stream BLAST results in the IMG classifier, and search
  centroids only for clusters with an IMG member.
- Fail a build on an empty rank inside a lineage, a modified git tree, or an
  existing output.
- Classification must run again before the next build (new `marker` column
  and `evidence_id`).

### Removed

- Legacy `silva-138-1_pr2-4-12` database layout, `--querydir` message, old
  `config/local.config` formats, and 12-column BLAST input.
- Nextflow profiles `debug`, `conda`, `local`, `slurm`, `test` and
  `config/environment.yml`.
- Pixi tasks `download-db`, `dag`, `run-verbose`, `report` and the `dev`
  environment.
- `database_manager.py resolve`, the single-job IMG search, and unused
  scripts, options, and figures.

## [1.2.2] - 2026-09-04

### Fixed

- Join only nearby, ordered fragments; apply the minimum length to supported
  bases; report model coverage and fragment coordinates.
- Bind resumed tasks to database contents and runtime code, and record both
  in run provenance.
- Require query coverage and candidate agreement for taxonomy calls; report
  reference labels and query assignments in separate fields.
- Keep ambiguous and unknown nearest tree references and limit tree calls to
  the supported route taxonomy.
- Use explicit BLASTN settings and honor the CPU limit in tree tasks.
- Include sample and model identifiers in extracted sequence identifiers.

### Changed

- Nonexact taxonomy calls require a runtime calibration in the database
  profile. Database v1.0.2 has none, so those queries report `Unclassified`
  and keep the reference evidence.
- Model coverage does not measure genome completeness or contamination.

## [1.2.1] - 2026-07-22

- Accept one `.fna`, `.fa`, or `.fasta` file as input; it writes to
  `results/<file-stem>` by default.
- Rename `--querydir` to `--query`, with `-q` as the short form in the Pixi
  wrapper.

## [1.2.0] - 2026-07-20

### Added

- `blast_top_hits.tsv` with ranked BLAST evidence, and `--top_hits` to set
  its size.
- Reference identifiers, reference versions, and the query sequence in
  `cmsearch_summary.tsv`.
- Optional `--tree_classification`: align the query with 100 references,
  infer an IQ-TREE 3 tree, and assign the lowest common taxonomy of the
  nearest named references. Evidence goes to `tree_nearest_neighbors.tsv`.

### Fixed

- Keep the BLAST assignment for a tree-mode query with fewer than three
  references.
- Run IQ-TREE with one thread for a stable order of near-equal distances.

## [1.1.1] - 2026-07-19

- Report the matched database sequence separately from the centroid that
  classifies an IMG cluster, with centroid names and taxonomy as new columns.
- Publish database v1.0.2 with the centroid evidence; sequences and member
  assignments are unchanged from v1.0.1.

## [1.1.0] - 2026-07-11

### Fixed

- Extract complete 1-based inclusive `cmsearch` intervals on both strands;
  earlier versions dropped one nucleotide.
- Keep FASTA records whose headers contain descriptions.
- Write valid empty summaries when no hit passes filtering.
- Accept absolute query, model, database, and output paths.

### Changed

- Replace the multi-script pipeline with typed hit parsing and four Nextflow
  processes.
- Select the highest-bit-score annotation and sort final tables
  deterministically.
- Publish only pipeline outputs; input assemblies are not copied into results.
- Validate identifiers, parameters, model lengths, and BLAST database files.

## [1.0.0] - 2025-06-06

- Migrate the pipeline from Snakemake to Nextflow.
- Add the RF00177 and RF01960 covariance models.
