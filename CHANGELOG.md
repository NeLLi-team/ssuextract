# Changelog

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
