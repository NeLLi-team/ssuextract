# Output files

| Path | Contents |
| --- | --- |
| `cmsearch_summary.tsv` | One row per extracted SSU region with coordinates and selected taxonomy. |
| `cmsearch_summary.tab` | Per-sample unique-contig counts by annotation category. |
| `blast_top_hits.tsv` | Ranked BLAST evidence with source identifiers, taxonomy metadata, and extracted query sequences. |
| `tree_nearest_neighbors.tsv` | References ordered by patristic distance from each query. Header only when tree mode is off. |
| `phylogeny/<sample>/<detected-model>/<query-key>/` | Per-query reference FASTA, covariance-model alignment, trimmed alignment, IQ-TREE files, assignment, neighbor table, and QC. Written in tree mode. |
| `extracted/*.fna` | Extracted SSU sequences for each sample/model pair. |
| `stats/*.hits.tsv` | Parsed covariance-model hit metadata. |
| `out/*.out` | Raw Infernal `cmsearch --tblout` output. |
| `m8/*.m8` | BLAST tabular output for each sample/model pair. |
| `m8/*.top_hits.tsv` | Per-model input rows for `blast_top_hits.tsv`. |
| `m8/merged.m8` | Deterministically merged BLAST output. |
| `run_provenance.json` | Database, manifest, calibration, and runtime-code content IDs, application version, and assignment policy. |
| `pipeline_info/` | Nextflow timeline, report, trace, and DAG. |

## Detailed taxonomy fields

| Column | Contents |
| --- | --- |
| `blast_sseqid` | Stable identifier of the top subject from the original detection-marker BLAST search. |
| `reference_identifiers` | Public source identifier for the selected sequence, prefixed by `IMG:`, `PR2:`, or `SILVA:`. Multiple exact-sequence records are separated by `|`. |
| `reference_versions` | Deduplicated reference releases for the selected sequence, prefixed by source. |
| `centroid_names` | Cluster centroid names associated with that sequence, separated by `|` when exact-sequence deduplication links more than one centroid. Empty for native SILVA and PR2 records. |
| `centroid_taxonomy` | Lowest common calibrated SILVA or PR2 taxonomy supported for those centroids before the cluster-member propagation limit is applied. |
| `centroid_taxonomy_source` | `SILVA`, `PR2`, or a joined source value for the centroid taxonomy. |
| `reference_source` | Source of the selected database sequence: `SILVA`, `PR2`, `IMG`, or an explicit joined value. |
| `taxonomy` | Supported taxonomy assigned to the extracted query. `Unclassified` means the evidence does not support a call. |
| `taxonomy_assignment_method` | Query decision rule, including exact native matches, calibrated calls, ambiguity, or the reason for abstention. |
| `blast_candidate_taxonomy` | Uncalibrated LCA of eligible near-best BLAST subjects. This is evidence, not the selected query call. |
| `reference_taxonomy`, `reference_taxonomy_source`, `reference_taxonomy_assignment_method` | Preferred taxonomy and method for the selected database subject, separate from the query call. |
| `blast_query_coverage` | Fraction of the query covered by the selected BLAST interval, from 0 to 1. |
| `assignment_candidate_count`, `assignment_unknown_count` | Number of subjects in the assignment set and number without a known lineage. |
| `assignment_candidates_truncated` | Whether the fetched results can omit subjects inside the assignment score window. |
| `query_sequence` | Extracted 16S rRNA gene or 18S rRNA gene sequence. |
| `taxonomy_mode` | `blast` or `tree`, identifying which method supplied `taxonomy`. |
| `blast_taxonomy`, `blast_taxonomy_source`, `blast_taxonomy_domain`, `blast_compartment`, `blast_taxonomy_assignment_method`, `blast_taxonomy_alternatives` | BLAST assignment from the original detection marker, retained independently of tree mode. |
| `tree_model`, `tree_marker` | Covariance model and marker selected after searching both marker indexes. |
| `tree_route_decision`, `tree_route_16s_votes`, `tree_route_18s_votes`, `tree_route_16s_best_bitscore`, `tree_route_18s_best_bitscore` | Hit counts and best scores used to select the alignment route. |
| `tree_taxonomy`, `tree_taxonomy_source`, `tree_taxonomy_domain`, `tree_compartment` | Query taxonomy supported by both the selected marker's BLAST call and direct-reference tree evidence. |
| `tree_assignment_method`, `tree_basis_neighbors` | Decision rule and number of references used for the neighbor LCA. |
| `tree_route_blast_taxonomy`, `tree_route_blast_taxonomy_source`, `tree_route_blast_taxonomy_domain`, `tree_route_blast_compartment`, `tree_route_blast_assignment_method` | Supported BLAST call on the marker selected for the tree. This marker can differ from the original detection route. |
| `tree_route_blast_candidate_taxonomy` | Uncalibrated candidate LCA on the tree-selected marker. |
| `tree_taxonomy_alternatives` | Conflicting lineage evidence retained by tree classification. |
| `tree_nearest_sseqid`, `tree_nearest_reference_identifiers`, `tree_nearest_distance` | Closest reference sequence and its patristic distance from the query. |
| `tree_query_edge_support` | SH-aLRT support on the query's adjacent internal branch when that branch is represented in the rooted output. |
| `tree_inference_model` | Nucleotide substitution model used by IQ-TREE. |

Reference columns, BLAST scores, and assignment-candidate counts describe the
original detection-marker search. They remain unchanged when tree routing
selects another marker. Use `tree_route_blast_*` for that route's BLAST decision
and the tree evidence files for its references. A skipped tree can therefore
have `taxonomy_mode=blast` and a final taxonomy that differs from
`blast_taxonomy`.

Exact-sequence deduplication can associate one database sequence with more than
one source cluster. `centroid_names` retains every centroid among the selected
taxonomy assignments; `centroid_taxonomy` is their lowest common taxonomy.
SILVA centroid names omit the semicolon-delimited lineage from the source
header. PR2 centroid names retain the original public PR2 header label.

## Ranked BLAST evidence

`blast_top_hits.tsv` reports the first `--top_hits` subjects for each extracted
query. It also reports the highest-ranked IMG, PR2, and SILVA subject among the
fetched BLAST candidates when that source occurs below the requested cutoff.
`hit_rank` is the rank among those candidates. `selection_reason` distinguishes
`overall_top_n`, `assignment_candidate`, `equal_best_assignment`, `best_IMG`, `best_PR2`, and
`best_SILVA` rows.

Each row joins the stable `blast_sseqid` to its public source identifier,
reference version, preferred taxonomy, centroid evidence, BLAST statistics, and
the extracted query sequence. The table can contain more than `--top_hits` rows
per query because all assignment candidates and the best subject from
another reference source are retained. The source-specific supplemental rows do
not change the selected taxonomy in `cmsearch_summary.tsv`.

An IMG subject without native or calibrated centroid taxonomy remains
`Unclassified`. SSUextract reports the best PR2 and SILVA evidence as separate
rows instead of assigning a lower-scoring lineage to that IMG subject.

## Tree-neighbor evidence

`tree_nearest_neighbors.tsv` contains every reference used in a tree. Rows are
ordered by `tree_distance` within each query. `tree_lineage_basis` identifies
the direct preferred taxonomy used for classification. Centroid labels remain
reference evidence. `used_for_assignment` marks the references included in the
LCA, including explicit ambiguity. Every reference tied by distance at the
configured neighbor boundary is included, so `tree_basis_neighbors` can exceed
`--tree_assignment_neighbors`.

Each per-query directory retains `references.tsv`, `references.fna`,
`cmalign_input.fna`, the raw and trimmed covariance-model alignments,
`alignment_qc.json`, the IQ-TREE report and tree, the assignment and neighbor
tables, and `tool_versions.txt`. `references.tsv` links the compact tree leaf
labels to stable database IDs and public source identifiers. The query directory
uses a deterministic `q_<hex>` key; `task.json` records the original query name
and the detected and tree-selected models.

If the selected marker has fewer than three reference subjects, SSUextract skips
the tree. It uses that marker route's BLAST decision with
`taxonomy_mode=blast`, and `tree_assignment_method` records
`tree_skipped_insufficient_references`. Other queries in the same run continue.

## Counting behavior

`cmsearch_summary.tsv` retains hits after overlapping RF00177 and RF01960
matches on the same strand have competed by Infernal E-value and bit score. In
`cmsearch_summary.tab`, one contig contributes at most once to each sample-level
annotation category.

Counts use the selected query assignment. A supported tree call can contribute
when the original marker's BLAST search has no subject.

Samples without accepted hits remain in `cmsearch_summary.tab`. Their detailed
summary contains a header and no data rows; per-model hit tables also contain
headers only.

## Coordinates and strand

Coordinates use 1-based inclusive contig positions. A simple record has one
Infernal hit interval. An assembled record spans its joined component intervals,
including any gap between them. Reverse-strand records contain the reverse
complement of the enclosing interval.

| Column | Contents |
| --- | --- |
| `length` | Length of the extracted sequence, including gaps between joined fragments. |
| `hit_length` | Union length of CM-hit sequence spans. Used by `--min_extract_length`. Insertions within one CM hit count toward this length. |
| `model_coverage` | Union length of the component model-coordinate intervals divided by model length (`CLEN`), from 0 to 1. |
| `fragment_count` | Number of CM hits in the extraction. |
| `component_coordinates` | Comma-separated, inclusive component intervals in gene order. |
| `is_assembled` | `True` when two or more partial hits form one extraction. |

These fields describe SSU detection. Model coverage includes positions between
hit endpoints; it does not count individually aligned model states. These fields
do not estimate genome completeness or contamination.

## Query identifiers and BLAST columns

Emitted query IDs use
`sample|model|contig|start-end|strand_<strand>|<sequence_type>`. FASTA, hit tables,
summaries, and merged BLAST rows use the same ID. Sample and model identifiers
distinguish identical contig names in separate inputs.

Runtime M8 output has 14 tab-separated columns: the 12 standard BLAST columns,
then `qlen` and `slen`. Scripts that require exactly 12 columns must account for
these two length fields. The annotation reader also accepts 12-column files,
but missing subject length cannot support an exact query assignment.
