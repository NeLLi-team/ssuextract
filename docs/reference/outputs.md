# Output files

| Path | Contents |
| --- | --- |
| `cmsearch_summary.tsv` | One row per extracted SSU region with coordinates and selected taxonomy. |
| `cmsearch_summary.tab` | Per-sample unique-contig counts by annotation category. |
| `extracted/*.fna` | Extracted SSU sequences for each sample/model pair. |
| `stats/*.hits.tsv` | Parsed covariance-model hit metadata. |
| `out/*.out` | Raw Infernal `cmsearch --tblout` output. |
| `m8/*.m8` | BLAST tabular output for each sample/model pair. |
| `m8/merged.m8` | Deterministically merged BLAST output. |
| `pipeline_info/` | Nextflow timeline, report, trace, and DAG. |

## Detailed taxonomy fields

| Column | Contents |
| --- | --- |
| `blast_sseqid` | Stable identifier of the selected database sequence. |
| `centroid_names` | Cluster centroid names associated with that sequence, separated by `|` when exact-sequence deduplication links more than one centroid. Empty for native SILVA and PR2 records. |
| `centroid_taxonomy` | Lowest common calibrated SILVA or PR2 taxonomy supported for those centroids before the cluster-member propagation limit is applied. |
| `centroid_taxonomy_source` | `SILVA`, `PR2`, or a joined source value for the centroid taxonomy. |
| `reference_source` | Source of the selected database sequence: `SILVA`, `PR2`, `IMG`, or an explicit joined value. |
| `taxonomy` | Taxonomy assigned to the selected sequence. IMG cluster-member assignments stop at domain. |
| `taxonomy_assignment_method` | Native, lowest-common-ancestor, ambiguity, or IMG cluster-assignment method. |

Exact-sequence deduplication can associate one database sequence with more than
one source cluster. `centroid_names` retains every centroid among the selected
taxonomy assignments; `centroid_taxonomy` is their lowest common taxonomy.
SILVA centroid names omit the semicolon-delimited lineage from the source
header. PR2 centroid names retain the original public PR2 header label.

## Counting behavior

`cmsearch_summary.tsv` retains hits after overlapping RF00177 and RF01960
matches on the same strand have competed by Infernal E-value and bit score. In
`cmsearch_summary.tab`, one contig contributes at most once to each sample-level
annotation category.

Samples without accepted hits remain in `cmsearch_summary.tab`. Their detailed
summary contains a header and no data rows; per-model hit tables also contain
headers only.

## Coordinates and strand

Hit coordinates are the 1-based inclusive values reported by Infernal. The
extracted sequence includes both endpoints. Reverse-strand records contain the
reverse complement of that interval.
