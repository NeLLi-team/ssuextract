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
