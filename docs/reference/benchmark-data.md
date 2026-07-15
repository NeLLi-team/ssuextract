# Benchmark data

`docs/data/pipeline_benchmark.tsv` records one warm-up and three measured
full-example runs for commit `620b1a6` and the 1.1.0 candidate. Each run used two
threads per search or BLAST task and a fresh Nextflow work directory. The source
timing files are preserved under `tasks/benchmark/` in the release worktree.

`docs/data/database_profile_benchmark.tsv` records the same warm-up and measured
trial structure for the legacy, curated, and IMG-enhanced database profiles.
All three profiles use the same assembly, models, two-thread search tasks, and
eight-CPU Nextflow task ceiling.

Elapsed time and maximum resident-set size come from GNU time around the
Nextflow command. The RSS value is the command's maximum reported value, not the
sum of memory used concurrently by every workflow process.

`docs/data/database_storage.tsv` records reproducible logical file bytes,
compressed archive bytes, and archive SHA-256 values. Every file under each
installed profile is counted by the same rule and partitioned into source
FASTA, BLAST indexes, Parquet metadata, or other files. The legacy row has no
release archive. Filesystem block allocation is excluded because it depends on
the installation host.

`scripts/prepare_benchmark_data.py` validates all 12 profile trials, both release
manifests, and both archives before writing the documentation tables atomically.

Warm-up trials remain in the source tables but are excluded from descriptive
medians. The notebook asserts profile order, row counts, positive measurements,
component totals, and archive hash format before drawing figures.
