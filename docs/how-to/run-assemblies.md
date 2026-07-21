# Run assembled genomes or metagenomes

Place `.fna`, `.fa`, or `.fasta` assemblies in one directory. The file basename
becomes the sample identifier.

## Run the directory

```bash
pixi run ssuextract --querydir data/my_dataset --outdir results/my_dataset --threads_per_job 4
```

SSUextract validates query, model, database, and output paths before scheduling
the workflow. Basenames may contain letters, numbers, `.`, `_`, and `-`.

## Report more reference hits

Set the number of ranked BLAST subjects written for each extracted query:

```bash
pixi run ssuextract --querydir data/my_dataset --outdir results/my_dataset --top_hits 10
```

Inspect the rank, selection reason, public reference identifier, taxonomy,
centroid taxonomy, identity, alignment length, and bit score:

```bash
cut -f1,4-5,7,9-10,17,19-20,28 results/my_dataset/blast_top_hits.tsv
```

The output contains the first 10 subjects, all equal-best assignment subjects,
and the highest-ranked IMG, PR2, and SILVA subjects among the fetched candidates.
The supplemental source rows do not change the assignment in
`cmsearch_summary.tsv`.

## Classify from tree neighbors

Enable tree classification for each extracted 16S rRNA gene or 18S rRNA gene:

```bash
pixi run ssuextract --querydir data/my_dataset --outdir results/my_dataset-tree --tree_classification --threads_per_job 4
```

Tree mode searches both marker indexes. The 100 best unique subjects determine
whether the 16S rRNA gene or 18S rRNA gene covariance model is used for
alignment. SSUextract aligns the query and references with `cmalign`, removes
covariance-model insert columns and match columns with more than 90% gaps, then
runs IQ-TREE 3 with `GTR+F+R4`, fast tree search, and 1,000 SH-aLRT replicates.

Inspect the selected tree taxonomy and the nearest reference branches:

```bash
cut -f1-3,29-53 results/my_dataset-tree/cmsearch_summary.tsv
cut -f1-11,13-20 results/my_dataset-tree/tree_nearest_neighbors.tsv
```

Per-query alignments, trees, model logs, reference FASTA files, and alignment
QC are under `results/my_dataset-tree/phylogeny/`. The BLAST result remains in
the `blast_taxonomy` fields and `blast_top_hits.tsv`.

If the selected marker yields fewer than three references, SSUextract keeps the
BLAST taxonomy for that query and records the skipped tree attempt. Other
queries continue.

## Bound local resource use

Set Nextflow task ceilings independently from the per-search thread count:

```bash
pixi run ssuextract --querydir data/my_dataset --outdir results/my_dataset --threads_per_job 4 --max_cpus 8 --max_memory 64.GB
```

`--threads_per_job` controls each `cmsearch`, BLAST, and `cmalign` task.
IQ-TREE uses one thread so near-tied branch distances are reproducible.
`--max_cpus` limits
the resources that Nextflow assigns to one task; it does not set the number of
concurrent tasks.

## Resume an interrupted run

Run the same command with Nextflow's resume flag:

```bash
pixi run ssuextract -resume --querydir data/my_dataset --outdir results/my_dataset
```

Keep the original `work/` directory and command parameters. Nextflow reuses
matching completed tasks and reruns missing work.
