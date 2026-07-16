# SSUextract

SSUextract searches assembled contigs for 16S and 18S rRNA. It extracts each
complete Infernal hit interval and assigns taxonomy with the reference database
for that marker.

[Run the bundled example](tutorials/first-run.md) ·
[Parameters](reference/cli.md) ·
[Output files](reference/outputs.md)

## Pipeline

[![SSUextract pipeline from assembly FASTA through Infernal detection, interval extraction, marker-specific BLAST, and result files](assets/figures/pipeline-architecture.svg)](assets/figures/pipeline-architecture.svg){ .docs-figure }

## Run the bundled example

```bash
git clone https://github.com/NeLLi-team/ssuextract.git
cd ssuextract
pixi install --frozen
pixi run setup
pixi run example
```

`pixi run setup` installs the default `curated` database profile. The example
writes extracted sequences, per-hit annotations, category counts, and Nextflow
reports under `results/smoke/`.

## Steps

1. Infernal searches each assembly with 16S and 18S covariance models.
2. SSUextract keeps accepted coordinates and extracts the complete interval on
   the reported strand.
3. Each sequence is searched against the 16S or 18S index selected by its model.
4. The workflow writes extracted FASTA files, one row per hit, category counts,
   and execution reports.

[Pipeline details](explanation/pipeline.md) ·
[Run assembled genomes or metagenomes](how-to/run-assemblies.md)
