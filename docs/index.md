# SSUextract

SSUextract searches assembled contigs for 16S rRNA genes and 18S rRNA genes.
It extracts each complete Infernal hit interval and assigns taxonomy with the
reference database for that marker.

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

`pixi run setup` lists the database profiles. A new checkout defaults to
`curated`; later runs use the saved selection. Before Nextflow starts, an
interactive run offers to install a newer database release when one is
available. The example writes extracted sequences, per-hit annotations,
category counts, and Nextflow reports under `results/smoke/`.

## Steps

1. Infernal searches each assembly with 16S rRNA gene and 18S rRNA gene
   covariance models.
2. Overlapping RF00177 and RF01960 hits on the same strand compete by Infernal
   E-value, then bit score.
3. SSUextract extracts each retained interval on the reported strand.
4. Each sequence is searched against the 16S rRNA gene or 18S rRNA gene index
   selected by its model.
5. The workflow writes extracted FASTA files, one row per hit, category counts,
   and execution reports.

[Pipeline details](explanation/pipeline.md) ·
[Run assembled genomes or metagenomes](how-to/run-assemblies.md)
