# Pipeline design

SSUextract detects candidate regions, extracts their sequences, selects a
marker-specific database, and reports taxonomy.

[![SSUextract pipeline from assembly FASTA through Infernal detection, interval extraction, marker-specific BLAST, and result files](../assets/figures/pipeline-architecture.svg)](../assets/figures/pipeline-architecture.svg){ .docs-figure }

*Click the figure to open the full-size SVG.*

## Detection and extraction

Infernal supplies 1-based inclusive coordinates. The extraction stage retains
both endpoints and reverse-complements negative-strand intervals. The typed hit
table carries sample, model, contig, coordinates, strand, and sequence identity
into annotation.

RF00177 and RF01960 belong to Rfam clan CL00111 and can recognize the same SSU
locus. For overlapping same-strand hits, SSUextract retains the lower Infernal
E-value. Equal E-values are resolved by the higher bit score. An exact tie stops
the run instead of selecting a model without supporting evidence.

Any overlapping pair of different models without an explicit scientific
competition rule stops the run; SSUextract does not silently retain both hits
when their relationship is unknown.

## Marker-specific annotation

RF00177 sequences use the 16S rRNA gene BLAST index. RF01960 sequences use the
18S rRNA gene index. The selected `curated` or `img` profile supplies both
indexes and their taxonomy tables.

## Results

The workflow writes extracted FASTA files, parsed hit tables, BLAST results, a
per-hit summary, category counts, and Nextflow execution reports. See
[output files](../reference/outputs.md) for paths and contents.
