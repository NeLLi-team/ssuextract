# Pipeline design

SSUextract detects candidate regions, extracts their sequences, selects a
marker-specific database, and reports taxonomy.

[![SSUextract pipeline from Infernal detection and marker-specific BLAST through default BLAST taxonomy or optional tree-neighbor classification](../assets/figures/pipeline-architecture.svg)](../assets/figures/pipeline-architecture.svg){ .docs-figure }

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

## Optional tree classification

`--tree_classification` adds one tree task per extracted gene. Both marker
indexes are searched before alignment. The marker represented by most of the
best 100 unique subjects is selected; best bit score and the accepted Infernal
model resolve ties.

The selected 100 reference sequences and query are aligned with RF00177 or
RF01960 using `cmalign`. Covariance-model insert columns are masked; match
columns with more than 90% gaps are then removed. IQ-TREE 3 runs a fast search
under `GTR+F+R4` with 1,000 SH-aLRT replicates and one inference thread.
SSUextract sorts references by
patristic distance from the query and assigns the lowest common taxonomy of the
five nearest named references. Equal-distance references at the fifth-neighbor
boundary are included in the same assignment.

A selected marker with fewer than three reference subjects cannot yield a tree.
That query retains its BLAST taxonomy and records the skipped tree attempt;
other extracted genes continue through tree classification.

## Results

The workflow writes extracted FASTA files, parsed hit tables, BLAST results, a
per-hit summary, category counts, and Nextflow execution reports. Tree mode also
writes reference sequences, alignments, trees, logs, QC, and a combined neighbor
table. See
[output files](../reference/outputs.md) for paths and contents.
