# Pipeline design

SSUextract separates detection, extraction, annotation, and reporting so each
stage owns one representation and one set of output files.

```mermaid
flowchart LR
    A[Assembly FASTA] --> B[Infernal search]
    B --> C[Typed hit table]
    C --> D[Complete interval extraction]
    D --> E{Marker route}
    E -->|RF00177| F[16S BLAST index]
    E -->|RF01960| G[18S BLAST index]
    F --> H[Taxonomy selection]
    G --> H
    H --> I[Detailed and category summaries]
```

Infernal supplies 1-based inclusive coordinates. The extraction stage retains
both endpoints and reverse-complements negative-strand intervals. The typed hit
table carries sample, model, contig, coordinates, strand, and sequence identity
into later stages without reconstructing metadata from filenames.

Marker routing keeps 16S queries out of the 18S database and 18S queries out of
the 16S database. Final reporting scans each BLAST result once and sorts emitted
rows for deterministic tabular output.

