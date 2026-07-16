# Database profiles

| Profile | Reference content | Additional content | Search indexes |
| --- | --- | --- | --- |
| `curated` | SILVA 138.2 NR99 and PR2 5.1.1 | None | Separate BLAST v5 indexes for 16S and 18S. |
| `img` | SILVA 138.2 NR99 and PR2 5.1.1 | IMG sequences from the eukcensus 16S/18S collections | Separate BLAST v5 indexes for 16S and 18S. |

RF00177 routes to the 16S index. RF01960 routes to the 18S index. The mapping is
stored in `config/model_markers.json`.

Database packages use a version independent of the SSUextract application.
The current database profile version is `1.0.0`.

## Runtime files

| Path | Contents |
| --- | --- |
| `manifest.json` | Profile identity, artifact paths, sizes, and SHA-256 digests. |
| `provenance.json` | Source versions, build contracts, software versions, and source-tree fingerprints. |
| `blast/16S.*` | BLAST v5 nucleotide index for the 16S marker. |
| `blast/18S.*` | BLAST v5 nucleotide index for the 18S marker. |
| `tables/sequences.parquet` | Content-addressed sequence identifiers, lengths, hashes, and marker membership. |
| `tables/preferred_taxonomy.parquet` | Selected taxonomy and explicit cross-domain conflict state. |
| `tables/source_records.parquet` | Normalized source record provenance. |
| `tables/taxonomy_assignments.parquet` | Native and derived taxonomy evidence. |
| `tables/img_location.parquet` | IMG taxon identifier and valid latitude/longitude values. |

Raw source FASTA files, source project descriptions, contacts, email addresses,
comments, cluster tables, and centroid identifiers are not distributed in a
runtime profile.

## Curated profile counts

The `curated` profile contains 609,298 unique exact sequences from 683,597
source records. The 16S index contains 416,021 sequences; the 18S index contains
193,282. Exact sequences assigned across domains remain marked as conflicts.
