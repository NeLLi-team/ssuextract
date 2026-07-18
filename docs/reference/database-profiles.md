# Database profiles

| Profile | Reference content | Additional content | Search indexes |
| --- | --- | --- | --- |
| `curated` | SILVA 138.2 NR99 and PR2 5.1.1 | None | Separate BLAST v5 indexes for 16S rRNA genes and 18S rRNA genes. |
| `img` | SILVA 138.2 NR99 and PR2 5.1.1 | IMG sequences from the eukcensus 16S rRNA gene and 18S rRNA gene collections | Separate BLAST v5 indexes for 16S rRNA genes and 18S rRNA genes. |

RF00177 routes to the 16S rRNA gene index. RF01960 routes to the 18S rRNA gene
index. The mapping is stored in `config/model_markers.json`.

Database packages use a version independent of the SSUextract application.
The version is stored in each profile's `manifest.json` and printed before a
pipeline run.

## Installation and update checks

`pixi run setup` lists each profile with the version and compressed download
size from `config/database_catalog.json`. The selected path and profile are
stored as `database_path` and `database_profile` in `config/local.config`.

The installer accepts a Zenodo release only when all of these values agree:

- the configured concept record and the latest release record;
- the exact four-file release inventory and Zenodo MD5 values;
- the release manifest, archive sizes, and archive SHA-256 values;
- the two entries in `SHA256SUMS`.

Archive extraction takes place in a staging directory. The profile becomes the
installed profile after its manifest files and BLAST indexes validate. A failed
replacement restores the installed profile.

Archive downloads use a 60-second network read timeout and stop after five
consecutive transient failures. A partial archive is keyed by profile, database
version, archive SHA-256 digest, and host. A later setup run resumes it only when
the server returns a matching HTTP byte range. A rejected range or failed
checksum permits one clean restart for that condition. The installer checks the
final byte count and SHA-256 digest before extraction. It removes the verified
download cache after using the archive.

Interactive terminals redraw one progress line with written bytes, transfer
rate, and estimated time remaining. Non-interactive logs retain each timed
update as a separate line.

`pixi run ssuextract` checks Zenodo before starting Nextflow. The check has a
five-second total timeout. A valid installed profile remains usable when the
check times out, Zenodo is unavailable, or the remote release contract fails
validation. Updates require `pixi run setup --update` or confirmation in
interactive setup.

## Runtime files

| Path | Contents |
| --- | --- |
| `manifest.json` | Profile identity, artifact paths, sizes, and SHA-256 digests. |
| `provenance.json` | Source versions, build contracts, software versions, and source-tree fingerprints. |
| `blast/16S.*` | BLAST v5 nucleotide index for 16S rRNA genes. |
| `blast/18S.*` | BLAST v5 nucleotide index for 18S rRNA genes. |
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
source records. The 16S rRNA gene index contains 416,021 sequences; the 18S
rRNA gene index contains 193,282. Exact sequences assigned across domains
remain marked as conflicts.
