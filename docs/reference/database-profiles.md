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
size from `config/database_catalog.json`. After the profile is ready, the
selected path and profile are stored as `database_path` and `database_profile`
in `config/local.config`. Commands read only these `key=value` lines and write
the file only when its content changes. `pixi run setup` accepts only
`--database_profile`, `--database_path`, `--update`, and `-h` or `--help`; any
other argument, a missing option value, or a profile name outside the catalog
exits with status 2 before any installation or Zenodo request.

## Validation

`database_manager.py validate --root DIR --profile NAME` prints the validation
error to standard error and sets the exit status:

| Status | Meaning |
| --- | --- |
| 0 | The profile validates. |
| 3 | The profile directory or its `manifest.json` is absent. |
| 4 | A manifest artifact is missing, empty, or differs in size or SHA-256 digest. |
| 2 | Any other failure, such as an invalid manifest, an artifact that is a symbolic link, a missing `blastdbcmd`, or a BLAST index check failure. |

`pixi run setup`, `pixi run ssuextract`, and `pixi run example` show the error
and install an absent profile (status 3). They replace a damaged profile
(status 4) by installing with `--force`: the old profile stays in place until
the new archive validates. They do not replace a damaged profile when
`--database_path` is given on the command line: the message names the directory
to remove before running setup again. Only `pixi run setup --update` with that
path replaces the profile, with a newer release. On every other failure they
stop with a nonzero status without a download. `database_manager.py install`
without `--force`
keeps an existing profile that validates, including one restored from an
interrupted replacement.

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

Interactive terminals redraw a width-bounded progress line with written bytes,
transfer rate, and estimated time remaining. Non-interactive logs include the
start, completion, and intermediate progress every 10 seconds.

`pixi run ssuextract` and `pixi run example` check Zenodo before starting
Nextflow. The check has a five-second total timeout. Interactive runs offer to
install an available release. Non-interactive runs and runs with an explicit
`--database_path` report the update and continue with the installed profile. A
valid installed profile remains usable when the check times out, Zenodo is
unavailable, or the remote release contract fails validation. Only
`pixi run setup --update` treats an unavailable check as a failure and exits
with a nonzero status.

## Runtime files

| Path | Contents |
| --- | --- |
| `manifest.json` | Profile identity, artifact paths, sizes, and SHA-256 digests. |
| `provenance.json` | Source versions, build contracts, software versions, and source-tree fingerprints. |
| `blast/16S.*` | BLAST v5 nucleotide index for 16S rRNA genes. |
| `blast/18S.*` | BLAST v5 nucleotide index for 18S rRNA genes. |
| `tables/sequences.parquet` | Content-addressed sequence identifiers, lengths, hashes, and marker membership. |
| `tables/preferred_taxonomy.parquet` | Selected sequence taxonomy, centroid names and taxonomy, and explicit cross-domain conflict state. |
| `tables/source_records.parquet` | Normalized source record provenance. |
| `tables/taxonomy_assignments.parquet` | Native and derived taxonomy evidence. |
| `tables/img_location.parquet` | IMG taxon identifier and valid latitude/longitude values. |

The optional manifest key `runtime_calibration` names a JSON file in the profile.
That file must also appear in the artifact list with its size and SHA-256 hash.
Runtime calibration uses schema 3, `calibration_use_case=runtime_query`,
`classification_policy=runtime_lca_v1`, and
`stratum_basis=predicted_candidate_source_set_and_lca_domain`. The root `status`
must be `calibrated` when rank caps exist and `failed` when none exist.
The positive integer `max_targets` must match the run's `--max_blast_targets`.
Its reference digest covers the BLAST files and preferred-taxonomy table.
Rank caps or identity rules apply to the predicted marker, source set, and domain.

The schema-3 `img_centroid` calibration has a separate purpose: it limits labels
assigned during IMG database construction. It cannot substitute for runtime
query calibration. Released v1.0.2 profiles lack runtime calibration; nonexact
query calls remain unclassified while reference evidence stays available.

Centroid calibration can supply optional `rank_rules` for each predicted source
and domain. Rules specify consecutive ranks from domain to the rank cap and
nondecreasing minimum candidate identities. Every source route must pass. The
classifier uses the lowest cap and highest required identity across those
routes. Applying corrected centroid calibration to stored IMG labels requires a
new profile build.

Raw source FASTA files, source project descriptions, contacts, email addresses,
comments, and cluster tables are not distributed in a runtime profile. The IMG
profile retains only the centroid names required to interpret cluster-derived
taxonomy.

## Curated profile counts

The `curated` profile contains 609,298 unique exact sequences from 683,597
source records. The 16S rRNA gene index contains 416,021 sequences; the 18S
rRNA gene index contains 193,282. Exact sequences assigned across domains
remain marked as conflicts.
