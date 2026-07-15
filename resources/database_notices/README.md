# SSUextract database profile

This archive contains one SSUextract database profile. `manifest.json` records
the SHA-256 digest and byte count of each file. `provenance.json` records source
versions, source checksums, sequence counts, transformations, and build tools.

The `blast/16S` and `blast/18S` prefixes are NCBI BLAST v5 nucleotide
databases. `tables/preferred_taxonomy.parquet` supplies the taxonomy used by
SSUextract. The remaining Parquet tables retain exact-sequence hashes, source
records, native and derived taxonomy assignments, and the privacy-filtered IMG
location allowlist.

SSUextract assigns SILVA taxonomy to bacterial and archaeal references. PR2
supplies eukaryotic taxonomy and host taxonomy for organellar records, with the
organelle stored in the `compartment` field. Exact sequences that have both a
prokaryotic SILVA assignment and a eukaryotic organellar PR2 assignment remain
explicitly ambiguous; `taxonomy_alternatives` retains both native paths.

See the files in `LICENSES/` before redistributing the profile.
