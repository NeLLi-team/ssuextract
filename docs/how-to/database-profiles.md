# Install and select database profiles

The `curated` profile is the default. The `img` profile uses the same current
SILVA and PR2 references and adds IMG-derived 16S and 18S sequences.

## Install the curated profile

```bash
pixi run setup -- --database_profile curated
```

## Install the IMG profile

```bash
pixi run setup -- --database_profile img
```

The installer resolves the archive from `config/database_catalog.json`, checks
its byte count and SHA-256 digest, validates the internal manifest, and publishes
the profile under the configured database root.

## Select a profile for one run

```bash
pixi run ssuextract -- \
  --querydir data/my_dataset \
  --outdir results/my_dataset-img \
  --database_profile img
```

## Use another database root

```bash
pixi run ssuextract -- \
  --database_path /path/to/ssuextract-databases \
  --database_profile curated
```

The root must contain the selected profile directory and a valid
`manifest.json`. See the [database reference](../reference/database-profiles.md)
for its contents and taxonomy contract.

