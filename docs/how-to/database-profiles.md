# Install and select database profiles

Run the interactive installer to view both profiles, their versions, and their
download sizes:

```bash
pixi run setup
```

Select a profile by number or name. The selection and database root are stored
in `config/local.config` for later runs.

## Install the curated profile

```bash
pixi run setup -- --database_profile curated
```

## Install the IMG profile

```bash
pixi run setup -- --database_profile img
```

The installer downloads the archive listed in `config/database_catalog.json`,
checks Zenodo for the latest database release, and verifies the release and
archive checksums before installation. The bundled catalog is used if Zenodo
cannot be reached during a first installation.

## Update an installed profile

Pipeline runs print the installed database version and check Zenodo for the
latest release. An available update does not interrupt the run.

Install the suggested update:

```bash
pixi run setup -- --database_profile curated --update
```

Without `--update`, interactive setup asks before replacing an installed
profile. The existing profile remains in place until the downloaded archive and
its contents pass validation.

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

The root must contain the selected profile directory and its `manifest.json`.
See [database profiles](../reference/database-profiles.md) for the installed
files and taxonomy sources.
