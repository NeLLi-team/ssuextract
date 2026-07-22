# Install and select database profiles

Run the interactive installer to view both profiles, their versions, and their
download sizes:

```bash
pixi run setup
```

Select a profile by number or name. If `config/local.config` contains a previous
selection, the prompt labels it as the default; press Enter to keep it. Bare
`pixi run ssuextract` and `pixi run example` commands use this saved profile.

## Install the curated profile

```bash
pixi run setup --database_profile curated
```

## Install the IMG profile

```bash
pixi run setup --database_profile img
```

The installer downloads the archive listed in `config/database_catalog.json`,
checks Zenodo for the latest database release, and verifies the release and
archive checksums before installation. The bundled catalog is used if Zenodo
cannot be reached during a first installation.

The terminal progress bar updates in place and fits the detected terminal
width. Non-interactive logs report the start, completion, and intermediate
progress every 10 seconds.

Setup stops after five consecutive transient download failures. If the retained
bytes fail verification or the server rejects their byte range, setup makes one
clean restart for that condition. Run the same setup command again after a
stopped transfer. The next run requests the remaining bytes and verifies the
complete archive size and SHA-256 digest before extracting it.

## Update an installed profile

Pipeline runs print the installed database version and check Zenodo for the
latest release. An interactive run asks before starting Nextflow when an update
is available:

```text
Database update available for img: v1.0.1 -> v1.0.2.
Update database now? [y/N]:
```

Enter `y` to download, validate, and install the release. Any other response
keeps the installed profile for that run. Non-interactive runs print the update
command and continue with the installed profile.

If the update fails, the installed profile remains in place and Nextflow does
not start. Rerun the command to resume the download or repeat the release check.

Install the suggested update:

```bash
pixi run setup --database_profile curated --update
```

Interactive `pixi run setup` without `--update` uses the same confirmation
prompt.

The existing profile remains in place until the downloaded archive and its
contents pass validation.

## Select a profile for one run

```bash
pixi run ssuextract --query data/my_dataset --outdir results/my_dataset-img --database_profile img
```

This option applies to one run and does not change the profile saved by setup.

## Use another database root

```bash
pixi run ssuextract --database_path /path/to/ssuextract-databases --database_profile curated
```

The root must contain the selected profile directory and its `manifest.json`.
Pipeline runs with an explicit `--database_path` report available releases but
do not offer to replace the custom database root. Use `pixi run setup` with the
same path and `--update` to replace it. The terminal prints this complete setup
command when an update is available.

See [database profiles](../reference/database-profiles.md) for the installed
files and taxonomy sources.
