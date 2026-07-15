# First SSUextract run

This tutorial installs SSUextract, runs one bundled assembly, and checks the
result. It uses the default `curated` database profile.

## Install the pinned environment

```bash
git clone https://github.com/NeLLi-team/ssuextract.git
cd ssuextract
pixi install --frozen
```

Pixi resolves the environment from `pixi.lock`. Confirm the entry point:

```bash
pixi run ssuextract -- --version
```

## Install the reference database

```bash
pixi run setup
```

Press Enter to accept `resources/database`, or enter another writable directory.
The installer verifies the downloaded archive checksum and every file listed in
the profile manifest before the profile becomes available.

## Run the bundled assembly

```bash
pixi run example
```

The run creates `results/smoke/`. Check the two principal reports:

```bash
test -s results/smoke/cmsearch_summary.tsv
test -s results/smoke/cmsearch_summary.tab
```

The first file contains one row per extracted region. The second contains
per-sample counts by annotation category.

## Inspect an extracted sequence

```bash
find results/smoke/extracted -name '*.fna' -size +0 -print
```

Each FASTA record spans the complete 1-based inclusive interval reported by
Infernal. Reverse-strand hits are reverse-complemented after interval extraction.

You have completed the standard SSUextract path. Continue with
[assembled genomes](../how-to/run-assemblies.md) or inspect the
[output reference](../reference/outputs.md).

