# Run the bundled example

This tutorial uses the default `curated` database profile.

## Install the pinned environment

```bash
git clone https://github.com/NeLLi-team/ssuextract.git
cd ssuextract
pixi install --frozen
```

Check the installed version:

```bash
pixi run ssuextract -- --version
```

## Install the reference database

```bash
pixi run setup
```

Press Enter to accept `resources/database`, or enter another writable directory.
The installer checks the archive and the files listed in its manifest.

## Run the bundled assembly

```bash
pixi run example
```

The run creates `results/smoke/`. Check the two summary files:

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

To use other FASTA files, continue with
[run assembled genomes or metagenomes](../how-to/run-assemblies.md). See the
[output reference](../reference/outputs.md) for every result path.
