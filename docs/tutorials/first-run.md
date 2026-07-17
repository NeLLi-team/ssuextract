# Run the bundled example

Select the `curated` database profile for this tutorial.

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

Setup lists both database profiles:

```text
Available database profiles:
  1) curated v1.0.1 (345.4 MiB) - PR2 5.1.1 and SILVA 138.2
  2) img v1.0.1 (828.8 MiB) - Curated profile plus IMG 16S rRNA gene and 18S rRNA gene sequences
Database profile [1, curated]:
```

Press Enter to select `curated`. Press Enter again to accept
`resources/database`, or enter another writable directory. The terminal reports
the download percentage, extraction, file validation, and the installed
database version.

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
