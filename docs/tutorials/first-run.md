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
pixi run ssuextract --version
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
Database profile (1-2) [default: 1 (curated)]:
```

The current default appears in brackets. A previous selection becomes the
default on later setup runs. Press Enter to keep it, or enter `1`, `2`,
`curated`, or `img`. For this tutorial, select `curated`. Press Enter again to
accept `resources/database`, or enter another writable directory.
During the download, one terminal line shows written bytes, transfer rate, and
estimated time remaining. It then reports extraction, file validation, and the
installed database version. If the download stops, rerun `pixi run setup`.
Setup continues the selected database download from the retained partial
archive.

## Run the bundled assemblies

```bash
pixi run example
```

The two assemblies contain 10 accepted loci: nine 16S rRNA gene annotations and
one 18S rRNA gene annotation. The command checks every annotation against the
selected database profile and stops if a marker, reference sequence, taxonomy,
or assignment method differs from the expected result.
The checked annotations are tied to the database version shown by setup. Update
SSUextract when installing a later database release so the example contract and
database remain in sync.

The run creates `results/smoke/`. Check the two summary files:

```bash
test -s results/smoke/cmsearch_summary.tsv
test -s results/smoke/cmsearch_summary.tab
```

The first file contains one row per extracted region. The second contains
per-sample counts by annotation category. To view the sample, marker model,
coordinates, reference source, and taxonomy for each annotation:

```bash
cut -f2,3,5,14,15 results/smoke/cmsearch_summary.tsv
```

## Inspect an extracted sequence

```bash
find results/smoke/extracted -name '*.fna' -size +0 -print
```

Each FASTA record spans the complete 1-based inclusive interval reported by
Infernal. Reverse-strand hits are reverse-complemented after interval extraction.

To use other FASTA files, continue with
[run assembled genomes or metagenomes](../how-to/run-assemblies.md). See the
[output reference](../reference/outputs.md) for every result path.
