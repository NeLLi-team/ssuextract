# SSUextract

SSUextract is a Nextflow pipeline that detects 16S and 18S ribosomal RNA in
assembled contigs, extracts complete hit intervals, and assigns taxonomy with
marker-specific reference databases.

## Quick start

```bash
git clone https://github.com/NeLLi-team/ssuextract.git
cd ssuextract
pixi install --frozen
pixi run setup
pixi run example
```

The [SSUextract documentation](https://nelli-team.github.io/ssuextract/) contains
the tutorial, database instructions, parameter reference, output schemas,
taxonomy policy, and measured benchmarks.

Run `pixi run ssuextract -- --help` for the command-line summary. SSUextract is
distributed under the [MIT license](LICENSE).
