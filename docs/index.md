<div class="ssue-hero" markdown>

# Find and classify 16S and 18S rRNA

SSUextract detects small-subunit ribosomal RNA in assembled contigs, extracts
the complete hit interval, and assigns taxonomy with a marker-specific reference
database.

[Run the tutorial](tutorials/first-run.md){ .md-button .md-button--primary }
[View the parameters](reference/cli.md){ .md-button }

</div>

<div class="metric-grid" markdown>
<div class="metric-card" markdown>
<strong>16S + 18S</strong>
<span>marker-specific search</span>
</div>
<div class="metric-card" markdown>
<strong>609,298</strong>
<span>unique curated sequences</span>
</div>
<div class="metric-card" markdown>
<strong>55–73 s</strong>
<span>median runtime, bundled example</span>
</div>
<div class="metric-card" markdown>
<strong>2.8 GiB</strong>
<span>median maximum reported RSS</span>
</div>
</div>

## Start with a complete example

```bash
git clone https://github.com/NeLLi-team/ssuextract.git
cd ssuextract
pixi install --frozen
pixi run setup
pixi run example
```

`pixi run setup` installs the default `curated` database profile. The example
writes extracted sequences, per-hit annotations, category counts, and Nextflow
reports under `results/smoke/`.

## What the pipeline does

1. Infernal searches each assembly with 16S and 18S covariance models.
2. SSUextract keeps accepted coordinates and extracts the complete interval on
   the reported strand.
3. Each sequence is searched against the 16S or 18S index selected by its model.
4. The workflow writes extracted FASTA files, one row per hit, category counts,
   and execution reports.

[See the pipeline diagram](explanation/pipeline.md) or
[run your assemblies](how-to/run-assemblies.md).
