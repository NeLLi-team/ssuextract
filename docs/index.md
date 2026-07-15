<div class="ssue-hero" markdown>

# Find and classify 16S and 18S rRNA

SSUextract detects small-subunit ribosomal RNA in assembled contigs, extracts
the complete Infernal hit interval, and assigns taxonomy against a database
routed to the marker type.

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
<strong>60.1% smaller</strong>
<span>curated installed footprint</span>
</div>
<div class="metric-card" markdown>
<strong>0.7%</strong>
<span>pipeline runtime difference</span>
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

`pixi run setup` installs the default `curated` database profile and records its
location in `config/local.config`. The example writes extracted FASTA records,
per-hit tables, taxonomy summaries, and Nextflow execution reports under
`results/smoke/`.

## Choose the page that matches the task

- Follow the [first-run tutorial](tutorials/first-run.md) to learn the complete
  analysis path.
- Use the [how-to guides](how-to/run-assemblies.md) for your assemblies or an
  alternate database profile.
- Consult the [reference](reference/cli.md) for parameters, file schemas, and
  database contracts.
- Read the [explanations](explanation/pipeline.md) for workflow, taxonomy, and
  benchmark rationale.

