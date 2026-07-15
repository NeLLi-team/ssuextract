# Run assembled genomes or metagenomes

Place `.fna`, `.fa`, or `.fasta` assemblies in one directory. The file basename
becomes the sample identifier.

## Run the directory

```bash
pixi run ssuextract -- \
  --querydir data/my_dataset \
  --outdir results/my_dataset \
  --threads_per_job 4
```

SSUextract validates query, model, database, and output paths before scheduling
the workflow. Basenames may contain letters, numbers, `.`, `_`, and `-`.

## Bound local resource use

Set Nextflow task ceilings independently from the per-search thread count:

```bash
pixi run ssuextract -- \
  --querydir data/my_dataset \
  --outdir results/my_dataset \
  --threads_per_job 4 \
  --max_cpus 8 \
  --max_memory 64.GB
```

`--threads_per_job` controls each `cmsearch` and BLAST task. `--max_cpus` limits
the resources that Nextflow assigns to one task; it does not set the number of
concurrent tasks.

## Resume an interrupted run

Run the same command with Nextflow's resume flag:

```bash
pixi run ssuextract -- \
  -resume \
  --querydir data/my_dataset \
  --outdir results/my_dataset
```

Keep the original `work/` directory and command parameters. Nextflow reuses
matching completed tasks and reruns missing work.

