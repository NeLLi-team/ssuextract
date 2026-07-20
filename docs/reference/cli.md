# Command-line parameters

Run `pixi run ssuextract --help` for the installed version's generated help.

| Parameter | Default | Description |
| --- | --- | --- |
| `--querydir` | `data/example` | Directory containing input FASTA files. |
| `--modeldir` | `resources/models` | Directory containing Infernal covariance models. |
| `--outdir` | `results/<querydir>` | Output directory. |
| `--database_path` | Configured path or `resources/database` | Root containing database profiles. |
| `--database_profile` | Saved profile or `curated` | Profile name: `curated` or `img`. The option overrides the saved profile for one run. |
| `--model_marker_map` | `config/model_markers.json` | Covariance-model to marker mapping. |
| `--max_blast_targets` | `500` | Candidate policy limit; one extra hit is fetched as a truncation sentinel. |
| `--min_extract_length` | `500` | Minimum accepted hit length in nucleotides; `0` disables the filter. |
| `--threads_per_job` | `2` | CPUs assigned to each `cmsearch` and BLAST task. |
| `--max_cpus` | `16` | Maximum CPUs assigned to one Nextflow task. |
| `--max_memory` | `128.GB` | Maximum memory assigned to one Nextflow task. |
| `--max_time` | `240.h` | Maximum wall time assigned to one Nextflow task. |

## Pixi commands

| Command | Effect |
| --- | --- |
| `pixi run setup` | Select, install, validate, or update a database profile. |
| `pixi run setup --database_profile img` | Install or inspect the IMG profile without the profile prompt. |
| `pixi run setup --database_profile curated --update` | Install the latest verified curated profile without an update prompt. |
| `pixi run example` | Run one bundled assembly. |
| `pixi run ssuextract` | Run the pipeline with supplied Nextflow arguments. |
| `pixi run test` | Run unit, integration, profile-routing, and version checks. |
| `pixi run dryrun` | Preview the Nextflow graph without executing tasks. |
