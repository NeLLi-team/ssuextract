# Command-line parameters

Run `pixi run ssuextract --help` for the installed version's generated help.

| Parameter | Default | Description |
| --- | --- | --- |
| `-q`, `--query` | `data/example` | One `.fna`, `.fa`, or `.fasta` file, or a directory containing those files. The short form is available through `pixi run ssuextract`. |
| `--modeldir` | `resources/models` | Directory containing Infernal covariance models. |
| `--outdir` | `results/<input-name>` | Output directory. The input name is the directory name or FASTA file stem. |
| `--database_path` | Configured path or `resources/database` | Root containing database profiles. |
| `--database_profile` | Saved profile or `curated` | Profile name: `curated` or `img`. The option overrides the saved profile for one run. |
| `--model_marker_map` | `config/model_markers.json` | Covariance-model to marker mapping. Marker values must be `16S` or `18S`. |
| `--max_blast_targets` | `500` | Candidate policy limit; one extra hit detects truncation. Must match `max_targets` when the profile has runtime calibration. |
| `--top_hits` | `5` | Number of ranked BLAST subjects reported per query. All assignment candidates, equal-best subjects, and the best fetched IMG, PR2, and SILVA subjects are retained below this cutoff. |
| `--tree_classification` | off | Use a query-plus-reference tree for taxonomy assignment. |
| `--tree_reference_count` | `100` | Unique BLAST reference sequences aligned with each query in tree mode. |
| `--tree_assignment_neighbors` | `5` | Neighbor count for direct-reference tree evidence; includes distance ties and nearby explicit ambiguity. |
| `--tree_trim_gap_fraction` | `0.9` | After masking covariance-model insert columns, remove match columns with a larger gap fraction. |
| `--min_extract_length` | `500` | Minimum union length of CM-hit sequence spans; gaps between joined fragments do not count. `0` disables the filter. |
| `--threads_per_job` | `2` | CPUs assigned to each Infernal, BLAST, or `cmalign` task. IQ-TREE uses one thread for reproducible neighbor ordering. |
| `--max_cpus` | `16` | Maximum CPUs assigned to one Nextflow task. |
| `--max_memory` | `128.GB` | Maximum memory assigned to one Nextflow task. |
| `--max_time` | `240.h` | Maximum wall time assigned to one Nextflow task. |

## Pixi commands

| Command | Effect |
| --- | --- |
| `pixi run setup` | Select, install, validate, or update a database profile. |
| `pixi run setup --database_profile img` | Install or inspect the IMG profile without the profile prompt. |
| `pixi run setup --database_profile curated --update` | Install the latest verified curated profile without an update prompt. |
| `pixi run example` | Run the bundled assemblies. |
| `pixi run ssuextract` | Run the pipeline with supplied Nextflow arguments. |
| `pixi run test` | Run unit, integration, profile-routing, and version checks. |
| `pixi run dryrun` | Preview the Nextflow graph without executing tasks. |
