# 📁 SSUextract Project Structure

This document describes the modern, state-of-the-art directory structure used in SSUextract.

## Directory Layout

```
ssuextract/
├── 📋 Snakefile                 # Main workflow definition
├── 📦 pixi.toml                 # Dependencies and task runner
├── 📄 LICENSE                   # Project license
├── 📄 PROJECT_STRUCTURE.md      # This file
│
├── 📁 config/                   # Configuration files
│   ├── default.yaml            # Default pipeline configuration
│   └── environment.yml         # Conda environment specification
│
├── 📁 scripts/                  # Pipeline scripts
│   ├── cmprocessing.py         # BLAST result processing
│   ├── get_cmsequences.py      # Sequence extraction
│   ├── get_cmstats.py          # Alignment statistics
│   ├── get_table.py            # Summary table generation
│   ├── rename_fnaheaders.py    # Header validation
│   └── cmsearchout_extract_by_position_size.py  # Legacy extraction
│
├── 📁 resources/               # Static resources
│   ├── models/                 # Covariance models
│   │   ├── RF00177.cm         # Bacterial/Archaeal SSU rRNA
│   │   └── RF01960.cm         # Eukaryotic SSU rRNA
│   └── database/              # Reference databases
│       ├── silva-138-1_pr2-4-12.fasta
│       ├── silva-138-1_pr2-4-12.nhr
│       ├── silva-138-1_pr2-4-12.nin
│       └── silva-138-1_pr2-4-12.nsq
│
├── 📁 data/                    # Input data
│   └── example/               # Example test data
│       ├── LKH462_P08_Rh.fna
│       └── LKH565_P11_Ci.fna
│
├── 📁 results/                 # Analysis results (auto-generated)
│   └── results_{dataset}/     # Dataset-specific results
│       ├── fna/              # Processed sequences
│       ├── out/              # cmsearch raw output
│       ├── stats/            # Alignment statistics
│       ├── extracted/        # Extracted SSU sequences
│       ├── m8/               # BLAST results
│       ├── cmsearch_summary.tab   # Category summary
│       └── cmsearch_summary.tsv   # Detailed results
│
├── 📁 docs/                   # Documentation
│   ├── README.md             # Main documentation
│   └── TUTORIAL.md           # Step-by-step tutorial
│
└── 📁 tests/                  # Unit tests (future)
    └── test_pipeline.py      # Pipeline tests
```

## Design Principles

### 🎯 **Separation of Concerns**
- **`config/`**: All configuration in one place
- **`scripts/`**: Executable code separate from workflow
- **`resources/`**: Static files that don't change
- **`data/`**: Input data separate from code
- **`results/`**: Output separate from everything else

### 📦 **Modern Package Management**
- **`pixi.toml`**: Reproducible dependency management
- **`config/environment.yml`**: Conda environment backup
- **Cross-platform compatibility** ensured

### 🔧 **Workflow Organization**
- **`Snakefile`**: Clean workflow definition
- **Modular scripts**: Each script has single responsibility
- **Configurable paths**: No hardcoded directories

### 📊 **Data Management**
- **Input/Output separation**: Clean data flow
- **Results versioning**: Each dataset gets own directory
- **Intermediate files**: Organized by processing step

## Configuration Files

### `config/default.yaml`
```yaml
# Directory containing covariance models (.cm files)
modeldir: "resources/models"

# Directory containing query sequences (.fna files)  
querydir: "data/example"

# Number of threads per job
threads_per_job: 2

# Minimum sequence length for extraction (bp)
min_extract_length: 30
```

### `config/environment.yml`
- Conda environment specification
- Bioinformatics tools (Infernal, BLAST)
- Python dependencies

## Script Organization

### Core Processing Scripts

#### `scripts/get_cmstats.py`
- **Purpose**: Extract alignment coordinates from cmsearch output
- **Input**: Raw cmsearch .out files
- **Output**: .seqmap files with coordinates

#### `scripts/get_cmsequences.py`
- **Purpose**: Extract sequences based on coordinates
- **Input**: .seqmap files + FASTA sequences
- **Output**: Extracted SSU sequences

#### `scripts/cmprocessing.py`
- **Purpose**: Process BLAST results into categories
- **Input**: BLAST m8 files
- **Output**: Category count table

#### `scripts/get_table.py`
- **Purpose**: Generate comprehensive summary table
- **Input**: All intermediate files
- **Output**: Final TSV with annotations

#### `scripts/rename_fnaheaders.py`
- **Purpose**: Validate and clean FASTA headers
- **Input**: Raw FASTA files
- **Output**: Cleaned FASTA files

## Output Structure

### Main Results Files

#### `cmsearch_summary.tsv`
Comprehensive table with columns:
- **name**: Full sequence identifier
- **sample**: Input file basename (sampleid)
- **model**: Covariance model used
- **length**: Sequence length
- **coordinates**: Genomic coordinates
- **strand**: DNA strand (+/-)
- **sequence_type**: Hit type (simple/complex)
- **contig_name**: Source contig
- **blast_sseqid**: Best BLAST hit with taxonomy
- **blast_pident**: Percent identity
- **blast_length**: Alignment length
- **blast_bitscore**: BLAST bit score
- **is_assembled**: Assembly status

#### `cmsearch_summary.tab`
Category count matrix:
```
               BacteriaSSU  ArchaeaSSU  EukaryotaSSU
Sample1        5            1           0
Sample2        3            0           2
```

## Usage Examples

### Basic Usage
```bash
# Run with default configuration
pixi run run

# Use custom config
pixi run snakemake --configfile config/my_config.yaml --cores 4
```

### Custom Data
```bash
# Create new dataset directory
mkdir data/my_genomes
cp /path/to/*.fna data/my_genomes/

# Create custom config
cat > config/my_config.yaml << EOF
modeldir: "resources/models"
querydir: "data/my_genomes"
threads_per_job: 8
min_extract_length: 50
EOF

# Run pipeline
pixi run snakemake --configfile config/my_config.yaml --cores 16
```

### Results will be in: `results_my_genomes/`

## Benefits of This Structure

### ✅ **Maintainability**
- Clear separation of code, config, data, and results
- Easy to find and modify components
- Follows bioinformatics best practices

### ✅ **Reproducibility**
- All dependencies in pixi.toml
- Configuration versioned with code
- Clear data provenance

### ✅ **Scalability**
- Easy to add new scripts or workflows
- Modular design supports extensions
- Clean interfaces between components

### ✅ **Usability**
- Intuitive directory names
- Self-documenting structure
- Easy onboarding for new users

---

This structure follows modern software engineering practices while being tailored for bioinformatics workflows and computational reproducibility.