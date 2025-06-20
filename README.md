# 🧬 SSUextract: Small Subunit rRNA Extraction Pipeline

[![GitHub release](https://img.shields.io/github/v/release/NeLLi-team/ssuextract?style=flat-square&color=green)](https://github.com/NeLLi-team/ssuextract/releases/latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/nextflow-≥22.10.1-brightgreen.svg)](https://nextflow.io)
[![pixi](https://img.shields.io/badge/pixi-supported-blue.svg)](https://pixi.sh)

A high-performance bioinformatics pipeline for extracting and annotating Small Subunit (SSU) rRNA sequences from genomic assemblies using covariance models.

> **🚀 Nextflow Pipeline**: Complete conversion from Snakemake to Nextflow for improved scalability and cloud compatibility! This is now the main branch (v1.0.0).

## 🎯 Overview

SSUextract identifies SSU rRNA sequences in your genomic data by:
- 🔍 Searching with Infernal's cmsearch using curated covariance models
- 🧮 Extracting high-quality SSU sequences based on configurable parameters
- 🏷️ Annotating sequences with taxonomic information via BLAST
- 📊 Generating comprehensive summary reports

## 🚀 Quick Start

### Prerequisites
- [pixi](https://pixi.sh) package manager
```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

### Installation

```bash
# Clone the repository
git clone https://github.com/NeLLi-team/ssuextract.git
cd ssuextract

# Install dependencies and download reference database
pixi run setup
```

### Basic Usage

```bash
# Run on example data
pixi run run

# Results will be in: results/example/
```

## 🔧 Pipeline Architecture

```mermaid
graph TD
    A[🗂️ Input FASTA files] --> B[📝 Header validation]
    B --> C[🔍 cmsearch<br/>RF00177 & RF01960]
    C --> D[📊 Extract coordinates]
    D --> E[✂️ Extract sequences]
    E --> F[🧬 BLAST annotation]
    F --> G[📈 Merge results]
    G --> H[📋 Summary table]
    
    style A fill:#e1f5fe
    style H fill:#c8e6c9
    style C fill:#fff3e0
    style F fill:#fff3e0
```

## 📁 Project Structure

```
ssuextract/
├── 📄 main.nf               # Nextflow pipeline definition
├── 📄 nextflow.config       # Nextflow configuration
├── 📦 pixi.toml             # Dependencies and tasks
├── 📂 config/               # Configuration files
│   └── default.yaml         # Pipeline configuration
├── 📂 scripts/              # Pipeline scripts
│   ├── cmprocessing.py      # BLAST result processing
│   ├── get_cmsequences.py   # Sequence extraction
│   ├── get_cmstats.py       # Alignment statistics
│   ├── get_table.py         # Summary table generation
│   └── rename_fnaheaders.py # Header validation
├── 📂 resources/            # Static resources
│   ├── models/              # Covariance models
│   │   ├── RF00177.cm      # Bacterial/Archaeal SSU
│   │   └── RF01960.cm      # Eukaryotic SSU
│   └── database/            # Reference database
├── 📂 data/                 # Input data
│   └── example/             # Example test data
└── 📊 results/              # Output directory
    └── {dataset}/           # Dataset-specific results
        ├── fna/             # Processed sequences
        ├── out/             # cmsearch outputs
        ├── stats/           # Alignment statistics
        ├── extracted/       # Extracted SSU sequences
        ├── m8/              # BLAST results
        └── *.tsv            # Summary tables
```

## ⚙️ Configuration

### Basic Configuration

Pipeline parameters can be set in `nextflow.config` or via command line:

```bash
# Run with custom parameters
nextflow run main.nf --querydir data/my_dataset --threads_per_job 4 --min_extract_length 50
```

### Custom Data

```bash
# Place your .fna files in data/your_dataset/
mkdir data/your_dataset
cp /path/to/*.fna data/your_dataset/

# Run pipeline on custom data
pixi run nextflow run main.nf --querydir data/your_dataset --threads_per_job 4
```

### Configuration Files

- **`nextflow.config`**: Main pipeline configuration
- **`config/base.config`**: Process-specific resource requirements
- **`config/environment.yml`**: Conda environment specification

Results will be in `results/your_dataset/`

## 📊 Output Files

### Main Output: `results/{dataset}/cmsearch_summary.tsv`

A comprehensive table containing:
- **name**: Sequence identifier
- **sample**: Source sample name (basename of input file)
- **model**: CM model used (RF00177/RF01960)
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

### Category Summary: `results/{dataset}/cmsearch_summary.tab`

Counts of SSU types per sample:
- BacteriaSSU
- ArchaeaSSU
- EukaryotaSSU
- MitochondriaSSU
- PlastidSSU
- And more specialized categories

## 🛠️ Available Commands

```bash
# Setup and installation
pixi run setup              # Install deps + download database

# Pipeline execution
pixi run run                # Run full pipeline (default config)
pixi run dryrun            # Preview what will be executed

# Custom execution  
pixi run nextflow run main.nf --querydir data/my_dataset

# Database management
pixi run download-db        # Download reference database

# Cleanup
pixi run clean             # Clean pipeline outputs
pixi run clean-results     # Remove all results directories

# Development
pixi run -e dev run-verbose # Verbose output
pixi run -e dev report      # Generate HTML report
pixi run -e dev dag         # Generate pipeline visualization
```

## 🔬 Pipeline Details

### Step 1: Input Validation
- Validates FASTA headers with `scripts/rename_fnaheaders.py`
- Ensures compatible formatting

### Step 2: Covariance Model Search
- Uses Infernal's cmsearch with `--anytrunc` flag
- Searches for both bacterial/archaeal (RF00177) and eukaryotic (RF01960) SSU

### Step 3: Hit Processing
- Extracts alignment coordinates with `scripts/get_cmstats.py`
- Handles truncated alignments
- Filters by minimum length (configurable: default 30bp)

### Step 4: Sequence Extraction
- Extracts sequences based on coordinates with `scripts/get_cmsequences.py`
- Handles reverse complement for minus strand hits

### Step 5: Taxonomic Annotation
- BLAST search against SILVA/PR2 database
- Processes results with `scripts/cmprocessing.py`

### Step 6: Report Generation
- Merges all results with `scripts/get_table.py`
- Creates summary statistics
- Generates final TSV report

## 🐛 Troubleshooting

### Empty BLAST results
- Check minimum length setting in config/default.yaml
- Verify extracted sequences meet length threshold
- Ensure database was downloaded correctly

### Pipeline locked error
```bash
pixi run nextflow clean -f
```

### Memory issues
- Reduce `threads_per_job` parameter
- Use profile configurations: `nextflow run main.nf -profile local`

### File path errors
- Ensure input files have `.fna`, `.fa`, or `.fasta` extensions
- Check that `querydir` path exists and contains sequence files

## 📝 Changelog

### v1.0.0 (Latest) - 2024-06-06
- 🚀 **Major Release: Migration to Nextflow** - Complete pipeline rewrite using Nextflow
- 🏗️ **Enhanced workflow structure** - Improved modularity and maintainability
- 📦 **Updated covariance models** - Added RF00177.cm and RF01960.cm resources
- 📚 **Comprehensive documentation** - Updated README with Nextflow-specific instructions
- 🛠️ **Improved script functionality** - Enhanced get_cmstats.py and other processing scripts
- 📁 **Repository restructure** - Legacy Snakemake version preserved in ssuextract-snk branch
- [View full release notes →](https://github.com/NeLLi-team/ssuextract/releases/tag/v1.0.0)

### Previous Versions
- **v0.9.0** - Last Snakemake version (now in ssuextract-snk branch)
- [View all releases →](https://github.com/NeLLi-team/ssuextract/releases)

## 📝 Citation

If you use SSUextract in your research, please cite:

```
SSUextract: A Nextflow pipeline for SSU rRNA extraction
[Your publication details here]
```

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- [Infernal](http://eddylab.org/infernal/) for covariance model searches
- [SILVA](https://www.arb-silva.de/) and [PR2](https://pr2-database.org/) databases
- [Nextflow](https://nextflow.io) workflow engine
- [pixi](https://pixi.sh) package manager

---

<p align="center">Made with ❤️ for the bioinformatics community</p>
