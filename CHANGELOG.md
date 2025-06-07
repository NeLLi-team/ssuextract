# Changelog

All notable changes to this project will be documented in this file.

## [2.0.0] - 2024-06-06

### 🚀 Major Release: Migration to Nextflow

This is a major release that transitions the ssuextract pipeline from Snakemake to Nextflow.

#### ✨ New Features
- **Complete pipeline rewrite using Nextflow** - Modern, scalable workflow management
- **Enhanced workflow structure** - Improved modularity and maintainability
- **Updated covariance models** - Added RF00177.cm and RF01960.cm resources
- **Comprehensive documentation** - Updated README with Nextflow-specific instructions
- **Improved script functionality** - Enhanced get_cmstats.py and other processing scripts

#### 🔄 Breaking Changes
- **Pipeline engine changed from Snakemake to Nextflow**
- **Configuration format updated** - Now uses Nextflow config syntax
- **Command line interface changed** - Use `nextflow run main.nf` instead of `snakemake`
- **Dependencies updated** - See pixi.toml for current requirements

#### 📁 Repository Structure Changes
- **Main branch** - Now contains the Nextflow implementation
- **ssuextract-snk branch** - Contains the legacy Snakemake implementation for reference

#### 🛠️ Migration Guide
For users migrating from the Snakemake version:
1. Update your environment using the new pixi.toml
2. Replace Snakemake commands with Nextflow equivalents
3. Update configuration files to use Nextflow syntax
4. Refer to the updated README.md for detailed usage instructions

#### 📋 Legacy Support
The previous Snakemake implementation is preserved in the `ssuextract-snk` branch and remains available for users who need to continue using the Snakemake version.

---

## [1.x.x] - Previous Versions

Previous versions used Snakemake as the workflow engine. For historical changes, please refer to the `ssuextract-snk` branch.
