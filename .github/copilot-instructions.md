# UCI Single-Cell RNA-seq Pipeline - Copilot Instructions

## Repository Overview

This repository contains a comprehensive Snakemake-based pipeline for single-cell RNA-seq analysis, specifically designed for processing 10x Genomics data. The pipeline is optimized for HPC environments using SLURM and includes both CellRanger processing and downstream analytical tools.

**Repository Type:** Bioinformatics workflow pipeline  
**Primary Language:** Snakemake with R and Python scripts  
**Size:** ~50 analysis scripts across 4 main workflow components  
**Target Environment:** SLURM-managed HPC clusters  
**Main Dependencies:** CellRanger 8.0.1, Snakemake, R/Bioconductor, Python

## Build and Validation Instructions

### Prerequisites
**ALWAYS ensure these are available before running any workflows:**

1. **Snakemake Installation:**
   ```bash
   # Install via conda (recommended)
   conda install -c conda-forge snakemake
   # OR install via pip
   pip install snakemake
   ```

2. **CellRanger Module (HPC environments):**
   ```bash
   module load cellranger/8.0.1
   ```

3. **R Dependencies:**
   ```bash
   # Run the package installation script
   Rscript scripts/install_pkgs.R outtxt="install_log.txt" annotation="Gencode" organism="Homo_sapiens" ncores="4"
   ```

### Workflow Validation

**ALWAYS run dry-run validation before executing workflows:**

1. **Single Sample Workflow:**
   ```bash
   # Edit config_cellranger.yaml first to set your sample_id and paths
   snakemake --snakefile Snakefile_cellranger --configfile config_cellranger.yaml --dryrun
   ```

2. **Multi-Sample Workflow:**
   ```bash
   # Edit config_cellranger_multi.yaml first to set sample list and paths
   # NOTE: This workflow includes SCENIC rules that may cause validation errors
   # due to undefined 'outputdir' variable
   snakemake --snakefile Snakefile_cellranger_multi --configfile config_cellranger_multi.yaml --dryrun
   ```

### Running Workflows

**Resource Requirements:** 128GB RAM, 8 cores, up to 48 hours

1. **Local Execution (testing only):**
   ```bash
   snakemake --snakefile Snakefile_cellranger --configfile config_cellranger.yaml --cores 8
   ```

2. **SLURM Cluster Execution (production):**
   ```bash
   # Single sample
   sbatch submit_cellranger.sh
   
   # Multi-sample
   snakemake --snakefile Snakefile_cellranger_multi --configfile config_cellranger_multi.yaml \
   --cluster "sbatch -A SBSANDME_LAB -p standard --mem=128G --time=48:00:00" --jobs 5
   ```

### Common Issues and Solutions

1. **Missing Input Files Error:**
   - Verify FASTQ file paths in configuration files
   - Ensure FASTQ files follow 10x naming convention: `{sample}_S1_L001_R1_001.fastq.gz`
   - Check that `fastqs` path in config points to correct directory

2. **Module Loading Issues:**
   - Always load CellRanger module before execution: `module load cellranger/8.0.1`
   - Unload module after completion: `module unload cellranger/8.0.1`

3. **Memory Issues:**
   - CellRanger requires minimum 128GB RAM - NEVER reduce below this
   - Increase `localmem` parameter in config if analysis fails

4. **Reference Genome Paths:**
   - Human: `/dfs8/commondata/cellranger/refdata-cellranger-GRCh38-3.0.0`
   - Mouse: `/dfs8/commondata/cellranger/refdata-cellranger-mm10-1.2.0`
   - Always verify paths exist before running

5. **Undefined Variables in Rules:**
   - Some rules reference `outputdir` variable that is not defined in current workflows
   - If using SCENIC/Numbat rules, define `outputdir` variable or remove the include statement
   - Multi-sample workflow includes `rules/scenic.smk` which will cause validation errors

## Project Architecture and Layout

### Root Directory Structure
```
├── Snakefile_cellranger          # Single sample CellRanger workflow
├── Snakefile_cellranger_multi    # Multi-sample CellRanger workflow  
├── config_cellranger.yaml        # Single sample configuration
├── config_cellranger_multi.yaml  # Multi-sample configuration
├── README_cellranger.md          # CellRanger workflow documentation
├── rules/                        # Snakemake rule definitions
└── scripts/                      # Analysis scripts (49 files)
```

### Workflow Components

1. **Primary Workflows (Root Level):**
   - `Snakefile_cellranger`: Processes single 10x sample through CellRanger count
   - `Snakefile_cellranger_multi`: Batch processes multiple samples with optional FastQC/MultiQC

2. **Rule Modules (rules/ directory):**
   - `scenic.smk`: Gene regulatory network analysis using SCENIC
   - `numbat.smk`: Copy number variation analysis in single cells
   - `pileup_and_phase.smk`: SNP pileup and phasing for CNV analysis
   - `jbrowse.smk`: Genome browser visualization setup

3. **Analysis Scripts (scripts/ directory - 49 files):**
   - **R scripts (42 files):** Downstream analysis, visualization, and reporting
   - **Python scripts (5 files):** Data processing and format conversion
   - **Shell scripts (2 files):** Utility and phasing operations

### Configuration Requirements

**ALWAYS update these paths before running:**

1. **In config_cellranger.yaml:**
   ```yaml
   sample_id: "YOUR_SAMPLE"           # Change this
   paths:
     fastqs: "PATH_TO_YOUR_FASTQS"    # Change this
     output: "OUTPUT_DIRECTORY"        # Change this
   ```

2. **In config_cellranger_multi.yaml:**
   ```yaml
   samples:                           # Change this list
     - "SAMPLE1"
     - "SAMPLE2"
   paths:
     fastqs: "PATH_TO_YOUR_FASTQS"    # Change this
   ```

### Expected Input/Output Structure

**Input FASTQ Structure:**
```
FastqFiles/
├── {sample}_merged/
│   ├── {sample}_S1_L001_R1_001.fastq.gz  # Cell barcode + UMI
│   ├── {sample}_S1_L001_R2_001.fastq.gz  # RNA sequence
│   └── {sample}_S1_L001_I1_001.fastq.gz  # Sample index
```

**Output Structure:**
```
{output_dir}/{sample}/
├── outs/
│   ├── web_summary.html                    # QC report
│   ├── filtered_feature_bc_matrix.h5       # Gene-cell matrix
│   ├── possorted_genome_bam.bam           # Aligned reads
│   └── molecule_info.h5                   # Per-molecule info
```

### Dependencies Not in Repository

**External Dependencies (must be available):**
- CellRanger software (`cellranger/8.0.1` module)
- Reference genomes in `/dfs8/commondata/cellranger/`
- R/Bioconductor packages (install via `scripts/install_pkgs.R`)
- SLURM cluster environment
- Email system for workflow notifications

**Missing Conda Environments:**
- Rules reference `../envs/environment_R.yaml` and `../envs/environment.yaml`
- These files are NOT in the repository - create them or remove conda directives

### Validation Steps for Changes

1. **Configuration Changes:**
   ```bash
   # Always validate syntax
   snakemake --snakefile Snakefile_cellranger --configfile config_cellranger.yaml --dryrun
   ```

2. **Rule Changes:**
   ```bash
   # Test individual rules
   snakemake --snakefile Snakefile_cellranger --configfile config_cellranger.yaml rulename --dryrun
   ```

3. **Script Changes:**
   ```bash
   # Test R scripts independently
   Rscript scripts/your_script.R arg1="value1" arg2="value2"
   ```

### Key Files for Understanding

**Essential reading order:**
1. `README_cellranger.md` - Workflow overview and usage
2. `config_cellranger.yaml` - Configuration structure
3. `Snakefile_cellranger` - Main workflow logic
4. `scripts/steps_to_single_cell_cnv_analysis.md` - Analysis approach
5. `rules/numbat.smk` - CNV analysis workflow
6. `scripts/pileup_and_phase.R` - Core analysis script

## Important Notes for Coding Agents

**TRUST THESE INSTRUCTIONS:** Only search for additional information if these instructions are incomplete or contain errors. This repository has been thoroughly analyzed.

**Critical Requirements:**
- NEVER reduce memory requirements below 128GB for CellRanger
- ALWAYS load CellRanger module before execution
- ALWAYS validate paths in configuration files before running
- ALWAYS run dry-run validation before executing workflows
- Email notifications require proper SLURM environment setup

**File Path Conventions:**
- Use absolute paths for input data and reference genomes
- Output paths can be relative to working directory
- FASTQ files must follow exact 10x naming conventions

**Workflow Execution Order:**
1. CellRanger processing (primary output: gene-cell matrices)
2. Optional downstream analyses (SCENIC, Numbat, etc.)
3. Visualization and reporting scripts

**Common Modification Points:**
- Sample lists and paths in YAML configuration files
- Resource requirements in Snakemake rules (threads, mem_mb)
- Analysis parameters in R scripts (passed as command-line arguments)
- Module loading commands for different HPC environments