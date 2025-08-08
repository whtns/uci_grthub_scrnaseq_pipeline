# CellRanger Snakemake Workflow

This Snakemake workflow has been converted from the original SLURM batch script for single-cell RNA-seq analysis using CellRanger.

## Files Created

### Single Sample Workflow:
- `Snakefile_cellranger` - Main workflow for single sample
- `config_cellranger.yaml` - Configuration file for single sample
- `cluster_cellranger.yaml` - SLURM cluster configuration
- `submit_cellranger.sh` - Script to submit single sample workflow

### Multi-Sample Workflow:
- `Snakefile_cellranger_multi` - Workflow for multiple samples
- `config_cellranger_multi.yaml` - Configuration for multiple samples

## Workflow Overview

The workflow performs the following steps:
1. **CellRanger Count**: Processes 10x Genomics single-cell RNA-seq data
   - Aligns reads to reference transcriptome
   - Generates gene-barcode matrices
   - Creates quality control metrics
   - Produces BAM files (optional)

## Prerequisites

1. **CellRanger**: Available as module `cellranger/8.0.1`
2. **Snakemake**: Install via conda
3. **Reference Genomes**: Available at `/dfs8/commondata/cellranger/`
   - Human: `refdata-cellranger-GRCh38-3.0.0`
   - Mouse: `refdata-cellranger-mm10-1.2.0`

## Usage

### Single Sample Analysis

1. **Edit configuration**:
   ```yaml
   # config_cellranger.yaml
   sample_id: "A172"  # Change to your sample ID
   ```

2. **Submit to SLURM**:
   ```bash
   sbatch submit_cellranger.sh
   ```

3. **Run locally** (for testing):
   ```bash
   snakemake --snakefile Snakefile_cellranger --configfile config_cellranger.yaml --cores 8
   ```

### Multi-Sample Analysis

1. **Edit configuration**:
   ```yaml
   # config_cellranger_multi.yaml
   samples:
     - "A172"
     - "A173"
     - "A174"  # Add your sample IDs
   ```

2. **Submit to SLURM**:
   ```bash
   snakemake --snakefile Snakefile_cellranger_multi --configfile config_cellranger_multi.yaml \
   --cluster "sbatch -A SBSANDME_LAB -p standard --mem=128G --time=48:00:00" --jobs 5
   ```

## Input Requirements

### FASTQ File Structure
CellRanger expects FASTQ files in a specific format:
```
FastqFiles/
├── A172_merged/
│   ├── A172_S1_L001_R1_001.fastq.gz
│   ├── A172_S1_L001_R2_001.fastq.gz
│   └── A172_S1_L001_I1_001.fastq.gz
```

### File Naming Convention
- R1: Read 1 (contains cell barcode and UMI)
- R2: Read 2 (contains actual RNA sequence)
- I1: Index read (sample index)

## Output Files

For each sample, CellRanger generates:

### Main Outputs:
- `web_summary.html` - Quality control report
- `filtered_feature_bc_matrix.h5` - Filtered gene-barcode matrix
- `possorted_genome_bam.bam` - Aligned reads (if `--create-bam=true`)
- `molecule_info.h5` - Per-molecule information

### Directory Structure:
```
{sample}/
├── outs/
│   ├── web_summary.html
│   ├── filtered_feature_bc_matrix.h5
│   ├── filtered_feature_bc_matrix/
│   ├── raw_feature_bc_matrix.h5
│   ├── raw_feature_bc_matrix/
│   ├── possorted_genome_bam.bam
│   ├── possorted_genome_bam.bam.bai
│   └── molecule_info.h5
```

## Resource Requirements

- **Memory**: 128GB (CellRanger requirement)
- **CPU**: 8 cores
- **Time**: Up to 48 hours for large datasets
- **Storage**: ~10-50GB per sample (depending on data size)

## Customization

### For Different Species:
Change the transcriptome reference in `config_cellranger.yaml`:
```yaml
references:
  transcriptome: "/dfs8/commondata/cellranger/refdata-cellranger-mm10-1.2.0"  # For mouse
```

### For Different Parameters:
Modify CellRanger parameters in the configuration:
```yaml
params:
  localcores: 16    # Increase for faster processing
  localmem: 256     # Increase for large datasets
```

### Additional CellRanger Options:
Add to the shell command in the Snakefile:
```bash
cellranger count --id={params.sample_id} \
--transcriptome={params.transcriptome} \
--fastqs={input.fastq_dir} \
--expect-cells=10000 \        # Expected cell count
--force-cells=5000 \          # Force cell count
--include-introns             # Include intronic reads
```

## Troubleshooting

1. **Memory Issues**: Increase `localmem` parameter
2. **FASTQ Format**: Ensure proper 10x FASTQ naming convention
3. **Reference Path**: Verify transcriptome path exists
4. **Disk Space**: Ensure sufficient space for output files

## Monitoring

- Check SLURM job status: `squeue -u $USER`
- View CellRanger progress: Check `*.out` files
- Monitor resource usage: `sacct -j JOBID --format=JobID,MaxRSS,Elapsed`
