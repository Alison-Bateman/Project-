#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --job-name=fastp_trim
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=slurm-fastp-%j.out

set -euo pipefail

echo "starting fastqc"
run_date=$(date +%m%d)
FASTQC_DIR=data/trimmed/fastqc

mkdir -p "$FASTQC_DIR" 

module load fastqc/0.12.1 
fastqc \
    --threads 8 \
    --outdir "$FASTQC_DIR" \
    data/trimmed/*"${run_date}"*.trimmed.fastq.gz
echo 'Finshed Fastqc'
#Run MultiQC on FastQC files
echo 'Starting MultiQC'
MULTIQC_DIR=data/trimmed/multiqc

mkdir -p "$MULTIQC_DIR"

apptainer exec \
  oras://community.wave.seqera.io/library/multiqc:1.31--e111a2e953475c51 \
  multiqc \
  data/trimmed \
  --outdir "$MULTIQC_DIR"
echo 'Finsihed MutliQC'
#Final Logging 
echo 'Completed Fastp script'