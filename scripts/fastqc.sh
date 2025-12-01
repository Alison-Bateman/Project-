#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-fastqc-%j.out
#SBATCH --job-name=fastqc
#SBATCH --time=01:00:00
#SBATCH --mem=4G

set -euo pipefail

# Load the OSC module for FastQC
module load fastqc/0.12.1

# Positional arguments: input files
fastq_file=$1

# Output base directory
OUT=/fs/ess/PAS2880/users/bateman139/project/data

# Initial logging
echo "Starting script fastqc.sh"
date
echo "Input FASTQ file:  $fastq_file"
echo "Output dir:         $OUT/fastqc"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$OUT/fastqc"

# Run FastQC
fastqc \
    --threads 4 \
    --outdir "$OUT/fastqc" \
    $fastq_file

# Final logging
echo "Finished fastqc.sh on $fastq_file" 
date
