#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=8
#SBATCH --time=30
#SBATCH --output=slurm-trimgalore-%j.out
#SBATCH --mail-type=FAIL
set -euo pipefail
   
# Constants
TRIMGALORE_CONTAINER=oras://community.wave.seqera.io/library/trim-galore:0.6.10--bc38c9238980c80e
   
# Copy the placeholder variables
R1="$1"
R2="$2"
outdir="$3"
   
# Initial logging
echo "# Starting script trimgalore.sh"
date
echo "# Input R1 FASTQ file:      $R1"
echo "# Input R2 FASTQ file:      $R2"
echo "# Output dir:               $outdir"
echo
   
# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs
   
# Run TrimGalore
apptainer exec "$TRIMGALORE_CONTAINER" trim_galore \
    --2colour 20 \
    --cores 8 \
    --paired \
    --fastqc \
    --output_dir "$outdir" \
    "$R1" \
    "$R2"

# Build the output file names:
sample_id=$(basename "$R1" _R1.fastq.gz)
R1_out_init="$outdir"/"$sample_id"_R1_val_1.fq.gz
R2_out_init="$outdir"/"$sample_id"_R2_val_2.fq.gz
R1_out_final="$outdir"/"$sample_id"_R1.fastq.gz
R2_out_final="$outdir"/"$sample_id"_R2.fastq.gz
  
# Rename the output files:
echo
echo "# Renaming the output files:"
mv -v "$R1_out_init" "$R1_out_final"
mv -v "$R2_out_init" "$R2_out_final"
  
# Final logging
echo
echo "# TrimGalore version:"
apptainer exec "$TRIMGALORE_CONTAINER" trim_galore -v
echo "# Successfully finished script trimgalore.sh"
date
