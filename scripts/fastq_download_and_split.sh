#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm-fastqc-%j.out

set -euo pipefail

# Positional arguments: 
RUNS_FILE=$1
OUT_DIR=$2

# Loop
RUNS=$(cat $RUNS_FILE)

for RUN in $RUNS; do
 
#Initial Logging
echo "Starting fastq_download&split.sh"
date
echo "NCBI Run:   $RUN"
echo "Output dir: $OUT_DIR"
fasterq-dump --version

mkdir -p $OUT_DIR

# Download the raw run from NCBI 
 # 1) Download the run using the SRA identifier code (Ex: SRR1555696)
prefetch -O $OUT_DIR $RUN  

# 2) Convert to FASTQ and split mates
echo "Running fasterq-dump on $RUN"
fasterq-dump $RUN --split-files --threads 16 -O $OUT_DIR

echo "$RUN converted fastq file with split mates"
# yields fastq/SRR1555696_1.fastq and fastq/SRR1555696_2.fastq

echo "Downloaded and split $RUN"

done 

echo 'compressing files'
# 3) Compress files 
gzip fastq/*.fastq
# yields fastq/SRR1555696_1.fastq.gz and fastq/SRR1555696_2.fastq.gz
echo "script complete"