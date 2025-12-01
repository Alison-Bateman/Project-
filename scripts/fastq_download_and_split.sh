#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-fastqc-%j.out

set -euo pipefail

# Set Var
RUN=$1
 
#Initial Logging
echo "Starting fastq_download&split.sh"
date
echo "NCBI Run:   $RUN"
echo "Output dir: data/fastq"
fasterq-dump --version

# Download the raw run from NCBI 
 # 1) Download the run using the SRA identifier code (Ex: SRR1555696)
prefetch -O data/fastq/SRA $RUN  

# 2) Convert to FASTQ and split mates
echo "Running fasterq-dump on $RUN"
fasterq-dump $RUN --split-files --threads 8 -O data/fastq

echo "$RUN converted fastq file with split mates"
# yields fastq/SRR1555696_1.fastq and fastq/SRR1555696_2.fastq

# 3) Compress files 
gzip data/fastq/$RUN*
# yields fastq/SRR1555696_1.fastq.gz and fastq/SRR1555696_2.fastq.gz
# 4) Validate file integrity
echo "Checking file integrity"
vdb-validate data/fastq/$RUN

echo "script complete"


