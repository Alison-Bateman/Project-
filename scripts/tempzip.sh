#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --job-name=fastp_trim
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=slurm-fastp-%j.out

set -euo pipefail

apptainer exec \
oras://community.wave.seqera.io/library/pbgzip:2016.08.04--c8e0d8eb135a301c \
  pbgzip -n 8 ref/figshare/*
  gzip -1 ref/gtf/*

echo "zipped"