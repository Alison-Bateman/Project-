#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --output=gff_convert_%j.out
#SBATCH --error=gff_convert_%j.err
#SBATCH --job-name=GFF_to_GTF

set -euo pipefail

# Positional arguments: 
GFF_FILE=$1
OUT_GTF=$2

# Make directory
mkdir -p ref/gtf

# initial Logging
echo 
echo "Converting $GFF_FILE to GTF file at $OUT_GTF"

# Run Command 
apptainer exec oras://community.wave.seqera.io/library/gffread:0.12.7--b08e770b84a4a126 \
gffread $GFF_FILE -T -o $OUT_GTF

#Final logging
echo "Completed File Conversion"