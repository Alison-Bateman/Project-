#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-multiqc-%j.out
set -euo pipefail

# Constants
MULTIQC_CONTAINER=oras://community.wave.seqera.io/library/multiqc:1.31--e111a2e953475c51

# Positional arguments: input files
fastqc_file=$1

# Output base directory
OUT=/fs/ess/PAS2880/users/bateman139/project/data

nopath1=${fastqc_file##*/} 

# Initial logging
echo "# Starting script multiqc.sh"
date
echo "# Input dir:                      $fastqc_file"
echo "# Output file:                    $OUT/multiqc"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$OUT/multiqc"

# Run MultiQC
apptainer exec "$MULTIQC_CONTAINER" multiqc \
    --outdir "$OUT/multiqc/$nopath1" \
    "$fastqc_file"

# Final logging
echo
echo "# Successfully finished script multiqc.sh with "$fastqc_file"
date


