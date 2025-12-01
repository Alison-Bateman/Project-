#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --output=gff_convert_%j.out
#SBATCH --error=gff_convert_%j.err
#SBATCH --job-name=GFF_to_GTF

set -euo pipefail
# Positional arguments: 

gff=ref/figshare/HelTub_1.0.gff
outfile=ref/gtf/HelTub_1.0.gtf
mkdir -p ref/gtf
apptainer exec oras://community.wave.seqera.io/library/gffread:0.12.7--b08e770b84a4a126 \
gffread $gff -T -o $outfile