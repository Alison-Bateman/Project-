#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --job-name=fastp_trim
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=slurm-fastp-%j.out

set -euo pipefail

# Set date so all outputs from this run are consistent (MMDD format)
echo "Starting script"

run_date=$(date +%m%d)

# Define output directories
TRIM_DIR=data/trimmed
FASTP_REPORT_DIR=data/trimmed/fastp_reports

# Create directories automatically
mkdir -p "$TRIM_DIR" "$FASTP_REPORT_DIR"

for R1 in data/fastq/*_1.fastq.gz; do
    # Define paired-end R2 read
    R2=${R1/_1.fastq.gz/_2.fastq.gz}

    # Extract sample name for file naming
    # Example: data/fastq/SRR1555696_1.fastq.gz → SRR1555696
    sample=${R1##*/}
    sample=${sample%_1.fastq.gz}

    # Run fastp inside Apptainer
    # -i/-I: input reads
    # -o/-O: trimmed outputs (date-stamped)
    # --cut_tail* : 3' quality trimming with 4-bp window at Q20
    # -n 5: discard reads with >5 Ns
    # -l 50: discard reads shorter than 50 bp after trimming
    # --detect_adapter_for_pe: safe adapter detection
    # -h/-j: HTML/JSON reports (date-stamped)
    # --thread 8: use 8 threads
    apptainer exec \
      oras://community.wave.seqera.io/library/fastp:1.0.1--a5a7772c43b5ebcb \
      fastp \
      -i "$R1" \
      -I "$R2" \
      -o "${TRIM_DIR}/${sample}_${run_date}_1.trimmed.fastq.gz" \
      -O "${TRIM_DIR}/${sample}_${run_date}_2.trimmed.fastq.gz" \
      --cut_tail \
      --cut_tail_window_size 4 \
      --cut_tail_mean_quality 20 \
      -n 5 \
      -l 50 \
      --detect_adapter_for_pe \
      -h "${FASTP_REPORT_DIR}/${sample}_${run_date}.html" \
      -j "${FASTP_REPORT_DIR}/${sample}_${run_date}.json" \
      --thread 8

done
echo "Completed trimming"
# Run FastQC on trimmed FASTQ files
echo "starting fastqc"
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