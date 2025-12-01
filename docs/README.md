# RNA-Seq analysis of _Culex pipiens_

- Author: <???>
- Date: <???>
- GitHub repository: <???>
- This was run at the <???> cluster of the Ohio Supercomputer Center (<www.osc.edu>)
- OSC working dir: <???>

## Prerequisites

- **Input files**:
  - The input files are:
    - Paired-end FASTQ files
    - Reference genomes files (assembly and annotation)

- **Software**:
  - The scripts will run programs using Singularity/Apptainer containers
    available online, with links contained inside the scripts
  - If running at OSC, no software installations are needed
  - If you're not at OSC, you need to have `apptainer` and `Slurm` installed.

- You don't need to create the output dirs because the scripts will do so.

## Set up: Get the input files

1. Download the input FASTQ files and the reference genome annotation GTF
   from a GitHub repo as follows:

   ```bash
   git clone https://github.com/jelmerp/garrigos-data
   ```

2. Download the reference genome assembly FASTQ from NCBI as follows:

   ```bash
   URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/801/865/GCF_016801865.2_TS_CPP_V2/GCF_016801865.2_TS_CPP_V2_genomic.fna.gz
   wget -P garrigos-data/ref "$URL"
   ```

3. Unzip the reference genome files:

   ```bash
   gunzip -v garrigos-data/ref/*
   ```

## Set up: Set file paths

```bash
# Inputs:
fastq_dir=<???>
ref_assembly=<???>
ref_annotation=<???>

# Base output dir:
outdir=results
```

## Step 1: Read QC with FastQC

The `fastqc.sh` script runs on 1 FASTQ file at a time and takes 2 arguments -
the input FASTQ file and output dir:

```bash
for fastq in "$fastq_dir"/*fastq.gz; do
    sbatch scripts/fastqc.sh <???> "$outdir"/fastqc
done
```

## Step 2: FastQC summaries with MultiQC

<???>

## Step 3: Read trimming with TrimGalore

<???>

## Step 4: Creating a reference genome index with STAR

<???>

## Step 5: Trimmed read alignment to the reference genome index with STAR

<???>