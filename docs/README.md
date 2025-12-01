---
title: "Protocol"
author: "Alison Bateman"
date: today 
date-fromat: long 
format: 
  html: 
    theme: cosmo 
    toc: true
---

# Project Outline

I'll be going through the RNA seq analysis workflow using FASTQ raw reads from Helianthus tuberosus (sunchoke) from a study in 2014, and a sunchoke reference genome from last year. Because I need to do RNA sequencing in the future with sunchoke, having scripts prepared for me to feed files into for analysis will help streamline my research.\
First, I will go through the workflow with practice data, to see if my data aligns with the data already avaliable to troubleshoots. This helps remove the "I think this might be right" with "I know this is right" before using the same or similar steps with the Helianthus data.

# **Preliminary steps:**

Download the following toolkits:

``` bash
{bash}
#Download Datasets command-line tools
echo "Downloading NCBI Datasets Command-line tools"
#Installed using guide at https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/
#Download datasets: 
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
#Download dataformat: 
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
#Make them executable: 
chmod +x datasets dataformat
echo "Command-line tool download complete"

# Download a prebuilt SRA Toolkit
echo "Starting Toolkit Download"
curl -L -o sratoolkit.tar.gz \
  https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-ubuntu64.tar.gz
# Unpack the archive
tar -xzf sratoolkit.tar.gz
tooldir=$(find . -maxdepth 1 -type d -name "sratoolkit.*")
export PATH="$PWD/$tooldir/bin:$PATH"

#To prevent having to reexport everytime use: 
export PATH="/fs/scratch/PAS2880/sratoolkit.3.2.1-ubuntu64/bin:$PATH"
echo "Successfully downloaded toolkit"
```

Download files, transfer them into vscode:

#### FASTQ files:

My FASTQ file set is from Jung et al., 2014 (Bioproject PRJNA258432). I will be downloading the RNA seq from the leaf sample first to troubleshoot through the steps and scripts, and then I will run the next 4 samples from the same study.

The files are interleaved, so they download as a single file containing both paired ends. We will be downloading the splitting the files in one script, the `fastq_download&split.sh` script.

The script is as follows:

``` bash
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
vdb-validate fastq/$RUN

# 5) Remove SRA files 
rm -r data/fastq/SRA 

echo "script complete"
```

For more than one run, input the list of SRA runs from NBCI into runs.txt, and input the following loop:

```         
RUNS=$(cat runs.txt)

for RUN in $RUNS; do
    echo "Submitting for $RUN"
    sbatch scripts/fastq_download_and_split.sh "$RUN"
done
```

to check on the job progress use:

```{bash}
squeue -u $USER -l
```

#### Reference Genome:

Because the reference genome and annotation data were not avaliable on NCBI, the FASTA (.fna) and GFF files were downloaded from figshare (https://figshare.com/articles/dataset/Annotated_reference_genome_of_Helianthus_tuberosus/22491205) using the following commands:

FASTA file:

```{bash}
wget \
  --content-disposition \
  --trust-server-names \
  --user-agent="Mozilla/5.0" \
  https://figshare.com/ndownloader/files/42806323 \
  -O ref/figshare/HelTub_1.0.fn.gz
ls -lh ref/figshare/HelTub_1.0.fn.gz
```

Unzip with:

```{bash}
gunzip -c ref/figshare/HelTub_1.0.fn.gz > ref/figshare/HelTub_1.0.fn
```

GFF file:

```{bash}
wget \
  --content-disposition \
  --trust-server-names \
  --user-agent="Mozilla/5.0" \
  https://figshare.com/ndownloader/files/42771826 \
  -O ref/figshare/HelTub_1.0.gff.gz
ls -lh ref/figshare/HelTub_1.0.gff.gz

```

Unzip with:

```{bash}
gunzip -c ref/figshare/HelTub_1.0.gff.gz > ref/figshare/HelTub_1.0.gff
```

In order to use the annotation data for alignment, it needs to be in a GTF file format. To convert the GFF to GTF file, use the following command:

(Learned of GFFreads from https://www.biostars.org/p/45791/#45800)

``` bash
{bash} 
gff=ref/figshare/HelTub_1.0.gff
outfile=ref/gtf/HelTub_1.0.gtf
mkdir -p ref/gtf
apptainer exec oras://community.wave.seqera.io/library/gffread:0.12.7--b08e770b84a4a126 \
gffread $gff -T -o $outfile
```

Validate with:

```{bash}
head ref/gtf/HelTub_1.0.gtf
```

Output:

```{bash}
Htub.Chr01.H5 AUGUSTUS transcript 11603 16467 0.04 + . transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144" Htub.Chr01.H5 AUGUSTUS exon 11603 12143 0.2 + . transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144"; Htub.Chr01.H5 AUGUSTUS exon 13567 13766 0.23 + . transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144";
```

After installing and converting all files, I had trouble zipping due to the file size, so I zipped all of the files using (optional if not indexing immediately):

```         
#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=8G

set -euo pipefail

apptainer exec \
oras://community.wave.seqera.io/library/pbgzip:2016.08.04--c8e0d8eb135a301c \
  pbgzip -n 8 ref/figshare/*
  gzip -1 ref/gtf/*

echo "zipped"
```

# **Trimming:**

1.  Preform QC on fastq files (R1 and R2) using FASTQC and MultiQC scripts to see the quality of the files and determine how they should be trimmed.

    Fastqc file script:

    ```{bash}
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
    ```

    Loop:

    ```{bash}
    for fastq_file in data/fastq/*fastq.gz; do
        echo "running fastqc on $fastq_file"
        bash scripts/fastqc.sh $fastq_file
    done
    ```

    Next we will unzip each of the fastqc .zip output files to retreive the summary files. We'll use these to create a summary file of all of the fastqc files.

    ```{bash}
    for fastqczip in data/fastqc/*.zip
    do
    unzip $fastqczip 
    done
    ```

    The unzipped files will go in the \~ directory. From there, find the summary of each read fastqc file and concattenate it into a summary file.

    ```{bash}
    cat */summary.txt >> data/fastqc/fastqc_summaries.txt
    ```

Remove the inflated files with `rm -r SRR*`.

#### MultiQC Basic Script (there is a faster method below this script)

```{bash}
#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-multiqc-%j.out
set -euo pipefail

# Constants
MULTIQC_CONTAINER=oras://community.wave.seqera.io/library/multiqc:1.31--e111a2e953475c51

# Positional arguments: input files
qc=$1

# Output base directory
OUT=/fs/ess/PAS2880/users/bateman139/project/data

nopath1=${qc##*/} 

# Initial logging
echo "# Starting script multiqc.sh"
date
echo "# Input dir:                      $qc"
echo "# Output file:                    $OUT/multiqc"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$OUT/multiqc"

# Run MultiQC
apptainer exec "$MULTIQC_CONTAINER" multiqc \
    --outdir "$OUT/multiqc/$nopath1" \
    "$qc"

# Final logging
echo
echo "# Successfully finished script multiqc.sh with "$nopath1"
date

```

```{bash}
for fastqc_file in data/fastqc/*fastqc.zip; do
    echo "running MultiQC on $fastq_file"
    bash scripts/multiqc.sh $fastqc_file
done
```

Another way to do it so that all of the reports are on the same HTML file is by inputting a directory instead of an argument:

```{bash}
mkdir data/multiqc/allfiles
MULTIQC_CONTAINER=oras://community.wave.seqera.io/library/multiqc:1.31--e111a2e953475c51
apptainer exec "$MULTIQC_CONTAINER" multiqc --outdir data/multiqc/allfiles /fs/ess/PAS2880/users/bateman139/project/data/fastqc 
```

This will input all of the fastqc files in the fastqc dir into MultiQC. MultiQC can have other inputs from other programs, such as from STAR, Cutadapt, FastQ Screen, and featurecounts. To have them all on the same report, all of the files need to be copied to one directory, and read there.

## **Analysis of the MultiQC:**

**As we learned about using AI in this class, I wanted to experiment using the AI summary function in MultiQC. For this, I made an OpenAI account, put money into the account, created an API key, and input that key into the MultiQC Toolbox "Explain with AI" interface. AI Provider was OpenAI, Model used was GPT-5, and annonyize samples was turned on.**

#### **Report AI Summary**

-   Marked duplicate levels in several samples:

    -   Severe: SRR1555737_1 (61.86%), SRR1555737_2 (62.54%)

    -   High: SRR1555736_1 (51.13%), SRR1555736_2 (51.63%)

    -   Others show moderately elevated duplication (30--40%)

-   Pronounced N-content issues:

    -   Sharp N spike at 34 bp in: SRR1555695_1, SRR1555696_1, SRR1555697_1, SRR1555736_1, SRR1555737_1 (up to \~7.9% Ns)

    -   The remaining samples show persistent low-level Ns across reads (\~0.3--0.6%)

<!-- -->

-   Overall

    -   Read counts, lengths, and base qualities are consistent across samples (100 bp; \~23--25M reads each). Adapter/overrepresented sequences are negligible (\<0.4%).

    -   GC content is stable (44--47%) with only minor shape warnings in a subset.

-   High duplication rates (library complexity concern)

    -   Severe:

        -   SRR1555737_1 (61.86% dups)

        -   SRR1555737_2 (62.54% dups)

    -   Elevated:

        -   SRR1555736_1 (51.13% dups)

        -   SRR1555736_2 (51.63% dups)

    -   Others are moderate (≈31--40% dups).

-   Cycle-specific N spike (instrument/cycle issue)

    -   A pronounced N-call spike at base \~34 (\~7.9% Ns) in:

        -   SRR1555695_1, SRR1555696_1, SRR1555697_1, SRR1555736_1, SRR1555737_1

    -   Remaining samples show low but consistent Ns across reads (\~0.3--0.6%), likely benign.

    -   Per-tile quality modules frequently warn/fail, aligning with a cycle-/tile-specific artifact.

-   Minor additional flags

    -   Slight GC-distribution warnings in SRR1555737_1 and SRR1555737_2 (likely secondary to duplication or library composition).

-   No other major FastQC issues (base quality profiles are good; adapter content is minimal).

-   Recommendations

    -   Address duplication

        -   Mark/remove duplicates for downstream analyses; interpret coverage accordingly.

        -   For SRR1555737_1 and SRR1555737_2, consider downsampling or excluding if analyses are sensitivity-limited by library complexity; assess library complexity (e.g., preseq) before deciding to resequence.

        -   For future runs, reduce PCR cycles and/or use UMIs to enable accurate deduplication.

    -   Mitigate cycle artifact

        -   Confirm which read (R1/R2) and lane show the base-34 issue; review run/tile metrics. If BCLs are available, investigate re-basecalling or masking that cycle.

        -   In current data, apply stringent per-base quality/"N"-aware filtering; many aligners tolerate isolated Ns, but consider soft-clipping/quality trimming around the affected cycle if it improves alignment metrics.

    -   Sanity checks

        -   Verify that GC-distribution warnings in SRR1555737_1 and SRR1555737_2 reflect true biology and not contamination (quick taxonomic screen if unexpected).

        -   Proceed without adapter trimming (contamination is negligible), but retain standard quality trimming at read ends.

\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--

I am going to ignore this for now, as I would like to figure out what to do myself first, but it's nice to keep in mind that its an opiton.

### Interpretation and Trimming:

In order to trim, we need to understand what the standard is for reads that will go through alignment. Can there be any fails? What about warnings?

The areas that had Fails or Warnings were:

-   Per Tile Sequence Quality (all read 1s are fails and all read 2s are warnings)

-   Per Base Sequence Content (Either fail or warning)

-   Per Sequence GC Content has a warning for SRR1555737 (1 and 2)

-   Per Base N Content (Warnings for all read 1s)

-   Sequence Duplication Levels (Warn or fail)

-   Over represented seqeunces (present in 737 1 and 2, and 736 1 )

The good things about our reads is that:

-   The average base and sequence length are optimal for alignment.

-   Don't need adapter trimming.

Based on this we need to:

-   Trim low-quality 3' tails (standard).

-   Keep reads over 50 bp after trimming for alignment.

-   Filter reads with high N content.

To do this I used fastp. I chose fastp as most of the information, tutorials and guides on Multiqc analysis and trimming used fastp, and I found it easier to pick up.

I used the script below to trim the reads and analyze them in fastqc and multiqc:

```{bash}
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
```

On the MultiQC report, the AI said the following:

**Report AI Summary**

-   High duplicate rates flagged in four libraries: SRR1555737_1129_2 (64.35%), SRR1555737_1129_1 (61.86%), SRR1555736_1129_2 (53.29%), and SRR1555736_1129_1 (51.07%); remaining FastQC libraries show moderate duplication (30--40%).

-   Otherwise, metrics are within expected ranges: per-base quality high, GC stable (\~43--47%), very low adapter content (fastp: 0.77--1.47%), and good yield (fastp passed: \~43--48M reads); no other major QC concerns detected.

<!-- -->

-   Base qualities are consistently high (Phred \>32 across read length), GC content is consistent (≈43--47%), and insert-size distributions look normal across all datasets.

-   fastp-processed samples show good yields and low contamination:

    -   Reads after filtering: \~43--48M; %PF ≈93--95%; adapter-trimmed reads ≈0.77--1.47% (low)

    -   Duplication before filtering is low for most, with a moderate outlier in SRR1555737 (20.4%)

-   FastQC-only cohort shows universally elevated duplication (10/10 samples, range 30.7--64.4%):

    -   Severe: SRR1555737_1129_1 (61.86%), SRR1555737_1129_2 (64.35%)

    -   High-moderate: SRR1555736_1129_2 (53.29%), SRR1555736_1129_1 (51.07%)

    -   Moderate: SRR1555695_1129_2 (40.70%), SRR1555695_1129_1 (38.80%), SRR1555697_1129_2 (37.58%), SRR1555697_1129_1 (35.96%), SRR1555696_1129_2 (32.53%), SRR1555696_1129_1 (30.73%)

    -   Duplicate reads exceed unique reads in this set, indicating low library complexity/over-amplification.

-   Overrepresented sequences are minimal (report notes \<1% in 10 samples), and adapter content is negligible.

-   Next steps:

    -   For severe-duplication samples:

        -   SRR1555737_1129_1, SRR1555737_1129_2: consider excluding from downstream analyses or re-sequencing; if retained, plan strict duplicate removal and assess unique coverage sufficiency.

    -   For high/moderate-duplication samples:

        -   Apply duplicate marking/removal in downstream alignment; evaluate library complexity (e.g., Preseq) to judge whether more sequencing would yield unique reads.

        -   Review library prep (PCR cycle count, input DNA/RNA quality/quantity) and pooling; avoid over-amplification in future runs.

    -   For fastp group:

        -   SRR1555737: monitor duplication post-alignment; otherwise proceed.

    -   General:

        -   Proceed to alignment and post-alignment QC (mapping rate, insert-size, duplication, coverage uniformity).

        -   Keep standard adapter/quality trimming; no additional trimming indicated by current adapter metrics.

\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--

On the report, there still appears to be a large N content jump at the 35 read postion, but because we filtered all N reads over 5, this must be a small segment under 5, that STAR should be able to handle.

# **Alignment:**

Index script (make sure fasta and gtf files are unzipped and validated):

```{#!/bin/bash}
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-star_index-%j.out
#SBATCH --mem=160G
#SBATCH --time=24:00:00
set -euo pipefail

# Constants
STAR_CONTAINER=oras://community.wave.seqera.io/library/star:2.7.11b--84fcc19fdfab53a4

# Copy the placeholder variables
fasta=ref/figshare/HelTub_1.0.fn
gtf=ref/gtf/HelTub_1.0.gtf
outdir=ref/index/logs

# Initial logging
echo "# Starting script star_index.sh"
date
echo "# Input assembly FASTA file:      $fasta"
echo "# Input annotation GTF file:      $gtf"
echo "# Output dir:                     $outdir"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p ref/index
mkdir -p ref/index/logs

# Run STAR
apptainer exec "$STAR_CONTAINER" STAR \
    --runMode genomeGenerate \
    --genomeFastaFiles "$fasta" \
    --sjdbGTFfile "$gtf" \
    --sjdbOverhang 99 \
    --genomeSAindexNbases 14 \
    --genomeDir "$outdir" \
    --runThreadN 16

# Explanation of key options:
#   --runMode genomeGenerate
#       Build a genome index instead of aligning reads.
#   --runThreadN 16
#       Match Slurm's --cpus-per-task (STAR uses 16 threads).
#   --sjdbOverhang
#       Set to max read length - 1. For ~100 bp reads, 99.
#   --genomeSAindexNbases
#       Controls suffix array depth (speed vs RAM).
#       Default ≈ min(14, log2(genome_length)/2 - 1).
#       For ~10.5 Gb genome, that gives ~15.5 → capped to 14.
#       14 can spike RAM on huge plant genomes, so I used to start 13
#       to reduce memory pressure with minimal performance cost, but it produced too 
#       many large files. 
# Final logging
echo
echo "# Used STAR version:"
apptainer exec "$STAR_CONTAINER" STAR --version
echo "# Successfully finished script star_index.sh"
date

```

1.  Align the reads using the created index.

    \--\> Figure out ideal parameters for my genome using --help and internet sources (my genome is very large and hexaploid, which may complicate the process).

Sample for STAR Align + count

```{bash}
STAR \
  --runThreadN 8 \
  --genomeDir /path/to/STAR_index \
  --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix sample1_ \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts
```

1.  Check to see if the alignment worked using STAR log file, MultiQC on the fastqc and star logs, samtools on the BAM and IGV for visuals.

***Star log file:***

```         
less sample1_Log.final.out
```

-   Uniquely mapped reads %

-   \% reads mapped to multiple loci

-   \% reads unmapped: too many mismatches / too short / other

-   Mismatch rate per base %

-   Splice stats (number of splices, canonical vs non-canonical, annotated vs novel if you used a GTF)

    \-\--\> need to know what realistic splice patterns look like, and what the standard is to consider "most reads to be mapping" or "reads mapping uniquely".

**Samtools on Bam (sanity check):**

```         
samtools flagstat sample1_Aligned.sortedByCoord.out.bam
```

-   Total reads

-   How many mapped

-   How many properly paired

-   How many duplicates, etc.

```         
samtools idxstats sample1_Aligned.sortedByCoord.out.bam
```

-   Per-chromosome read counts

**Multiqc (optional for summary):**

```         
multiqc .
```

-   Mapping rates

-   Insert sizes

-   Per-sample comparison

**IGV (Optional for visualization)**

https://github.com/igvteam/igv-reports

Install:

```         
pip install igv-reports
```

-   Allows visualization of alignments using FASTA, BAM+BAI and GTF.

Can help determine if :

-   splices look clean

-   there are missing chunks in the genome

-   reads are landing on genes

-   The trimming is poor

-   rRNA or contamination

-   over-collasped repeats / homeologs

    \--\> Need to read the github page to figure out the commands to use. This can be a great resource but I have no idea how to use it at this moment.

**Generate Count Table**

Use Salmon, featurecounts or `STAR -quantMode GeneCounts`

Salmon:

1.  Need t make a transcript.fa using a tool like `gffread` inputting the genome.fa and the transcripts.fa

    ```{gffread annotation.gtf -g genome.fa -w transcripts.fa}
    ```

2.  Build salmon index

    Sample:

    ```         
    salmon index \
      -t transcripts.fa \
      -i salmon_index \
      -k 31
    ```

3.  Quantify each sample from the STAR transcriptom BAM

    Sample:

    ```         
    salmon quant \
      -i salmon_index \
      -l A \
      -a sample_Aligned.toTranscriptome.out.bam \
      -p 8 \
      -o salmon_quant/sample
    ```

    This outputs a quant.sf file.

4.  Make the gene count table (TSV) from salmon outputs in R using tximport.

    Sample R input:

    ```{R}
    samples <- c("sample1")  # add more as you get more samples
    files <- paste0(samples, "_ReadsPerGene.out.tab")
    names(files) <- samples

    count_list <- lapply(files, function(f) {
      tab <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
      # col 1 = gene, col 2 = unstranded counts (use 3 or 4 if stranded)
      counts <- tab[, c(1, 2)]
      colnames(counts) <- c("gene_id", sub("_ReadsPerGene.out.tab", "", f))
      counts
    })

    counts_merged <- Reduce(function(x, y) merge(x, y, by = "gene_id"),
                            count_list)

    # remove STAR summary rows
    counts_merged <- counts_merged[!grepl("^N_", counts_merged$gene_id), ]

    write.table(counts_merged,
                file = "gene_counts_matrix.tsv",
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
    ```

STAR:

Do alignment and counts at the same time using `--quantMode GeneCounts` .

Feature Counts:

**Differential Expression Analysis\
\
Functional Enrichment Analysis**

## Step 0: Set up
