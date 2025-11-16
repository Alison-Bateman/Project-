For my final project, I will set up a reproducible RNA-seq analysis pipeline for Helianthus tuberosus (sunchoke) that will also be directly useful for my future research. I chose this project as my research project on breeding sunchoke will heavily involve sequencing, and building and documenting a RNA-seq pipeline will give me practical experience with the tools and workflows I’ll need, and will provide a template that I can adapt when I receive my own sequencing data.

A reference genome assembly for sunchoke (HelTub_1.0, GCA_030545235.1; Wang et al., 2024) was recently published to NCBI, with all files being available on Figshare (annotation files aren’t available on NBCI), and I plan to use this genome together with available sunchoke RNA-seq data to create a pipeline that I can later apply in my project	.
The RNA-seq data I will use come from Jung et al. (2014), specifically the H. tuberosus run SRR1555696, which I will download as raw FASTQ files from NCBI SRA:
https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR1555696&display=download

For the reference genome, I will use HelTub_1.0 (GCA_030545235.1) from Wang et al. (2024), obtained from Figshare, including the FASTA and GFF annotation files:
https://figshare.com/articles/dataset/Annotated_reference_genome_of_Helianthus_tuberosus/22491205
To transfer these files to my directory in the OSC, I’ll be using FileZilla (although I’m still having issues with it) for the FASTA and GFF files downloaded from Figshare, and the SRA program using prefetch, and fasterq-dump to separate the interleaved files into two separate read files.

- The main inputs to my pipeline will therefore be:
- The raw paired-end FASTQ files from SRR1555696
- The reference genome FASTA
- The corresponding genome annotation (GFF)

The main output I aim to produce is a TSV/CSV gene count table suitable for downstream analysis in R (e.g., with DESeq2), along with QC and alignment outputs created in the process.
 I will maintain the project in a version-controlled Git repository in OSC under /fs/ess/PAS2880/users/bateman139/project, and push the repository to GitHub. I plan to organize the project with a clear directomkdir ry structure, for example:

data/fastq/ – raw FASTQ files
data/trimmed/ – trimmed reads
ref/ – genome FASTA and annotation
scripts/ – bash scripts for each pipeline step
results/qc/ – FastQC and MultiQC reports
results/alignment/ – BAM files and STAR logs
results/counts/ – count matrices and related files
results/figures/ – plots generated in R

All commands will live in bash scripts under version control so that I can re-run the full pipeline or individual steps.
For tools, I plan to use:

SRA Toolkit (prefetch, fasterq-dump, vdb-validate) to download and validate the raw reads.
Bioperl to convert the GFF file to GTF for tools like STAR.
FastQC and MultiQC for quality control of raw and trimmed reads.
fastp for adapter and quality trimming while preserving read length appropriate for alignment.
STAR for alignment, with parameters that allow multi-mapping reads because the sunchoke genome is hexaploid and contains many homeologs.
Either featureCounts (with settings for multi-mapping reads) or Salmon after creating transcriptome for quantification.

The expected outputs include:
- QC: FastQC and MultiQC HTML reports for both raw and trimmed reads.
- Alignment: Sorted, indexed BAM files and STAR log files.
- Quantification: A gene count matrix  in TSV/CSV format.
- Downstream analysis: If time permits, I will perform differential expression analysis and Functional enrichment analysis in R using packages such as DESeq2 and clusterprofiler respectively. I would use RMarkdown or Quarto to summarize the analysis workflow and key findings, and save it under results/analysis.

Because this project involves many new tools and a significant amount of troubleshooting, I may not fully complete the downstream differential expression analysis within the course timeline. However, I expect to at least complete the pipeline through QC, trimming, alignment, and generation of a count table.

I am still uncertain about some technical details, particularly the optimal trimming parameters for this dataset, the best alignment settings for a hexaploid genome, and the most appropriate way to handle multi-mapping reads during quantification. I plan to rely on what we have covered in class, tool documentation, and online resources to figure out these issues as I go.

Jung WY, Lee SS, Kim CW, Kim H-S, Min SR, Moon JS, et al. (2014) RNA-Seq Analysis and De Novo Transcriptome Assembly of Jerusalem Artichoke (Helianthus tuberosus Linne). PLoS ONE 9(11): e111982. https://doi.org/10.1371/journal.pone.0111982

Sen Wang, Anqi Wang, Rong Chen, Dong Xu, Hengchao Wang, Fan Jiang, Hangwei Liu, Wanqiang Qian, Wei Fan, Haplotype-resolved chromosome-level genome of hexaploid Jerusalem artichoke provides insights into its origin, evolution, and inulin metabolism, Plant Communications, Volume 5, Issue 3, 2024,
100767, ISSN 2590-3462, https://doi.org/10.1016/j.xplc.2023.100767.




