\ ## 
module load miniconda3/24.1.2-py310
conda create -n ncbi-datasets-cli
/users/PAS28804/bateman139/.conda/envs/ncbi-datasets-cli
conda install -c conda-forge ncbi-datasets-cli

## I had to create Conda env in scratch folder becuase I got an error saying that Disk size exceeded
conda create -y -p '/fs/scratch/PAS2880/conv-gff-to-gtf' 
conda activate '/fs/scratch/PAS2880/conv-gff-to-gtf' 
## Downloaded packages always go inside ~/.conda/pkgs by default, so change it 
## to scratch directory becuase you were getting an error that quota size exceeded
export CONDA_PKGS_DIRS=/fs/scratch/PAS2880/conda_pkgs
mkdir -p $CONDA_PKGS_DIRS
conda install -c conda-forge -c bioconda -c defaults perl-bioperl perl-lwp-protocol-https perl-yaml

## learn about the command using 
bp_genbank2gff3 --help
## Convert the gbff file to gff file
bp_genbank2gff3 GCA_030545235.1_HelTub_1.0_genomic.gbff.gz > GCA_030545235.1.gff3

## The gff file has been moved to Alison_help folder
## gff file contains the commnets for metadata than gff record
head Alison_help/GCA_030545235.1_HelTub_1.0_genomic.gbff.gz.gff -n 5000 > wrapped.gff
grep -v '^#' wrapped.gff > no_comments.gff
## Check the unique values in your third column
awk '{print $3}' no_comments.gff | sort -u

## Print the rows that contain chromosme in third column
awk '$3=="chromosome"' Alison_help/GCA_030545235.1_HelTub_1.0_genomic.gbff.gz.gff > chromosome_detail.gff

grep 'chromosome' Alison_help/GCA_030545235.1_HelTub_1.0_genomic.gbff.gz.gff
## The issue is it only output top level features- Chromosome and assembly gap.
## If you want to wrap each row separately use this command
fold -w 500 Alison_help/GCA_030545235.1_HelTub_1.0_genomic.gbff.gz.gff > wrapped.gff


zgrep "mRNA" Alison_help/GCA_030545235.1_HelTub_1.0_genomic.gbff.gz
grep "mRNA" Alison_help/GCA_030545235.1_HelTub_1.0_genomic.gbff.gz.gff | head 
zcat Alison_help/GCA_030545235.1_HelTub_1.0_genomic.gbff.gz | less

grep -P "\tgene\t" Alison_help/GCA_030545235.1.gff3 > genes.gff3


## Try to unzip the file and see if the output changes
gunzip Alison_help/GCA_030545235.1_HelTub_1.0_genomic.gbff.gz Alison_help/GCA_030545235.1_HelTub_1.0_genomic.gbff

