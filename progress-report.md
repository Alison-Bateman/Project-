Progess report: 
There is a lot more information in my Quarto file about all of the scripts and commands, as well as the planned scripts and commands. 

Scripts to be used in order: 
1. Tool kit download
 -Downloads NBCI Command-line tools and SRA toolkit to download and split interleaved fastq RNA-seq files from NBCI. 
Status: Complete 
2. Fastq file download and split
- one script and one loop for multiple files
Status: Complete 
3. FASTA and GFF file download from figshare
-Still needs explanation text because I don't completely understand the command
Status: Almost Complete 
4. FASTQC 
Status: Complete, but can be better- see note below 
5. MultiQC 
Status: Complete but I could combine the MultiQC and FASTQC so that the FASTQC loops and once finished the MultiQC runs on the dir.
6. FastP trimming and FASTQC and MultiQC script: 
Status: Very complete 
7. Index script
-Ran into issues with reaching quota 
Status: in Progress. 
8. Alignment script
- I have a rough draft of it/idea of it 
Status: Not started 
9. Script to generate count table:
Status: have not started and am still not sure if I'll use feature counts, salmon, or thwe built in STAR function
Status: Not Started 
10. R Analysis Scripts: 
Status: Have not started or thought about.

To-do List: 
1. Run the index script (12/1)
2. Create the alignment script (11/30)
3. figure out how to analyze the alignment results. (12/2)
4. Create the counts script to make a TSV table (12/2)
5. give that table to R for Analysis (12/2)
6. Figure out how to do Differential Expression Analysis and Functional Enrichment Analysis (12/3)
7.  Interpret results. (12/3)
8. Start presentation (12/3)
9. Finish presentation (12/6)
10. practice presentation (12/7+8)
11. Give presentation and submit project (12/9) 


Feedback request: 
1. I tried using sources online to figure out how to trim my reads for alignment. I'm not sure if I did so correctly. I just trimmed the ends and filtered reads with N greater than 5, and removed reads smaller than 50 bp. Would this be the correct way to do it? 

2. What information is my project_protocol missing? It's still in it's early stages but I would still like to know what it can improve on.  

3. For the alignment and thereafter, I'm not sure what I'm doing. I've been finding tutorials and guides online to help as I go, but I may need more help in the near future. If you have any advice that you can offer, that would be amazing. 
