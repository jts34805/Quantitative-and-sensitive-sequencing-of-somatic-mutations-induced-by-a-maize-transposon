This directory includes all code to map and process the raw fastq files into the processed data files available from GEO. This code will produce the 3 "Supplementary files" from GSE279993 and the 2 supplementary files from GSE296286. 

For these scripts, the read 1 and read 2 fastqs will need to be in separate files. We recommend downloading fastqs using SRA-Toolkit and fasterq-dump to get the separated, paired end, fastq files. 

Starting from raw fastqs, the scripts should be run in the following order:
1. Data_Mapping_Pipeline.sh – this script will map fastq files to the maize genome, handle UMI information, and output a txt table that can be loaded into R
2. Data_Processing_1_Read_Filtering_and_Molecular_Counting.R – this script will load the txt table into R, filter reads with low mapping quality, lack of match to the transposon sequence, filter reads lacking match to expected adapter barcode sequence, and associate reads with initial molecules (molecular counting) using UMIs
3. Data_Processing_2_Final_Data_Processing.R – this script connects molecules sequencing out of the left and right TE border and outputs a matrix of TE-spanning molecules at any given site and sample
4. Family_Genotyping_For_GSE296286.R – this script makes genotype calls for the offspring and parents from entry GSE296286

The following additional files are used internal to these scripts:
1. adapter2.1.fa - used internal to Data_Mapping_Pipeline.sh script to trim reads that sequence the reverse complement of adapter (read 2)
2. blastW22.out - blast file generated from MU TIRs (from Lisch 2002 Mu transposon review), used in Data_Processing_2_Final_Data_Processing to match transposon spanning molecules that map to Mu elements in the reference genome
3. historicals.txt used in Family_Genotyping_for_GSE296286 to call historical Mu insertions in the inherited dataset file "Transmitted Mu insertions from bulk pollen outcross experiments.csv"
4. blacklist.txt used in Family_Genotyping_for_GSE296286 to remove insertions mapping to blacklisted genome regions. 



