The directory "Data_Processing_for_GSE279993_and_GSE296286" includes all code to map and process the raw fastq files into the processed data files available from GEO. This code will produce the 3 "Supplementary files" from GSE279993 and the 2 supplementary files from  GSE296286. Starting from raw fastqs, the scripts should be run in the following order:

1. Data_Processing_for_GSE279993_and_GSE296286/Data_Mapping_Pipeline
2. Data_Processing_for_GSE279993_and_GSE296286/Data_Processing_1_Read_Filtering_and_Molecular_Counting
3. Data_Processing_for_GSE279993_and_GSE296286/Data_Processing_2_Final_Data_Processing
4. Code-to-reproduce-figures.R

This file describes the scripts and files in the Data_Processing_for_GSE279993_and_GSE296286, which convert the raw fastq.gz files available on GEO into the processed data files, that are also available on GEO. If you are interested in analyzing the exact data we published in these GEO series, we recommend donwloading the processed data files available on GEO under GSE 279993 and GSE296286. For MuSeq2 data generated elsewhere, we recommend following the scripts in the Data_Processing_for_GSE279993_and_GSE296286 directory

Raw MuSeq2 fastq.gz files can be mapped to the W22 V2 chromosome scaffold genome build and converted into .txt files consisting of columns chromosome, position, position reverse, Validation sequence, UMI, Barcode, read length, mapping quality, and ligation base using the script titled "Data_Mapping_Pipeline"

The .txt output files from "Data_Mapping_Pipeline" script can then be converted into an R list object using the "Data_Processing_1_Read_Filter_and_Molecular_Counting.R" script. This R list object, named "CleanData" is a list of molecular counts per genome site for each MuSeq2 library. 

"CleanData" can then be converted into an R matrix named "MuCounts" using the "Data_Processing_2_Final_Data_Processing.R" script in the "MuInsertionQuantification_R" directory. This process included converting the molecular count data to Mu insertions, and converging the list object into one matrix that reports the numbers of molecules detected per genome location per sample.

This R matrix can be used to quantify Mu Insertion Allele frequencies using scripts "Code-to-reproduce-figures" in the base directory. Alternatively, the processed data file "GSE279993_MuSeq2_data_for_maize_tissues.csv.gz" can be loaded directly from GSE279993, or "GSE296286_MuSeq2_count_data_for_bulk_pollen_outcross_experiments.csv.gz" from GSE296286.

file "adapter2.1.fa" is used by fastp commands in the "Data_Mapping_Pipeline" script to trim reads that sequence the reverse complement of read 2 in read 1 (usually occurs with short reads that are not informative). 

File "blastW22.out" is used to identify the ancestral (present in the W22 reference genome) Mu sequences, relevant for the "Data_Processing_2_CleanDatatoMuCountsMatrix script". "blastW22.out" was generated using the Mu TIR sequence in Lisch 2002 review of Mutator transposons and blasting this sequence (NCBI blast) against the W22 chromosome scaffold version 2 genome build.








