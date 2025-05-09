This directory includes all code to map and process the raw fastq files into the processed data files available from GEO. This code will produce the 3 "Supplementary files" from GSE279993 and the 2 supplementary files from  GSE296286. Starting from raw fastqs, the scripts should be run in the following order:

1. Data_Mapping_Pipeline
2. Data_Processing_1_Read_Filtering_and_Molecular_Counting
3. Data_Processing_2_Final_Data_Processing
4. Code-to-reproduce-figures.R
