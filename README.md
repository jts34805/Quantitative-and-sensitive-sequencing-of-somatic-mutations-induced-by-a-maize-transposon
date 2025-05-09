This directory includes all code to map and process the raw fastq files into the processed data files available from GEO. This code will produce the 3 "Supplementary files" from GSE279993 and the 2 supplementary files from  GSE296286. Starting from raw fastqs, the scripts should be run in the following order:

1. Data_Processing_for_GSE279993_and_GSE296286/Data_Mapping_Pipeline
2. Data_Processing_for_GSE279993_and_GSE296286/Data_Processing_1_Read_Filtering_and_Molecular_Counting
3. Data_Processing_for_GSE279993_and_GSE296286/Data_Processing_2_Final_Data_Processing
4. Code-to-reproduce-figures.R

More detailed information for processing raw fastq.gz files into processed data files is available in the README.md file in the "Data_Processing_for_GSE279993_and_GSE296286" directory 


Readme: Code for "Quantitative and sensitive sequencing of somatic mutations induced by a maize transposon"

This repository is for code related to the analysis, and mapping/quantification of Museq2 data associated with the paper "Quantitative and sensitive sequencing of somatic mutations induced by a maize transposon"
The subdirectory "Data_Processing_for_GSE279993_and_GSE296286" includes scripts to map and process raw fastq files (available from GEO) to a matrix of molecule counts for each transposon insertion site in each sample (also available from GEO). Alternatively, the processed data files can be downloaded from the "Supplementary files" section of the respective GEO entries and these files can be used directly into the "Code to reproduce figures.R".

"Code to reproduce figures.R" contains all code to analyze the processed data files (which can be downloaded from GEO) and produce all main text figures in the paper.
