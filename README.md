# Code for "Quantitative and sensitive sequencing of somatic mutations induced by a maize transposon"
This repository is for code related to the analysis, and mapping/quantitation of Museq2 data associated with the paper "Quantitative and sensitive sequencing of somatic mutations induced by a maize transposon"


To navigate this repository, raw fastqs should first be processed into txt files compatible with R using the script "Data_Mapping_Pipeline" in the "DataMappingPipeline-forGSE279993-and-GSE296286" directory. 

The processed .txt files from Data_Mapping_Pipeline script should then be converted into an R list object using the "Data_Processing_1_Read_Filter_and_Molecular_Counting.R" script in the "MuInsertionQuantification_R" directory. 
This R list object, named "CleanData" can then be converted into an R matrix named "MuCounts" using the "Data_Processing_2_CleanDataToMuCountsMatrix" script in the "MuInsertionQuantification" directory. 
This R matrix can be used to quantify Mu Insertion Allele frequencies using scripts "MuCounts_Analysis" in the "MuInsertionQuantification_R" directory. 
File "blastW22.out" is used to identify the ancestral (present in the W22 reference genome) Mu sequences, relevant for the "Data_Processing_2_CleanDatatoMuCountsMatrix script". "blastW22.out" was generated using the Mu TIR sequence in Lisch 2002 review of Mutator transposons and blasting this sequence (NCBI blast) against the W22 chromosome scaffold version 2 genome build. 
