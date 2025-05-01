# Mu_AFS
This repository is for code related to the analysis, and mapping/quantitation of Museq2 data associated with the paper "Tissue development shapes the abundance of somatic mutations"

To navigate this repository, raw fastqs should first be processed using the script in the "DataMappingPipeline" directory. The output txt files should be then converted into an R list object using the "Data_Processing_1" script in the "MuInsertionQuantification" directory. This R list object, named "CleanData" can then be converted into an R matrix named "MuCounts" using the "Data_Processing_2" script in the "MuInsertionQuantification" directory. this R matrix can be used to quantify Mu Insertion Allele frequencies 
