##this code is for converting the mapped MuSeq2 reads to a matrix that can be analyzed in R. 
library(plyr) #this package is used to concatenate insertion locations based on UMI's
f = dir() #set working directory to directory where mapped MuSeq2 txt files are located.
if ('CleanData.rda' %in% f) { load('CleanData.rda') } else { CleanData = list() } #making "CleanData" can be memory intensive. Sometimes, this loop needs to be run in different sessions to reduce memory pressure. This "if" statement prevent overwriting CleanData if it has already been initialized. 
QCTable=matrix(data=NA, nrow=length(f), ncol=9) #initialize a matrix to store the QC metrics
colnames(QCTable)=c("PercentMatchValidation","PercentMatchAdapter","CorrectBarcode","IncorrectBarcode1","IncorrectBarcode2","IncorrectBarcode3","PercentIncorrect1","PercentIncorrect2","PercentIncorrect3")
rownames(QCTable)=gsub("_.+","",f) #name the rows of the table to match the txt files that are being processed
AdapterBarcodeMatching=read.csv("AdapterBarcodeTable")#create a csv consisting of two columns. First column is the sample name that matches the txt files, while the second column is the 10 bp adapter barcode sequence (see "MuSeq2 adapter preparation" section, table S2.) This process is optional if MuSeq2 libraries were made without the sample specific barcode. if so, use "Data_Processing_1B" scripts
rownames(AdapterBarcodeMatching)=AdapterBarcodeMatching[,1]
f = f[grepl('.txt', f)] #use grepl to grab any files in the directory that end in ".txt" which are all the mapped MuSeq2 txt files. 
for (f2 in f) { #begin loop processing through all txt files in the directory
  A = read.table(f2)
  colnames(A) = c("Chromosome","Position","PositionRev","Reverse","Validation","UMI","MAPQ","Length", "Barcode", "LigationBase") 
  A$Position[A$Reverse < 0] = A$PositionRev[A$Reverse < 0] - A$Reverse[A$Reverse < 0] - 1  # Shift position to edge of TIR rather than left-most base
  valmatch=table(A$Validation)
  QCTable[gsub("_.+","",f2),1]=(valmatch["TATCTC"]/sum(valmatch)*100) #fill in column 1 of the QCTable with the % of reads that have an exact match to the validation sequence "TATCTC"
  A = A[(A$MAPQ >= 10) & (A$Validation %in% c('TATCTC','TGTCTC')), -c(3,5,7:8,10)]# Require mapping quality greater than 10
  barcodes=table(A$Barcode)
  QCTable[gsub("_.+","",f2),2]=(barcodes[AdapterBarcodeMatching[gsub("_.+","",f2),2]]/sum(barcodes)*100) #calcs percent match to correct adapter
  QCTable[gsub("_.+","",f2),3]=names(barcodes[AdapterBarcodeMatching[gsub("_.+","",f2),2]]) #adds name of correct barcode to the table
  barcodes1=barcodes[names(barcodes)!=AdapterBarcodeMatching[gsub("_.+","",f2),2]] #removes the correct barcode from table
  QCTable[gsub("_.+","",f2),4:6]=names(barcodes1[order(barcodes1,decreasing=TRUE)][1:3]) #gets name of the 3 most common incorrect barcodes aNd adds to table
  QCTable[gsub("_.+","",f2),7:9]=(barcodes1[order(barcodes1,decreasing=TRUE)][1:3]/sum(barcodes)*100) #calculated the %for each incorrect barcode, important to include the correct barcode in this calc
  A=A[A$Barcode==AdapterBarcodeMatching[gsub("_.+","",f2),2],] #Filters A to only include samples with the correct barcode
  A = ddply(A, colnames(A), nrow) #use ddply to concatenate rows that are identical, adds a count column for each row that was concatenated
  colnames(A)[6] = 'Reads' 
  A$Direction = sign(A[,3])  
  A = A[order(A[,1],A[,2],-A$Direction),c(1,2,7,6,4)]
  CleanData[[f2]]=A #adds the newly processed file to the CleanData list. 
  save(CleanData, file = 'CleanData.rda') #saves CleanData list so that code saves progress, in case of excessive memory pressure. loop begins again with the next txt file in the directory.
}
