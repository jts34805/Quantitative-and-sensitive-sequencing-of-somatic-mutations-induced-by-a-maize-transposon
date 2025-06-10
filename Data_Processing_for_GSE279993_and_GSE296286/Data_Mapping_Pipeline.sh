mkdir rawdata
mkdir processed
mkdir filtered
mkdir mapped
mkdir dedup
mkdir txtfiles ##making directories that will be used in the following loop. Read 1 and Read 2 of fastq.gz files should be loaded into the rawdata directory. 
for file in rawdata/*R1* ##set up a loop to iterate over all files in the rawdata directory, which should contain all fastqs from a Museq2 run  
do
  	file2="${file:8:-11}" #trim the first 8 and last 11 characters of the file name for the loop. Can be different values depending on the name of the fastq file. 
        fastp -w 8 -i $file -o processed/"$file2""R1A.fastq.gz" -f 23 -A -G -Q -L ###read 1 begins with the primer sequence that is 23 bp long. this step removes that. output file is put in the processed directory with the "R1A.fastq.gz" text appended to the end. 

        fastp -i processed/"$file2""R1A.fastq.gz" -I rawdata/"$file2"".R2.fastq.gz" -o processed/"$file2""R1B.fastq.gz" -O processed/"$file2""R2B.fastq.gz" -U --umi_loc read1 --umi_len 6 --umi_prefix TIR -h html/"$file2"".html" -A -G -Q -L ###the next 6 bp of read 1 is the edge of the TIR not apart of the primer sequence. It should match to "TATCTC" validation sequence. This line moves these 6 bp to the header so that later steps can filter based on matching to the expected sequence. 

       	fastp -i processed/"$file2""R1B.fastq.gz" -I processed/"$file2""R2B.fastq.gz" -o processed/"$file2""R1C.fastq.gz" -O processed/"$file2""R2C.fastq.gz" -U --umi_loc read2 --umi_len 8 --umi_prefix UMI -A -G -Q -L #Read 2 contains the UMI barcode for molecule counting. this step moves the 8 bp UMI to the header to enable molecule counting later on. 
        fastp -i processed/"$file2""R1C.fastq.gz" -I processed/"$file2""R2C.fastq.gz" -o filtered/"$file2""R1.fastq.gz" -O filtered/"$file2""R2.fastq.gz" -U --umi_loc read2 --umi_len 11 --umi_prefix BC --length_required 40 --trim_poly_x --cut_tail --adapter_fasta adapter2.1.fa #Read 2 also contains the 10 bp sample specific barcode, used for filtering library cross talk. This segment is similarly moved to the header of each read2 (including the "T" overhang). We also trim reads at this step, and require a minimum length of 40 as well as other routine fastp filters on our reads. we provide the adapter.fa file so that fastp can trim reads that contain a read through of our adapter sequence in read1.(see adapter2.1.fa in Data-Processing-for-GSE279993-and-GSE296286)




bowtie2 --threads 8 -x /scratch/jts34805/W22Build/W22chrscaff --phred33 -X 1000 --no-mixed --no-discordant -1 filtered/"$file2""R1.fastq.gz" -2 filtered/"$file2""R2.fastq.gz" | samtools view -@ 8 -b -o mapped/"$file2"".bam" #Using Bowtie2, the W22 V2 scaffold genome version, and Samtools, we convert our filtered reads to a sam file and subsequently convert to a bam file


        samtools sort -@ 8 mapped/"$file2"".bam" -o mapped/"$file2"".bam" #routine sorting of bam files
        samtools index -@ 8 mapped/"$file2"".bam" #routine indexing of bam files


        umi_tools group -I mapped/"$file2"".bam" --paired --chimeric-pairs discard --unpaired-reads discard --output-bam -S dedup/"$file2"".Dedup.bam" #using umi_tools, we deduplicate reads, remove chimeric pairs and unpaired reads.


        samtools view -f 64 -F 4 dedup/"$file2"".Dedup.bam" | awk -F" " '{split($1, a, "_"); print $3, $4, $8, $9, substr(a[2],1,6), substr(a[3],1,8), $5, length($10),substr(a[4],1,10), substr(a[4],11,11)}' > txtfiles/"$file2"".W22.txt")}' #using samtools view to convert the bam file to a sam file which is piped into the awk command to convert our sam files to a txt file which is a csv like file that can easily be loaded into R. "split" is used to split the header of each sam file entry so that we can get out the TIR, UMI, and Barcode, respectively. In the end, the txt file consists of the chromosome, position, position reverse, TIR Validation sequence, UMI, mapping quality, read length, Barcode, and ligation base (T overhang)

done 
