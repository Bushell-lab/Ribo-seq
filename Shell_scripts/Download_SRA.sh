#!/usr/bin/env bash

###This script will download <.sra> files and convert them to <.fastq> files. You will need to ensure you have the SRA toolkit correctly installed

#make a list of all the SRA numbers that you want to download (these normally begin SRR)
SRAs='SRR00000001 SRR00000002 SRR00000003 SRR00000004' #this is just a template and needs to be edited

#set the directory where SRA toolkit will download the <.sra> files. This will have been user defined when you installed SRA toolkit
SRA_dir='~/Downloads/SRA_downloads/sra' #this is just a template and needs to be edited

#set the directory for where the <.fastq> files will be saved. This needs to be the same as the the common_variables.sh script.
fastq_files='/home/local/BICR/jwaldron/data/R11/external_sequencing_data/Author_year/fastq_files' #this is just a template and needs to be edited

#make a temporary directory to store temporary files for fasterq-dump
temp_dir=${fastq_files}/temp
mkdir $temp_dir

#download data with a for loop
for SRA in $SRAs
do
prefetch $SRA #this downloads the file as a <.sra>
fasterq-dump ${SRA_dir}/${SRA}.sra -O $fastq_files -t $temp_dir #this converts the <.sra> to <.fastq>
rm ${SRA_dir}/${SRA}.sra #this deletes the <.sra> file which is no longer needed
done

#delete the temporary directory
rmdir $temp_dir

#rename the files and compress them to save space

mv ${fastq_files}/SRR00000001.fastq ${fastq_files}/RPFs_1.fastq
gzip ${fastq_files}/RPFs_1.fastq

mv ${fastq_files}/SRR00000002.fastq ${fastq_files}/RPFs_2.fastq
gzip ${fastq_files}/RPFs_2.fastq

mv ${fastq_files}/SRR00000003.fastq ${fastq_files}/RPFs_3.fastq
gzip ${fastq_files}/RPFs_3.fastq

mv ${fastq_files}/SRR00000004.fastq ${fastq_files}/RPFs_4.fastq
gzip ${fastq_files}/RPFs_4.fastq
