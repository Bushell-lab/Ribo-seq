#!/usr/bin/env bash

#read in variables
source common_variables.sh

###This script will download <.sra> files and convert them to <.fastq> files. You will need to ensure you have the SRA toolkit correctly installed

#make a list of all the SRA numbers that you want to download (these normally begin SRR)
SRAs='SRR00000001 SRR00000002 SRR00000003 SRR00000004' #this is just a template and needs to be edited

#set the directory where SRA toolkit will download the <.sra> files. This will have been user defined when you installed SRA toolkit
SRA_dir='~/Downloads/SRA_downloads/sra' #this is just a template and needs to be edited

#make a temporary directory to store temporary files for fasterq-dump
temp_dir=${fastq_dir}/temp
mkdir $temp_dir

#download data with a for loop
for SRA in $SRAs
do
prefetch $SRA #this downloads the file as a <.sra>
fasterq-dump ${SRA_dir}/${SRA}.sra -O $fastq_dir -t $temp_dir #this converts the <.sra> to <.fastq>
rm ${SRA_dir}/${SRA}.sra #this deletes the <.sra> file which is no longer needed
done

#delete the temporary directory
rmdir $temp_dir

#rename the files and compress them to save space

mv ${fastq_dir}/SRR00000001.fastq ${fastq_dir}/RPFs_1.fastq
gzip ${fastq_dir}/RPFs_1.fastq

mv ${fastq_dir}/SRR00000002.fastq ${fastq_dir}/RPFs_2.fastq
gzip ${fastq_dir}/RPFs_2.fastq

mv ${fastq_dir}/SRR00000003.fastq ${fastq_dir}/RPFs_3.fastq
gzip ${fastq_dir}/RPFs_3.fastq

mv ${fastq_dir}/SRR00000004.fastq ${fastq_dir}/RPFs_4.fastq
gzip ${fastq_dir}/RPFs_4.fastq
