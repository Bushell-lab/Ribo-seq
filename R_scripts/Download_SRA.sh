#!/usr/bin/env bash

###This script downloads all the RPF and Total RNA-seq <.fastq> files from the Steinberger et al. 2020 paper. Paths will need to be changed for each user and SRA and filenames adapted for alternative datasets

#make a list of all the SRA numbers that you want to download
SRAs='SRR11915236 SRR11915237 SRR11915238 SRR11915239 SRR11915240 SRR11915241 SRR11915242 SRR11915243 SRR11915244 SRR11915245 SRR11915246 SRR11915247'

#set the directory where SRA toolkit will download the <.sra> files
SRA_dir='~/Downloads/SRA_downloads/sra'

#set the parent directory for where the <.fastq> files will be saved
parent_dir='/home/local/BICR/jwaldron/data/R11/external_sequencing_data/Steinberger_2020'

#make a directory to store <.fastq> files
mkdir ${parent_dir}/fastq_files

#make a temporary directory to store temporary files for fasterq-dump
mkdir ${parent_dir}/temp

#download data with a for loop
for SRA in $SRAs
do
prefetch $SRA #this downloads the file as a <.sra>
fasterq-dump ${SRA_dir}/${SRA}.sra -O ${parent_dir}/fastq_files -t ${parent_dir}/temp #this converts the <.sra> to <.fastq>
rm ~/Downloads/SRA_downloads/sra/${SRA}.sra #this deletes the <.sra> file which is no longer needed
done
wait

#delete the temporary directory
rmdir ${parent_dir}/temp

#rename the files and compress them to save space

mv ${parent_dir}/fastq_files/SRR11915236.fastq ${parent_dir}/fastq_files/Ctrl_1_RPFs.fastq
gzip ${parent_dir}/fastq_files/Ctrl_1_RPFs.fastq

mv ${parent_dir}/fastq_files/SRR11915237.fastq ${parent_dir}/fastq_files/Ctrl_2_RPFs.fastq
gzip ${parent_dir}/fastq_files/Ctrl_2_RPFs.fastq

mv ${parent_dir}/fastq_files/SRR11915238.fastq ${parent_dir}/fastq_files/Ctrl_3_RPFs.fastq
gzip ${parent_dir}/fastq_files/Ctrl_3_RPFs.fastq

mv ${parent_dir}/fastq_files/SRR11915239.fastq ${parent_dir}/fastq_files/Hipp_1_RPFs.fastq
gzip ${parent_dir}/fastq_files/Hipp_1_RPFs.fastq

mv ${parent_dir}/fastq_files/SRR11915240.fastq ${parent_dir}/fastq_files/Hipp_2_RPFs.fastq
gzip ${parent_dir}/fastq_files/Hipp_2_RPFs.fastq

mv ${parent_dir}/fastq_files/SRR11915241.fastq ${parent_dir}/fastq_files/Hipp_3_RPFs.fastq
gzip ${parent_dir}/fastq_files/Hipp_3_RPFs.fastq

mv ${parent_dir}/fastq_files/SRR11915242.fastq ${parent_dir}/fastq_files/Ctrl_1_Totals.fastq
gzip ${parent_dir}/fastq_files/Ctrl_1_Totals.fastq

mv ${parent_dir}/fastq_files/SRR11915243.fastq ${parent_dir}/fastq_files/Ctrl_2_Totals.fastq
gzip ${parent_dir}/fastq_files/Ctrl_2_Totals.fastq

mv ${parent_dir}/fastq_files/SRR11915244.fastq ${parent_dir}/fastq_files/Ctrl_3_Totals.fastq
gzip ${parent_dir}/fastq_files/Ctrl_3_Totals.fastq

mv ${parent_dir}/fastq_files/SRR11915245.fastq ${parent_dir}/fastq_files/Hipp_1_Totals.fastq
gzip ${parent_dir}/fastq_files/Hipp_1_Totals.fastq

mv ${parent_dir}/fastq_files/SRR11915246.fastq ${parent_dir}/fastq_files/Hipp_2_Totals.fastq
gzip ${parent_dir}/fastq_files/Hipp_2_Totals.fastq

mv ${parent_dir}/fastq_files/SRR11915247.fastq ${parent_dir}/fastq_files/Hipp_3_Totals.fastq
gzip ${parent_dir}/fastq_files/Hipp_3_Totals.fastq
