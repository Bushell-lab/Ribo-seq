#!/usr/bin/env bash

#This script uses cd-hit-dup to remove PCR duplicates based on the UMIs
#It will output a fastq file with only unique reads
#It then runs fastQC on output to check it is as expected
#-i specifies the fastq input file
#-o specifies the fastq output file
#-e specifies the number of mismatches allowed. By setting to 0 means only unique reads will be kept

#For more info on cd-hit-dup see https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#cdhitdup

#read in variables
source common_variables.sh

#activate cd-hit environment
conda activate cdhit

#read deduplication
for filename in $RPF_filenames
do
cd-hit-dup -i $fastq_dir/${filename}_cutadapt.fastq -o $fastq_dir/${filename}_cdhitdup.fastq -e 0
done

#deactivate cd-hit environment and activate fastQC environment
conda deactivate
conda activate fastQC

#run fastqc on cd-hit-dup output
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}_cdhitdup.fastq --outdir=$fastqc_dir &
done
wait

#deactivate fastQC environment
conda deactivate
