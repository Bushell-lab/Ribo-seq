#!/usr/bin/env bash

#This script runs fastQC on all your raw fastq files and outputs them in the fastQC directory

#read in variables
source common_variables.sh

#activate fastQC conda environment
conda activate fastQC

#run fastQC
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}.fastq --outdir=$fastqc_dir &
done
wait

#deactivate fastQC conda environment
conda deactivate
