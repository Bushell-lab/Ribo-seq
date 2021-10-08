#!/usr/bin/env bash

#This script runs fastQC on all your raw fastq files and outputs them in the fastQC directory

#read in variables
source common_variables.sh

#run fastQC
for filename in $filenames
do
fastqc $fastq_dir/${filename}.fastq --outdir=$fastqc_dir &
done
wait
