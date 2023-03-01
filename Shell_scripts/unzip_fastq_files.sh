#!/usr/bin/env bash

#read in variables
source common_variables.sh

#unzip files
for filename in $RPF_filenames
do
gunzip $fastq_dir/${filename}.fastq.gz &
done
wait
