#!/usr/bin/env bash

#read in variables
source common_variables.sh

#rename the totals files
mv ${fastq_dir}/SRR00001.fastq ${fastq_dir}/RPFs_1.fastq
mv ${fastq_dir}/SRR00002.fastq ${fastq_dir}/RPFs_2.fastq

