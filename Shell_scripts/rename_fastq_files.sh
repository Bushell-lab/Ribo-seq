#!/usr/bin/env bash

#read in variables
source common_variables.sh

#rename the totals files
mv ${fastq_files}/SRR00001.fastq ${fastq_files}/RPFs_1.fastq
mv ${fastq_files}/SRR00002.fastq ${fastq_files}/RPFs_2.fastq

