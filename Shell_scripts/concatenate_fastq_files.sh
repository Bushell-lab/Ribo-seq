#!/usr/bin/env bash

#read in variables
source common_variables.sh

#concatenate seperate fastq files into one
cat ${fastq_dir}/RPFs_1a.fastq ${fastq_dir}/RPFs_1b.fastq > ${fastq_dir}/RPFs_1.fastq
cat ${fastq_dir}/RPFs_2a.fastq ${fastq_dir}/RPFs_2b.fastq > ${fastq_dir}/RPFs_2.fastq


