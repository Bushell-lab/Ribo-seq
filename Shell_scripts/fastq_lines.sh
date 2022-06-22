#!/usr/bin/env bash

#This script will calculate the lines in every fastq file (lines /4 = reads)

#read in variables
source common_variables.sh

wc -l $fastq_dir/*.fastq > $log_dir/fastq_lines.txt
