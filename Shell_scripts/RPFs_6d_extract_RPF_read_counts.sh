#!/usr/bin/env bash

#This script extracts the read counts at each stage of the pipeline from the log files. It is paired with the RPF_read_counts.R, which will will then make plots for this QC.

#read in variables
source common_variables.sh

for filename in $RPF_filenames
do
extract_read_counts.py ${filename} RPFs -log_dir $log_dir &
done
wait
