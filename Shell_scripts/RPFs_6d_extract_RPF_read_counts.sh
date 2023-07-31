#!/usr/bin/env bash

#read in variables
source common_variables.sh

for filename in $RPF_filenames
do
extract_read_counts.py ${filename} RPFs -log_dir $log_dir &
done
wait
