#!/usr/bin/env bash

#read in variables
source common_variables.sh

for filename in $RPF_filenames
do
summing_UTR5_counts.py ${filename}_pc_final.counts $region_lengths -in_dir $counts_dir -out_dir $UTR5_counts_dir &
done
wait
