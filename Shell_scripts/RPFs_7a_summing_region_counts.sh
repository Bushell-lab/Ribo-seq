#!/usr/bin/env bash

#read in variables
source common_variables.sh

#set offset
offset=15

#run summing_region_counts.py script
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
summing_region_counts.py ${filename}_pc_L${length}_Off0.counts $offset $region_lengths -in_dir $counts_dir -out_dir $region_counts_dir &
done
done
wait
