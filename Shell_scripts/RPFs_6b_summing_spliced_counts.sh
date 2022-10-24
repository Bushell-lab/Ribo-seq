#!/usr/bin/env bash

#read in variables
source common_variables.sh

#set number of nt to splice
n=50

#run summing_spliced_counts.py script
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
summing_spliced_counts.py ${filename}_pc_L${length}_Off0.counts $n $region_lengths -in_dir $counts_dir -out_dir $spliced_counts_dir &
done
done
wait
