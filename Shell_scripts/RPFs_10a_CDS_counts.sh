#!/usr/bin/env bash

#read in variables
source common_variables.sh

#set output directory
cds_counts_dir=$parent_dir/'Analysis/CDS_counts'

#run summing_CDS_counts.py script
for filename in $filenames
do
summing_CDS_counts.py ${filename}_pc_all_L29-33_Off12-13.counts $region_lengths -remove_end_codons -in_dir $counts_dir -out_dir $cds_counts_dir &
done
wait