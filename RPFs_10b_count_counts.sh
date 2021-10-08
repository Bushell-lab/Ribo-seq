#!/usr/bin/env bash

#read in variables
source common_variables.sh

#set output directory
codon_counts_dir=$parent_dir/'Analysis/codon_counts'

#run summing_CDS_counts.py script
for filename in $filenames
do
count_codon_occupancy.py ${filename}_pc_L27-29_Off0.counts $pc_fasta $region_lengths -in_dir $Count_dir -out_dir $codon_counts_dir &
done
wait
