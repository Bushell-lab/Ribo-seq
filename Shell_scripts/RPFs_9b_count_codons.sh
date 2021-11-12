#!/usr/bin/env bash

#read in variables
source common_variables.sh

#run summing_CDS_counts.py script
#use the best alignments so that all transcripts can be used, without any reads counting twice

for filename in $RPF_filenames
do
count_codon_occupancy.py ${filename}_pc_best_final.counts $pc_fasta $region_lengths -in_dir $counts_dir -out_dir $codon_counts_dir &
done
wait
