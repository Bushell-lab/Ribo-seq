#!/usr/bin/env bash

#read in variables
source common_variables.sh

#run summing_CDS_counts.py script to count all reads for each CDS for use with DEseq2
#the -remove_end_codons doesn't count any reads which correspond to A-site occupation within the first 20 or last 10 codons

for filename in $RPF_filenames
do
summing_CDS_counts.py ${filename}_pc_final.counts $region_lengths -remove_end_codons -in_dir $counts_dir -out_dir $cds_counts_dir &
done
wait
