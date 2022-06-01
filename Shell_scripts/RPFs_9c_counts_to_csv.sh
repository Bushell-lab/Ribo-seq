#!/usr/bin/env bash

#read in variables
source common_variables.sh

for filename in $AK_RPF_filenames
do
counts_to_csv.py ${filename}_pc_final.counts $AK_most_abundant_fasta -all_transcripts -in_dir $counts_dir -out_dir $csv_counts_dir &
done
wait

for filename in $APC_RPF_filenames
do
counts_to_csv.py ${filename}_pc_final.counts $APC_most_abundant_fasta -all_transcripts -in_dir $counts_dir -out_dir $csv_counts_dir &
done
wait