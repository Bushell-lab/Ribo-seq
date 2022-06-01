#!/usr/bin/env bash

#read in variables
source common_variables.sh

for filename in $RPF_filenames
do
counts_to_csv.py ${filename}_pc_final.counts $most_abundant_fasta -all_transcripts -in_dir $counts_dir -out_dir $csv_counts_dir &
done
wait
