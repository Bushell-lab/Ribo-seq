#!/usr/bin/env bash

#read in variables
source common_variables.sh

for filename in $filenames
do
for length in $(seq 25 35)
do
counting_script.py -bam $BAM_dir/${filename}_pc_best_sorted.bam -fasta $pc_fasta -len $length -out_file ${filename}_pc_best_L${length}_Off0.counts -out_dir $counts_dir &
done
done
wait
