#!/usr/bin/env bash

#read in variables
source common_variables.sh

#run the counting_script.py with a range of read lengths (adjust below if required, currently set to 25-35)
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
counting_script.py -bam $BAM_dir/${filename}_pc_sorted.bam -fasta $most_abundant_fasta -len $length -out_file ${filename}_pc_L${length}_Off0.counts -out_dir $counts_dir &
done
done
wait
