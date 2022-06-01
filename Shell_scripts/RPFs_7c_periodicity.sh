#!/usr/bin/env bash

#read in variables
source common_variables.sh

#set offset
offset=15

#run periodicity.py script
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
periodicity.py ${filename}_pc_L${length}_Off0.counts $region_lengths -offset $offset -in_dir $counts_dir -out_dir $periodicity_dir &
done
done
wait
