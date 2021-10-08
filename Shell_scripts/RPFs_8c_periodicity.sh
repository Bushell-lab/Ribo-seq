#!/usr/bin/env bash

#read in variables
source common_variables.sh

#set output directory
periodicity_dir=$parent_dir/'Analysis/periodicity'

#set offset
offset=15

#run periodicity.py script
for filename in $filenames
do
for length in $(seq 25 35)
do
periodicity.py ${filename}_pc_best_L${length}_Off0.counts $region_lengths -offset $offset -in_dir $counts_dir -out_dir $periodicity_dir &
done
done
wait
