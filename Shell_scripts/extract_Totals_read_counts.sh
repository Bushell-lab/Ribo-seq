#!/usr/bin/env bash

#read in variables
source common_variables.sh

#monosomes
for filename in $Totals_filenames
do
extract_read_counts.py ${filename} Totals -log_dir $log_dir &
done
wait