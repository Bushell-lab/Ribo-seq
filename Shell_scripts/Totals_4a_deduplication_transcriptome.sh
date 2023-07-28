#!/usr/bin/env bash

#This script uses UMItools to deduplicate the protein coding <.bam> file

#read in variables
source common_variables.sh

#run UMI tools deduplication function
for filename in $Totals_filenames
do
umi_tools dedup -I $BAM_dir/${filename}_pc_sorted.bam -S $BAM_dir/${filename}_pc_deduplicated.bam --output-stats=$log_dir/${filename}_deduplication 1> $log_dir/${filename}_deduplication_log.txt &
done
wait
