#!/usr/bin/env bash

#This script uses UMItools to deduplicate the protein coding <.bam> file
#It then uses samtools to sort the <.bam> files and make <.bai> index files.
#When using the sorted <.bam> file as input to the counting_script.py script, ensure the corresponding <.bai> (index file) is in the same directory

#read in variables
source common_variables.sh

#run UMI tools deduplication function
for filename in $RPF_filenames
do
umi_tools dedup -I $BAM_dir/${filename}_pc_sorted.bam -S $BAM_dir/${filename}_pc_deduplicated.bam --output-stats=$log_dir/${filename}_deduplication 1> $log_dir/${filename}_deduplication_log.txt &
done
wait

#sort bam
for filename in $RPF_filenames
do
samtools sort $BAM_dir/${filename}_pc_deduplicated.bam -o $BAM_dir/${filename}_pc_deduplicated_sorted.bam -@ threadN -m 1G
done

#index bam
for filename in $RPF_filenames
do
samtools index $BAM_dir/${filename}_pc_deduplicated_sorted.bam $BAM_dir/${filename}_pc_deduplicated_sorted.bai &
done
wait
