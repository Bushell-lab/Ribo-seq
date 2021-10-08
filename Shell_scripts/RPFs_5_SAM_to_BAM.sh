#!/usr/bin/env bash

#This script uses SAMtools to convert <.sam> files to <.bam> files, sort the <.bam> files and make an <.bai> index files.
#When using the sorted <.bam> file as input to the counting_script.py script, ensure the corresponding <.bai> (index file) is in the same directory

#read in variables
source common_variables.sh

#convert sam to sam 
for filename in $RPF_filenames
do
samtools view -bS $SAM_dir/${filename}_pc_all.sam > $BAM_dir/${filename}_pc_all.bam &
samtools view -bS $SAM_dir/${filename}_pc_best.sam > $BAM_dir/${filename}_pc_best.bam &
done
wait

#sort bam
for filename in $RPF_filenames
do
samtools sort $BAM_dir/${filename}_pc_all.bam -o $BAM_dir/${filename}_pc_all_sorted.bam -@ 6 -m 1G
samtools sort $BAM_dir/${filename}_pc_best.bam -o $BAM_dir/${filename}_pc_best_sorted.bam -@ 6 -m 1G
done

#index bam
for filename in $RPF_filenames
do
samtools index $BAM_dir/${filename}_pc_all_sorted.bam $BAM_dir/${filename}_pc_all_sorted.bai &
samtools index $BAM_dir/${filename}_pc_best_sorted.bam $BAM_dir/${filename}_pc_best_sorted.bai &
done
wait
