#!/usr/bin/env bash

#read in variables
source common_variables.sh

#set lengths and offsets. These are the final ones which you will use for any DE analysis.
#Use the plots generated by the periodicity and offset scripts to decide what values to set these at

lengths='27,28,29,30,31,32'
offsets='12,12,12,12,13,13'

#run counting script to generate <.counts> files
#use all alignments so that the most abundant transcript per gene (based on paired total RNA-seq data) can be selected downstream

for filename in $RPF_filenames
do
counting_script.py -bam $BAM_dir/${filename}_pc_all_sorted.bam -fasta $pc_fasta -len $lengths -offset $offsets -out_file ${filename}_pc_all_final_CDS.counts -out_dir $counts_dir &
done
wait

#run summing_CDS_counts.py script to count all reads for each CDS for use with DEseq2
#the -remove_end_codons doesn't count any reads which correspond to A-site occupation within the first 20 or last 10 codons

for filename in $RPF_filenames
do
summing_CDS_counts.py ${filename}_pc_all_final_CDS.counts $region_lengths -remove_end_codons -in_dir $counts_dir -out_dir $cds_counts_dir &
done
wait

