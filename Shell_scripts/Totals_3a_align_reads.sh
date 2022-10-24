#!/usr/bin/env bash

#This script uses bbmap to align reads.
#in specificies the input <.fastq> file.
#out specificies the output <.SAM> file.
#ref specifies the <.fasta> file to use as a reference. bbmap will use this to make an index. As this is much quicker than other alignment programs, we use the nodisk option so that this isn't written to file
#ambigous specifies how to treat multimapped reads. We use all so that RSEM can calculated relative isoform expression downstream
#2> stores the text that is printed to the screen as a log

#read in variables
source common_variables.sh

#Align to protein coding transcriptome
for filename in $Totals_filenames
do
bbmap.sh in=$fastq_dir/${filename}_UMI_clipped.fastq out=$SAM_dir/${filename}_pc.sam ref=$pc_fasta ambiguous=all nodisk threads=$threadN 2> $log_dir/${filename}_pc_log.txt
done
