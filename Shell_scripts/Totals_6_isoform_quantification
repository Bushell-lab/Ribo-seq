#!/usr/bin/env bash

#This script uses RSEM to calculate isoform and gene level expression from deduplicated BAM file.
#--strandedness describes the strand of the genome that the sequencing reads should align to. For the CORALL kit this is forward, but for a lot of standed Illumina Tru-seq kits this will be reverse. If this is not known then it is best to try both and the alignment logs should tell you which is correct
#--fragment-length-mean 300 --fragment-length-sd 100 sets the mean and standard deviation of the library fragment size. These do not need to be exact but best estimates. These are good starting values to use for the CORALL kit

#read in variables
source common_variables.sh

#Align to protein coding transcriptome
for filename in $Totals_filenames
do
rsem-calculate-expression --strandedness forward --fragment-length-mean 300 --fragment-length-sd 100 --alignments $BAM_dir/${filename}_pc_deduplicated.bam $rsem_index $rsem_dir/${filename} &
done
wait
