#!/usr/bin/env bash

#This script uses RSEM to map reads to the protein coding transcriptome. It is strongly recommended to use the filtered transcriptome unless there is a good reason otherwise
#--strandedness describes the strand of the genome that the sequencing reads should align to. For the CORALL kit this is forward, but for a lot of standed Illumina Tru-seq kits this will be reverse. If this is not known then it is best to try both and the alignment logs should tell you which is correct
#--bowtie2 sets bowtie2 as the aligner
#--fragment-length-mean 300 --fragment-length-sd 100 sets the mean and standard deviation of the library fragment size. These do not need to be exact but best estimates. These are good starting values to use for the CORALL kit
#if reads are paired end, then the following line needs to be added and also inputs for --fragment-length-mean and --fragment-length-sd don't need to be set as these will be calculated from the data
#--paired-end $fastq_dir/${filename}_UMI_clipped_R1.fastq $fastq_dir/${filename}_UMI_clipped_R2.fastq
#Note that the options for extracting UMIs with UMI tools in the previous script will also have had to been altered for paired end reads

#read in variables
source common_variables.sh

#Align to protein coding transcriptome
for filename in $Totals_filenames
do
rsem-calculate-expression --strandedness forward --bowtie2 --fragment-length-mean 300 --fragment-length-sd 100 $fastq_dir/${filename}_UMI_clipped.fastq $rsem_index $rsem_dir/${filename} &
done
wait

