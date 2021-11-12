#!/usr/bin/env bash

#This script uses cutadpat to remove UMIs. For the nextflex library prep kit that we use for RPFs, these are 4nt at either end of the read.
#Edit the -u options if this is different for your specific library prep. Positive values remove bases from the 5' end of the read and negative values remove bases from the 3' end
#The standard McGlincy and Ingolia 2017 method, you need to remove 2 bases from the 5' end and 5 bases from the 3' end (-u 2 -u -5)

#read in variables
source common_variables.sh

#run cutadapt with -u 4 -u -4 to remove 4nt from either end of all reads
for filename in $RPF_filenames
do
cutadapt $fastq_dir/${filename}_cdhitdup.fastq -u 4 -u -4 -o $fastq_dir/${filename}_UMIremoved.fastq &
done
wait

#run fastqc on UMIremoved output
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}_UMIremoved.fastq --outdir=$fastqc_dir &
done
wait
