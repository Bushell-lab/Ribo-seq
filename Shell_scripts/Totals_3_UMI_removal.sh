#!/usr/bin/env bash

#This script uses cutadpat to remove UMIs. For the LEXOGEN CORALL Total RNA-Seq Library Prep Kit that we use for total RNA-seq, these are 12nt at the 5' end of the read.
#Edit the -u options if this is different for your specific library prep. Positive values remove bases from the 5' end of the read and negative values remove bases from the 3' end

#read in variables
source common_variables.sh

#run cutadapt with -u 12 to remove 12nt from the 5' end of all reads
for filename in $Totals_filenames
do
cutadapt $fastq_dir/${filename}_cdhitdup.fastq -u 12 -o $fastq_dir/${filename}_UMIremoved.fastq &
done
wait

#run fastqc on UMIremoved output
for filename in $Totals_filenames
do
fastqc $fastq_dir/${filename}_UMIremoved.fastq --outdir=$fastqc_dir &
done
wait
