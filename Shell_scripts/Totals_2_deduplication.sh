#!/usr/bin/env bash

#This script uses UMItools to extract the UMIs from the reads and add them to the read name
#The script is written for libraries which contain 12nt UMIs at the 5'end of the read. If this is not the case for your libraries you will need to change the following part of the command
#"--bc-pattern=NNNNNNNNNNNN"
#For more info on UMItools see https://github.com/CGATOxford/UMI-tools and https://umi-tools.readthedocs.io/en/latest/

#It will output a new fastq file with the suffix _UMI_clipped
#It then runs fastQC on output to check it is as expected

#read in variables
source common_variables.sh

#read deduplication
for filename in $Totals_filenames
do
umi_tools extract -I $fastq_dir/${filename}_cutadapt.fastq -S $fastq_dir/${filename}_UMI_clipped.fastq --bc-pattern=NNNNNNNNNNNN --log=$log_dir/${filename}_extracted_UMIs.log &
done
wait

#run fastqc on output
for filename in $Totals_filenames
do
fastqc $fastq_dir/${filename}_UMI_clipped.fastq --outdir=$fastqc_dir &
done
wait
