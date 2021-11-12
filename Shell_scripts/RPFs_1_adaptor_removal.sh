#!/usr/bin/env bash

#This script runs cutadapt on all your raw fastq files. For more info on settings see https://cutadapt.readthedocs.io/en/stable/guide.html

#You need to make sure the 3' adaptor sequence in the variables.txt file is correct. The -a option specifies this here

#The --nextseq-trim=20 option will trim bases from the 3' end if the quality score is below 20. This is only for sequencing data generated on the
#Next-seq (which is what we have here at the Beatson). If using external data set that has been sequenced on a Hi-seq, replace this command with -q 20
#It should state which sequencing platform was used on the GEO page
#The following text from the cutadpat manual explains why this is

#Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq.
#In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. However, dark cycles also occur when sequencing “falls off” the end of the fragment.
#The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.
#Since the regular quality-trimming algorithm cannot deal with this situation, you need to use the --nextseq-trim option:
#This works like regular quality trimming (where one would use -q 20 instead), except that the qualities of G bases are ignored.

#-m and -M specify the minimum and Maximum read lengths following adaptor removal and base trimming. RPFs should be roughly 30nt but can vary.
#However if the libraries also contain UMIs then you need to add this on too. UMIs from the nextflex kit that we use are 4nt on each end of the read
#therefore a 30nt RPF should be 38nt with UMIs following adaptor removal. Using -m 30 -M 50 means after UMI removal you will be left with read lengths 22-42
#which is suitable for this situation but will need to be modified to suit specific needs if using external data which may have different UMI lengths or may not contain UMIs

#1> causes all the text that is normally printed to the screen to be saved in a log file in your logs directory for each sample

#once cutadapt is complete, fastQC is ran on the output fastq files to check they are as expected

#read in variables
source common_variables.sh

#run cutadapt
for filename in $RPF_filenames
do
cutadapt $fastq_dir/${filename}.fastq -a $RPF_adaptor --nextseq-trim=20 -m 30 -M 50 -o $fastq_dir/${filename}_cutadapt.fastq 1> $log_dir/${filename}_cutadapt_log.txt &
done
wait

#run fastqc on cutadapt output
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}_cutadapt.fastq --outdir=$fastqc_dir &
done
wait
