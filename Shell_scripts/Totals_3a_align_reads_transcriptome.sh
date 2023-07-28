#!/usr/bin/env bash

#This script uses bowtie2 to align reads.
#-S specifies the <.sam> output file name
#-U specificies the input <.fastq> file.
#-x specifies the index to use as a reference. This is the same one as used for rsem to calculate isoform expression
#--threads specifies the number of threads
#"--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1" are the arguments recommended if using rsem downstream

#it then uses samtools to sort and index the <.bam> file

#read in variables
source common_variables.sh

#Align to protein coding transcriptome
for filename in $Totals_filenames
do
bowtie2 -S $SAM_dir/${filename}_pc.sam -U $fastq_dir/${filename}_UMI_clipped.fastq -x $rsem_index --threads $threadN --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 2> $log_dir/${filename}_pc_log.txt
done

#convert sam to bam 
for filename in $Totals_filenames
do
samtools view -b $SAM_dir/${filename}_pc.sam > $BAM_dir/${filename}_pc.bam &
done
wait

#sort bam
for filename in $Totals_filenames
do
samtools sort $BAM_dir/${filename}_pc.bam -o $BAM_dir/${filename}_pc_sorted.bam -@ $threadN -m 1G
done

#index bam
for filename in $Totals_filenames
do
samtools index $BAM_dir/${filename}_pc_sorted.bam $BAM_dir/${filename}_pc_sorted.bai &
done
wait
