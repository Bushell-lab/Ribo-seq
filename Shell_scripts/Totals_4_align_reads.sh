#!/usr/bin/env bash

#This script first uses bbmap to map to and remove reads corresponding to rRNAs, tRNAs and mitocondrial mRNAs.
#in specificies the input <.fastq> file.
#out specificies the output <.SAM> file.
#ref specifies the <.fasta> file to use as a reference. bbmap will use this to make an index. As this is much quicker than other alignment programs, we use the nodisk option so that this isn't written to file
#outm and outu specificies filenames to write <.fastq> files for all reads that either align or do not align respectively
#ambigous specifies how to treat multimapped reads. 
#2> stores the text that is printed to the screen as a log

#RSEM is then used to map the non rRNA/tRNA/mito mRNA reads to the same protein coding transcriptome as used for the RPFs.

#read in variables
source common_variables.sh

#Align to rRNA
for filename in $Totals_filenames
do
bbmap.sh in=$fastq_dir/${filename}_UMIremoved.fastq out=$SAM_dir/${filename}_rRNA.sam ref=$rRNA_fasta outm=$fastq_dir/${filename}_rRNA.fastq outu=$fastq_dir/${filename}_non_rRNA.fastq ambiguous=best nodisk 2> $log_dir/${filename}_rRNA_log.txt
done

#Align to protein coding transcriptomes
for filename in $Totals_filenames
do
rsem-calculate-expression --strandedness forward --bowtie2 --fragment-length-mean 300 --fragment-length-sd 100 $fastq_dir/${filename}_non_rRNA.fastq $rsem_index $rsem_dir/${filename} &
done
wait

