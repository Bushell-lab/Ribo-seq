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

#Align to genome
for filename in $Totals_filenames
do
STAR --readFilesIn $fastq_dir/${filename}_non_rRNA.fastq --runThreadN 12 --genomeDir $STAR_index --outFilterMultimapNmax 5 --outFilterMismatchNmax 5 --outSAMprimaryFlag AllBestScore --alignEndsType EndToEnd --outSAMtype BAM Unsorted --outSAMunmapped None --sjdbGTFfile $STAR_GTF --outFileNamePrefix $STAR_dir/${filename}
done
