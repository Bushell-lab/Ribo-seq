#!/usr/bin/env bash

#This script uses STAR to align reads to a genome, which may be neccessary to check for specific knockouts, where specific exons have been removed etc
#It is worth creating a new conda environment just for STAR
#for more info on STAR see https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

#"--outFilterMultimapNmax" sets maximum number of loci the read is allowed to map to (default is 20, but we set to 5 to be more stringent)
#"--outFilterMismatchNmax" sets the maximum number of mismatches (default is 10, but we set to 5 to be more stringent)
#"--outSAMprimaryFlag AllBestScore" will output all alignments with the best score as primary alignments, rather than default behaviour which will only mark one alignment as primary
#"--alignEndsType EndToEnd" forces end-to-end read alignment, do not soft-clip
#"--outSAMtype BAM SortedByCoordinate" outputs a sorted <.bam> file, circumventing the requirement to run samtools afterwards
#"--outSAMunmapped None" means unmapped reads are not output (set as default)

#read in variables
source common_variables.sh

#Align to genome
for filename in $Totals_filenames
do
STAR --readFilesIn $fastq_dir/${filename}_UMI_clipped.fastq --runThreadN $threadN --genomeDir $STAR_index --outFilterMultimapNmax 5 --outFilterMismatchNmax 5 --outSAMprimaryFlag AllBestScore --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --sjdbGTFfile $STAR_GTF --outFileNamePrefix $STAR_dir/${filename}
done
