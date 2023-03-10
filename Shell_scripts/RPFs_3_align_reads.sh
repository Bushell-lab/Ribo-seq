#!/usr/bin/env bash

#This script uses bbmap to align reads.
#in specificies the input <.fastq> file.
#out specificies the output <.SAM> file.
#ref specifies the <.fasta> file to use as a reference. bbmap will use this to make an index. As this is much quicker than other alignment programs, we use the nodisk option so that this isn't written to file
#outm and outu specificies filenames to write <.fastq> files for all reads that either align or do not align respectively
#ambigous specifies how to treat multimapped reads. We use best (keeps the highest scored alignment).
#trimreaddescription=t removes any white space and following text from the fastq files when writing the <.bam> files. This is important for using UMItools downstream if any samples have come from more than one sequencing run
#2> stores the text that is printed to the screen as a log
#We first align to rRNA and then tRNA, spliting the reads into aligned and unaligned <.fastq> files (without outputting a <.sam> file).
#All non-rRNA and non-tRNA reads are then aligned to a protein coding transcriptome (prefarably a <.fasta> containing the most abundant transcripts created in the total RNA-seq pipeline).

#fastQC is then ran on all the output <.fastq> files

#samtools is then used to sort and index the <.bam> file

#read in variables
source common_variables.sh

#Align to rRNA
for filename in $RPF_filenames
do
bbmap.sh in=$fastq_dir/${filename}_UMI_clipped.fastq ref=$rRNA_fasta outm=$fastq_dir/${filename}_rRNA.fastq outu=$fastq_dir/${filename}_non_rRNA.fastq ambiguous=best nodisk threads=$threadN 2> $log_dir/${filename}_rRNA_log.txt
done

#Align to tRNA fasta
for filename in $RPF_filenames
do
bbmap.sh in=$fastq_dir/${filename}_non_rRNA.fastq ref=$tRNA_fasta outm=$fastq_dir/${filename}_tRNA.fastq outu=$fastq_dir/${filename}_non_rRNA_tRNA.fastq ambiguous=best nodisk threads=$threadN 2> $log_dir/${filename}_tRNA_log.txt
done

#Align to protein coding transcriptome
for filename in $RPF_filenames
do
bbmap.sh in=$fastq_dir/${filename}_non_rRNA_tRNA.fastq out=$BAM_dir/${filename}_pc.bam ref=$most_abundant_fasta outm=$fastq_dir/${filename}_pc.fastq outu=$fastq_dir/${filename}_unaligned.fastq ambiguous=best nodisk trimreaddescription=t threads=$threadN 2> $log_dir/${filename}_pc_log.txt
done

#sort bam
for filename in $RPF_filenames
do
samtools sort $BAM_dir/${filename}_pc.bam -o $BAM_dir/${filename}_pc_sorted.bam -@ $threadN -m 1G
done

#index bam
for filename in $RPF_filenames
do
samtools index $BAM_dir/${filename}_pc_sorted.bam $BAM_dir/${filename}_pc_sorted.bai &
done
wait

#run fastqc on mapped reads
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}_rRNA.fastq --outdir=$fastqc_dir &
fastqc $fastq_dir/${filename}_tRNA.fastq --outdir=$fastqc_dir &
fastqc $fastq_dir/${filename}_pc.fastq --outdir=$fastqc_dir &
fastqc $fastq_dir/${filename}_unaligned.fastq --outdir=$fastqc_dir &
done
wait
