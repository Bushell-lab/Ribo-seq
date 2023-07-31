#!/usr/bin/env bash

#Read in variables
source common_variables.sh

#Run cutadapt
for filename in $Totals_filenames
do
cutadapt $fastq_dir/${filename}.fastq -a $Totals_adaptor --nextseq-trim=20 -m 30 --cores=0 -o $fastq_dir/${filename}_cutadapt.fastq 1> $log_dir/${filename}_cutadapt_log.txt
done

#Extract UMIs
for filename in $Totals_filenames
do
umi_tools extract -I $fastq_dir/${filename}_cutadapt.fastq -S $fastq_dir/${filename}_UMI_clipped.fastq --bc-pattern=NNNNNNNNNNNN --log=$log_dir/${filename}_extracted_UMIs.log &
done
wait

#Align to protein coding transcriptome
for filename in $Totals_filenames
do
bowtie2 -S $SAM_dir/${filename}_pc.sam -U $fastq_dir/${filename}_UMI_clipped.fastq -x $rsem_index --threads $threadN --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 2> $log_dir/${filename}_pc_log.txt
done

#Convert sam to bam 
for filename in $Totals_filenames
do
samtools view -b $SAM_dir/${filename}_pc.sam > $BAM_dir/${filename}_pc.bam &
done
wait

#Sort bam
for filename in $Totals_filenames
do
samtools sort $BAM_dir/${filename}_pc.bam -o $BAM_dir/${filename}_pc_sorted.bam -@ $threadN -m 1G
done

#Index bam
for filename in $Totals_filenames
do
samtools index $BAM_dir/${filename}_pc_sorted.bam $BAM_dir/${filename}_pc_sorted.bai &
done
wait

#Run UMI tools deduplication function
for filename in $Totals_filenames
do
umi_tools dedup -I $BAM_dir/${filename}_pc_sorted.bam -S $BAM_dir/${filename}_pc_deduplicated.bam --output-stats=$log_dir/${filename}_deduplication 1> $log_dir/${filename}_deduplication_log.txt &
done
wait

#Run RSEM to quantify gene and isoform level expression
for filename in $Totals_filenames
do
rsem-calculate-expression --strandedness forward --fragment-length-mean 300 --fragment-length-sd 100 --alignments $BAM_dir/${filename}_pc_deduplicated.bam $rsem_index $rsem_dir/${filename} &
done
wait
