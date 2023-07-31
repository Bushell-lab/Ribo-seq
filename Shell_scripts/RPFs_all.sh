#!/usr/bin/env bash

#read in variables
source common_variables.sh

#run cutadapt
for filename in $RPF_filenames
do
cutadapt $fastq_dir/${filename}.fastq -a $RPF_adaptor --nextseq-trim=20 -m 30 -M 50 --cores=0 -o $fastq_dir/${filename}_cutadapt.fastq 1> $log_dir/${filename}_cutadapt_log.txt
done

#extract UMIs
for filename in $RPF_filenames
do
umi_tools extract -I $fastq_dir/${filename}_cutadapt.fastq --extract-method=regex --bc-pattern='^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$' -S $fastq_dir/${filename}_UMI_clipped.fastq --log=$log_dir/${filename}_extracted_UMIs.log &
done
wait

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

#run UMI tools deduplication function
for filename in $RPF_filenames
do
umi_tools dedup -I $BAM_dir/${filename}_pc_sorted.bam -S $BAM_dir/${filename}_pc_deduplicated.bam --output-stats=$log_dir/${filename}_deduplication 1> $log_dir/${filename}_deduplication_log.txt &
done
wait

#sort bam
for filename in $RPF_filenames
do
samtools sort $BAM_dir/${filename}_pc_deduplicated.bam -o $BAM_dir/${filename}_pc_deduplicated_sorted.bam -@ $threadN -m 1G
done

#index bam
for filename in $RPF_filenames
do
samtools index $BAM_dir/${filename}_pc_deduplicated_sorted.bam $BAM_dir/${filename}_pc_deduplicated_sorted.bai &
done
wait

#make an fai (fasta index) file from the fasta using samtools. This is required for the counting script and needs to exist before running counting_script.py
samtools faidx $most_abundant_fasta

#run the counting_script.py with a range of read lengths (adjust below if required, currently set to 25-35)
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
counting_script.py -bam $BAM_dir/${filename}_pc_deduplicated_sorted.bam -fasta $most_abundant_fasta -len $length -out_file ${filename}_pc_L${length}_Off0.counts -out_dir $counts_dir &
done
done
wait

#set offset
offset=15

#run summing_region_counts.py script
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
summing_region_counts.py ${filename}_pc_L${length}_Off0.counts $offset $region_lengths -in_dir $counts_dir -out_dir $region_counts_dir &
done
done
wait

#set number of nt to splice
n=50

#run summing_spliced_counts.py script
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
summing_spliced_counts.py ${filename}_pc_L${length}_Off0.counts $n $region_lengths -in_dir $counts_dir -out_dir $spliced_counts_dir &
done
done
wait

#run periodicity.py script
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
periodicity.py ${filename}_pc_L${length}_Off0.counts $region_lengths -offset $offset -in_dir $counts_dir -out_dir $periodicity_dir &
done
done
wait

#Extract the read counts from the log files for each sample
for filename in $RPF_filenames
do
extract_read_counts.py ${filename} RPFs -log_dir $log_dir &
done
wait
