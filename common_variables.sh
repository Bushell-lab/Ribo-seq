#!/usr/bin/env bash

#filenames
filenames='R11_test_RPF_1 R11_test_RPF_2'

#adaptors
adaptor='TGGAATTCTCGGGTGCCAAGG'

#paths
parent_dir='/home/local/BICR/jwaldron/data/JWALDRON/Scripts/python/Ribosome_profiling/test_data'

fastq_dir=${parent_dir}/fastq_files
fastqc_dir=${parent_dir}/fastQC_files
SAM_dir=${parent_dir}/SAM_files
BAM_dir=${parent_dir}/BAM_files
log_dir=${parent_dir}/logs
counts_dir=${parent_dir}/Counts_files

#Fastas
fasta_dir='/home/local/BICR/jwaldron/data/JWALDRON/Indexes/mouse'

rRNA_fasta=${fasta_dir}/rRNA/mouse_rRNA.fa
tRNA_fasta=${fasta_dir}/tRNA/mouse_tRNA.fa
mito_fasta=${fasta_dir}/GENCODE/vM27/gencode.vM27.mito_transcripts.fa
pc_fasta=${fasta_dir}/GENCODE/vM27/gencode.vM27.pc_transcripts_filtered.fa

#fasta info
region_lengths=${fasta_dir}/GENCODE/vM27/gencode.vM27.pc_transcripts_region_lengths.csv

