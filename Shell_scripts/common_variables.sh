#!/usr/bin/env bash

###filenames
#These are the filenames for all RPF and Total RNA-seq samples (without the <.fastq> or any alternative extension)

RPF_filenames='RPF_1 RPF_2 RPF_3'
Totals_filenames='Totals_1 Totals_2 Totals_3'

###adaptor
#This is the sequence of the 3' adaptor that was used in the library prep. Common sequences are below, unhash the correct one if present, or if not enter it as a variable

#adaptor='TGGAATTCTCGGGTGCCAAGG' #this is the adaptor used in the nextflex small RNA library kit
#adaptor='CTGTAGGCACCATCAAT' #this is the adaptor that seems to have been more commonly used in older ribosome-footprinting studies such as Wolfe 2014
#adaptor='AGATCGGAAGAGCAC' #this is the one stated in the McGlincy and Ingolia 2017 methods paper
#adaptor=''

###paths
#This is the path to the parent directory that contains all the data and where all the processed data will be saved

parent_dir='/Path/to/data'

#The below directories are where all the processed data will be saved. These all need to be created prior to starting the analysis

fastq_dir=${parent_dir}/fastq_files
fastqc_dir=${parent_dir}/fastQC_files
SAM_dir=${parent_dir}/SAM_files
BAM_dir=${parent_dir}/BAM_files
log_dir=${parent_dir}/logs
counts_dir=${parent_dir}/Counts_files
rsem_dir=${parent_dir}/rsem

#The below directories are where all the csv files that are used as input into R will be saved

region_counts_dir=$parent_dir/'Analysis/region_counts'
spliced_counts_dir=$parent_dir/'Analysis/spliced_counts'
periodicity_dir=$parent_dir/'Analysis/periodicity'
cds_counts_dir=$parent_dir/'Analysis/CDS_counts'
codon_counts_dir=$parent_dir/'Analysis/codon_counts'

#Fastas
fasta_dir='/Path/to/FASTAs'

rRNA_fasta=${fasta_dir}/rRNA/human_rRNA.fa
tRNA_fasta=${fasta_dir}/tRNA/human_tRNA.fa
mito_fasta=${fasta_dir}/GENCODE/v38/gencode.v38.mito_transcripts.fa
pc_fasta=${fasta_dir}/GENCODE/v38/gencode.v38.pc_transcripts_filtered.fa

###fasta info
#The below needs to point to a <.csv> file that contains the following information for all transcripts within the protein coding FASTA
#transcript_ID,5'UTR length,CDS length,3'UTR length
#Running the Filter_GENCODE_FASTA.py script will generate this file as one of its outputs

region_lengths=${fasta_dir}/GENCODE/vM27/gencode.v38.pc_transcripts_region_lengths.csv

