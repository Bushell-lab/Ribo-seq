#!/usr/bin/env bash

###filenames
#These are the filenames for all RPF and Total RNA-seq samples (without the <.fastq> or any alternative extension)

RPF_filenames='RPF_1 RPF_2 RPF_3'
Totals_filenames='Totals_1 Totals_2 Totals_3'

###adaptor
#This is the sequence of the 3' adaptor that was used in the library prep. Common sequences are below, unhash the correct one if present, or if not enter it as a variable

#RPF adaptors
RPF_adaptor='TGGAATTCTCGGGTGCCAAGG' #this is the adaptor used in the nextflex small RNA library kit
#RPF_adaptor='CTGTAGGCACCATCAAT' #this is the adaptor that seems to have been more commonly used in older ribosome-footprinting studies such as Wolfe 2014
#RPF_adaptor='AGATCGGAAGAGCAC' #this is the one stated in the McGlincy and Ingolia 2017 methods paper
#RPF_adaptor=''

#Totals adaptors
Totals_adaptor='AGATCGGAAGAG' #this is the adaptor used in the LEXOGEN CORALL Total RNA-Seq Library Prep Kit

###paths

parent_dir='/Path/to/data' #This is the path to the parent directory that contains all the data and where all the processed data will be saved

#The following directories are where all the processed data will be saved. These all need to be created prior to starting the analysis

fastq_dir=${parent_dir}/fastq_files
fastqc_dir=${parent_dir}/fastQC_files
SAM_dir=${parent_dir}/SAM_files
BAM_dir=${parent_dir}/BAM_files
log_dir=${parent_dir}/logs
counts_dir=${parent_dir}/Counts_files
csv_counts_dir=${parent_dir}/Counts_files/csv_files

rsem_dir=${parent_dir}/rsem

#The following directories are where all the csv files that are used as input into R will be saved

analysis_dir=${parent_dir}/Analysis

region_counts_dir=${analysis_dir}/region_counts
spliced_counts_dir=${analysis_dir}/spliced_counts
periodicity_dir=${analysis_dir}/periodicity
cds_counts_dir=${analysis_dir}/CDS_counts
UTR5_counts_dir=$analysis_dir/UTR5_counts
codon_counts_dir=${analysis_dir}/codon_counts
most_abundant_transcripts_dir=${analysis_dir}/most_abundant_transcripts
DESeq2_dir=${analysis_dir}/DESeq2_output

#The following directories are where all the plots generated in R will be saved
plots_dir=${parent_dir}/plots
summed_counts_plots_dir=${plots_dir}/summed_counts
periodicity_plots_dir=${plots_dir}/periodicity
offset_plots_dir=${plots_dir}/offset
heatmaps_plots_dir=${plots_dir}/heatmaps
DE_analysis_dir=${plots_dir}/DE_analysis
PCA_dir=${plots_dir}/PCAs
Interactive_scatters_dir=${plots_dir}/Interactive_scatters
fgsea_dir=${plots_dir}/fgsea
fgsea_scatters_dir=${plots_dir}/fgsea/scatters
fgsea_interactive_scatters_dir=${plots_dir}/fgsea/Interactive_scatters

#Fastas
fasta_dir='/Path/to/FASTAs'

rRNA_fasta=${fasta_dir}/rRNA/human_rRNA.fa
tRNA_fasta=${fasta_dir}/tRNA/human_mature_tRNA.fa
mito_fasta=${fasta_dir}/GENCODE/v38/filtered/gencode.v38.mito_transcripts.fa
pc_fasta=${fasta_dir}/GENCODE/v38/filtered/gencode.v38.pc_transcripts_filtered.fa
rsem_index=${fasta_dir}/GENCODE/v38/filtered/rsem_bowtie2_index/gencode.v38.pc_transcripts_filtered
STAR_index=${fasta_dir}/GENCODE/v38/original/STAR_index
STAR_GTF=${fasta_dir}/GENCODE/v38/original/gencode.v38.annotation.gtf
most_abundant_fasta=$most_abundant_transcripts_dir/most_abundant_transcripts.fa #this needs to be created for each specific project

###fasta info
#The below needs to point to a <.csv> file that contains the following information for all transcripts within the protein coding FASTA
#transcript_ID,5'UTR length,CDS length,3'UTR length
#Running the Filter_GENCODE_FASTA.py script will generate this file as one of its outputs

region_lengths=${fasta_dir}/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_region_lengths.csv

