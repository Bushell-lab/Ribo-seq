#!/usr/bin/env bash

###filenames
#These are the filenames for all RPF and Total RNA-seq samples (without the <.fastq> or any alternative extension)

RPF_filenames='RPF_1 RPF_2 RPF_3'
Totals_filenames='Totals_1 Totals_2 Totals_3'

###set the sumber of threads available to use
#It is recomended to use one or two less than what is available and also consider whether any else is being run at the same time
#Some of the packages used do not support multi-threading and so the for loops run in parallel, so that all files are run at the same time. For this reason do not run on more files than the number of cores available to use
threadN=16

###adaptors
#RPF adaptors
#This is the sequence of the 3' adaptor that was used in the library prep. Common sequences are below, unhash the correct one if present, or if not enter it as a variable

RPF_adaptor='TGGAATTCTCGGGTGCCAAGG' #this is the adaptor used in the nextflex small RNA library kit
#RPF_adaptor='CTGTAGGCACCATCAAT' #this is the adaptor that seems to have been more commonly used in older ribosome-footprinting studies such as Wolfe 2014
#RPF_adaptor='AGATCGGAAGAGCAC' #this is the one stated in the McGlincy and Ingolia 2017 methods paper
#RPF_adaptor=''

#Totals adaptors
Totals_adaptor='AGATCGGAAGAG' #this is the adaptor used in the LEXOGEN CORALL Total RNA-Seq Library Prep Kit

###paths
parent_dir='/Path/to/data' #This is the path to the parent directory that contains all the data and where all the processed data will be saved

#The following directories are where all the processed data will be saved. These all need to be created prior to starting the analysis

#set the directory where the raw bcl data is. the directory that contains the raw sequencing data in bcl format. This is what you get from a sequencing run and needs to be demulitplexed to write the <.fastq> files.
#If you have more than one bcl directory (you will get one for each sequencing run), then hash one out and write a new one below, each time you re-run the demultiplex.sh script script, so that this acts as a log for all the bcl directories associated with this project
bcl_dir='Path/to/bcl/data'

fastq_dir=${parent_dir}/fastq_files
fastqc_dir=${parent_dir}/fastQC_files
SAM_dir=${parent_dir}/SAM_files
BAM_dir=${parent_dir}/BAM_files
log_dir=${parent_dir}/logs
counts_dir=${parent_dir}/Counts_files
csv_counts_dir=${parent_dir}/Counts_files/csv_files
csv_R_objects=${parent_dir}/Counts_files/R_objects

STAR_dir=${parent_dir}/STAR
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
reads_summary_dir=${analysis_dir}/reads_summary
fgsea_dir=${analysis_dir}/fgsea

#The following directories are where all the plots generated in R will be saved
plots_dir=${parent_dir}/plots

summed_counts_plots_dir=${plots_dir}/summed_counts
periodicity_plots_dir=${plots_dir}/periodicity
offset_plots_dir=${plots_dir}/offset
heatmaps_plots_dir=${plots_dir}/heatmaps
DE_analysis_dir=${plots_dir}/DE_analysis
PCA_dir=${plots_dir}/PCAs
Interactive_scatters_dir=${plots_dir}/Interactive_scatters
fgsea_plots_dir=${plots_dir}/fgsea
fgsea_scatters_dir=${plots_dir}/fgsea/scatters
fgsea_interactive_scatters_dir=${plots_dir}/fgsea/Interactive_scatters
read_counts_summary_dir=${plots_dir}/read_counts_summary
binned_plots_dir=${plots_dir}/binned_plots
single_transcript_binned_plots_dir=${plots_dir}/binned_plots/single_transcripts
normalisation_binned_plots_dir=${plots_dir}/binned_plots/normalisation


#Fastas
fasta_dir='/Path/to/FASTAs'

rRNA_fasta=${fasta_dir}/rRNA/human_rRNA.fa #this needs to point to a fasta file containing rRNA sequences for the correct species
tRNA_fasta=${fasta_dir}/tRNA/human_mature_tRNA.fa #this needs to point to a fasta file containing tRNA sequences for the correct species
pc_fasta=${fasta_dir}/GENCODE/v38/filtered/gencode.v38.pc_transcripts_filtered.fa #this needs to point to a protein coding fasta. See GitHub README file for more information on what is most recommended
rsem_index=${fasta_dir}/GENCODE/v38/filtered/rsem_bowtie2_index/gencode.v38.pc_transcripts_filtered #this needs to point to a index that has been generated for alignment, that is also compatible for RSEM usage. Bowtie2 is recommended for this
STAR_index=${fasta_dir}/GENCODE/v38/original/STAR_index #This needs to point to a STAR genome index that needs to have been previously created
STAR_GTF=${fasta_dir}/GENCODE/v38/original/gencode.v38.annotation.gtf #This needs to point to the GTF file used to create the STAR index
most_abundant_fasta=$most_abundant_transcripts_dir/most_abundant_transcripts.fa #this needs to be created for each specific project (see GitHub README file for more information)

###fasta info
#The below needs to point to a <.csv> file that contains the following information for all transcripts within the protein coding FASTA
#transcript_ID,5'UTR length,CDS length,3'UTR length
#Running the Filter_GENCODE_FASTA.py script will generate this file as one of its outputs

region_lengths=${fasta_dir}/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_region_lengths.csv

