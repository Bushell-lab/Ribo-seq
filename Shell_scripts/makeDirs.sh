#!/usr/bin/env bash

#read in variables
source common_variables.sh

#make directories
mkdir $fastq_dir
mkdir $fastqc_dir
mkdir $SAM_dir
mkdir $BAM_dir
mkdir $log_dir
mkdir $counts_dir
mkdir $rsem_dir
mkdir $STAR_dir

mkdir $analysis_dir
mkdir $region_counts_dir
mkdir $spliced_counts_dir
mkdir $periodicity_dir
mkdir $cds_counts_dir
mkdir $UTR5_counts_dir
mkdir $csv_counts_dir
mkdir $csv_R_objects
mkdir $codon_counts_dir
mkdir $DESeq2_dir
mkdir $most_abundant_transcripts_dir
mkdir $reads_summary_dir
mkdir $fgsea_dir

mkdir $plots_dir
mkdir $summed_counts_plots_dir
mkdir $periodicity_plots_dir
mkdir $offset_plots_dir
mkdir $heatmaps_plots_dir
mkdir $DE_analysis_dir
mkdir $PCA_dir
mkdir $Interactive_scatters_dir
mkdir $fgsea_plots_dir
mkdir $fgsea_scatters_dir
mkdir $fgsea_interactive_scatters_dir
mkdir $read_counts_summary_dir
mkdir $binned_plots_dir
mkdir $single_transcript_binned_plots_dir
mkdir $normalisation_binned_plots_dir


