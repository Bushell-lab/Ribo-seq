#!/usr/bin/env bash

#read in variables
source common_variables.sh

#uses the flat text file containing the most abundant transcripts per gene created with calculate_most_abundant_transcript.R to filter the protein coding fasta

filter_FASTA.py $pc_fasta $most_abundant_transcripts_dir/most_abundant_transcripts.txt $most_abundant_fasta


