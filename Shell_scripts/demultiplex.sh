#!/usr/bin/env bash

#This script will demulitplex the bcl data and write a <.fastq> file for each sample. Make sure the bcl_dir and fastq_dir are specified in the common_variables.sh script.
#All raw sequencing data needs to be saved in the R11/raw_sequencing_data. Create your own directory in here and place each directory for each sequencing run (this is the bcl_dir) in that directory.
#bcl2fastq is used to convert bcl files to fastq. This can be downloaded with conda.
#You will need to complete a sampleSheet.csv file for each run, which specifies which samples have which barcodes. This needs to be in bcl_dir.
#-p denotes how many threads to use
#use the --no-lane-splitting option so that you get one fastq file per sample and not four (there are four seperate lanes in the NextSeq, but unless you suspect there have been any sequencing issues these can be combined
#use --barcode-mismatches 0 so that only the extact barcodes are used.

#read in variables
source common_variables.sh

#set the directory where the raw bcl data is.
#If you have more than one bcl directory (you will get one for each sequencing run), then hash one out and write a new one below, each time you re-run this script, so that this acts as a log for all the bcl directories associated with this project

bcl_dir='Path/to/raw_seq_data' #This is the path to the directory that contains the raw sequencing data in bcl format (this is what you get from a sequencing run and needs to be demulitplexed to write the <.fastq> files)

#run blc2fastq
bcl2fastq -p 12 --no-lane-splitting --runfolder-dir $bcl_dir --output-dir $fastq_dir --barcode-mismatches 0
