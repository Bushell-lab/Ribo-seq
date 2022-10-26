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

#run blc2fastq
bcl2fastq -p $threadN --no-lane-splitting --runfolder-dir $bcl_dir --output-dir $fastq_dir --barcode-mismatches 0
