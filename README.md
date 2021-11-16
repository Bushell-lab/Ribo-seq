# Ribo-seq
This analysis pipeline processes sequencing data from both Ribosome Protected Footprints (RPFs) and the associated total RNA-seq for Ribosome-footprinting data.

**See the RiboSeq.png figure for an overview of the pipeline.**

**Please ensure you read the following instructions carefully before running this pipeline.**

## Dependencies
Make sure you have all dependencies installed by following the instructrions [here](https://github.com/Bushell-lab/Ribo-seq/tree/main/Installation)

## Pipeline
This pipeline uses custom shell scripts to call external programs and custom python scripts, to process the data. The processed data can then be used as input into the custom R scripts to either generate library QC plots, perform differential expression (DE) analysis with DEseq2 or to determine codon elongation rates or ribosome pause sites.

### Shell scripts
The shell scripts <.sh> are designed to serve as a template for processing the data but require the user to modify them so that they are specific to each experiment. This is because different library prep techniques are often used, such as the use of different **adaptor sequences** or the use of **UMIs**. It is therefore essential that the user knows how the libraries were prepared before starting the analysis. If this is their own data then this should already be known, but for external data sets this is often not as straight forward as expected. Also, it can be unclear whether the data uploaded to GEO is the raw unprocessed <.fastq> files or whether initial steps such as adaptor removal or de-duplication and UMI removal have already been carried out. This is why each processing step is carried out seperately and why the output from these steps is checked with fastQC, so it is essential that the user checks these files by visual inspection after each step, so that the user can be certain that the step has processed the data as expected. Each shell script has annotation at the top describing what the script does and which parts might need editing. **It is therefore strongly recommended that the user opens up each shell script and reads and understands them fully before running them**

### R scripts
The R scripts will read in the processed data generated from the custom python scripts and generate plots and perform DE analysis. These shouldn't need to be edited as the final processed data should be in a standard format, although the user is free to do what they wish with these and change the styles of the plots or add further analyses/plots should they wish. The common_variables.R (see below) script will need to be edited to add the filenames and path to the parent directory, as well as the read lengths that they wish to inspect with the library QC plots. The common_variables.R script needs to be in the current working directory when running the other <.R> scripts.

### Python scripts
The <.py> python scripts should not need to be edited. These can be used for multiple projects and so it is recommended that these are placed in a seperate directory. If you set the $PATH to this directory, they can be called from any directory and therefore be used for all Ribo-seq analyses.

## Setting up the project
- Before running any scripts, create a new directory for the experiment. This will be the parent directory which will contain all raw and processed data for the experiment as well as any plots generated.
- Then create a folder within this directory to place all the shell <.sh> and R <.R> scripts from this GitHub repository. Ensure the $PATH is set to the directory cotaining all the <.py> scripts from this repository.
- There is a common_variables.sh and common_variables.R script that will both need to be edited before running any of the other scripts. The filenames for both the RPF and Totals <.fastq> files (without the <.fastq> extension) will need to be added as well as the path to the parent directory and the adaptor sequences. The path to the FASTA files and RSEM index which will be used for mapping also needs to be added. These should be common between all projects from the same species so should be stored in a seperate directory from the project directory.
- A region lengths <.csv> file that contains ***transcript_ID, 5'UTR length, CDS length, 3'UTR length*** in that stated order without a header line, for all transcripts within the protein coding FASTA, also needs to be created and the path to this file should also be stated in the common_variables.sh script.
- Once the common_variables.sh script has been completed, run the makeDirs.sh to set up the file structure to store all raw and processed data within the parent directory. Alternativly you can create these directories manually without the command line.
- **It is highly recommended that this data structure is followed as the scripts are designed to output the data in these locations and this makes it much easier for other people to understand what has been done and improves traceability. The filenames are also automatically generated within each script and should contain all important information. Again it is highly recommended that this is not altered for the same reasoning.**
- Once the directories have been set up, place all the raw <.fastq> files in the fastq_files directory. These will be the input into the *RPFs_0_QC.sh* and *RPFs_1_adaptor_removal.sh* scripts so check that that extensiones match. It is fine if these files a <.gz> compressed as both fastQC and cutadapt can use compressed files as input, but again make sure that the shell scripts have the .gz extension included.

## Processing RPFs
**Ensure you activate the RiboSeq conda environment before running the RPF shell scripts, with the following command**
```console
conda activate RiboSeq
```
### Sequencing QC
Before processing any data it is important to use fastQC to see what the structure of the sequencing reads is.

**RPFs_0_QC.sh** will run fastQC on all RPF <.fastq> files and output the fastQC files into the fastQC directory.

The output will tell you the number of reads for each <.fastq> file as well as some basic QC on the reads.

A good indication of whether the <.fastq> files have already been processed or not is the sequence length distribution. If no processing has been done, then all reads should be the same length, which will be the number of cycles used when sequenced. For example, if 75 cycles were selected when setting up the sequencing run, then all reads would be 75 bases long, even if the library fragment length was much shorter or much longer than this. Therefore, if for example with a standrad RPF library with 4nt UMIs on either end of the RPF, the fragment length will be roughly 30nt (RPF length) plus 8nt (UMIs) plus the length of the 3' adaptor. If the adaptors had already been removed prior to uploading the <.fastq> files to GEO, then the sequence length distribution will be a range of values, peaking at roughly 38. If the peak was closer to 30nt then it could be presumed that the UMIs had also been removed. The adaptor content will also give a good indication of this. For example, in the above example, if adaptors hadn't been removed, you should expect to see adaptor contamination coming up in the reads from roughly 38nts into the reads.

### Remove adaptors
The 3' adaptor used in the library prep will be sequenced immediately after the fragment (and UMIs if used). These therefore needs to be removed so that they do not affect alignment. The ***RPFs_1_adaptor_removal.sh*** script uses cutadapt for this, which removes this sequence (specified in the common_variables.sh script) and any sequence downstream of this. It also trims low quality bases from the 3' end of the read below a certain quality score (user defined, we use q20) and removes reads that are shorter or longer than user defined values. For RPFs (~30) with 4nt UMIs at each end, we filter reads so that they are 30-50nt. If UMIs have been used then you need to set the minimum read length to 30nt as otherwise cd-hit-dup has issues de-duplicating the reads. If UMIs haven't been used, you will need to change these settings to 20-40.

After cutadapt has finished, fastQC is run on the output <.fastq> files. **Visual inspection of these fastQC files is essential to check that cutadapt has done what you think it has**

### De-duplication and UMI removal
If UMIs have been used in the library prep, reads that are PCR duplicates can be removed from the <.fastq> file, ensuring that all remaining reads originated from unique RPFs. The ***RPFs_2_deduplication.sh*** script uses cd-hit-dup to make a new <.fastq> file containing only unique reads.

After cd-hit-dup has finished, fastQC is run on the output <.fastq> files. **Visual inspection of these fastQC files is essential to check that cd-hit-dup has done what you think it has.**

The UMIs can now be removed. The ***RPFs_3_UMI_removal.sh*** script uses cutadapt to remove a set number of bases from the 5' and 3' end of all reads. This needs to be set to match the structure of the UMIs used. For the nextflex library prep kit that we use for RPFs, these are 4nt at either end of the read.

After cutadpat has finished, fastQC is run on the output <.fastq> files. **Visual inspection of these fastQC files is essential to check that cutadapt has done what you think it has.**

**If the library prep did not include UMIs then *RPFs_2_deduplication.sh* and *RPFs_3_UMI_removal.sh* should be skipped. If this is the case you need to edit the names of the input <.fastq> files in the** ***RPFs_4_align_reads.sh*** **script to the names of the output <.fastq> files from the** ***RPFs_1_adaptor_removal.sh*** **script.**

### Read alignment
The processing of the reads up to this point should have removed any sequences introduced during the library prep, and any PCR duplicates. This means the reads should reflect the exact sequences of the extracted RNA and so can now be aligned to a transcriptome.

**It is very important to give some consideration to what transcriptome you use and how to handle multimapped reads. It is strongly recommended that you use the gencode protein coding transcriptome that has been filtered to include**
- only Havana protein coding transcripts
- that have both 5' and 3'UTRs
- which the CDS is equally divisble by 3, starts with a start codon and finishes with a stop codon

The ***RPFs_4_align_reads.sh*** script uses bbmap to align reads first to the rRNAs, tRNAs and mitochondrial mRNAs (will be filtered from the above fasta due to lack of UTRs). Each alignment will create two new <.fastq> files containing the reads that did and did not align, as well as a <.SAM> file containing the alignments. The reads that didn't align to either of these transcriptomes are then aligned to the protein coding transcriptome, firstly keeping all alignments (including multi-mapped reads) and then keeping only the highest scoring alignment for any mulit-mapped reads.

fastQC is then used to inspect the QC and read length distribution of these different alignments. You would expect to see a nice peak of read lengths 28-30nt for the protein coding aligned reads but a wider distribition of reads for the rRNA (this should reflect the size you cut on the RNA extraction gel).

### SAM to BAM
Sequence alignments are written to <.SAM> files. These can be compressed to <.BAM> files which take up less space and can also be sorted and indexed. The ***RPFs_5_SAM_to_BAM.sh*** script will do exactly this.

### Count reads
In order to do any downstream analysis, we need to know how many reads aligned to which mRNAs at which positions. Also for library QC it is important to be able to distinguish between different read lengths, as certain read lengths may be filtered to remove those reads that are less likely to be true RPFs.

The ***counting_script.py*** script was adpated from the [RiboPlot package](https://pythonhosted.org/riboplot/ribocount.html). This script creates <.counts> files, which are plain text files in the following structure;

Transcript_1

0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 150 20 34 85 34 58 75 22 27 85 53 24 85.....................................................

Transcript_2

0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 103 10 37 83 24 57 45 28 7 89 43 26 55.....................................................

Transcript_3

0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 87 20 34 85 34 48 79 12 19 75 51 14 95.....................................................


where each two lines represents one transcript, with the first of each two lines containg the transcript ID and the second of each two lines containing tab seperated values of the read counts that start at that position within the transcript. The number of values for each transcript should therefore match the length of that transcript.

The ***RPFs_6_Extract_counts_all_lengths.sh*** uses the ***counting_script.py*** script to generate a <.counts> file for each sample for each read length defined in the for loop and stores all thes files in the Counts_files directory. This uses the sorted <.BAM> file as input and also needs the associated index <.BAI> file to be in the same directory. **The script uses the best alignments so every read only counts once but that all transcripts can be considered.** These are the input files for the downstream analysis.

### Library QC
The ***RPFs_7a_summing_region_counts.sh; RPFs_7b_summing_spliced_counts.sh and RPFs_7c_periodicity.sh*** scripts utlise the custom python scripts, reading in the counts files generated above and creating <.csv> files that the ***region_counts.R; heatmaps.R; offset_plots.R and periodicity.R*** scripts use to generate the library QC plots. From these plots you should be able to determine whether the RPF libraries have the properties that would argue they are truelly RPFs. These are;
- read length distribution peaking at 28-30nt
- strong periodicity
- strong enrichment of reads within the CDS and depletion of reads within the 3'UTR

From these plots you can then determine what read lengths you want to include in your downstream analysis for DE and codon level analyses.

The offset plots should also allow you to determine what to use for the offset. This is the value that you use in the ***RPFs_8_Extract_final_counts.sh*** so that the counts in the final <.counts> files are referring to the position within each transcript which is the first nt of the codon positioned within the P-site of the ribosome which was protecting that RPF, rather than the start of the read. This can be determined from the position at which the first peak of reads is observed just upstream of the start codon, as these reads correspond to RPFs which were protected by ribosomes with the P-site situated at the start codon. This is typically 12-13nt, but it is likely that different read lengths will require slightly different offsets.

**Once you are happy that the data has been processed properly you should delete the following intermediary files that are no longer required**
- <cutadapt.fastq> files generated from *RPFs_1_adaptor_removal.sh*
- <cdhitdup.fastq> files generated from *RPFs_2_deduplication.sh*
- <UMI_removed.fastq> files generated from *RPFs_3_UMI_removal.sh*
- protein coding <.sam> files which have already been converted to <.bam> files (keep the rRNA/tRNA/mito <.sam> files in case you use these for anything)

**Do not delete the raw <.fastq> files**

### Extract final counts
Once you know what read lengths and offsets to use, you can use these values with the ***RPFs_8_Extract_final_counts.sh*** script to create a final <.counts> file that contains only the specified read lengths with the specified offsets applied. The script makes a <.counts> file for the best and all alignments seperately. **For DE expression it is then best to use all alignments, and then select the most abundant transcript per gene (based on the total RNA-seq data) during the DEseq2 analysis. For codon level analyses which is done transcriptome wide, it is best to use the best mappings.**

### Summing CDS counts
The ***RPFs_9a_CDS_counts.sh*** uses the ***summing_CDS_counts.py*** to sum all the read counts that are within the CDS. **This will be the input into DESeq2.**

The *summing_CDS_counts.py* has an option to remove the first 20 and last 10 codons, which is recommended (and set as default in *RPFs_9a_CDS_counts.sh* script) to avoid biases at the start and stop codons, essentially meaning that only activly elongating ribosomes are counted.

There is also the option to only include reads that are in frame. However, although periodicty indicates that the majority of the reads are truely RPFs, it doesn't neccessarily mean that reads that are not in frame are not RPFs and the majority of the reads in the CDS will most likely be RPFs. **It is therefore recommended to include reads in all frames for DE analysis. For codon level analyses, only reads in frame should be used as it is not possible to determine codon level resolution with high confidence for reads not in frame.**

### Counting codon occupancy
The ***RPFs_9b_count_codons.sh*** uses the ***count_codon_occupancy.py*** to determine which codon was positioned at the A,P and E-site plus two codons either side, for every RPF read and sum them all together. The ***codon_occupancy.R*** script then takes this data and uses it to measure relative elongation rates for each codon based on the number of RPFs where that codon was at the A-site compared to the number of RPFs where that codon was at either of the 7 sites described above. This therefore accounts for differing mRNA abundances and initiation rates transcriptome-wide.


# Common troubleshooting
### remove \r end lines
The end of line character for windows is \r but for linux and mac it is \n. Sometimes, when you open a script on your PC in a text editor it will automatically add both \n and \r to the end of any new lines created. However, as the shell scripts are intended to be run on a linux/mac platform, this will cause issues and will return the following error message.

***/usr/bin/env: ‘bash\r’: No such file or directory***

To check if this is the case, in notepad++ select View->Show symbol->Show all characters to see hidden characters. If \r characters have been added to the end of lines, use find and replace (with regular expressions ticked) to remove them all, leaving just \n characters in their place
### check the path to directories is right
The path to the parent directory needs to be set in both the common_variables.sh and the common_variables.R scripts. Although these should point to the same directory, the path will be slightly different as the path in the shell script needs to be the linux path and the path for the R script needs to be the PC path.

- To find the linux path, go to that directory in the terminal and use pwd to see what the full path is and then copy this into the shell script
- To find the PC path, open R studio by doubleclicking on the common_variables.R script and use the getwd() function to find the current working directory. The parent directory will be a couple of directories up from this
