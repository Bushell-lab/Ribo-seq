# Ribo-seq
Analysis pipeline for Ribosome footprinting data. **Please ensure you read the following instructions carefully before running this pipeline.**

## Dependencies
Make sure you have all dependencies installed by following the instructrions [here](https://github.com/Bushell-lab/Ribo-seq/tree/main/Installation)

## Pipeline
This pipeline uses custom shell scripts to call several external programs as well as custom python scripts to process the data. The processed data can then be used as input into the custom R scripts to either generate library QC plots, perform differential expression (DE) analysis with DEseq2 or to determine codon elongation rates or ribosome pause sites

### Shell scripts
The shell scripts <.sh> are designed to serve as a template for processing the data but require the user to modify them so that they are specific to each experiment. This is because different library prep techniques are often used, such as the use of different **adaptor sequences** or the use of **UMIs**. It is therefore essential that the user knows how the libraries were prepared before starting the analysis. If this is there own data then this should already be known, but for external data sets this is often not as straight forward as expected. Also, it can be unclear whether the data uploaded to GEO is the raw unprocessed <.fastq> files or whether initial steps such as adaptor removal or de-duplication and UMI removal have already been carried out. This is why each processing step is carried out seperately and why the output from these steps is checked with fastQC so it is essential that the user checks these files by viusal inspection after each step, so that the user can be certain that the step has processed the data as expected. Each shell script has annotation at the top describing what the script does and which parts might need editing. **It is therefore strongly recommended that the user opens up each shell script and reads and understands them fully before running them**

### R scripts
The R scripts will read in the processed data generated from the custom python scripts and generate plots and perform DE analysis. These shouldn't need to be edited as the final processed data should be in a standard format, although the user is free to do what they wish with these and change the styles of the plots or add further analyses/plots should they wish. The common_variables.R (see below) script will need to be edited to add the filenames and path to the parent directory as well as the read lengths that they wish to inspect with the library QC plots.

### Python scripts
The python scripts <.py> should not need to be edited. These can be used for multiple projects and so it is recommended that these are placed in a seperate directory. If you set the $PATH to this directory, they can be called from any directory and therefore be used for all Ribo-seq analyses. 

## Setting up the project
- Before running any scripts, create a new directory for the experiment. This will be the parent directory which will contain all raw and processed data for the experiment as well as any plots generated.
- Then create a folder within this directory to place all the shell <.sh> and R <.R> scripts from this GitHub repository. Ensure the $PATH is set to the directory cotaining all the <.py> scripts from this repository
- There is a common_variables.sh and common_variables.R script that will both need to be edited before running any of the other scripts. The filenames for both the RPF and Totals <.fastq> files (without the <.fastq> extension) will need to be added as well as the path to the parent directory. The path to the FASTA files and RSEM index which will be used for mapping also need to be added. These should be common between all projects from the same species so should be stored in a seperate directory from the project directory.
- A region lengths <.csv> file that contains ***transcript_ID, 5'UTR length, CDS length, 3'UTR length*** in that stated order without a header line, for all transcripts within the protein coding FASTA, also needs to be created and the path to this file should also be stated in the common_variables.sh script.
- Once the common_variables.sh script has been completed, run the makeDirs.sh to set up the file structure to store all raw and processed data within the parent directory. Alternativly you can create these directories manually without the command line.
- **It is highly recommended that this data structure is followed as the scripts are designed to output the data in these locations and this makes it much easier for other people to understand what has been done and improves traceability. The filenames are also automatically generated within each script and should contain all important information. Again it is highly recommended that this is not altered for the same reasoning.**

## Processing RPFs
### Sequencing Quality Checks
Before processing any data it is important to use fastQC to see what the structure of the sequencing reads is.

**RPFs_0_fastQC.sh** will run fastQC on all RPF <.fastq> files and output the fastQC files into the fastQC directory. It will do this in parallel to use the maximum number of cores to process all files at the same time, greatly improving the speed.

The output will tell you the number of reads for each <.fastq> file as well as some basic QC on the reads.

A good indication of whether the <.fastq> files have already been processed or not is the sequence length distribution. If no processing has been done, then all reads should be the same length, which will be the number of cycles used when sequenced. For example, if 75 cycles were selected when setting up the sequencing run, then all reads would be 75 bases long, even if the library fragment length was much shorter or much longer than this. Therefore, if for example with a standrad RPF library with 4nt UMIs on either end of the RPF, the fragment length will be roughly 30nt (RPF length) plus 8nt (UMIs) plus the length of the 3' adaptor. If the adaptors had already been removed prior to uploading the <.fastq> files to GEO, then the sequence length distribution will be a range of values, peaking at roughly 38. If the peak was closer to 30nt then it could be presumed that the UMIs had also been removed. The adaptor content will also give a good indication of this, as in the above example, if adaptors hadn't been removed, you should expect to see adaptor contamination coming up in the reads from roughly 38nts into the reads.


# Common troubleshooting
### remove \r end lines
### check the path to directories is right. Go to that directory in the terminal and use pwd to see what the full path is and then copy this into the scripts
