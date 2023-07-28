# Dependencies
Use conda to install the following dependencies

Install conda by following the instuctions here:
[conda installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html)

For more information on conda, see [here](https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307) and [Conda cheat sheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)

**Make a seperate environment for processing the RPFs and for processing the total RNA-seq data**

## RiboSeq environment
The RiboSeq environment is for processing RPFs and requires the following programs to be installed
#### fastQC [manual](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
#### cutadapt [manual](https://cutadapt.readthedocs.io/en/stable/guide.html)
#### UMItools [manual](https://umi-tools.readthedocs.io/en/latest/QUICK_START.html)
#### bbmap [manual](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
#### SAMtools (this has to be v1.9) [manual](http://www.htslib.org/doc/samtools.html)

It also requires the python packaes [pysam](https://github.com/pysam-developers/pysam) and [biopython](https://biopython.org/)

This environment can be created as follows
**Do not copy and paste multiple lines into the terminal at once as some commands will require you to follow prompts, i.e. typing y to proceed with installation**

```console
conda create --name RiboSeq
conda activate RiboSeq
conda install -c bioconda fastqc
conda install -c bioconda cutadapt
conda install -c bioconda umi_tools
conda install -c bioconda bbmap
conda install -c bioconda samtools=1.9
conda install -c bioconda pysam
conda install -c anaconda biopython
conda deactivate
```

## RNAseq environment
The RNAseq environment is for processing totals and requires the following programs to be installed. This environment can also be used to process standard RNA-seq data
#### RSEM [manual](https://deweylab.github.io/RSEM/README.html)
RSEM also requires either bowtie, bowtie2 or STAR to align the reads. We use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

Installing RSEM with conda will also install samtools, however you need to force it to install version 1.9 to avoid getting an error while loading shared libraries: libcrypto.so.1.0.0

There was an issue when installing bowtie2 with conda that required downgrading the tbb package [see here](https://www.biostars.org/p/494922/)

```console
conda create --name RNAseq
conda activate RNAseq
conda install -c bioconda fastqc
conda install -c bioconda cutadapt
conda install -c bioconda umi_tools
conda install -c bioconda rsem
conda install -c bioconda samtools=1.9 --force-reinstall
conda install -c bioconda bowtie2
conda install tbb=2020.2
conda install -c bioconda bbmap
conda install -c anaconda biopython
conda deactivate
```



