# Dependencies
Use conda to install the following dependencies

Install conda by following the instuctions here:
[conda installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html)

For more information on conda, see [here](https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307) and [Conda cheat sheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)

**Do not copy and paste multiple lines into the terminal at once as some commands will require you to follow prompts, i.e. typing y to proceed with installation**

**Make a seperate environment for processing the RPFs and for processing the total RNA-seq data**

## RiboSeq environment for processing RPFs
The RiboSeq environment requires the following programs to be installed
#### fastQC [manual](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
#### cutadapt [manual](https://cutadapt.readthedocs.io/en/stable/guide.html)
#### cd-hit-dup [manual](https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#cdhitdup)
#### bbmap [manual](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
#### SAMtools [manual](http://www.htslib.org/doc/samtools.html)
```console
$ conda create --name fastQC
$ conda activate fastQC
$ conda install -c bioconda fastqc
$ conda deactivate
```


```console
conda create --name cutadapt
conda activate cutadapt
conda install -c bioconda cutadapt
conda deactivate
```


```console
conda create --name SAMtools
conda activate SAMtools
conda install -c bioconda samtools
conda deactivate
```

### RSEM [manual](https://deweylab.github.io/RSEM/README.html)
RSEM also requires either bowtie, bowtie2 or STAR to align the reads. We use bowtie2 which also needs to be installed within this environment. There was an issue when installing bowtie2 with conda that required downgrading the tbb package [see here](https://www.biostars.org/p/494922/)
```console
conda create --name rsem
conda activate rsem
conda install -c bioconda rsem
conda install -c bioconda bowtie2
conda install tbb=2020.2
conda deactivate
```


```console
conda create --name cdhit
conda activate cdhit
conda install -c bioconda cd-hit-auxtools
conda deactivate
```


```console
conda create --name bbmap
conda activate bbmap
conda install -c bioconda bbmap
conda deactivate
```


