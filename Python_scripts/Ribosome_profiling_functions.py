#contains the commonly used functions for ribosome profiling analysis

from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#functions
def read_counts(counts_fyle):
    '''Reads in RPF counts'''
    adict = {}
    with open(counts_fyle,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,counts = [n.strip() for n in next_n_lines]
            adict[transcript] = counts.split('\t')
    return adict

def read_region_lengths(csv_fyle):
    '''reads in region lengths csv file and saves as a dictionary'''
    adict = {}
    with open(csv_fyle, 'r') as f:
        for line in f:
            ENST = line.strip().split(',')[0]
            UTR5_len = line.strip().split(',')[1]
            cds_len = line.strip().split(',')[2]
            UTR3_len = line.strip().split(',')[3]
            adict[ENST] = (int(UTR5_len), int(cds_len), int(UTR3_len))
    return adict

def read_gene_IDs(csv_fyle):
    '''reads in transcript to gene ID csv file and saves as a dictionary'''
    adict = {}
    with open(csv_fyle, 'r') as f:
        for line in f:
            ENST = line.strip().split(',')[0]
            ENSG = line.strip().split(',')[1]
            gene_sym = line.strip().split(',')[2]
            adict[ENST] = (ENSG, gene_sym)
    return adict

def read_in_fasta(fasta_fyle):
    '''Reads in a fasta file to a dictionary'''
    fasta_dict = {}
    fasta_sequences,fasta_dict = SeqIO.parse(open(fasta_fyle),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict
    


