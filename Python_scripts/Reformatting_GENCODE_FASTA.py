#!/usr/bin/env python

#Imports
import re
import argparse
from Bio import SeqIO

#Functions
def read_in_fasta(afasta):
    '''Reads in a fasta file to a dictionary'''
    fasta_dict = {}
    fasta_sequences = SeqIO.parse(open(afasta),'fasta')
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq).upper()
    return fasta_dict

def reformat_fasta(fasta_transcripts_in, fasta_translations_in):
    '''Uses transcript and translation fasta dicts to return a five dicts, one for each fasta with just the transcript or protein ID as key, one for region lengths and one for gene and protein IDs'''
        
    fasta_transcripts_out, fasta_translations_out, region_lens, gene_IDs, protein_IDs = {},{},{},{},{}
    
    for k,v in fasta_transcripts_in.items():
        
        #extract transcript/gene IDs
        transcript_ID = k.split('|')[0]
        gene_ID = k.split('|')[1]
        gene_sym = k.split('|')[5]
        
        #extract 5'UTR length
        UTR5_search = re.findall(r"\|\UTR5:\d+\-\d+\|",k)
        if UTR5_search != []:
            UTR5_len = int(UTR5_search[0].strip('|').split('-')[1])
        else:
            UTR5_len = 0
        
        #extract CDS length
        cds_search = re.findall(r"\|\CDS:\d+\-\d+\|",k)
        cds_range = cds_search[0].split(':')[1].strip('|')
        cds_len = int(cds_range.split('-')[1]) - int(cds_range.split('-')[0]) + 1
        
        #extract 3'UTR length
        UTR3_search = re.findall(r"\|\UTR3:\d+\-\d+\|",k)
        if UTR3_search != []:
            UTR3_range = UTR3_search[0].split(':')[1].strip('|')
            UTR3_len = int(UTR3_range.split('-')[1]) - int(UTR3_range.split('-')[0]) + 1
        else:
            UTR3_len = 0
        
        total_length = int(re.findall(r"\|\d+\|",k)[0].strip('|'))
        if total_length == UTR5_len + cds_len + UTR3_len: #sanity check that nothing has gone wrong
            #save sequence to fasta dict
            fasta_transcripts_out[transcript_ID] = v
            #save region lengths and gene/transcript IDs to dicts
            region_lens[transcript_ID] = (str(UTR5_len), str(cds_len), str(UTR3_len))
            gene_IDs[transcript_ID] = [gene_ID, gene_sym]
    
    #extract protein IDs
    for k,v in fasta_translations_in.items():
        protein_ID = k.split('|')[0]
        transcript_ID = k.split('|')[1]
        protein_IDs[protein_ID] = [transcript_ID] + gene_IDs[transcript_ID]
        
        #save sequence to fasta dict
        fasta_translations_out[protein_ID] = v
                                    
    return (fasta_transcripts_out, fasta_translations_out, region_lens, gene_IDs, protein_IDs)
    
def write_fasta(adict, outfyle, LW=80):
    '''takes a fasta dictionary and writes a fasta'''
    with open(outfyle, 'w') as g:
        for k,v in adict.items():
            g.write('>' + k + '\n')
            for i in range(0, len(v), LW):
                g.write(v[i:i+LW] + '\n')

def write_csv(adict, outfyle):
    '''takes a dict and writes to csv'''
    with open(outfyle, 'w') as g:
        for k,v in adict.items():
            g.write(k + ',' + ','.join(v) + '\n')
    
#Main Function
def main():
    parser = argparse.ArgumentParser(description='Takes a GENCODE FASTA file containing transcript sequences and a \
    GENCODE FASTA file containing amino acid sequences and reformats them so that the header line is just the transcript or protein ID \
    and exporting the UTR and CDS region lengths and gene/transcript/protein IDs to seperate csv files')
    parser.add_argument('FASTA_transcripts', type=str, help='GENCODE FASTA file containing transcript sequences')
    parser.add_argument('FASTA_translations', type=str, help='GENCODE FASTA file containing amino acid sequences')
    args = parser.parse_args()
    
    #read in FASTA
    original_fasta_transcripts = read_in_fasta(args.FASTA_transcripts)
    original_fasta_translations = read_in_fasta(args.FASTA_translations)
    
    #reformat FASTAs
    final_fasta_transcripts, final_fasta_translations, region_lens, gene_IDs, protein_IDs = reformat_fasta(original_fasta_transcripts, original_fasta_translations)
    
    #write FASTAs
    fasta_fylename = args.FASTA_transcripts.replace('.fa', '_reformatted.fa')
    write_fasta(final_fasta_transcripts, fasta_fylename)
    
    fasta_fylename = args.FASTA_translations.replace('.fa', '_reformatted.fa')
    write_fasta(final_fasta_translations, fasta_fylename)
    
    #write CSVs
    write_csv(region_lens, args.FASTA_transcripts.replace('.fa', '_region_lengths.csv'))
    write_csv(gene_IDs, args.FASTA_transcripts.replace('.fa', '_gene_IDs.csv'))
    write_csv(protein_IDs, args.FASTA_transcripts.replace('.fa', '_protein_IDs.csv'))

if __name__ == '__main__': 
    main()
 
