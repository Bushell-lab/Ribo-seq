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

def extract_HAVANA_IDs(GTFfyle):
    '''Reads a gtf file line by line. Extracts the transcript ID from each feature line if transcript is annotated by HAVANA as being protein coding'''
    HAVANA_IDs_dict = {}
    with open(GTFfyle,'r') as f:
        for line in f:
            if line.startswith('#'):
                pass
            else:
                source = line.split('\t')[1]
                if source == "HAVANA":
                    attribute = line.split('\t')[8]
                    pc_transcript = re.findall(r'transcript_type\ \"protein_coding\"', attribute)
                    if pc_transcript != []:
                        transcript = re.findall(r'\"E\w+T\d+\.\d+\"', attribute)
                        if transcript != []:
                            transcriptID = transcript[0].strip('"')
                            if transcriptID in HAVANA_IDs_dict:
                                pass
                            else:
                                HAVANA_IDs_dict[transcriptID] = None
    return HAVANA_IDs_dict

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

def filter_fasta(fasta_dict, HAVANA_IDs, region_lens):
    '''filters fasta'''
        
    filtered_dict = {}
    
    for transcript,seq in fasta_dict.items():
        if transcript[-5:] != 'PAR_Y': #removes any PAR_Y transcripts
            if transcript in HAVANA_IDs: #Ensures the transcript has been manually annotated by HAVANA
                
                UTR5_len = int(region_lens[transcript][0])
                CDS_len = int(region_lens[transcript][1])
                UTR3_len = int(region_lens[transcript][2])
                
                if UTR5_len > 0 and UTR3_len > 0: #ensures the transcript has both a 5' and 3'UTR
                    if CDS_len % 3 == 0: #ensures the CDS is equally divisible by 3
                        cds_seq = seq[UTR5_len:-UTR3_len]
                        start_codon = cds_seq[:3]
                        stop_codon = cds_seq[-3:]
                        if start_codon == "ATG" or start_codon == "CTG" or start_codon == "TTG" or start_codon == "GTG": #ensures the CDS starts with an nUG start codon
                            if stop_codon == "TAA" or stop_codon == "TGA" or stop_codon == "TAG": #ensures the CDS ends with a stop codon
                                filtered_dict[transcript] = seq
    return (filtered_dict)

def write_fasta(adict, outfyle, LW=80):
    '''takes a fasta dict and writes fasta'''
    with open(outfyle, 'w') as g:
        for k,v in adict.items():
            g.write('>' + k + '\n')
            for i in range(0, len(v), LW):
                g.write(v[i:i+LW] + '\n')

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Takes a fasta file, a GTF file and a region lengths file and filters FASTA to include only those transcripts\
    that have been manually annotated by HAVANA, that have both a 5\' and 3\'UTR and whose CDS is equally divisible by 3 and starts with an ATG and ends with a stop codon')
    parser.add_argument('FASTA', type=str, help='GENCODE FASTA file')
    parser.add_argument('GTF', type=str, help='GENCODE GTF file')
    parser.add_argument('region_lengths_fyle', type=str, help='region lengths file. Needs to be csv file in the following format; Transcript ID,5\'UTR length,CDS length,3\'UTR length')
    args = parser.parse_args()
    
    #read in FASTA
    original_fasta = read_in_fasta(args.FASTA)
    
    #extract HAVANA protein coding transcripts IDs
    HAVANA_IDs = extract_HAVANA_IDs(args.GTF)
    
    #read in the region lengths
    region_lengths = read_region_lengths(args.region_lengths_fyle)
    
    #filter FASTA
    filtered_dict = filter_fasta(original_fasta, HAVANA_IDs, region_lengths)
    
    #write FASTA
    fasta_fylename = args.FASTA.replace('_reformatted', '_filtered')
    write_fasta(filtered_dict, fasta_fylename)

if __name__ == '__main__': 
    main()
 
