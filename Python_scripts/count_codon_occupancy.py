#!/usr/bin/env python

#Imports
import argparse
from Ribosome_profiling_functions import read_counts
from Ribosome_profiling_functions import read_region_lengths
from Ribosome_profiling_functions import read_in_fasta

#Functions
def splice_cds(in_dict, region_lengths_dict):
    '''takes a dictionary and splices the CDS from the values'''
    out_dict = {}
    for k, v in in_dict.items():
        UTR5_len = region_lengths_dict[k][0]
        cds_len = region_lengths_dict[k][1]
        cds_end = UTR5_len + cds_len
        
        out_dict[k] = v[UTR5_len:cds_end]
    return out_dict

def count_codons(counts_dict, seqs_dict, outfyle):
    ''''''
    minus_4_dict = {"AAA":0,"AAC":0,"AAG":0,"AAT":0,"ACA":0,"ACC":0,"ACG":0,"ACT":0,"ATT":0,"ATC":0,"ATA":0,"ATG":0,"AGT":0,"AGC":0,"AGA":0,"AGG":0,"CAA":0,"CAC":0,"CAG":0,"CAT":0,"CCA":0,"CCC":0,"CCG":0,"CCT":0,"CTT":0,"CTC":0,"CTA":0,"CTG":0,"CGT":0,"CGC":0,"CGA":0,"CGG":0,"GAA":0,"GAC":0,"GAG":0,"GAT":0,"GCA":0,"GCC":0,"GCG":0,"GCT":0,"GTT":0,"GTC":0,"GTA":0,"GTG":0,"GGT":0,"GGC":0,"GGA":0,"GGG":0,"TAA":0,"TAC":0,"TAG":0,"TAT":0,"TCA":0,"TCC":0,"TCG":0,"TCT":0,"TTT":0,"TTC":0,"TTA":0,"TTG":0,"TGT":0,"TGC":0,"TGA":0,"TGG":0}
    minus_3_dict = {"AAA":0,"AAC":0,"AAG":0,"AAT":0,"ACA":0,"ACC":0,"ACG":0,"ACT":0,"ATT":0,"ATC":0,"ATA":0,"ATG":0,"AGT":0,"AGC":0,"AGA":0,"AGG":0,"CAA":0,"CAC":0,"CAG":0,"CAT":0,"CCA":0,"CCC":0,"CCG":0,"CCT":0,"CTT":0,"CTC":0,"CTA":0,"CTG":0,"CGT":0,"CGC":0,"CGA":0,"CGG":0,"GAA":0,"GAC":0,"GAG":0,"GAT":0,"GCA":0,"GCC":0,"GCG":0,"GCT":0,"GTT":0,"GTC":0,"GTA":0,"GTG":0,"GGT":0,"GGC":0,"GGA":0,"GGG":0,"TAA":0,"TAC":0,"TAG":0,"TAT":0,"TCA":0,"TCC":0,"TCG":0,"TCT":0,"TTT":0,"TTC":0,"TTA":0,"TTG":0,"TGT":0,"TGC":0,"TGA":0,"TGG":0}
    E_site_dict = {"AAA":0,"AAC":0,"AAG":0,"AAT":0,"ACA":0,"ACC":0,"ACG":0,"ACT":0,"ATT":0,"ATC":0,"ATA":0,"ATG":0,"AGT":0,"AGC":0,"AGA":0,"AGG":0,"CAA":0,"CAC":0,"CAG":0,"CAT":0,"CCA":0,"CCC":0,"CCG":0,"CCT":0,"CTT":0,"CTC":0,"CTA":0,"CTG":0,"CGT":0,"CGC":0,"CGA":0,"CGG":0,"GAA":0,"GAC":0,"GAG":0,"GAT":0,"GCA":0,"GCC":0,"GCG":0,"GCT":0,"GTT":0,"GTC":0,"GTA":0,"GTG":0,"GGT":0,"GGC":0,"GGA":0,"GGG":0,"TAA":0,"TAC":0,"TAG":0,"TAT":0,"TCA":0,"TCC":0,"TCG":0,"TCT":0,"TTT":0,"TTC":0,"TTA":0,"TTG":0,"TGT":0,"TGC":0,"TGA":0,"TGG":0}
    P_site_dict = {"AAA":0,"AAC":0,"AAG":0,"AAT":0,"ACA":0,"ACC":0,"ACG":0,"ACT":0,"ATT":0,"ATC":0,"ATA":0,"ATG":0,"AGT":0,"AGC":0,"AGA":0,"AGG":0,"CAA":0,"CAC":0,"CAG":0,"CAT":0,"CCA":0,"CCC":0,"CCG":0,"CCT":0,"CTT":0,"CTC":0,"CTA":0,"CTG":0,"CGT":0,"CGC":0,"CGA":0,"CGG":0,"GAA":0,"GAC":0,"GAG":0,"GAT":0,"GCA":0,"GCC":0,"GCG":0,"GCT":0,"GTT":0,"GTC":0,"GTA":0,"GTG":0,"GGT":0,"GGC":0,"GGA":0,"GGG":0,"TAA":0,"TAC":0,"TAG":0,"TAT":0,"TCA":0,"TCC":0,"TCG":0,"TCT":0,"TTT":0,"TTC":0,"TTA":0,"TTG":0,"TGT":0,"TGC":0,"TGA":0,"TGG":0}
    A_site_dict = {"AAA":0,"AAC":0,"AAG":0,"AAT":0,"ACA":0,"ACC":0,"ACG":0,"ACT":0,"ATT":0,"ATC":0,"ATA":0,"ATG":0,"AGT":0,"AGC":0,"AGA":0,"AGG":0,"CAA":0,"CAC":0,"CAG":0,"CAT":0,"CCA":0,"CCC":0,"CCG":0,"CCT":0,"CTT":0,"CTC":0,"CTA":0,"CTG":0,"CGT":0,"CGC":0,"CGA":0,"CGG":0,"GAA":0,"GAC":0,"GAG":0,"GAT":0,"GCA":0,"GCC":0,"GCG":0,"GCT":0,"GTT":0,"GTC":0,"GTA":0,"GTG":0,"GGT":0,"GGC":0,"GGA":0,"GGG":0,"TAA":0,"TAC":0,"TAG":0,"TAT":0,"TCA":0,"TCC":0,"TCG":0,"TCT":0,"TTT":0,"TTC":0,"TTA":0,"TTG":0,"TGT":0,"TGC":0,"TGA":0,"TGG":0}
    plus_1_dict = {"AAA":0,"AAC":0,"AAG":0,"AAT":0,"ACA":0,"ACC":0,"ACG":0,"ACT":0,"ATT":0,"ATC":0,"ATA":0,"ATG":0,"AGT":0,"AGC":0,"AGA":0,"AGG":0,"CAA":0,"CAC":0,"CAG":0,"CAT":0,"CCA":0,"CCC":0,"CCG":0,"CCT":0,"CTT":0,"CTC":0,"CTA":0,"CTG":0,"CGT":0,"CGC":0,"CGA":0,"CGG":0,"GAA":0,"GAC":0,"GAG":0,"GAT":0,"GCA":0,"GCC":0,"GCG":0,"GCT":0,"GTT":0,"GTC":0,"GTA":0,"GTG":0,"GGT":0,"GGC":0,"GGA":0,"GGG":0,"TAA":0,"TAC":0,"TAG":0,"TAT":0,"TCA":0,"TCC":0,"TCG":0,"TCT":0,"TTT":0,"TTC":0,"TTA":0,"TTG":0,"TGT":0,"TGC":0,"TGA":0,"TGG":0}
    plus_2_dict = {"AAA":0,"AAC":0,"AAG":0,"AAT":0,"ACA":0,"ACC":0,"ACG":0,"ACT":0,"ATT":0,"ATC":0,"ATA":0,"ATG":0,"AGT":0,"AGC":0,"AGA":0,"AGG":0,"CAA":0,"CAC":0,"CAG":0,"CAT":0,"CCA":0,"CCC":0,"CCG":0,"CCT":0,"CTT":0,"CTC":0,"CTA":0,"CTG":0,"CGT":0,"CGC":0,"CGA":0,"CGG":0,"GAA":0,"GAC":0,"GAG":0,"GAT":0,"GCA":0,"GCC":0,"GCG":0,"GCT":0,"GTT":0,"GTC":0,"GTA":0,"GTG":0,"GGT":0,"GGC":0,"GGA":0,"GGG":0,"TAA":0,"TAC":0,"TAG":0,"TAT":0,"TCA":0,"TCC":0,"TCG":0,"TCT":0,"TTT":0,"TTC":0,"TTA":0,"TTG":0,"TGT":0,"TGC":0,"TGA":0,"TGG":0}
    total_dict = {"AAA":0,"AAC":0,"AAG":0,"AAT":0,"ACA":0,"ACC":0,"ACG":0,"ACT":0,"ATT":0,"ATC":0,"ATA":0,"ATG":0,"AGT":0,"AGC":0,"AGA":0,"AGG":0,"CAA":0,"CAC":0,"CAG":0,"CAT":0,"CCA":0,"CCC":0,"CCG":0,"CCT":0,"CTT":0,"CTC":0,"CTA":0,"CTG":0,"CGT":0,"CGC":0,"CGA":0,"CGG":0,"GAA":0,"GAC":0,"GAG":0,"GAT":0,"GCA":0,"GCC":0,"GCG":0,"GCT":0,"GTT":0,"GTC":0,"GTA":0,"GTG":0,"GGT":0,"GGC":0,"GGA":0,"GGG":0,"TAA":0,"TAC":0,"TAG":0,"TAT":0,"TCA":0,"TCC":0,"TCG":0,"TCT":0,"TTT":0,"TTC":0,"TTA":0,"TTG":0,"TGT":0,"TGC":0,"TGA":0,"TGG":0}
    
    for transcript, counts in counts_dict.items():
        seq = seqs_dict[transcript]
        
        for i in range(60, (len(counts) - (12 + 30)), 3): #this will create a range to iterate through each codon in the CDS. Starting at 60 removes the first 20 codons. Ending at -(12 + 30) removes the final 10 codons
            RPF_seq = seq[(i - 9):(i + 12)] #this accounts for the previously applied offset, meaning that the each count refers to the codon in the P-site
            RPF_counts = int(counts[i])
            
            minus_4_codon = (RPF_seq[:3])
            minus_4_dict[minus_4_codon] += RPF_counts
            total_dict[minus_4_codon] += RPF_counts
            
            minus_3_codon = (RPF_seq[3:6])
            minus_3_dict[minus_3_codon] += RPF_counts
            total_dict[minus_3_codon] += RPF_counts
            
            E_site_codon = (RPF_seq[6:9])
            E_site_dict[E_site_codon] += RPF_counts
            total_dict[E_site_codon] += RPF_counts
            
            P_site_codon = (RPF_seq[9:12])
            P_site_dict[P_site_codon] += RPF_counts
            total_dict[P_site_codon] += RPF_counts
            
            A_site_codon = (RPF_seq[12:15])
            A_site_dict[A_site_codon] += RPF_counts
            total_dict[A_site_codon] += RPF_counts
            
            plus_1_codon = (RPF_seq[15:18])
            plus_1_dict[plus_1_codon] += RPF_counts
            total_dict[plus_1_codon] += RPF_counts
            
            plus_2_codon = (RPF_seq[18:21])
            plus_2_dict[plus_2_codon] += RPF_counts
            total_dict[plus_2_codon] += RPF_counts
            
    with open(outfyle, 'w') as g:
        g.write('codon,-4,-3,E,P,A,1,2,total\n')
        for codon in total_dict.keys():
            g.write(','.join([codon, str(minus_4_dict[codon]), str(minus_3_dict[codon]), str(E_site_dict[codon]), str(P_site_dict[codon]), str(A_site_dict[codon]), str(plus_1_dict[codon]), str(plus_2_dict[codon]), str(total_dict[codon])]) + '\n')

def filter_dictonary(in_dict,filter_dict):
    '''Removes entries from a dictionary'''
    for k in in_dict.keys():
        if k not in filter_dict:
            del in_dict[k]

def read_in_target_transcripts(txt_file):
    '''Read in an transcript ID file'''
    transcript_dict = {}
    with open(txt_file,'r') as f:
        for line in f:
            transcript_dict[line.strip()] = None
    return transcript_dict
    
    

def main():
    parser = argparse.ArgumentParser(description='This script will, for each codon, count the number of reads that correspond to that codon being at each position within the ribosome.\nThe <.counts> input must have already been offset so that the position of each count refers to the start of the P-site for that read')
    parser.add_argument('infyle', type=str, help='counts file to pull values from')
    parser.add_argument('fasta', type=str, help='FASTA file. Must be the same file that was used to align reads')
    parser.add_argument('region_lengths_fyle', type=str, help='region lengths file. Needs to be csv file in the following format; Transcript ID,5\'UTR length,CDS length,3\'UTR length')
    parser.add_argument('-transcripts', type=str, default=None, help='List of transcript IDs. Carry out analysis on these transcripts only, default is all transcripts')
    parser.add_argument('-in_dir', type=str, default=None, help='Change the input folder, default is current directory')
    parser.add_argument('-out_dir', type=str, default=None, help='Change the output folder, default is current directory')
    args = parser.parse_args()
    
    #Read in the counts file
    if args.in_dir == None:
        counts = read_counts(args.infyle)
    else:
        counts = read_counts(args.in_dir + '/' + args.infyle)
    
    #filter transcripts by list if option selected
    if args.transcripts != None:
        restrict_dict = read_in_target_transcripts(args.transcripts)
        filter_dictonary(counts, restrict_dict)
        
    #read in fasta
    seqs = read_in_fasta(args.fasta)
    
    #read in the region lengths
    region_lengths = read_region_lengths(args.region_lengths_fyle)
        
    #splice cds counts and seq
    cds_seq = splice_cds(seqs, region_lengths)
    cds_counts = splice_cds(counts, region_lengths)
    
    #count codons and write to file
    if args.out_dir == None:
        if args.transcripts == None:
            fylename = args.infyle.replace('.counts', '_counts.csv')
        else:
            fylename = args.infyle.replace('.counts', '_counts_') + args.transcripts.replace('txt', 'csv')
    else:
        if args.transcripts == None:
            fylename = args.out_dir + '/' + args.infyle.replace('.counts', '_counts.csv')
        else:
            fylename = args.out_dir + '/' + args.infyle.replace('.counts', '_counts_') + args.transcripts.split('/')[-1].replace('txt', 'csv')
        
    count_codons(cds_counts, cds_seq, fylename)
    
if __name__ == '__main__':
    main()


