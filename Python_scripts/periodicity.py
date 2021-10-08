#!/usr/bin/env python

#Imports
import argparse
from Ribosome_profiling_functions import read_counts
from Ribosome_profiling_functions import read_region_lengths

#Functions
def calculate_periodicity(in_dict, offset, region_lengths_dict, outfyle):
    '''takes a counts dict and sums all the counts within each CDS frame for each transcript, using a user defined offset'''
    f0_dict, f1_dict, f2_dict = {},{},{}
    
    for k, v in in_dict.items():
        UTR5_len = region_lengths_dict[k][0]
        cds_len = region_lengths_dict[k][1]
        cds_end = UTR5_len + cds_len
                
        if UTR5_len >= offset:
            cds_counts = v[(UTR5_len - offset):(cds_end - offset)]
            f0_counts = cds_counts[::3]
            f1_counts = cds_counts[1::3]
            f2_counts = cds_counts[2::3]
                        
            f0_dict[k] = sum(map(int,f0_counts))
            f1_dict[k] = sum(map(int,f1_counts))
            f2_dict[k] = sum(map(int,f2_counts))
    
    with open(outfyle, 'w') as g:
        g.write('transcript,frame,counts\n')
        for k,v in f0_dict.items():
            f0_out = ','.join((k, 'f0', str(v)))
            f1_out = ','.join((k, 'f1', str(f1_dict[k])))
            f2_out = ','.join((k, 'f2', str(f2_dict[k])))
            g.write('\n'.join((f0_out, f1_out, f2_out)) + '\n')
    
def main():
    parser = argparse.ArgumentParser(description='Takes a RPF counts file and sums all counts at each position n nt either side of\nstart and stop codons and the first and last n nt of the 5\' and 3\'UTRs\n(only for transcripts with sufficient lengths to cover entire regions being spliced)')
    parser.add_argument('infyle', type=str, help='counts file to pull values from')
    parser.add_argument('region_lengths_fyle', type=str, help='region lengths file. Needs to be csv file in the following format; Transcript ID,5\'UTR length,CDS length,3\'UTR length')
    parser.add_argument('-offset', type=int, default=0, help='number of nt to offset reads by')
    parser.add_argument('-in_dir', type=str, default=None, help='Change the input folder, default is current directory')
    parser.add_argument('-out_dir', type=str, default=None, help='Change the output folder, default is current directory')
    args = parser.parse_args()
    
    #Read in the counts file
    if args.in_dir == None:
        input_counts = read_counts(args.infyle)
    else:
        input_counts = read_counts(args.in_dir + '/' + args.infyle)
    
    #read in the region lengths
    region_lengths = read_region_lengths(args.region_lengths_fyle)
    
    #generate output filename
    if args.out_dir == None:
        fylename = args.infyle.replace('.counts', '_periodicity.csv')
    else:
        fylename = args.out_dir + '/' + args.infyle.replace('.counts', '_periodicity.csv')
        
    #calculate periodicity and write out
    calculate_periodicity(input_counts, args.offset, region_lengths, fylename)
    
if __name__ == '__main__':
    main()
