#!/usr/bin/env python

#Imports
import argparse
from Ribosome_profiling_functions import read_counts
from Ribosome_profiling_functions import read_region_lengths

#Functions
def sum_CDS_counts(in_dict, frame, region_lengths_dict, remove_end_codons, outfyle):
    '''takes a counts dict and sums all the counts within each region, using a user defined offset'''
    summed_CDS_counts = {}
    
    for k, v in in_dict.items():
        UTR5_len = region_lengths_dict[k][0]
        cds_len = region_lengths_dict[k][1]
        cds_end = UTR5_len + cds_len
                
        cds_counts = v[UTR5_len:cds_end]
        if remove_end_codons == True:
            cds_counts = cds_counts[60:-30] #removes 20 codons from the 5' end and 10 codons from the 3' end
        
        if frame == None:
            summed_CDS_counts[k] = sum(map(int,cds_counts))
        else:
            framed_cds_counts = cds_counts[frame::3]
            summed_CDS_counts[k] = sum(map(int,framed_cds_counts))
    
    with open(outfyle, 'w') as g:
        g.write('transcript,CDS_counts\n')
        for k,v in summed_CDS_counts.items():
            g.write(k + "," + str(v) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Takes a RPF counts file and sums all counts CDS')
    parser.add_argument('infyle', type=str, help='counts file to pull values from')
    parser.add_argument('region_lengths_fyle', type=str, help='region lengths file. Needs to be csv file in the following format; Transcript ID,5\'UTR length,CDS length,3\'UTR length')
    parser.add_argument('-frame', type=int, default=None, help='Specificy frame 0/1/2. Use reads only within this frame. Default is to use all reads within CDS')
    parser.add_argument('-remove_end_codons', action="store_true", default=False, help = 'removes the first 20 and last 10 codons from the CDS (recommended')
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
        if args.frame == None:
            fylename = args.infyle.replace('.counts', '_CDS_counts_all_frames.csv')
        else:
            fylename = args.infyle.replace('.counts', '_CDS_counts_frame_') + str(args.frame) + '.csv'
    else:
        if args.frame == None:
            fylename = args.out_dir + '/' + args.infyle.replace('.counts', '_CDS_counts_all_frames.csv')
        else:
            fylename = args.out_dir + '/' + args.infyle.replace('.counts', '_CDS_counts_frame_') + str(args.frame) + '.csv'
    
    #write summed counts for each region to file
    sum_CDS_counts(input_counts, args.frame, region_lengths, args.remove_end_codons, fylename)

if __name__ == '__main__':
    main()
