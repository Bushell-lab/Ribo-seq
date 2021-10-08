#!/usr/bin/env python

#Imports
import argparse
from Ribosome_profiling_functions import read_counts
from Ribosome_profiling_functions import read_region_lengths

#Functions
def sum_region_counts(in_dict, offset, region_lengths_dict, outfyle):
    '''takes a counts dict and sums all the counts within each region, using a user defined offset'''
    summed_UTR5_counts, summed_CDS_counts, summed_UTR3_counts = {},{},{}
    
    for k, v in in_dict.items():
        UTR5_len = region_lengths_dict[k][0]
        cds_len = region_lengths_dict[k][1]
        cds_end = UTR5_len + cds_len
        UTR3_len = region_lengths_dict[k][2]
        
        if UTR5_len >= offset:
            UTR5_counts = v[:(UTR5_len - offset)]
            summed_UTR5_counts[k] = sum(map(int,UTR5_counts))
        else:
            summed_UTR5_counts[k] = 'NA'
        
        cds_counts = v[(UTR5_len - offset):(cds_end - offset)]
        summed_CDS_counts[k] = sum(map(int,cds_counts))
        
        UTR3_counts = v[(cds_end - offset):]
        summed_UTR3_counts[k] = sum(map(int,UTR3_counts))
    
    with open(outfyle, 'w') as g:
        g.write('transcript,region,counts\n')
        for k,v in summed_UTR5_counts.items():
            UTR5_out = ','.join((k, 'UTR5', str(v)))
            cds_out = ','.join((k, 'CDS', str(summed_CDS_counts[k])))
            UTR3_out = ','.join((k, 'UTR3', str(summed_UTR3_counts[k])))
            g.write('\n'.join((UTR5_out, cds_out, UTR3_out)) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Takes a RPF counts file and sums all counts within each 5\'UTR, CDS and 3\'UTR')
    parser.add_argument('infyle', type=str, help='counts file to pull values from')
    parser.add_argument('offset', type=int, help='number of nt to offset reads by')
    parser.add_argument('region_lengths_fyle', type=str, help='region lengths file. Needs to be csv file in the following format; Transcript ID,5\'UTR length,CDS length,3\'UTR length')
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
    
    #generate file name
    if args.out_dir == None:
        fylename = args.infyle.replace('.counts', '_region_counts.csv')
    else:
        fylename = args.out_dir + '/' + args.infyle.replace('.counts', '_region_counts.csv')
        
    #calculate and write region counts
    sum_region_counts(input_counts, args.offset, region_lengths, fylename)

if __name__ == '__main__':
    main()
