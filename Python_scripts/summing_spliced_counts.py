#!/usr/bin/env python

#Imports
import argparse
from Ribosome_profiling_functions import read_counts
from Ribosome_profiling_functions import read_region_lengths

#Functions
def splice_counts(in_dict, region_lengths_dict, n, polyA_lib):
    '''takes a counts dict and for all transcripts with a 5'UTR, CDS and 3'UTR greater than 2 * n (user defined integer),
    it splices the first n nt of the 5'UTR, n nt either side of start and stop codons and last n nt of the 3'UTR'''
    fp_end_counts, start_site_counts, stop_site_counts, tp_end_counts = {},{},{},{}
    for k, v in in_dict.items():
        UTR5_len = region_lengths_dict[k][0]
        cds_len = region_lengths_dict[k][1]
        UTR3_len = region_lengths_dict[k][2]
        cds_end = UTR5_len + cds_len
        
        if UTR5_len >= (2 * n) and cds_len >= (2 * n) and UTR3_len >= (2 * n):
            fp_end_counts[k] = v[:n]
            start_site_counts[k] = v[(UTR5_len - n):(UTR5_len + n)]
            stop_site_counts[k] = v[(cds_end - n):(cds_end + n)]
            if polyA_lib == True:
                tp_end_counts[k] = v[-(n + 100):-100]
            else:
                tp_end_counts[k] = v[-n:]
    return (fp_end_counts, start_site_counts, stop_site_counts, tp_end_counts)

def sum_counts(adict, n):
    '''sums all counts at each position within a dict and returns as a list'''
    summed_counts = [0]*n
    for counts in adict.values():
        summed_counts = [x + y for x, y in zip(summed_counts, map(int, counts))]
    return summed_counts

def write_junction_counts(alist, outfyle, n):
    '''Takes a list of summed junction counts and writes out as a csv with positions'''
    with open(outfyle, 'w') as g:
        g.write('position,counts\n')
        zipped = zip(range(-n, n), alist)
        for i in zipped:
            g.write(','.join(map(str,i)) + '\n')

def write_end_counts(alist, outfyle, n):
    '''Takes a list of summed end counts and writes out as a csv with positions'''
    with open(outfyle, 'w') as g:
        g.write('position,counts\n')
        zipped = zip(range(1, n+1), alist)
        for i in zipped:
            g.write(','.join(map(str,i)) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Takes a RPF counts file and sums all counts at each position n nt either side of\nstart and stop codons and the first and last n nt of the 5\' and 3\'UTRs\n(only for transcripts with sufficient lengths to cover entire regions being spliced)')
    parser.add_argument('infyle', type=str, help='counts file to pull values from')
    parser.add_argument('n', type=int, help='n number of nt to splice')
    parser.add_argument('region_lengths_fyle', type=str, help='region lengths file. Needs to be csv file in the following format; Transcript ID,5\'UTR length,CDS length,3\'UTR length')
    parser.add_argument('-polyA_lib', action="store_true", default=False, help = 'Removes the 3\' most 100nt due to bias in poly(A) library prep')
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
    
    #splice counts
    fp_end_counts, start_site_counts, stop_site_counts, tp_end_counts = splice_counts(input_counts, region_lengths, args.n, args.polyA_lib)
    
    #sum spliced counts
    summed_start_site_counts = sum_counts(start_site_counts, (args.n * 2))
    summed_stop_site_counts = sum_counts(stop_site_counts, (args.n * 2))
    
    summed_fp_end_counts = sum_counts(fp_end_counts, args.n)
    summed_tp_end_counts = sum_counts(tp_end_counts, args.n)
    
    #write summed spliced counts to file
    if args.out_dir == None:
        start_site_fylename = args.infyle.replace('.counts', '_start_site.csv')
        stop_site_fylename = args.infyle.replace('.counts', '_stop_site.csv')
        fp_end_fylename = args.infyle.replace('.counts', '_UTR5_start.csv')
        tp_end_fylename = args.infyle.replace('.counts', '_UTR3_end.csv')
    else:
        start_site_fylename = args.out_dir + '/' + args.infyle.replace('.counts', '_start_site.csv')
        stop_site_fylename = args.out_dir + '/' + args.infyle.replace('.counts', '_stop_site.csv')
        fp_end_fylename = args.out_dir + '/' + args.infyle.replace('.counts', '_UTR5_start.csv')
        tp_end_fylename = args.out_dir + '/' + args.infyle.replace('.counts', '_UTR3_end.csv')
    
    write_junction_counts(summed_start_site_counts, start_site_fylename, args.n)
    write_junction_counts(summed_stop_site_counts, stop_site_fylename, args.n)
    
    write_end_counts(summed_fp_end_counts, fp_end_fylename, args.n)
    write_end_counts(summed_tp_end_counts, tp_end_fylename, args.n)

if __name__ == '__main__':
    main()
