#!/usr/bin/env python

#Imports
import argparse

#Functions
def extract_cutadapt_counts(fylename):
    '''takes a cutadapt log file and extracts the read counts'''
    logfyle=fylename + "_cutadapt_log.txt"
    with open(logfyle,'r') as f:
        for line in f:
            if line.startswith('Total reads processed'):
                in_counts = line.split(':')[1].strip().replace(',','')
            if line.startswith('Reads written (passing filters):'):
                out_counts = line.split(':')[1].split('(')[0].strip().replace(',','')
    return (in_counts, out_counts)
    
def extract_bbmap_counts(fylename,RNA_molecule):
    '''takes a bbmap log and extracts mapped and umapped read counts'''
    logfyle=fylename + "_" + RNA_molecule + "_log.txt"
    with open(logfyle,'r') as f:
        for line in f:
            if line.startswith('Reads Used:'):
                in_counts = line.split('\t')[1]
            if line.startswith('mapped:'):
                out_counts = line.split('\t')[2].strip()
    return (in_counts, out_counts)
    
def extract_UMI_clipped_counts(fylename):
    '''takes a UMItools extracted UMI log files and extracts the read counts'''
    logfyle=fylename + "_extracted_UMIs.log"
    with open(logfyle,'r') as f:
        for line in f:
            if 'INFO Input Reads:' in line:
                in_counts = line.split(':')[-1].strip()
            if 'INFO regex matches read1:' in line:
                out_counts = line.split(':')[-1].strip()
    return (in_counts, out_counts)
    
def extract_deduplication_counts(fylename):
    '''takes a UMItools deduplication log files and extracts the read counts'''
    logfyle=fylename + "_deduplication_log.txt"
    with open(logfyle,'r') as f:
        for line in f:
            if 'INFO Reads: Input Reads:' in line:
                in_counts = line.split(':')[-1].strip()
            if 'INFO Number of reads out:' in line:
                out_counts = line.split(':')[-1].strip()
    return (in_counts, out_counts)
    
def main():
    parser = argparse.ArgumentParser(description='Takes a filename and extracts the read counts for each step of the library')
    parser.add_argument('infyle', type=str, help='filename to pull values from the log files')
    parser.add_argument('-log_dir', type=str, default=None, help='directory containing log files')
    args = parser.parse_args()
    
    cutadapt_in, cutadapt_out = extract_cutadapt_counts(args.log_dir + '/' + args.infyle)
    
    UMI_clipped_in, UMI_clipped_out = extract_UMI_clipped_counts(args.log_dir + '/' + args.infyle)
    
    rRNA_in, rRNA_out = extract_bbmap_counts(args.log_dir + '/' + args.infyle, "rRNA")
    tRNA_in, tRNA_out = extract_bbmap_counts(args.log_dir + '/' + args.infyle, "tRNA")
    pc_in, pc_out = extract_bbmap_counts(args.log_dir + '/' + args.infyle, "pc")
    
    deduplication_in, deduplication_out = extract_deduplication_counts(args.log_dir + '/' + args.infyle)
    
    #check everything adds up
    if int(cutadapt_out) != int(UMI_clipped_in):
        print("warning: cutadapt out does not equal UMI clipped in for sample: " + args.infyle)
    if int(UMI_clipped_out) != int(rRNA_in):
        print("warning: UMI clipped out does not equal rRNA in for sample: " + args.infyle)
    if (int(rRNA_in) - int(rRNA_out)) != int(tRNA_in):
        print("warning: non-rRNA reads does not equal tRNA in for sample: " + args.infyle)
    if (int(tRNA_in) - int(tRNA_out)) != pc_in:
        print("warning: non-rRNA_tRNA reads does not equal pc in for sample: " + args.infyle)
    if int(pc_out) != int(deduplication_in):
        print("warning: pc reads does not input to depulication for sample: " + args.infyle)
    
    #write output
    outfyle = args.log_dir + '/' + args.infyle + "_read_counts.csv"
    with open(outfyle,'w') as g:
        g.write("cutadapt_in,cutadapt_out,UMI_clipped_in,UMI_clipped_out,rRNA_in,rRNA_out,tRNA_in,tRNA_out,pc_in,pc_out,deduplication_in,deduplication_out\n")
        g.write((',').join((cutadapt_in,cutadapt_out,UMI_clipped_in,UMI_clipped_out,rRNA_in,rRNA_out,tRNA_in,tRNA_out,pc_in,pc_out,deduplication_in,deduplication_out)))

if __name__ == '__main__':
    main()
