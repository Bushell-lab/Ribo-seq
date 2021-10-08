#!/usr/bin/env python

#Imports
import sys
import pysam
import os
import argparse
import logging
import collections

#Functions
def lengths_offsets(value):
    """Split the given comma separated value to multiple integer values. """
    values = []
    for item in value.split(','):
        item = int(item)
        values.append(item)

    return values


def get_RPF_counts(pysamobj, transcript_name, transcript_length, read_lengths, read_offsets):
    read_counts = collections.defaultdict(int)
    for _ in range(1, transcript_length + 1):
        read_counts[_] = 0

    total_reads = 0
    for record in pysamobj.fetch(transcript_name):
        query_length = record.query_length
        position_ref = record.pos + 1
        for index, read_length in enumerate(read_lengths):
            position = position_ref
            if read_length == 0 or read_length == query_length:
                position += read_offsets[index]
                total_reads += 1
                read_counts[position] += 1
            else:
                continue

    return (read_counts, total_reads)


# create logger for the entire program
log = logging.getLogger('riboplot')
log.setLevel(logging.DEBUG)

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s - %(module)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
ch.setFormatter(formatter)
log.addHandler(ch)

def create_parser():
    """Argument parser. """
    parser = argparse.ArgumentParser(
        prog='countingScript.py', description='Output read counts for all transcripts')

    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-bam', '--BAMinput', help='Ribo-Seq alignment file in BAM format', required=True)
    required.add_argument('-fasta', '--transcriptome_fasta', help='FASTA format file of the transcriptome', required=True)

    # optional arguments
    parser.add_argument('-len', '--read_lengths', help='Read lengths to consider (default: %(default)s). '
                        'Multiple read lengths should be separated by commas. If multiple read lengths '
                        'are specified, corresponding read offsets should also be specified. If you do '
                        'not wish to apply an offset, please input 0 for the corresponding read length',
                        default='0', type=lengths_offsets)
    parser.add_argument('-offset', '--read_offsets', help='Read offsets (default: %(default)s). '
                        'Multiple read offsets should be separated by commas',
                        default='0', type=lengths_offsets)

    parser.add_argument('-out_dir', '--outpath', help='Files are saved in this directory', default='output')
    parser.add_argument('-out_file', '--outfile', help='name to save file', default='outputfile')

    return parser

def main(args):
    (BAMinput,fasta_file,read_lengths,read_offsets,outpath,outfile) = \
         (args.BAMinput,args.transcriptome_fasta,args.read_lengths,args.read_offsets,args.outpath,args.outfile)
    
    log.debug('Supplied arguments\n{}'.format(
        '\n'.join(['{:<20}: {}'.format(k, v) for k, v in vars(args).items()])))

    #also need BAM index ****
    #should add in QC steps here to check inputs are valid
    log.info('check 1')

    with pysam.AlignmentFile(BAMinput, 'rb') as b, pysam.FastaFile(fasta_file) as f:
    
        count=0
    
    #create output directories
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        os.chdir(outpath)
        log.info('directory created')

        log.info('Getting RPF counts for all transcripts in FASTA')
        with open(outfile,'aw') as outputfile:
            for transcript in f.references:
                transcript_sequence=f.fetch(transcript)
                ribo_counts, ribo_reads = get_RPF_counts(pysamobj=b,transcript_name=transcript,transcript_length=len(transcript_sequence),read_lengths=read_lengths,read_offsets=read_offsets)
                if ribo_reads>1:
                    outputfile.write(str(transcript)+'\n'+'\t'.join(map(str,ribo_counts.values()))+'\n')

                
        outputfile.close()
        log.info('RPF counting complete')

def run():
    """Run program"""
    parsed = create_parser()
    args = parsed.parse_args()
    log.info('args')
    print(args)
    main(args)

if __name__ == '__main__':
    run()
