#!/usr/bin/env python

"""
Filter FASTQ reads according to quality and length
"""

from __future__ import print_function
import argparse, gzip, HTSeq, sys
import multiprocessing as mp

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input',
                        dest = 'fastq_file',
                        required = True)
    
    parser.add_argument('-o', '--output', 
                        dest = 'output', 
                        required = False, 
                        default = "-")
    
    parser.add_argument('-l', '--max-length',
                        dest = 'max_length',
                        default = 26, 
                        type = int)
    
    parser.add_argument('-L', '--long-read-file',
                        dest = 'long_read_file',
                        default = "", 
                        type = str)
    
    parser.add_argument('-q', '--mean-quality',
                        dest = 'mean_quality',
                        default = 30, 
                        type = int)
    
    parser.add_argument('-e', '--quality-encoding',
                        dest = 'quality_encoding',
                        default = "phred")
    
    parser.add_argument('-r', '--report-file',
                        dest = 'report_file',
                        required = False,
                        default = 'filter_report.txt')
    
    parser.add_argument('-p', '--min_base_qual',
                        dest = 'min_base_qual',
                        default = 20,
                        required = False,
                        type = int)

    args = parser.parse_args()

    ## handle compressed files 
    if args.fastq_file.endswith('.gz'):
        f = gzip.open(args.fastq_file, 'rt')

    elif args.fastq_file == '-':
        f = sys.stdin
    
    else:
        f = open(args.fastq_file, 'r')

    fastq_file = HTSeq.FastqReader(f, args.quality_encoding)

    if args.output == '-':
        outfile = sys.stdout
    
    elif args.output.endswith('.gz'):
        outfile = gzip.open(args.output, 'wt')
    else:
        outfile = open(args.output, 'w')

    if not args.long_read_file == '':
        
        ## get file handler to save discarded reads from length filter
        if args.long_read_file.lower().endswith(".gz"):
            long_read_file = gzip.open(args.long_read_file, 'wt')
        
        else:
            long_read_file = open(args.long_read_file, 'w')

    low_base_qual_read_count = 0
    low_mean_base_qual_read  = 0
    long_read_count = 0
    processed_read_count = 0

    for read in fastq_file:
        processed_read_count += 1
        read_len = len(read.qual)
        read_mean_qual = sum(read.qual)/float(read_len)
        read_lowqual_bases = sum(read.qual < args.min_base_qual)
        
        if read_lowqual_bases < 2:
            if read_mean_qual >= args.mean_quality:
                if read_len <= args.max_length:
                    read.write_to_fastq_file(outfile)

                else:
                    long_read_count += 1
                    if not args.long_read_file == '':
                        read.write_to_fastq_file(long_read_file)

            else:
                low_mean_base_qual_read += 1

        else:
            low_base_qual_read_count += 1

    if not args.long_read_file == '':
        long_read_file.close()
    
    outfile.close()

    with open(args.report_file, "w") as report_file:
        print('''{} processed reads'''.format(processed_read_count), 
              file = report_file)

        print('''{} reads discarded because > '''\
              '''{} bases had quality < '''\
              '''{}'''.format(low_base_qual_read_count,
                              2, 
                              args.min_base_qual), 
              file = report_file)

        print('''{} reads discarded because mean quality '''\
              '''was < {}'''.format(low_mean_base_qual_read, 
                                    args.mean_quality), 
              file = report_file)

        print('''{} reads were > {} nt long'''.format(long_read_count, 
                                                 args.max_length), 
              file = report_file)


