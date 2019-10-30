#!/usr/bin/env python

"""
Filter FASTQ reads according to quality and length
"""

from __future__ import print_function
import argparse, gzip, HTSeq, sys
import multiprocessing as mp
from multiprocessing import Pool, Value

def filter_n_split_read(intuple):

    read = intuple[0]
    min_base_qual = intuple[1] 
    mean_quality = intuple[2]
    max_length = intuple[3]
    
    global processed_read_count
    global low_mean_base_qual_read
    global long_read_count
    global low_base_qual_read_count

    with processed_read_count.get_lock():
        processed_read_count.value += 1

    htread = HTSeq.SequenceWithQualities(read[0].encode(), 
                                         read[1],
                                         read[2].encode())

    read_len = len(htread.qual)
    read_mean_qual = sum(htread.qual)/float(read_len)
    read_lowqual_bases = sum(htread.qual < min_base_qual)

    ## format FASTQ entry to print
    read = ('@' + read[1], read[0], '+', read[2])

    ## default no filter passing code
    filter_code = -1
   
    if read_lowqual_bases < 2:
        if read_mean_qual >= mean_quality:
            if read_len <= max_length:
                ## set reads-passing-filters file code
                filter_code = 0

            else:
                with long_read_count.get_lock():
                    long_read_count.value += 1
                ## set lenght-discarded read file code
                filter_code = 1

        else:
            with low_mean_base_qual_read.get_lock():
                low_mean_base_qual_read.value += 1

    else:
        with low_base_qual_read_count.get_lock():
            low_base_qual_read_count.value += 1

    ## return filter code and FASTQ formatted read
    return (filter_code, read)


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

    parser.add_argument('-j', '--jobs',
                        dest = 'jobs', 
                        default = 1,
                        required = False,
                        type = int,
                        help = 'The number of parallel jobs to execute')

    args = parser.parse_args()

    ## handle compressed files 
    if args.fastq_file.endswith('.gz'):
        f = gzip.open(args.fastq_file, 'rt')

    elif args.fastq_file == '-':
        f = sys.stdin
    
    else:
        f = open(args.fastq_file, 'r')

    fastq_file = HTSeq.FastqReader(f, 
                                   args.quality_encoding,
                                   raw_iterator = True)
    
    ## open file to save reads passing the filters
    if args.output == '-':
        outfile = sys.stdout
    
    elif args.output.endswith('.gz'):
        outfile = gzip.open(args.output, 'wt', compresslevel = 6)

    else:
        outfile = open(args.output, 'w')
    
    ## if enabled, open the file to save length-discarded reads
    long_read_file = None
    if not args.long_read_file == '':        
        ## get file handler
        if args.long_read_file.lower().endswith(".gz"):
            long_read_file = gzip.open(args.long_read_file, 'wt', 
                                       compresslevel = 6)
        
        else:
            long_read_file = open(args.long_read_file, 'w')

    ## initialize (global) counter variables
    low_base_qual_read_count = Value('i', 0)
    low_mean_base_qual_read  = Value('i', 0)
    long_read_count = Value('i', 0)
    processed_read_count = Value('i', 0)

    pool = mp.Pool(processes = args.jobs)

    ## loop input and pass it to multiple processes
    ## while reading results as soon as they are computed
    for res in pool.imap(filter_n_split_read, 
                             ((read,
                               args.min_base_qual, 
                               args.mean_quality, 
                               args.max_length) for read in fastq_file), 
                             chunksize = 10000):
        if res[0] == 0:
            outfile.write('\n'.join(res[1]) + '\n')

        elif res[0] == 1 and long_read_file:
            long_read_file.write('\n'.join(res[1]) + '\n')
    
    pool.close()
    pool.join()

    if long_read_file:
        long_read_file.close()
    
    outfile.close()

    with open(args.report_file, "w") as report_file:
        print('''{} processed reads'''.format(processed_read_count.value), 
              file = report_file)

        print('''{} reads discarded because > '''\
              '''{} bases had quality < '''\
              '''{}'''.format(low_base_qual_read_count.value,
                              2, 
                              args.min_base_qual), 
              file = report_file)

        print('''{} reads discarded because mean quality '''\
              '''was < {}'''.format(low_mean_base_qual_read.value, 
                                    args.mean_quality), 
              file = report_file)

        print('''{} reads were > {} nt long'''.format(long_read_count.value, 
                                                 args.max_length), 
              file = report_file)


