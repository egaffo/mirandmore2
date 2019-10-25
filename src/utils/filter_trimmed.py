#!/usr/bin/env python

"""
Filter FASTQ reads according to quality and length
"""

from __future__ import print_function
import argparse, gzip, HTSeq, sys
import multiprocessing as mp
from multiprocessing import Pool, Value

### https://gist.github.com/ngcrawford/2237170
#from itertools import izip_longest
#
#def process_chunk(d):
#    """Replace this with your own function
#    that processes data one line at a
#    time"""
#
#    d = d.strip() + ' processed'
#    return d 
#
#def grouper(n, iterable, padvalue=None):
#    """grouper(3, 'abcdefg', 'x') -->
#    ('a','b','c'), ('d','e','f'), ('g','x','x')"""
#
#    return izip_longest(*[iter(iterable)]*n, fillvalue=padvalue)
#
#def main():
#p = multiprocessing.Pool(4)
#
#    # Use 'grouper' to split test data into
#    # groups you can process without using a
#    # ton of RAM. You'll probably want to 
#    # increase the chunk size considerably
#    # to something like 1000 lines per core.
#
#    # The idea is that you replace 'test_data'
#    # with a file-handle
#    # e.g., testdata = open(file.txt,'rU')
#
#    # And, you'd write to a file instead of
#    # printing to the stout
#
#    for chunk in grouper(10, test_data):
#       results = p.map(process_chunk, chunk)
#      for r in results:
#                  print r
#                  # replace
#                  # with
#                  # outfile.write()

def mirna_reads_writer(q):
    ## assume a outfile file handler has been instantiated
    while 1:
        read = q.get()
        if str(read) == 'kill':
            break
        outfile.write('\n'.join(read) + '\n')
        outfile.flush()

def long_reads_writer(q):
    ## assume a long_read_file file handler has been instantiated
    while 1:
        read = q.get()
        if str(read) == 'kill':
            break
        if long_read_file:
            long_read_file.write('\n'.join(read) + '\n')
            long_read_file.flush()

def filter_n_split_read(read, min_base_qual, mean_quality, max_length, 
                        outfile_queue, long_read_file_queue):
    
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
   
    if read_lowqual_bases < 2:
        if read_mean_qual >= mean_quality:
            if read_len <= max_length:
                ## save to reads-passing-filters file
                outfile_queue.put(read)

            else:
                with long_read_count.get_lock():
                    long_read_count.value += 1
                ## save to lenght-discarded read file
                long_read_file_queue.put(read)

        else:
            with low_mean_base_qual_read.get_lock():
                low_mean_base_qual_read.value += 1

    else:
        with low_base_qual_read_count.get_lock():
            low_base_qual_read_count.value += 1

    return None


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

    if args.output == '-':
        outfile = sys.stdout
    
    elif args.output.endswith('.gz'):
        outfile = gzip.open(args.output, 'wt')
    else:
        outfile = open(args.output, 'w')

    long_read_file = None
    if not args.long_read_file == '':        
        ## get file handler to save discarded reads from length filter
        if args.long_read_file.lower().endswith(".gz"):
            long_read_file = gzip.open(args.long_read_file, 'wt')
        
        else:
            long_read_file = open(args.long_read_file, 'w')

    ## initialize global variables
    low_base_qual_read_count = Value('i', 0)
    low_mean_base_qual_read  = Value('i', 0)
    long_read_count = Value('i', 0)
    processed_read_count = Value('i', 0)

    manager = mp.Manager()
    
    pool = mp.Pool(processes = args.jobs)

    outfile_queue = manager.Queue()
    watcher_1 = pool.apply_async(mirna_reads_writer, (outfile_queue,))

    long_read_file_queue = manager.Queue()
    watcher_2 = pool.apply_async(long_reads_writer, (long_read_file_queue,))

    #res = pool.starmap(filter_n_split_read, 
    res = pool.starmap_async(filter_n_split_read, 
                             [(read,
                               args.min_base_qual, 
                               args.mean_quality, 
                               args.max_length,
                               outfile_queue,
                               long_read_file_queue) for read in fastq_file], 
                             chunksize = 10000)
    
    outfile_queue.put('kill')
    long_read_file_queue.put('kill')
    pool.close()
    pool.join()

    if not args.long_read_file == '':
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


