#!/usr/bin/env python

"""
filtering is hardcoded but the pipeline in part or in toto
can be built programmatically using reduce.

For example:
chain = lambda b,a: b.__pipe__(a) if hasattr(b, '__pipe__') else b(a)

chain is a function version of >> stream operator

range(10) >> reduce(chain, [stream.filter(lambda x: x%2 ==0), stream.filter(lambda x: x > 5),stream.filter(lambda x : x <= 8 )]) >> list
[0, 2, 4, 6, 8]
"""

from __future__ import print_function
import argparse
import stream
from  stream import filter
from deepfilter import *
import gzip

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',dest='input',required=True)
    parser.add_argument('-o','--output',dest='output', required = False, 
                        default = "-")
    parser.add_argument('-l','--max-length',dest='max_length',default=26, type=int)
    parser.add_argument('-L','--long-read-file',dest='long_read_file',
                        default="", type=str)
    parser.add_argument('-q','--mean-quality',dest='mean_quality',default=30, type=int)
    parser.add_argument('-e','--quality-encoding',dest='quality_encoding',default="phred")
    parser.add_argument('-r','--report-file',dest='report_file',required=False,
                        default = 'filter_report.txt')

    args = parser.parse_args()

    total = Counter()
    passed = Counter()

    if not args.long_read_file == '':
        ## get file handler to save discarded reads from length filter
        if args.long_read_file.lower().endswith(".gz"):
            long_read_file = gzip.open(args.long_read_file, 'wb')
        else:
            long_read_file = open(args.long_read_file, 'w')
        
        lenfilter      = MaxLenFilter(args.max_length, long_read_file)
    else:
        lenfilter      = MaxLenFilter(args.max_length)

    meanqualfilter = MeanQualityFilter(args.mean_quality)
    posfilter      = MorePositionsLessThan(2,20)
    
    writer = FastqWriter(args.output)

    dsopen(args.input,args.quality_encoding) >> total >> filter(meanqualfilter) >> filter(posfilter) >> filter(lenfilter) >> passed  >> writer
    
    filters = [meanqualfilter, posfilter, lenfilter]
    
    if not args.long_read_file == '':
        long_read_file.close()

    with open(args.report_file, "w") as report_file:
        for filter_ in filters:
            print("%d reads discarded because %s" % (filter_.discarded,filter_.message), 
                  file = report_file)
    
        print("%d passed" % passed.n, file = report_file)
        print("%d total"  % total.n, file = report_file)
        print("%2.2f%% keeped" % ( (passed.n/float(total.n))*100 ), 
                file = report_file)


