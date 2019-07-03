#!/usr/bin/env python

import HTSeq 
import os
import sys
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
import cPickle as pickle
from rna import build_pre_to_mature_table, PreExactResultSet #MIRANDMORE_HOME, 
from coroutines import  SamPump, SamCountFilter, MultipleHitsGenomicFilter, Meter, Filter 
import argparse
from pipeline import Pipeline
import pdb


def count_field(results,i):
    return len(set(map(itemgetter(i),results)))

class Processor(Filter):
    def __init__(self, base_name, mature_table, pre_exact_result_set):
        Filter.__init__(self)
        self.table = build_pre_to_mature_table(mature_table)
        self.base_name = base_name
        self.pre_exact_result_set = pre_exact_result_set

    def __call__(self):
        while True:
            results = []
            item = yield
            name, alnmts = item
            #print "processing %s" % name
            for alnmt in alnmts:
                if alnmt.aligned and alnmt.iv.strand=="+":
                    n_of_mismatches = alnmt.optional_field("NM") 
                    sequence    = alnmt.read.seq
                    pre_name    = alnmt.iv.chrom
                    start       = alnmt.iv.start+1
                    end         = alnmt.iv.end
                    
                    base_data = (sequence, pre_name, start, end)
                    if n_of_mismatches == 0:
                        results.append(base_data)
            if results:
                n_pre = count_field(results,1)
                for result in results:
                    self.pre_exact_result_set.add(result + (1/float(n_pre),))

    def dump(self):
        #exact_out = open(self.base_name+"_exact.blob","wb")
        exact_out = open(self.base_name,"wb")
        pickle.dump(self.pre_exact_result_set,exact_out)
                

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--mature-table",dest="mature_table",required=True)
    parser.add_argument("-t","--threshold",dest="threshold",type=int,required=True)
    parser.add_argument("-g","--genomic-hits",dest="genomic_hits",required=True)
    parser.add_argument("-i","--input",dest="input",required=True)
    parser.add_argument("-o","--output",dest="outfile",required=True)
    parser.add_argument("-c", "--min_count", dest = "MIN_COUNT", required = True,
                        help = "", type = int)

    args = parser.parse_args()
    options = args
    MIN_COUNT = args.MIN_COUNT

    base_file = args.input #args[0]
    base_name = options.outfile #base_file.split(".")[0]
    sam_reader = SamPump(base_file)
    sam_count_filter = SamCountFilter(MIN_COUNT)
    multiple_hits_filter = MultipleHitsGenomicFilter(options.threshold, options.genomic_hits)
    processor = Processor(base_name, options.mature_table, PreExactResultSet())
    #meter = Meter()
    processor.enqueue('dump')
    #pipe = Pipeline(sam_reader, meter, sam_count_filter, multiple_hits_filter, processor)
    pipe = Pipeline(sam_reader, sam_count_filter, multiple_hits_filter, processor)
    pipe.run()

if __name__ == '__main__':
    main()
