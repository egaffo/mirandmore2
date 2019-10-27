#!/usr/bin/env python

## This program populates precursor objects by scanning read alignments
## to the miRNA precursors and the annotation from the (known) miRs.
## Only alignments with no mismatches (exact alignments) are considered.
## In addition, alignments with more than the max number (given as a 
## parameter in previous scripts for the multiple hits filter) of mappings
## in the whole genome outside miRNA precursors will be discarded.
## Results in form of PreExactResulSet object are saved (serialized) 
## in a bytecode file (the legacy sample_exact.blob).

import HTSeq, os, sys, pickle, argparse
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
from rna import build_pre_to_mature_table, PreExactResultSet
from coroutines import SamPump, SamCountFilter, MultipleHitsGenomicFilter, Filter
from pipeline import Pipeline

def count_field(results, i):
    return len(set(map(itemgetter(i), results)))

class Processor(Filter):

    def __init__(self, outfile, mature_table, pre_exact_result_set):
        Filter.__init__(self)
        self.table = build_pre_to_mature_table(mature_table)
        self.outfile = outfile
        self.pre_exact_result_set = pre_exact_result_set

    def __call__(self):
        ## consume the input passed from the pipeline
        while True:
            results = []
            ## get the alignments to a precursor from the input 
            item = yield
            name, alnmts = item
            for alnmt in alnmts:
                if alnmt.aligned and alnmt.iv.strand == "+":
                    n_of_mismatches = alnmt.optional_field("NM") 
                    sequence    = alnmt.read.seq
                    pre_name    = alnmt.iv.chrom
                    start       = alnmt.iv.start+1
                    end         = alnmt.iv.end
                    
                    base_data = (sequence, pre_name, start, end)
                    if n_of_mismatches == 0:
                        results.append(base_data)
            if results:
                n_pre = count_field(results, 1)
                for result in results:
                    ## the expression computed for miRs is rescaled
                    ## by the number of homolog precursors so the 
                    ## expression counts will sum up to the actual
                    ## number of reads (not alignments).
                    ## This scale will probably be drop in future
                    ## TODO: do not divide by the number of multiple
                    ## precursor alignments
                    ##self.pre_exact_result_set.add(result + (1,))
                    self.pre_exact_result_set.add(result + (1/float(n_pre),))

    def dump(self):
        exact_out = open(self.outfile,"wb")
        pickle.dump(self.pre_exact_result_set, exact_out)
                

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mature-table", 
                        dest = "mature_table",
                        help = 'known mature miRNA annotation in table format'\
                               ' as output by gff2maturetable.py',
                        required = True)
    parser.add_argument("-t", "--threshold", 
                        dest = "threshold",
                        type = int,
                        help = 'Max allowed number of multiple genomic '\
                               'alignments of a read outside miRNA precursors',
                        required = True)
    parser.add_argument("-g", "--genomic-hits", 
                        dest = "genomic_hits",
                        help = 'A serialized object of the multiple genomic '\
                               'alignments, as output by genomic_hits_to_blob.py',
                        required = True)
    parser.add_argument("-i", "--input", 
                        dest = "input",
                        help = 'Read alignments to the miRNA precursors in '\
                               'SAM format. Gzipped SAM is supported',
                        required = True)
    parser.add_argument("-o", "--output", 
                        dest = "outfile",
                        help = 'File name of serialized PreExactResutSet output',
                        required = True)
    parser.add_argument("-c", "--min_count", 
                        dest = "MIN_COUNT", 
                        required = True,
                        help = 'Minimum number of alignments required to '\
                               'cover a miRNA for being considered as '\
                               'expressed (i.e. saved in the serialized '\
                               'output).', 
                        type = int)

    args = parser.parse_args()
    
    sam_reader = SamPump(args.input)
    
    sam_count_filter = SamCountFilter(args.MIN_COUNT)
    
    multiple_hits_filter = MultipleHitsGenomicFilter(args.threshold, 
                                                     args.genomic_hits)
    
    processor = Processor(args.outfile, 
                          args.mature_table, 
                          PreExactResultSet())
    processor.enqueue('dump')
    
    pipe = Pipeline(sam_reader, 
                    sam_count_filter, 
                    multiple_hits_filter, 
                    processor)
    
    pipe.run()

if __name__ == '__main__':
    main()
