#!/usr/bin/env python

## This program populates precursor objects by scanning read alignments
## to the miRNA precursors and the annotation from the (known) miRs.
## Only alignments with no mismatches (exact alignments) are considered.
## In addition, alignments with more than the max number (given as a 
## parameter in previous scripts for the multiple hits filter) of mappings
## in the whole genome outside miRNA precursors will be discarded.
## Results in form of PreExactResulSet object are saved (serialized) 
## in a bytecode file (the legacy sample_exact.blob).

import HTSeq, os, sys, pickle, argparse, gzip
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
from rna import build_pre_to_mature_table, PreExactResultSet
from coroutines import SamPump, SamCountFilter, MultipleHitsGenomicFilter, Filter
from pipeline import Pipeline

import multiprocessing as mp
from multiprocessing import Pool, Value
from collections import Counter
from HTSeq._HTSeq import *
import pickle

## ----------- (old) serial code ----------------
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
                #n_pre = count_field(results, 1)
                for result in results:
                    ## the expression computed for miRs is rescaled
                    ## by the number of homolog precursors so the 
                    ## expression counts will sum up to the actual
                    ## number of reads (not alignments).
                    #self.pre_exact_result_set.add(result + (1/float(n_pre),))
                    ## NEW: do not divide by the number of multiple
                    ## precursor alignments
                    self.pre_exact_result_set.add(result + (1,))

    def dump(self):
        exact_out = open(self.outfile,"wb")
        pickle.dump(self.pre_exact_result_set, exact_out)
        
## --------------- parallel code --------------
class ReadTag():
    def __init__(self, sequence, read_count = 0, alignment_position = None):
        self.tag = sequence
        self.read_count = read_count
        self.alignment_positions = Counter()
        self.alignment_positions[alignment_position] += 1

    def update(self, size_increase, alignment_position):
        self.read_count += size_increase
        self.alignment_positions[alignment_position] += 1

class ReadTags():
    def __init__(self):
        self.sequenceTags = {}

    def addTag(self, tag_tuple):
        try:
            sequence = tag_tuple[0]
            alignment_position = tag_tuple[1]
            try:
                self.sequenceTags[sequence].update(1, alignment_position)
            except KeyError:
                self.sequenceTags[sequence] = ReadTag(sequence, 1,
                                                      alignment_position)
        except TypeError:
            pass

    def getTag(self, tag):
        if tag in self.sequenceTags:
            return self.sequenceTags[tag]
        else:
            return None

    def getAllTags(self):
        return self.sequenceTags.keys()

def samline_to_tuple(line):
    try:
        alnmt = SAM_Alignment.from_SAM_line(line)
    except ValueError:
        return None
    if alnmt.aligned and alnmt.iv.strand == "+":
        n_of_mismatches = alnmt.optional_field("NM") 
        sequence    = alnmt.read.seq.decode()
        pre_name    = alnmt.iv.chrom
        start       = alnmt.iv.start+1
        end         = alnmt.iv.end
        
        if n_of_mismatches == 0:
            return (sequence, (pre_name, start, end))
        else:
            return None

class SAMReader():
    def __init__(self, infilename):
        if infilename.endswith('.gz'):
            self.infile = gzip.open(infilename, 'rt')
        elif infilename == '-':
            self.infile = sys.stdin
        else:
            self.infile = open(infilename, 'r')

    def __iter__(self):
        for line in self.infile:
            yield line

class Processor2():
    def __init__(self, pre_exact_result_set, filename, MIN_COUNT,
                 genomic_filter, jobs):
        self.pre_exact_result_set = pre_exact_result_set
        self.filename = filename
        self.MIN_COUNT = MIN_COUNT
        self.genomic_filter = pickle.load(open(genomic_filter[0], 'rb'))
        self.GENMULTI_THR = genomic_filter[1]
        self.genomic_filter_dist = Counter()
        self.genomic_filter_discarded = Counter()
        self.read_tag_filter = Counter()
        self.jobs = jobs

    def parallel_process(self):
    
        sam_file = SAMReader(self.filename)
        
        pool = mp.Pool(processes = self.jobs)
        
        readTags = ReadTags()
        for result in pool.imap(samline_to_tuple, 
                                (alnmt for alnmt in sam_file), 
                                chunksize = 50000):

            readTags.addTag(result)

        pool.close()
        pool.join()

        for sequence in (tag for tag in readTags.getAllTags() if tag):
            rt = readTags.getTag(sequence)
            ## apply min seq tag count filter
            if rt.read_count >= self.MIN_COUNT:
                self.read_tag_filter['passed'] += rt.read_count
                ## apply multi gen hit filter
                pre = set()
                for aln in rt.alignment_positions.keys():
                    try:
                        pre.add(aln[0])
                    except TypeError:
                        continue

                delta = self.genomic_filter.get(sequence, 0) -\
                        len(pre)
                if delta <= self.GENMULTI_THR:
                    ## update gen filter obj for summary statistics
                    self.genomic_filter_dist[delta] += rt.read_count
                    ## populate the PreExactResultSet
                    for alignment in rt.alignment_positions:
                        base_data = (rt.tag.encode(), alignment[0], alignment[1],
                                     alignment[2])
                        self.pre_exact_result_set.add(base_data + (rt.read_count,))
                else:
                    self.genomic_filter_discarded[delta] += rt.read_count
            else:
                self.read_tag_filter['discarded'] += rt.read_count
        

    def write_filters_stats(self, filename):

        def format_dict_as_table(d, sep = '\t'):
            return '\n'.join([sep.join((str(l), str(r))) for l,r in d.items()])

        with open(filename, 'w') as out:
            print(self.read_tag_filter['discarded'], 
                  '\treads discarded because their sequence was represented by <', 
                  str(self.MIN_COUNT), 'reads', 
                  file = out)
            print(self.read_tag_filter['passed'], 
                  '\treads kept because their sequence was represented by <', 
                  str(self.MIN_COUNT), 'reads, of these', 
                  file = out)

            print(str(sum(self.genomic_filter_discarded.values())),
                  '\treads discarded because aligned in >',
                  str(self.GENMULTI_THR),
                  'genomic loci outside miRNA precursors',
                  file = out)
            print(str(sum(self.genomic_filter_dist.values())), 
                  '\treads kept because aligned in >', 
                  str(self.GENMULTI_THR), 
                  'genomic loci outside miRNA precursors',
                  file = out)

            print('non_mir_loci\treads', file = out)
            print(format_dict_as_table(self.genomic_filter_dist),
                  file = out)
            print(format_dict_as_table(self.genomic_filter_discarded),
                  file = out)


    def serialize(self, outfile):
        with open(outfile, "wb") as exact_out:
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
    parser.add_argument("-j", "--jobs",
                        dest = 'jobs',
                        default = 1,
                        required = False,
                        type = int,
                        help = 'The number of parallel jobs')
    parser.add_argument("-L", "--legacy",
                        dest = 'legacy',
                        default = False,
                        required = False,
                        type = bool,
                        help = 'Run the legacy serial code')

    args = parser.parse_args()
    

    if args.legacy:
        ## run the (old) serial code
        multiple_hits_filter = MultipleHitsGenomicFilter(args.threshold, 
                                                         args.genomic_hits)
    
        multiple_hits_summary_file = os.path.join(os.path.dirname(args.outfile),
                                                      "multiple_hits_summary.txt")

        multiple_hits_filter.enqueue("write_summary",
                                     multiple_hits_summary_file)
                                     
        sam_reader = SamPump(args.input)
        
        sam_count_filter = SamCountFilter(args.MIN_COUNT)
        
                
        processor = Processor(args.outfile, 
                              args.mature_table, 
                              PreExactResultSet())
        processor.enqueue('dump')
        
        pipe = Pipeline(sam_reader, 
                        sam_count_filter, 
                        multiple_hits_filter, 
                        processor)
        
        pipe.run()
    
    else:
        ## run parallelized code
        processor = Processor2(PreExactResultSet(),
                               args.input,
                               args.MIN_COUNT,
                               (args.genomic_hits, args.threshold),
                               args.jobs)
        
        processor.parallel_process()
        
        processor.serialize(args.outfile)

        ## save statistics on multi hits gen filter
        processor.write_filters_stats(os.path.join(os.path.dirname(args.outfile),
                                           "filters_stats.txt"))


if __name__ == '__main__':
    main()
