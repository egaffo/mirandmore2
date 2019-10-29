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
import pprint
from HTSeq._HTSeq import *

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
                n_pre = count_field(results, 1)
                for result in results:
                    ## the expression computed for miRs is rescaled
                    ## by the number of homolog precursors so the 
                    ## expression counts will sum up to the actual
                    ## number of reads (not alignments).
                    ## This scale will probably be drop in future
                    ## TODO: do not divide by the number of multiple
                    ## precursor alignments
                    self.pre_exact_result_set.add(result + (1,))
                    #self.pre_exact_result_set.add(result + (1/float(n_pre),))

    def dump(self):
        exact_out = open(self.outfile,"wb")
        pickle.dump(self.pre_exact_result_set, exact_out)
        #[pprint.pprint(dict(per.data)) for per in self.pre_exact_result_set]
        
## --------------- parallel code --------------
class ReadTag():
    def __init__(self, sequence, read_count = 0, alignment_position = None):
        self.tag = sequence
        self.read_count = read_count
        self.alignment_positions = Counter()
        self.alignment_positions[alignment_position] += 1
        #if alignment_position:
        #    self.alignment_positions.add(alignment_position)

    def update(self, size_increase, alignment_position):
        self.read_count += size_increase
        #if alignment_position:
        #    self.alignment_positions.add(alignment_position)
        self.alignment_positions[alignment_position] += 1

class ReadTags():
    def __init__(self):
        self.sequenceTags = {}

    def addTag(self, tag_tuple):
        try:
            ## tag_tuple == (sequence, alignment_position)
            sequence = tag_tuple[0]
            alignment_position = tag_tuple[1]
            #if sequence in self.sequenceTags:
            #    self.sequenceTags[sequence].update(1, alignment_position)
            #else:
            #    self.sequenceTags[sequence] = ReadTag(sequence, 1,
            #                                          alignment_position)
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
    #alnmt = SAMLine(line.rstrip().split('\t'))
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

#def sum_cigar(CIGAR):
#    ## computes the reference bases consumed 
#    ## by the CIGAR operators
#    ## CIGAR: CIGAR string. The CIGAR operations are given in the following table (set
#    ## ‘*’ if unavailable):
#    ## Op   BAM Description                             Consumes    query   reference
#    ## M    0   alignment match (can be a sequence match or mismatch)  yes  yes
#    ## I    1   insertion to the reference                             yes  no
#    ## D    2   deletion from the reference                            no   yes
#    ## N    3   skipped region from the reference                      no   yes
#    ## S    4   soft clipping (clipped sequences present in SEQ)       yes  no
#    ## H    5   hard clipping (clipped sequences NOT present in SEQ)   no   no
#    ## P    6   padding (silent deletion from padded reference)        no   no
#    ## =    7   sequence match                                         yes  yes
#    ## X    8   sequence mismatch                                      yes  yes
#
#    count = 0
#    digits = []
#    for char in CIGAR:
#        if char.isdigit():
#            digits.append(char)
#        if char in ('M', 'D', 'N', '=', 'X'):
#            count += int(''.join(digits))
#            digits = []
#
#    return count
#
#class SAMiv():
#    def __init__(self, chrom, start, CIGAR, flag):
#        self.chrom = chrom
#        self.start = int(start) - 1
#        # sum up int values only for M D N = X CIGAR operations        
#        self.end = self.start + sum_cigar(CIGAR)
#        self.strand = '-' if int(flag) & 0x0010 else '+'
#
#def optFieldTodict(optfields):
#    ## converts 'NM:i:2 XA:i:2' like optional fields into
#    ## a dictionary. The value type identifier is dropped.
#    ## TODO: type convert the values accordig to the type identifier
#    return dict([[x.split(':')[i] for i in [0, 2]] for x in optfields])
#
#class SAMRead():
#    def __init__(self, seq):
#        self.seq = seq
#
#class SAMLine():
#    def __init__(self, samline):
#        self.aligned = (int(samline[1]) & 0x0004) == 0
#        self.iv = SAMiv(samline[2], samline[3], samline[5], samline[1])
#        self.read = SAMRead(samline[9])
#        self.optional_fields = optFieldTodict(samline[11:])
#
#    def optional_field(self, field):
#        try:
#            f = self.optional_fields[field]
#            try:
#                f = int(f)
#            except ValueError:
#                pass
#        except KeyError:
#            f = 0
#            
#        #if field in self.optional_fields:
#        #    f = self.optional_fields[field]
#        #else:
#        #    ## TODO: just a mock value good for NM. Not valid of other fields
#        #    f = None
#        #if field == 'NM':
#        #    f = int(f) if f else 0
#        return f

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
            #if line.startswith("@"):
            #    continue
            #yield SAMLine(line.rstrip().split('\t'))
            yield line

    #def __iter__(self):
    #    return self
    #
    #def __next__(self):
    #    line = next(self.infile).decode().split('\t')
    #    while len(line) < 11:
    #        line = next(self.infile).decode().split('\t')
    #    
    #    return SAMLine(line)

class Processor2():
    def __init__(self, pre_exact_result_set, filename, MIN_COUNT,
                 genomic_filter, jobs):
        self.pre_exact_result_set = pre_exact_result_set
        self.filename = filename
        self.MIN_COUNT = MIN_COUNT
        #self.genomic_filter = genomic_filter
        self.genomic_filter = True
        self.jobs = jobs

    def parallel_process(self):
    
        #sam_file = HTSeq.SAM_Reader(self.filename)
        sam_file = SAMReader(self.filename)
        
        pool = mp.Pool(processes = self.jobs)
        
        readTags = ReadTags()
        for result in pool.imap(samline_to_tuple, 
                                (alnmt for alnmt in sam_file), 
                                chunksize = 50000):

            readTags.addTag(result)

        ## using map, this code  does not work...
        #map(readTags.addTag, pool.imap(samline_to_tuple,
        #                               (alnmt for alnmt in sam_file),
        #                               chunksize = 30000))
        
        pool.close()
        pool.join()

        for sequence in readTags.getAllTags():
            rt = readTags.getTag(sequence)
            ##TODO: set proper genomic hits filter
            if rt.read_count >= self.MIN_COUNT and self.genomic_filter:
                for alignment in rt.alignment_positions:
                    #self.pre_exact_result_set.add((rt.tag.encode(), 
                    #                               alignment[0],
                    #                               alignment[1],
                    #                               alignment[2]) +
                    #                              (rt.read_count,))
                    base_data = (rt.tag.encode(), alignment[0], alignment[1],
                                 alignment[2])
                    self.pre_exact_result_set.add(base_data + (rt.read_count,))



    def serialize(self, outfile):
        with open(outfile, "wb") as exact_out:
            pickle.dump(self.pre_exact_result_set, exact_out)
        #[pprint.pprint(dict(per.data)) for per in self.pre_exact_result_set]


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

    args = parser.parse_args()
    

    #import HTSeq
    #sam_reader = HTSeq.SAM_Reader(args.input)
    #c = 0
    #for l in sam_reader:
    #    if l.aligned:
    #        print(l)
    #        c += 1
    #        if c > 5:
    #            #exit(0)
    #            break

    #c = 0
    #for l in SAMReader(args.input):
    #     if l.aligned:
    #        print(l.iv.chrom, l.iv.start, l.iv.end, l.iv.strand)
    #        c += 1
    #        if c > 5:
    #            exit(0)

    if args.jobs < 2:

        ## run the (old) serial code
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
    
    else:
        ## run parallelized code
        genomic_filter = True ##TODO

        processor = Processor2(PreExactResultSet(),
                               args.input,
                               args.MIN_COUNT,
                               genomic_filter,
                               args.jobs)
        
        processor.parallel_process()
        
        processor.serialize(args.outfile)


if __name__ == '__main__':
    main()
