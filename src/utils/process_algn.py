#!/usr/bin/env python
import HTSeq
import sys
import os
import csv
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
#import cPickle as pickle
import pickle
import re
from rna import MatureResultSet, PreResultSet, Mature, build_pre_to_mature_table#, MIRANDMORE_HOME
#from config import ALLOWED_OVERHANG, ALLOWED_OVERHANG_MORNA, EXACT, SHORTER_OR_LONGER, MIS_1, MIS_2, FIVE_PRIME, THREE_PRIME, MIN_COUNT
from coroutines import  SamPump, SamCountFilter, MultipleHitsGenomicFilter, Meter, Filter
from pipeline import Pipeline
import argparse
import pdb


def overlap(alnmt, mature):
    return abs(alnmt.iv.start+1 - mature.start) <= ALLOWED_OVERHANG and \
           abs(alnmt.iv.end - mature.end) <= ALLOWED_OVERHANG

#def is_morna(mature):
#    return mature.name.find("moR")!=-1

def extract(results, category):
    try:
        idx = list(map(itemgetter(0),results)).index(category)
    except ValueError:
        idx = -1
    return idx

#def is_exact_morna(alnmt, mature):
#    mismatches = alnmt.optional_field("NM")
#    morna_len  = alnmt.iv.end - alnmt.iv.start
#
#    five_prime_morna_exact = mature.order == "5p" and \
#                             abs(alnmt.iv.end - mature.start) <= ALLOWED_OVERHANG_MORNA and \
#                             alnmt.iv.chrom == mature.pre and \
#                             mismatches == 0 and \
#                             morna_len >= MIN_MORNA_LEN
#
#    three_prime_morna_exact = mature.order == "3p" and \
#                              abs(alnmt.iv.start + 1 - mature.end) <= ALLOWED_OVERHANG_MORNA and \
#                              alnmt.iv.chrom == mature.pre and \
#                              mismatches == 0 and \
#                              morna_len >= MIN_MORNA_LEN
#
#    if five_prime_morna_exact or three_prime_morna_exact:
#        return True
#
#    return False
#
#def is_mismatch_1_morna(alnmt, mature):
#    mismatches = alnmt.optional_field("NM")
#    morna_len  = alnmt.iv.end - alnmt.iv.start
#
#    five_prime_morna_exact = mature.order == "5p" and \
#                             abs(alnmt.iv.end - mature.start) <= ALLOWED_OVERHANG_MORNA and \
#                             alnmt.iv.chrom == mature.pre and \
#                             mismatches == 1 and \
#                             morna_len >= MIN_MORNA_LEN
#    three_prime_morna_exact = mature.order == "3p" and \
#                              abs(alnmt.iv.start + 1 - mature.end) <= ALLOWED_OVERHANG_MORNA and \
#                              alnmt.iv.chrom == mature.pre and \
#                              mismatches == 1 and \
#                              morna_len >= MIN_MORNA_LEN
#
#    if five_prime_morna_exact or three_prime_morna_exact:
#        return True
#
#    return False


def is_exact_mirna(alnmt, mature):
    if alnmt.iv.start + 1 == mature.start and \
       alnmt.iv.end == mature.end and \
       alnmt.iv.chrom == mature.pre and \
       alnmt.optional_field("NM") == 0:
        return True
    return False

def is_shorter_or_longer_mirna(alnmt, mature):
    if overlap(alnmt,mature) and \
       alnmt.iv.chrom == mature.pre and \
       alnmt.optional_field("NM") == 0:
        return True
    return False

def is_mismatch_1_mirna(alnmt, mature):
    if alnmt.optional_field("NM") == 1 and \
       overlap(alnmt,mature) and \
       alnmt.iv.chrom == mature.pre :
        return True
    return False


def is_mismatch_2_mirna(alnmt, mature):
    md = alnmt.optional_field("MD")
    ## two mismatches are allowed only on 3p end
    three_prime_addiction_re = "[0-9]+[ACTG]0[ACTG]0"
    if alnmt.optional_field("NM") == 2 and \
       overlap(alnmt,mature) and \
       alnmt.iv.chrom == mature.pre and \
       re.match(three_prime_addiction_re,md):
        return True
    return False


def count_field(results, i):
    return len(set(map(itemgetter(i), results)))


class Processor(Filter):

    MD_RE = re.compile("(\d+)[ACTG]*\d*")

    def __init__(self, base_name,mature_table, mature_result_set, pre_result_set):
        Filter.__init__(self)
        self.table     = build_pre_to_mature_table(mature_table)
        self.base_name = base_name
        self.mature_result_set = mature_result_set
        self.pre_result_set    = pre_result_set
        
    def dump_stats(self, count_filter):
        #pdb.set_trace()
        output = open(self.base_name + "_stats.txt","w")
        output.write("%s;%d\n" % ("reads kept", count_filter.keeped()))
        output.write("rejected because number of reads per tag < %d;%d\n" \
                     % (MIN_COUNT, count_filter.rejected()))
        output.write("hairpins with at least one read mapped;%d\n" \
                     % (len(self.pre_result_set.keys())))
        output.close()

    def __call__(self):
        
        CATEGORY = 0
        MATURE   = 1
        SEQUENCE = 2
        PRE = 3
        
        while True:
            results=[]
            item = yield
            name, alnmts = item
            for alnmt in alnmts:
                ## we use only alignments in the forward strand, since
                ## alignments on the reverse strand might be reads coming
                ## from the complementary sequence on the opposite arm
                ## of the precursor
                if alnmt.aligned and alnmt.iv.strand == "+":
                    sequence    = alnmt.read.seq
                    pre_name    = alnmt.iv.chrom
                    start       = alnmt.iv.start + 1
                    end         = alnmt.iv.end
                    ## get number of mismatches in the alignment
                    n_of_mismatches = alnmt.optional_field("NM")

                    ## set a bool var to check whether the alignment match any sRNA
                    was_assigned = False
                    ## restrict the loop to sRNA defined on the reference precursor
                    matures = self.table[alnmt.iv.chrom]
                    for mature in matures:
                        mature_name = mature.name
                        variant_end = ""
                        base_data = (mature_name, sequence, pre_name, start, end)

                        previous_assignments = []
                        if was_assigned:
                            for result in results:
                                if result[3] == pre_name:
                                    previous_assignments.append(result)

                        if  is_mismatch_1_mirna(alnmt, mature):
                            if was_assigned:
                                sys.stderr.write("WARNING: alignment " + \
                                                 str(alnmt) + \
                                                 " was already assigned to " + \
                                                 str(previous_assignments) +\
                                                 " but it also fits " + \
                                                 str(mature) + \
                                                 " with mismatch_1\n")
                            before_mismatch = Processor.MD_RE.match(alnmt.optional_field("MD")).groups()[0]
                            if int(before_mismatch) < alnmt.iv.length/2:
                                variant_end = "5p"
                            else:
                                variant_end = "3p"

                            results.append( (MIS_1,) + base_data + (variant_end,))
                            was_assigned = True

                        elif is_mismatch_2_mirna(alnmt, mature):
                            if was_assigned:
                                sys.stderr.write("WARNING: alignment " + \
                                                 str(alnmt) + \
                                                 " was already assigned to " + \
                                                 str(previous_assignments) +\
                                                 " but it also fits " + \
                                                 str(mature) + \
                                                 " with mismatch_2\n")

                            variant_end = "3p"
                            results.append( (MIS_2,) + base_data + (variant_end,))
                            was_assigned = True

                        elif is_exact_mirna(alnmt, mature):
                            if was_assigned:
                                sys.stderr.write("WARNING: alignment " + \
                                                 str(alnmt) + \
                                                 " was already assigned to " + \
                                                 str(previous_assignments) +\
                                                 " but it also fits " + \
                                                 str(mature) + \
                                                 " exactly \n")

                            results.append( (EXACT,) + base_data + (variant_end,))
                            was_assigned = True

                        elif is_shorter_or_longer_mirna(alnmt, mature):
                            if was_assigned:
                                sys.stderr.write("WARNING: alignment " + \
                                                 str(alnmt) + \
                                                 " was already assigned to " + \
                                                 str(previous_assignments) +\
                                                 " but it also fits " + \
                                                 str(mature) + \
                                                 " with shorter_or_longer\n")

                            if start != mature.start:
                                variant_end = "5p"
                            if end != mature.end:
                                variant_end = "3p"
                            if start != mature.start and end != mature.end:
                                variant_end = "both"
                            
                            #pdb.set_trace()
                            results.append( (SHORTER_OR_LONGER,) + base_data + (variant_end,))
                            was_assigned = True

                        ## output also moRNA variants
                        #elif is_exact_morna(alnmt, mature):
                        #    ## TODO: fix moR name (let miRs, etc.)
                        #    #sys.stderr.write("Exact moR in " + str(mature) + "\n")
                        #    mor_name = mature_name.replace("miR", "moR")
                        #    base_data = (mor_name, sequence, pre_name, start, end)
                        #    results.append( (EXACT,) + base_data + (variant_end,))
                        #    was_assigned = True
                        
                        #elif is_mismatch_1_morna(alnmt, mature):
                        #    ## TODO: fix moR name (let miRs, etc.)
                        #    mor_name = mature_name.replace("miR", "moR")
                        #    base_data = (mor_name, sequence, pre_name, start, end)
                        #    before_mismatch = Processor.MD_RE.match(alnmt.optional_field("MD")).groups()[0]
                        #    if int(before_mismatch) < alnmt.iv.length/2:
                        #        variant_end = "5p"
                        #    else:
                        #        variant_end = "3p"
                        
                        #    results.append( (MIS_1,) + base_data + (variant_end,))
                        #    #sys.stderr.write("Mismatch_1 moR in " + str(mature) + "\n")
                        #    was_assigned = True

                if alnmt.iv.strand == "+" and not was_assigned:
                    ## TODO: set a logging at DEBUG level. Mind that the
                    ## following output can be as large as the alignment file
                    #sys.stderr.write("WARNING: sequence {0} unassigned.\n".format(str(alnmt)))
                    pass

                if alnmt.iv.strand != "+":
                    ## TODO: set a logging at DEBUG level. However, this check is
                    ## unnecessary by the logic of alignmets. See the above comment
                    ## on this variable.
                    #sys.stderr.write("WARNING: alignment strand not handled "\
                    #                 "{0}: {1}\n".format(str(alnmt),
                    #                                     str(alnmt.iv.strand)))
                    pass


            if results:
                categories = list(set(map(itemgetter(0), results)))
                if len(categories) > 1:
                    ## (? is the code below interpreted right ?)
                    ## give a priority to equally scored alignments:
                    ## exact alignments are preferred over all the others,
                    ## and lower number of mismatches is preferred
                    exact = extract(results, EXACT)
                    if exact >= 0:
                        results = [results[exact]]
                    
                    shorter_or_longer = extract(results, SHORTER_OR_LONGER)
                    if shorter_or_longer >= 0:
                        results =  [results[shorter_or_longer]]

                    mis_1 = extract(results, MIS_1)
                    if mis_1 >= 0:
                        results = [results[mis_1]]

                    mis_2 = extract(results, MIS_2)
                    if mis_2 >= 0:
                        results = [results[mis_2]]

                n_matures = count_field(results, MATURE)
                n_pre = count_field(results, PRE)
                categories = list(set(map(itemgetter(0), results)))
                if len(categories) > 1:
                    sys.stderr.write("ERROR process_algn.py: "\
                                     "multiple isoform categories "\
                                     "for mature {0}: {1} -> {2}".format(set(map(itemgetter(1), results)),
                                                                         str(categories),
                                                                         str(results)))
                    #exit(1)
                for result in results:
                    self.mature_result_set.add(result + (1/float(n_matures),))
                    self.pre_result_set.add(result + (1/float(n_pre),))

    def dump(self):
        matures_out = open(self.base_name + "_matures.blob", "wb")
        pre_out     = open(self.base_name + "_pre.blob", "wb")
        pickle.dump(self.mature_result_set, matures_out)
        pickle.dump(self.pre_result_set, pre_out)
        matures_out.close()
        pre_out.close()




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mature-table", dest="mature_table",required=True)
    parser.add_argument("-t", "--threshold", dest="threshold",type=int,required=True)
    parser.add_argument("-g", "--genomic-hits", dest="genomic_hits",required=True)
    parser.add_argument("-i", "--input", dest="input",required=True)
    parser.add_argument('-b', '--outfilebasename', dest = 'outbase', required = False,
                        help = 'A prefix for the output files name.')

    parser.add_argument("-a", "--allowed_overhang", dest = "ALLOWED_OVERHANG", required = True,
                        help = "", type = int)
    parser.add_argument("-e", "--exact", dest = "EXACT", required = True,
                        help = "", type = int)
    parser.add_argument("-l", "--shorter_or_longer", dest = "SHORTER_OR_LONGER", required = True,
                        help = "", type = int)
    parser.add_argument("-s", "--mis1", dest = "MIS_1", required = True,
                        help = "", type = int)
    parser.add_argument("-S", "--mis2",dest = "MIS_2", required = True,
                        help = "", type = int)
    parser.add_argument("-f", "--five_prime", dest = "FIVE_PRIME", required = True,
                        help = "", type = int)
    parser.add_argument("-r", "--three_prime", dest = "THREE_PRIME", required = True,
                        help = "", type = int)
    parser.add_argument("-c", "--min_count", dest = "MIN_COUNT", required = True,
                        help = "", type = int)

    ## These parameters about moRNAs are actually not used by this script's functions
    parser.add_argument("-o", "--morna_allowed_overhang", dest = "ALLOWED_OVERHANG_MORNA",
                        required = True, help = "", type = int)
    parser.add_argument("-n", "--min_morna_len", dest = "MIN_MORNA_LEN", required = True,
                        help = "", type = int)

    args = parser.parse_args()

    ## Setting these variables as global is an ugly hack to avoid the import of
    ## a config file and let these variable be visible to all the methods here defined
    global ALLOWED_OVERHANG
    global EXACT
    global SHORTER_OR_LONGER
    global MIS_1
    global MIS_2
    global FIVE_PRIME
    global THREE_PRIME
    global MIN_COUNT
    global ALLOWED_OVERHANG_MORNA
    global MIN_MORNA_LEN

    ALLOWED_OVERHANG = args.ALLOWED_OVERHANG
    EXACT = args.EXACT
    SHORTER_OR_LONGER = args.SHORTER_OR_LONGER
    MIS_1 = args.MIS_1
    MIS_2 = args.MIS_2
    FIVE_PRIME = args.FIVE_PRIME
    THREE_PRIME = args.THREE_PRIME
    MIN_COUNT = args.MIN_COUNT
    ALLOWED_OVERHANG_MORNA = args.ALLOWED_OVERHANG_MORNA
    MIN_MORNA_LEN = args.MIN_MORNA_LEN

    base_file = args.input
    base_name = base_file.split(".")[0]
    if args.outbase:
        base_name = args.outbase
    sam_reader = SamPump(base_file)
    processor = Processor(base_name, args.mature_table, MatureResultSet(), PreResultSet())
    sam_count_filter = SamCountFilter(MIN_COUNT)
    multiple_hits_filter = MultipleHitsGenomicFilter(args.threshold, args.genomic_hits)
    #meter = Meter()
    #processor = Processor(base_name, args.mature_table, MatureResultSet(), PreResultSet())
    #processor = Processor(base_name, args.mature_table, MatureResultSet(), PreResultSet())
    #sam_count_filter = SamCountFilter(it,MIN_COUNT)
    #multiple_hits_filter = MultipleHitsGenomicFilter(sam_count_filter,args.threshold,args.genomic_hits)
    #for name, alnmts in multiple_hits_filter:
        #processor.process(alnmts)
    #processor.dump()
    #dump_stats(base_name,sam_count_filter,processor)
    #multiple_hits_summary = multiple_hits_filter.summary()
    #multiple_hits_summary.write_to_file(base_name+"_multiple_hits_summary.txt")
    processor.enqueue("dump")
    processor.enqueue("dump_stats", sam_count_filter)
    multiple_hits_filter.enqueue("write_summary", base_name + "_multiple_hits_summary.txt")
    #pipe = Pipeline(sam_reader, meter, sam_count_filter,  multiple_hits_filter, processor)
    pipe = Pipeline(sam_reader, sam_count_filter,  multiple_hits_filter, processor)
    pipe.run()

if __name__ == "__main__":
    main()
