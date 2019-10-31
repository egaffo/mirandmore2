#!/usr/bin/env python
import HTSeq, sys, os, csv, pickle, re, argparse, gzip
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
from rna import MatureResultSet, PreResultSet, Mature, build_pre_to_mature_table
from coroutines import SamPump, SamCountFilter, MultipleHitsGenomicFilter, Filter
from pipeline import Pipeline

import multiprocessing as mp
from multiprocessing import Pool, Value
import pprint

def overlap(alnmt, mature, ALLOWED_OVERHANG):
    return abs(alnmt.iv.start+1 - mature.start) <= ALLOWED_OVERHANG and \
           abs(alnmt.iv.end - mature.end) <= ALLOWED_OVERHANG

def extract(results, category):
    try:
        idx = list(map(itemgetter(0),results)).index(category)
    except ValueError:
        idx = -1
    return idx

def is_exact_mirna(alnmt, mature, ALLOWED_OVERHANG):
    if alnmt.iv.start + 1 == mature.start and \
       alnmt.iv.end == mature.end and \
       alnmt.iv.chrom == mature.pre and \
       alnmt.optional_field("NM") == 0:
        return True
    return False

def is_shorter_or_longer_mirna(alnmt, mature, ALLOWED_OVERHANG):
    if overlap(alnmt, mature, ALLOWED_OVERHANG) and \
       alnmt.iv.chrom == mature.pre and \
       alnmt.optional_field("NM") == 0:
        return True
    return False

def is_mismatch_1_mirna(alnmt, mature, ALLOWED_OVERHANG):
    if alnmt.optional_field("NM") == 1 and \
       overlap(alnmt, mature, ALLOWED_OVERHANG) and \
       alnmt.iv.chrom == mature.pre :
        return True
    return False


def is_mismatch_2_mirna(alnmt, mature, ALLOWED_OVERHANG):
    md = alnmt.optional_field("MD")
    ## two mismatches are allowed only on 3p end
    three_prime_addiction_re = "[0-9]+[ACTG]0[ACTG]0"
    if alnmt.optional_field("NM") == 2 and \
       overlap(alnmt, mature, ALLOWED_OVERHANG) and \
       alnmt.iv.chrom == mature.pre and \
       re.match(three_prime_addiction_re,md):
        return True
    return False


def count_field(results, i):
    return len(set(map(itemgetter(i), results)))

##------------- legacy serial code -----------
class Processor(Filter):

    MD_RE = re.compile("(\d+)[ACTG]*\d*")

    def __init__(self, base_name,mature_table, mature_result_set, pre_result_set):
        Filter.__init__(self)
        self.table     = build_pre_to_mature_table(mature_table)
        self.base_name = base_name
        self.mature_result_set = mature_result_set
        self.pre_result_set    = pre_result_set
        self.unassigned_seq = {}
        
    def dump_stats(self, count_filter):
        output = open(self.base_name + "_stats.txt","w")
        output.write("%s;%d\n" % ("reads kept", count_filter.keeped()))
        output.write("rejected because number of reads per tag < %d;%d\n" \
                     % (MIN_COUNT, count_filter.rejected()))
        output.write("hairpins with at least one read mapped;%d\n" \
                     % (len(self.pre_result_set.keys())))
        output.close()

    def dump_unassigned_log(self):
        ## print one line per sequence in a table format where
        ## 1st column is the sequence name, 
        ## 2nd column is the count of multimaps on precursors
        ## 3rd column is the positions of alignments, with mismatches
        ##                   formatted as '|' separated list of
        ##                   precursor:start:end:number-of-mismatches

        with open(self.base_name + '_unassigned_seqs.log', 'w') as out:
            header = '\t'.join(['sequence', 'count', 'alignments_and_mismatches'])
            out.write(header + '\n')

            if bool(self.unassigned_seq):
                for seq in self.unassigned_seq.keys():
                    positions = []
                    for pos in self.unassigned_seq[seq]['positions']:
                        nm = self.unassigned_seq[seq]['positions'][pos]
                        pos_with_nm = [str(tag) for tag in pos] + [str(nm)]
                        positions.append(':'.join(pos_with_nm))

                    line = '\t'.join([str(seq.decode()), 
                                      str(self.unassigned_seq[seq]['count']),
                                      '|'.join(positions)])

                    out.write(line + '\n')


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
                            
                            results.append( (SHORTER_OR_LONGER,) + base_data + (variant_end,))
                            was_assigned = True


                if alnmt.iv.strand == "+" and not was_assigned:
                    
                    ## keep track of unassigned sequences/alignmets
                    position = (pre_name, start, end)
                    
                    if sequence not in self.unassigned_seq:
                        self.unassigned_seq[sequence] = {'count': 1,
                                                         'positions': {position: n_of_mismatches}}

                    else:
                        self.unassigned_seq[sequence]['count'] += 1

                        if position not in self.unassigned_seq[sequence]['positions']:
                            self.unassigned_seq[sequence]['positions'][position] = n_of_mismatches 

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

##-------------- parallel code ----------------

def smartOpen(infilename):
    if infilename.endswith('.gz'):
        return gzip.open(infilename, 'rt')
    elif infilename == '-':
        return sys.stdin
    else:
        return open(infilename, 'r')


def assign_mature(intuple):
    
    line = intuple[0]
    table = intuple[1] ## the mature table with all precursor elements
    ALLOWED_OVERHANG = intuple[2] 
    #MIN_COUNT = intuple[3]

    ## constants
    UNASSIGNED        = -1
    EXACT             = 0 
    SHORTER_OR_LONGER = 1 
    MIS_1             = 2
    MIS_2             = 3 
    FIVE_PRIME        = 5 
    THREE_PRIME       = 6 

    MD_RE = re.compile("(\d+)[ACTG]*\d*")

    try:
        alnmt = HTSeq.SAM_Alignment.from_SAM_line(line)
    except ValueError:
        ## error could raise for SAM header line
        ## and/or not well formatted SAM line
        return None
    
    ## set a bool var to check whether the alignment match any 
    ## (and if > 1) sRNA(s) in the precursor
    was_assigned = False
    result = None

    ## we use only alignments in the forward strand, since
    ## alignments on the reverse strand might be reads coming
    ## from the complementary sequence on the opposite arm
    ## of the precursor. Bear in mind that reference are the
    ## (extended) precursor sequences
    if alnmt.aligned and alnmt.iv.strand == "+":
        sequence    = alnmt.read.seq
        pre_name    = alnmt.iv.chrom
        start       = alnmt.iv.start + 1
        end         = alnmt.iv.end
        ## get number of mismatches in the alignment
        n_of_mismatches = alnmt.optional_field("NM")
        
        ## prepare default result
        base_data = (None, sequence, pre_name, start, end)
        result = (was_assigned,
                  ((UNASSIGNED, ) + base_data + ('', )))

        ## loop through sRNA defined on the reference precursor
        matures = table[pre_name]
        
        ## 'results' will store assignments in the precursor
        ## we expect only one result per alignment.
        ## If multiple results per one alignment, then chose
        ## according to (from most preferred to last):
        ## exact, shorter_or_longer, mismatch 1, mismatch 2
        results = []
        
        ## IDEA for future improvement:
        ## we could chose which sRNA by checking the 
        ## best overlapping among 'matures' defined in the
        ## precursor. Yet, we need to define 'best overlap'
        for mature in matures:
            mature_name = mature.name
            variant_end = ""
            base_data = (mature_name, sequence, pre_name, start, end)
   
            previous_assignments = []
            if was_assigned:
                for r in results:
                    if r[3] == pre_name:
                        previous_assignments.append(r)
    
            if is_mismatch_1_mirna(alnmt, mature, ALLOWED_OVERHANG):
                if was_assigned:
                    sys.stderr.write("WARNING: alignment " + \
                                     str(alnmt) + \
                                     " was already assigned to " + \
                                     str(previous_assignments) +\
                                     " but it also fits " + \
                                     str(mature) + \
                                     " with mismatch_1\n")
                before_mismatch = MD_RE.match(alnmt.optional_field("MD")).groups()[0]
                if int(before_mismatch) < alnmt.iv.length/2:
                    variant_end = "5p"
                else:
                    variant_end = "3p"

                was_assigned = True
                part_result = ( (MIS_1,) + base_data + (variant_end,) )
                results.append(part_result)
                result = (was_assigned, 
                          part_result)
    
            elif is_mismatch_2_mirna(alnmt, mature, ALLOWED_OVERHANG):
                if was_assigned:
                    sys.stderr.write("WARNING: alignment " + \
                                     str(alnmt) + \
                                     " was already assigned to " + \
                                     str(previous_assignments) +\
                                     " but it also fits " + \
                                     str(mature) + \
                                     " with mismatch_2\n")
    
                variant_end = "3p"

                was_assigned = True
                part_result = ((MIS_2,) + base_data + (variant_end,))
                results.append(part_result)
                result = (was_assigned, part_result)

            elif is_exact_mirna(alnmt, mature, ALLOWED_OVERHANG):
                if was_assigned:
                    sys.stderr.write("WARNING: alignment " + \
                                     str(alnmt) + \
                                     " was already assigned to " + \
                                     str(previous_assignments) +\
                                     " but it also fits " + \
                                     str(mature) + \
                                     " exactly \n")
    
                was_assigned = True
                part_result = ((EXACT,) + base_data + (variant_end,))
                results.append(part_result)
                result = (was_assigned,
                          part_result)
   
            elif is_shorter_or_longer_mirna(alnmt, mature, ALLOWED_OVERHANG):
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
                
                was_assigned = True
                part_result = ((SHORTER_OR_LONGER,) + base_data + (variant_end,))
                results.append(part_result)
                result = (was_assigned,
                          part_result)

            ## check if the alignment was assigned to >1 sRNAs in the
            ## precursor, with different isoform types. Then, chose the
            ## 'best' assignment as actual return value.
            ## However, this shold not happen!
            if results and was_assigned:
                categories = list(set(map(itemgetter(0), results)))
                if len(categories) > 1:
                    ## give a priority to equally scored alignments:
                    ## exact alignments are preferred over all the others,
                    ## and lower number of mismatches is preferred
                    for cat in (EXACT, SHORTER_OR_LONGER, MIS_1, MIS_2):
                        if extract(results, cat) >= 0:
                            result = results[extract(results, cat)]
                            break
            
    return result

def fasta_to_dict(fasta_filename):
    ## parse the result of fastq_to_unique_fasta.py and 
    ## return a dictionary with the amount of reads
    ## for each sequence tag
    n_seqtags = defaultdict(int)
    for s in HTSeq.FastaReader(fasta_filename):
        n_seqtags[s.seq] = int(s.name.split('_')[-1].replace('x', ''))

    return n_seqtags        
    
class PProcessor():

    def __init__(self, SAM_pre_aln, outfiles_base_name, mature_table, 
                 ALLOWED_OVERHANG, MIN_COUNT, 
                 genomic_hits,
                 GENHIT_THR, unique_seq_fasta, jobs):

        self.table  = build_pre_to_mature_table(mature_table)
        self.infile = SAM_pre_aln
        self.base_name = outfiles_base_name
        self.ALLOWED_OVERHANG = ALLOWED_OVERHANG
        self.MIN_COUNT = MIN_COUNT
        self.genomic_hits = pickle.load(open(genomic_hits, 'rb'))
        self.GENHIT_THR = GENHIT_THR
        self.jobs = jobs
        #self.mature_result_set = MatureResultSet()
        self.pre_result_set    = PreResultSet()
        self.CATEGORY = 0
        self.MATURE   = 1
        self.SEQUENCE = 2
        self.PRE = 3
        ## 3 level nested defaultdict
        self.unassigned_seq = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        self.stats = defaultdict(int)
        self.seqtag_tab = fasta_to_dict(unique_seq_fasta)

        #print(pprint.pprint(self.genomic_hits))

    def dump_stats(self):
        output = open(self.base_name + "_stats.txt","w")
        output.write("%s : %d\n" % ("alignments kept", self.stats['assigned']))
        output.write('{} : rejected because number of '\
                     'reads per tag < {}\n'.format(self.stats['seqtag_filter'], 
                                                   self.MIN_COUNT))
        output.write('{} : rejected because number of genomic alignments '\
                     'out of miRNA precursors > {}\n'.format(self.stats['genmulti_filter'], 
                                                             self.GENHIT_THR))
        #output.write("hairpins with at least one read mapped;%d\n" \
        #             % (len(self.pre_result_set.keys())))
        output.close()

    def dump_unassigned_log(self):
        ## print one line per sequence in a table format where
        ## 1st column is the sequence name, 
        ## 2nd column is the count of multimaps on precursors
        ## 3rd column is the positions of alignments, with mismatches
        ##                   formatted as '|' separated list of
        ##                   precursor:start:end:number-of-mismatches

        with open(self.base_name + '_unassigned_seqs.log', 'w') as out:
            header = '\t'.join(['sequence', 'count', 'alignment_and_count'])
            out.write(header + '\n')

            for seq in self.unassigned_seq.keys():
                positions = []
                for pos, count in self.unassigned_seq[seq]['position'].items():
                    #nm = self.unassigned_seq[seq]['position'][pos]
                    #pos_with_nm = [str(tag) for tag in pos] + [str(nm)]
                    #positions.append(':'.join(pos_with_nm))
                    positions.append(':'.join([pos, str(count)]))

                line = '\t'.join([str(seq.decode()), 
                                  str(self.unassigned_seq[seq]['count']['count']),
                                  '|'.join(positions)])

                out.write(line + '\n')

    def dump(self):
        #matures_out = open(self.base_name + "_matures.blob", "wb")
        pre_out     = open(self.base_name + "_pre.blob", "wb")
        #pickle.dump(self.mature_result_set, matures_out)
        pickle.dump(self.pre_result_set, pre_out)
        #matures_out.close()
        pre_out.close()

    def parallel_process(self):
       
        sam_file = smartOpen(self.infile)
        pool = mp.Pool(processes = self.jobs)
        
        assignments = []
        for assignment in pool.imap(assign_mature, 
                                ((alnmt, 
                                  self.table,
                                  self.ALLOWED_OVERHANG) for alnmt in sam_file), 
                                chunksize = 10000):
        
            if assignment:
                sequence = assignment[1][2]
                position = (assignment[1][3], 
                            assignment[1][4],
                            assignment[1][5])

                ## check if alignment was assigned
                if assignment[0]:
                    ## apply MINCOUNT per seqtag and filters
                    if self.seqtag_tab[sequence] >= self.MIN_COUNT:
                        assignments.append(assignment[1])
                    else:
                        self.stats['seqtag_filter'] += 1
                else:
                    ## alignment was not assigned
                    pos_str = ':'.join([str(tag) for tag in position])
                    self.unassigned_seq[sequence]['count']['count'] += 1
                    self.unassigned_seq[sequence]['position'][pos_str] += 1

        pool.close()
        pool.join()

        ## get sequences and precursors
        seq_to_pre = defaultdict(set)
        for seq, pre in map(itemgetter(2, 3), assignments):
            seq_to_pre[seq].add(pre)
        
        for s in assignments:
            sequence = s[2]
            n_pre = len(seq_to_pre[sequence])
            n_genhits = self.genomic_hits.get(sequence.decode(), 0)
            if n_genhits - n_pre <= self.GENHIT_THR:
                self.stats['assigned'] += 1
                ## populate result set
                self.pre_result_set.add(s + (1,))
            else:
                self.stats['genmulti_filter'] += 1

            #n_matures = count_field(results, MATURE)
            #n_pre = count_field(results, PRE)
            #categories = list(set(map(itemgetter(0), results)))
            #if len(categories) > 1:
            #    sys.stderr.write("ERROR process_algn.py: "\
            #                     "multiple isoform categories "\
            #                     "for mature {0}: {1} -> {2}".format(set(map(itemgetter(1), results)),
            #                                                         str(categories),
            #                                                         str(results)))
            
            #    for result in results:
            #        self.mature_result_set.add(result + (1/float(n_matures),))
            #        self.pre_result_set.add(result + (1/float(n_pre),))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mature-table", 
                        dest = "mature_table",
                        required = True,
                        help = 'The annotated table of (extended) miRNA precursors\' '\
                               'small RNAs as predited by mir_discretizer.py')
    parser.add_argument("-t", "--threshold", 
                        dest = "GENHIT_THR",
                        type = int,
                        default = 10,
                        required = True,
                        help = 'The allowed number of genomic multiple mappings '\
                               'of a read/sequence outside miRNA precursors')
    parser.add_argument("-g", "--genomic-hits", 
                        dest = "genomic_hits",
                        required = True,
                        help = 'The serialized object of genomic hits')
    parser.add_argument("-i", "--input", 
                        dest = "input",
                        required = True,
                        help = 'Read alignemnts to miRNA precursors file in '\
                               'SAM format')
    parser.add_argument('-b', '--outfilebasename', 
                        dest = 'outbase', 
                        required = False,
                        help = 'A prefix for the output files name.')
    parser.add_argument("-a", "--allowed_overhang", 
                        dest = "ALLOWED_OVERHANG", 
                        required = True,
                        help = "", 
                        type = int)
    parser.add_argument("-c", "--min_count", 
                        dest = "MIN_COUNT", 
                        required = True,
                        help = "The minimum read count allowed for a '\
                               'sequence tag to be considered",
                        type = int)
    parser.add_argument("-u", "--unique_seq_fasta", 
                        dest = "unique_seq_fasta", 
                        required = True,
                        help = 'A FASTA file with unique sequences and reads '\
                               'per sequence in header (must end with _xN, '\
                               'where N is the read count), as output by '\
                               'fastq_to_uniques_fasta.py',
                        type = str)
    parser.add_argument('-j', '--jobs',
                        dest = 'jobs',
                        default = 1,
                        required = False,
                        type = int,
                        help = 'Number of parallel processes')
    parser.add_argument('-L', '--legacy',
                        dest = 'legacy',
                        default = False,
                        required = False,
                        type = bool,
                        help = 'Run legacy code (this is not guaranteed to work)')

    args = parser.parse_args()

    base_name = os.path.basename(args.input).replace('.gz', '').replace('.sam', '')
    if args.outbase:
        base_name = args.outbase
        if os.path.isdir(base_name):
            base_name = os.path.join(base_name, 
                                     os.path.basename(args.input).replace('.gz',
                                     '').replace('.sam', ''))

    if args.legacy:

        sam_reader = SamPump(base_file)
        processor = Processor(base_name, args.mature_table, MatureResultSet(), PreResultSet())
        sam_count_filter = SamCountFilter(args.MIN_COUNT)
        multiple_hits_filter = MultipleHitsGenomicFilter(args.GENHIT_THR, args.genomic_hits)
        processor.enqueue("dump")
        processor.enqueue("dump_stats", sam_count_filter)
        multiple_hits_filter.enqueue("write_summary", base_name + "_multiple_hits_summary.txt")
        processor.enqueue('dump_unassigned_log')
        pipe = Pipeline(sam_reader, sam_count_filter,  multiple_hits_filter, processor)
        pipe.run()

    else:

        processor = PProcessor(args.input, 
                               base_name,
                               args.mature_table,
                               args.ALLOWED_OVERHANG,
                               args.MIN_COUNT,
                               args.genomic_hits,
                               args.GENHIT_THR,
                               args.unique_seq_fasta,
                               args.jobs)
        processor.parallel_process()
        processor.dump()
        processor.dump_stats()
        processor.dump_unassigned_log()



if __name__ == "__main__":
    main()

