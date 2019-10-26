#!/usr/bin/env python

import os
import sys
import pdb
import rna
import argparse
import pickle
import itertools
from RNAMachine import rna_generator
from rna import AssignmentError, TwoMatureAssignmentError
# from asq.initiators import query
from functools import partial
#from IPython.Shell import IPShellEmbed
#ipshell = IPShellEmbed([]) 


class Processor:

    ### legacy constructor (deprecated)
    #def __init__(self,base_name,pre,exact,mature_table,annotations):
    #    self.table = rna.build_pre_to_mature_table(mature_table)
    #    self.pre_map = rna.build_pre_to_mature_map(mature_table)
    #    self.base_name = base_name
    #    self.data = {}
    #    self.garbage = []
    #    self.oracle = rna.Oracle(mature_table,annotations)
    #    self.pre = pre         
    #    self.exact = exact 
    #    self.process()
    #    self.post_process_pre_blob()

    #def __init__(self, base_name, pre, exact, mature_table, annotations,
    def __init__(self, base_name, exact, mature_table, annotations,
                 MIN_MORNA_LEN, SPECIES, NEW_RNA_OBJECT_THRESHOLD, SISTER_OVERHANG_LEN, 
                 SISTER_MATCHES_THRESHOLD, ALLOWED_OVERHANG):
        ## the (extended) precursors annotations
        self.table = rna.build_pre_to_mature_table(mature_table)
        ## the (extended) precursors annotations as a dictionary (?)
        self.pre_map = rna.build_pre_to_mature_map(mature_table)
        self.base_name = base_name ## sample name to name output files
        self.data = {}
        #self.garbage = [] ## not used
        self.oracle = rna.Oracle(mature_table, annotations)
    	#self.pre = pre ## the alignments to hairpins as processed by process_algn.py
        self.exact = exact ## the exact alignemnts as processed by dump_exact_blob.py
        
        ## parameters for rules
        self.NEW_RNA_OBJECT_THRESHOLD= NEW_RNA_OBJECT_THRESHOLD
        self.MIN_MORNA_LEN           = MIN_MORNA_LEN
        self.SPECIES                 = SPECIES
        self.SISTER_OVERHANG_LEN     = SISTER_OVERHANG_LEN
        self.SISTER_MATCHES_THRESHOLD= SISTER_MATCHES_THRESHOLD
        self.ALLOWED_OVERHANG        = ALLOWED_OVERHANG

        ## the actual initializator of object's fields 
        self.process()
        
        ## write the python-friendly output file
        #self.post_process_pre_blob()

    def process(self):
        with_assignment_issues = 0

        ## load exact alignments (as processed by dump_exact_blob.py)
        with open(self.exact, 'rb') as b:
            blob = pickle.load(b)

        ## open the assignment log file
        with open(self.base_name + "_assign.log", "w") as log:
            ## for each exact alignment
            for item in blob.items():
                name, pre_blob = item
                presummary = rna.PreSummary(name, self.table, self.oracle, 
                                            self.MIN_MORNA_LEN, self.SPECIES, 
                                            self.SISTER_OVERHANG_LEN, 
                                            self.SISTER_MATCHES_THRESHOLD,
                                            self.ALLOWED_OVERHANG)
                try:
                    presummary.populate(rna_generator(pre_blob, self.NEW_RNA_OBJECT_THRESHOLD))
                except AssignmentError as e:
                    #log.write(e.message + "\n")
                    log.write(str(e) + "\n")
                    with_assignment_issues += 1
                except TwoMatureAssignmentError as e:
                    #log.write(e.message + "\n")
                    log.write(str(e) + "\n")
                    with_assignment_issues += 1
                finally:
                    presummary.clean()
                
                self.data[name] = presummary
            log.write("%d rna objects with assignment issues\n" % with_assignment_issues)
        
        ## write the precursors_summary data in a python-friendly file format
        with open(self.base_name + "_pre_summary.blob", "wb") as presummary_blob:
            pickle.dump(self.data, presummary_blob)

#    def is_not_exact(self,name,entry):
#        seq, start, end = entry[0]
#        for mature in self.table[name]:
#            if mature.start == start and mature.end == end:
#                return False
#        return True
#
#
#    def post_process_pre_blob(self):
#        with open(self.pre) as f:
#            b = pickle.load(f)
#        
#        mir_attrs = ["five_prime_mir", "three_prime_mir"]
#        for name, presummary in self.data.items():
#            for mir_attr in mir_attrs:
#                mir = getattr(presummary, mir_attr)
#                if mir and mir.name in self.pre_map[name]:
#                    b[name][mir.name].shorter_or_longer = {}
#                    f = partial(self.is_not_exact, name)
#                    reads = filter(f, mir.assigned_reads)
#                    b[name][mir.name].shorter_or_longer = dict(reads)
#
#        with open(self.base_name + "_pre_processed.blob", "w") as out:
#            pickle.dump(b,out)
        
        
    def number_of_matures(self,pre):
        return len(self.table[pre])


    def pad_list(self,l,m):
        n = len(l)
        delta = m-n
        if delta > 0:
            for i in xrange(0,delta):
                l.append(0)
        return l


    #def write_excel_table(self):
    #    with open(self.base_name+"_mir_table_excel.txt","w") as f:
    #        header = ["pre","mor-5p","mir-5p","loop","mir-3p","mor-3p","n_known_matures"]
    #        f.write(";".join(header)+"\n")
    #        for pre, presummary in self.data.items():
    #            items=[presummary.stdName]
    #            for attr in ["five_prime_mor", "five_prime_mir", "loop", 
    #                         "three_prime_mir","three_prime_mor"]:
    #                obj = getattr(presummary, attr)
    #                if obj:
    #                    items.append(str(obj.count))
    #                else:
    #                    items.append("0")
    #            items.append(str(presummary.n_matures))
    #            f.write(";".join(items)+"\n")
                
                
    def write_normalized_table(self):
        with open(self.base_name + "_mir_table.txt", "w") as f:
            header=["pre", "name", "category", "start", "end", "expression",
                    "new", "n_known_matures"]
            f.write(";".join(header).replace('"','') + "\n")
            for pre, presummary in self.data.items():
                for attr in ["five_prime_mor", "five_prime_mir", "loop", 
                             "three_prime_mir","three_prime_mor"]:
                    obj = getattr(presummary, attr)
                    if obj:
                        items = [presummary.stdName, obj.name, attr, str(obj.start), 
                                 str(obj.end), str(obj.count), str(obj.new).upper(), 
                                 str(presummary.n_matures)]
                        f.write(";".join(items) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('base_name')
    #parser.add_argument('-p','--pre',dest="pre", required = True)
    parser.add_argument('-e','--exact',dest="exact", required = True)
    parser.add_argument('-m','--mature-table',dest="mature_table", required = True)
    parser.add_argument('-a','--annotations',dest='annotations', required = True)
    
    parser.add_argument('-w', '--new_rna_object_threshold', dest = 'NEW_RNA_OBJECT_THRESHOLD', 
                        type = float, required = True)
    parser.add_argument('-l', '--min_morna_length', dest = 'MIN_MORNA_LEN', 
                        type = int, required = True)
    parser.add_argument('-s', '--species_prefix', dest = 'SPECIES', 
                        type = str, required = True)
    parser.add_argument('-o', '--sister_overhang_len', dest = 'SISTER_OVERHANG_LEN', 
                        type = int, required = True)
    parser.add_argument('-t', '--sister_matches_threshold', dest = 'SISTER_MATCHES_THRESHOLD', 
                        type = int, required = True)
    parser.add_argument('-g', '--allowed_overhang', dest = 'ALLOWED_OVERHANG', 
                        type = int, required = True,
                        help = 'The allowed overhang (either '\
                               'on the left and/or right) of alignment when '\
                               'assigning reads to precursor elements')

    args = parser.parse_args()

    NEW_RNA_OBJECT_THRESHOLD= args.NEW_RNA_OBJECT_THRESHOLD
    MIN_MORNA_LEN           = args.MIN_MORNA_LEN
    SPECIES                 = args.SPECIES
    SISTER_OVERHANG_LEN     = args.SISTER_OVERHANG_LEN
    SISTER_MATCHES_THRESHOLD= args.SISTER_MATCHES_THRESHOLD
    ALLOWED_OVERHANG        = args.ALLOWED_OVERHANG

    #processor = Processor(args.base_name, args.pre, args.exact, 
    #                      args.mature_table, args.annotations)
    #processor = Processor(args.base_name, args.pre, args.exact, 
    processor = Processor(args.base_name, args.exact, 
                          args.mature_table, args.annotations,
                          MIN_MORNA_LEN, SPECIES, NEW_RNA_OBJECT_THRESHOLD, 
                          SISTER_OVERHANG_LEN, SISTER_MATCHES_THRESHOLD,
                          ALLOWED_OVERHANG)
    processor.write_normalized_table()
    #processor.write_excel_table()


if __name__ == '__main__':
    main() 
