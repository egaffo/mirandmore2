#!/usr/bin/env python

import rna
import argparse
import HTSeq
import pdb
import os
from operator import itemgetter
from itertools import groupby

def basename(name):
    return ".".join(name.split(".")[:-1])

def clean_results(results):
    cleaned_results=[]
    for name, mirs  in groupby(results,itemgetter(0)):
        mirs = list(mirs)
        expression = sum(map(itemgetter(2),mirs))
        seq = mirs[0][1]
        result=[name,seq,expression]
        cleaned_results.append(result)
    return cleaned_results

def save_results(results,name):
    with open(name,"w") as out:
        out.write("name;sequence;expression\n")
        for result in results:
            name,seq,expression = result
            out.write("%s;%s;%.1f\n" % (name,seq,expression))
    

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract sequence of novel miRNAs and moRNAs')
    parser.add_argument('-i','--input',dest="input",required=True)
    parser.add_argument('-t','--threshold',dest="threshold",required=True,type=int)
    parser.add_argument('-m','--mirna',dest="mirna",required=True)
    parser.add_argument('-n','--morna',dest="morna",required=True)
    parser.add_argument('-x', '--hairpin_extended_to_fold', 
                        dest = 'HAIRPIN_EXTENDED_TO_FOLD', required = True)
    args = parser.parse_args()

    hairpins_fasta = open(args.HAIRPIN_EXTENDED_TO_FOLD,"r")
    
    seq_dict = {}
    for hairpin in HTSeq.FastaReader(hairpins_fasta):
        seq_dict[hairpin.name]= hairpin.seq
        

    mir_results = []
    mor_results  = []
    base = os.path.basename(args.input)
    base = basename(base)
    
    with open(args.input) as table:
        table.readline()
        for line in table:
            line = line.strip()
            pre,name,category,start,end,expression,new,n_known_matures = line.split(";")
            start = int(start)
            end  = int(end)
            expression = float(expression)
            if new == "TRUE":
                new = True
            else:
                new = False
            if category in ["five_prime_mir","three_prime_mir"] and new and expression > args.threshold:
                seq = seq_dict[pre+"-ext"][start-1:end].upper()
                result = [name,seq,expression]
                mir_results.append(result)
            if category in ["five_prime_mor","three_prime_mor"] and expression > args.threshold:
                seq = seq_dict[pre+"-ext"][start-1:end].upper()
                result = [name,seq,expression]
                mor_results.append(result)
    mir_results = sorted(mir_results,key=itemgetter(0))
    mor_results = sorted(mor_results,key=itemgetter(0))
           
        # pdb.set_trace()
    mir_results = clean_results(mir_results)
    mir_results = sorted(mir_results,key=itemgetter(2),reverse=True)
    save_results(mir_results,args.mirna)

    mor_results = clean_results(mor_results)
    mor_results = sorted(mor_results,key=itemgetter(2),reverse=True)
    save_results(mor_results,args.morna)
