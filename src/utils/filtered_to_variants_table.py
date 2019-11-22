#!/usr/bin/env python

import argparse
import sys
import pickle
from rna import PreResultSet, Mature, build_pre_to_mature_table #, MIRANDMORE_HOME
import pdb

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",dest="input",required=True)
    parser.add_argument("-o","--output",dest="output",required=True)
    args = parser.parse_args()

    CATEGORIES = ["exact","shorter_or_longer","mismatch_1","mismatch_2"]

    with open(args.input, 'rb') as f:
        pre_result_set = pickle.load(f)

    with open(args.output,"w") as out:
        header=["pre","mature","seq","start","end","category","type","count"]
        out.write(";".join(header)+"\n")
        for pre_name, pre in pre_result_set.items():
            for mature_name, mature in pre.items():
                for category in CATEGORIES:
                    _d = getattr(mature,category)
                    #try:
                    for seq, count in _d.items():
                        type_ = seq[3]
                        sequence = seq[0].decode()
                        start = str(seq[1])
                        end   = str(seq[2])
                        line = [str(pre_name), 
                                str(mature_name),
                                str(sequence),
                                str(start),
                                str(end),
                                str(category),
                                str(type_),
                                str(count)]
                        out.write(";".join(line)+"\n")
                    #except:
                        #pdb.set_trace()

if __name__ == '__main__':
    main()
