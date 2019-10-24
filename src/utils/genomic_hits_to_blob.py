#!/usr/bin/env python

import argparse
from itertools import groupby
from rna import reverse_complement
#import cPickle as pickle
import pickle

def align_generator(filename):
    with open(filename,"r") as f:
        for line in f:
            line=line.strip()
            items = line.split("\t")
            strand = items[1]
            seq    = items[4]
            if strand == "-":
                seq = reverse_complement(seq)
            yield seq

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",dest="input",required=True)
    parser.add_argument("-o","--output",dest="output",required=True)
    args = parser.parse_args()

    the_align_generator = align_generator(args.input)
    d_ = {}
    for k,l in groupby(the_align_generator):
        d_[k]=len(list(l))

    with open(args.output,"wb") as out:
        pickle.dump(d_,out)


if __name__ == '__main__':
    main()
