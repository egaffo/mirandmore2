#!/usr/bin/env python

'''
Script which collapses a fastq read file into a fasta file with unique and counted reads.
The ids of the fasta are in a mirdeep2 compatible format:
>seq_1_x2 1 = unique sequence number	2 = number of reads with that sequence

input:
 * -i --input		fastq file
 * -b --basename	Three letter prefix which identifies a sample
output:
 * -o --output		fasta file
'''

import argparse
import HTSeq

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',dest="input",required=True)
parser.add_argument('-b','--basename',dest="basename",required=True)
parser.add_argument('-o','--output',dest="output",required=True)
args = parser.parse_args()

uniques={}
for s in HTSeq.FastqReader(args.input):
    if s.seq in uniques:       
        uniques[s.seq]+=1
    else:
        uniques[s.seq]=1

progressive_num=1
with open(args.output,"w") as out:
 for sequence in uniques.keys():
     name = "{0}_{1}_x{2}".format(args.basename, progressive_num, uniques[sequence]) 
     fasta_read = HTSeq.Sequence(name=name,seq=sequence)
     progressive_num+=1
     fasta_read.write_to_fasta_file(out)

