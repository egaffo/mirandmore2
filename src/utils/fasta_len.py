#!/usr/bin/env python

from __future__ import print_function
import sys, gzip
from Bio import SeqIO

'''Reads a (multiple) FASTA file and gives the length of the sequence(s)'''
'''Usage: fasta_len.py mysequences.fa > mysequences.genome'''
ext = str(sys.argv[1]).rpartition('.')[-1]
if ext.lower() == 'gz':
    FastaFile = gzip.open(sys.argv[1], 'rb')
else:
    FastaFile = open(sys.argv[1], 'r')

for rec in SeqIO.parse(FastaFile, 'fasta'):
    name = rec.id
    seq = rec.seq
    seqLen = len(rec)
    print(name + '\t' + str(seqLen))

FastaFile.close()
