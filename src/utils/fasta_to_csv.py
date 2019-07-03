#!/usr/bin/env python 

import sys
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help="input file",required=True)
parser.add_argument('-o','--output',help="output file",required=True)
args = parser.parse_args()

i = 0
records = []
for line in open(args.input):
    if (i % 2 == 0):
        name = line.strip().replace(">","")
    else:
        seq = line.strip()
        records.append({'name':name,
                        'seq':seq.upper()})
    i+=1


writer = csv.DictWriter(open(args.output,"w"),fieldnames=['name','seq'])
writer.writeheader()
for record in records:
    writer.writerow(record)


