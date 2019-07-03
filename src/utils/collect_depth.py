#!/usr/bin/env python 

import re
import os
import sys
import csv
import argparse

def extract_number(line):
    return re.search("([0-9]+)",line).groups()[0]

def extract_sample_name(filepath):
    return os.path.basename(filepath).replace("_trimming_report.txt","")

def process_single_file(filepath):
    the_file = open(filepath) 
    the_file.readline()
    the_file.readline()
    the_file.readline()
    input_reads = extract_number(the_file.readline())
    output_reads = extract_number(the_file.readline())
    sample_name = extract_sample_name(filepath)
    the_file.close()
    return {'sample':sample_name,'input_reads':input_reads,'output_reads':output_reads}
     



parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help="input file",required=True)
parser.add_argument('-o','--output',help="output file",required=True)
args = parser.parse_args()

the_files = args.input.split()

out = open(args.output,'w')
fieldnames = ['sample','input_reads','output_reads']
writer = csv.DictWriter(out,fieldnames=fieldnames)
writer.writeheader()
for filepath in the_files:
    writer.writerow(process_single_file(filepath))

