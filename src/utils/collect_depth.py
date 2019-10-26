#!/usr/bin/env python 

import re
import os
import sys
import csv
import argparse

def extract_number(line):
    return re.search("([0-9]+)",line).groups()[0]

def extract_sample_name(filepath):    
    return os.path.basename(filepath).replace("_trimming_report.txt",
                                       "").replace('_cutadapt.log',
                                                   '')

## cutadapt log
#This is cutadapt 2.5 with Python 3.5.2
#Command line parameters: -j 8 -a AGATCGGAAGAGCACACGTCT --no-indels -n 10 etc.
#Processing reads on 8 cores in single-end mode ...
#Finished in 69.53 s (4 us/read; 13.70 M reads/minute).
#
#=== Summary ===
#
#Total reads processed:              15,870,553
#Reads with adapters:                15,630,006 (98.5%)
#Reads that were too short:             150,855 (1.0%)
#Reads that were too long:            6,953,440 (43.8%)
#Reads written (passing filters):     8,766,258 (55.2%)
#
#Total basepairs processed:   793,527,650 bp
#Total written (filtered):    244,256,148 bp (30.8%)
#
#=== Adapter 1 ===
#
#etc.

## fastx clipper report
#Clipping Adapter: AGATCGGAAGAGCACACGTCT
#Min. Length: 15
#Non-Clipped reads - discarded.
#Input: 15870553 reads.
#Output: 15606767 reads
#discarded 141407 too-short reads.
#discarded 11996 adapter-only reads.
#discarded 43714 non-clipped reads.
#discarded 66669 N reads


def process_single_file(filepath):
    
    with open(filepath, 'r') as the_file:
        if 'cutadapt' in os.path.basename(filepath):
            ## process cutadapt log
            for line in the_file:
                if 'Total reads processed' in line:
                    input_reads = extract_number(line.replace(',', ''))
                if 'Reads written' in line:
                    output_reads = extract_number(line.replace(',', ''))

        elif 'trimming_report' in os.path.basename(filepath):
            ## parse fastx clipper report
            the_file.readline()
            the_file.readline()
            the_file.readline()
            input_reads = extract_number(the_file.readline())
            output_reads = extract_number(the_file.readline())

    sample_name = extract_sample_name(filepath)

    return {'sample': sample_name,
            'input_reads': input_reads,
            'output_reads': output_reads}
     
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="input file",required=True)
    parser.add_argument('-o','--output',help="output file",required=True)
    args = parser.parse_args()
    
    the_files = args.input.split()
    
    with open(args.output,'w') as out:
        fieldnames = ['sample','input_reads','output_reads']
        writer = csv.DictWriter(out, fieldnames = fieldnames)
        writer.writeheader()
        for filepath in the_files:
            writer.writerow(process_single_file(filepath))
        
