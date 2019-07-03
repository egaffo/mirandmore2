#!/usr/bin/env python 

'''
Script which generates a gff3 files with the predicted precursors and matures
from a output.mrd file (miRDeep2_core_algorithm.pl output) and a precursors.coords
file (output of excise_precursors_iterative_final.pl)
##parse_mrd_to_gff.py
input:
 * -i --input	       output of miRDeep2_core_algorithm.pl
 * -c --coordinates    file with the genomic coordinates of the excised precursors
 * -p --prefix         prefix to use for matures and precursors names
output:
 * -o --output     gff3 file with the predicted precursors and matures
'''

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help="input .mrd file",required=True)
parser.add_argument('-c','--coordinates',help="input coordinates file",required=True)
parser.add_argument('-o','--output',help="output gff3 file",required=True)
parser.add_argument('-p','--prefix',help="naming prefix",required=True)
args = parser.parse_args()

coords_in = open(args.coordinates)

def get_branch(pre_start ,pre_end, strand, mature_start):
	length = pre_end - pre_start + 1
	if mature_start - pre_start + 1 < length / 2.0:
		if strand == '+':
			branch = '-5p'
		else:
			branch = '-3p'
	else:
		if strand == '+':
			branch = '-3p'
		else:
			branch = '-5p'
	return branch

def find_coords(coords_file,pre_id):
    '''
         Finds the chromosome, start, end and strand in the precursors.coords file
         format: >name strand start end (1 based but the start of a chromosome is 0)
    '''     
    coords_file.seek(0)                         # restart reading from the beginning of the file
    chrom, start, end, strand = "XXX", "XXX", "XXX", "XXX"
    for line in coords_file:
        line=line.rstrip()
        ls = line.split()
        if ls[0]==('>'+pre_id):       
            ls = line.split()
            chrom = ls[0][:ls[0].rfind('_')][1:]
            strand = ls[1]
            start = str(max(int(ls[2]),1))      # 1 based 
            end = ls[3]                         # 1 based 
            break
    return  chrom, start, end, strand
           
def parse_mature(line):
    '''
        parses the 'exp' output.mrd line to extract the offsets of the mature
        start and end from the start of the precursor:
        format:     exp     ffMMllSSff
    ''' 
    line = line.split()[1]          
    start_offset = line.find('M')   # get the 0 based index of the first occurrence of 'M'
    end_offset = line.rfind('M')    # get the 0 based index of the last occurrence of 'M'
    return start_offset, end_offset   

def make_gff_lines(basename, pre_id, chrom, pre_start ,pre_end, start_offset,end_offset, score, strand, number):
    pre_name = basename + '-mir-' + pre_id
    pre_attributes = "ID=PRE" + str(number) + ";Name=" + pre_name
    pre_list = [chrom, 'mm','miRNA_primary_transcript', str(pre_start), str(pre_end), score, strand,'.',pre_attributes ]
    pre_line = '\t'.join(pre_list) + '\n' 
    mature_start =  int(pre_start) + start_offset          # 1 based start + 0 based start_offset = 1 based start
    mature_end = int(pre_start) + end_offset               # 1 based start + 0 based end_offset = 1 based end
    mature_name = basename + '-miR-' + pre_id + get_branch(int(pre_start), int(pre_end), strand, mature_start)
    mature_attributes = "ID=MIR" + str(number) + ";Name=" + mature_name + ";Derives_from=PRE" + str(number)
    mature_list = [chrom, 'mm','miRNA', str(mature_start), str(mature_end), score, strand,'.',mature_attributes ]
    mature_line = '\t'.join(mature_list) + '\n' 
    return pre_line, mature_line


with open(args.output,"w") as out_file: 
    progressive_number = 1
    for line in open(args.input):
        line=line.rstrip()
        if line.startswith('>'):                #  format:  >chr1_1
            pre_id = line[1:]
        elif line.startswith('score total'): 
            score = line.split()[2]             #  format: score total    12
        elif line.startswith('exp'): 
            start_offset, end_offset = parse_mature(line)
            chrom, pre_start, pre_end, strand = find_coords(coords_in,pre_id)
            pre_line, mature_line = make_gff_lines(args.prefix, pre_id, chrom, pre_start ,pre_end,
                                                    start_offset, end_offset, score, strand, progressive_number)
            out_file.writelines([pre_line,mature_line])
            progressive_number+=1

coords_in.close()

