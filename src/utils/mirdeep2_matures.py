#!/usr/bin/env python
import re
import argparse

ID_re_str = r"ID=([^;\s]+)"		# extracts the ID from the attributes field of a .gff3 file
ID_re = re.compile(ID_re_str)
Name_re_str = r"Name=([^;\s]+)" # extracts the name from the attributes field of a .gff3 file
Name_re = re.compile(Name_re_str)

def write_line(name, ID, chrom, strand, start, end, score, overlapping, table):
    overlapping_str = ";".join(overlapping)
    line_list = [name, ID, chrom, strand, start, end, score, overlapping_str]
    line = "\t".join(line_list) + "\n"
    table.write(line)

def get_overlap(entry):
	overlap_match = Name_re.search(entry[17])
	if overlap_match:
		overlap = overlap_match.group(1)
		return overlap
	return None

def parse_ID(line):
	entry = line.strip().split()
	ID = ID_re.search(entry[8]).group(1)
	return ID

def parse_line(line):
	entry = line.strip().split() 
	chrom = entry[0]
	start = entry[3]    # 1 based (read from a .gff3)
	end = entry[4]      # 1 based (read from a .gff3)
	strand = entry[6]
	score = entry[5]    # miRDeep2 score         
	name = Name_re.search(entry[8]).group(1)
	overlap = get_overlap(entry)
	ID = ID_re.search(entry[8]).group(1)
	return name, ID, chrom, strand, start, end, score, overlap

def main(in_file,table):
	#header = "name\tID\tchr\tstrand\tstart\tend\tscore\tknown-pre_overlaps\n"
	#table.write(header)
	ID = ""
	overlapping = []
	line = in_file.readline()
	if not line:
		return
	while 1:
		name, old_ID, chrom, strand, start, end, score, overlap = parse_line(line)
		if overlap:
			overlapping.append(overlap)
		line = in_file.readline()
		if not line:
			write_line(name, old_ID, chrom, strand, start, end, score, overlapping, table)
			return
		ID = parse_ID(line)
		while ID == old_ID:	# in case the mature overlaps with more than one precursor
			name, ID, chrom, strand, start, end, score, overlap = parse_line(line)
			overlapping.append(overlap)
			line = in_file.readline()
			if not line:
				write_line(name, ID, chrom, strand, start, end, score, overlapping, table)
				return
			ID = parse_ID(line)
		write_line(name, old_ID, chrom, strand, start, end, score, overlapping, table)
		old_ID = ID
		overlapping = []
         
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates a table with miRDeep2 matures')
    parser.add_argument('-i','--input',dest="input",required=True)
    parser.add_argument('-o','--output',dest="output",required=True)
    args = parser.parse_args()
    with open(args.output, 'w') as table:
        with open(args.input) as in_file:
            main(in_file,table)
