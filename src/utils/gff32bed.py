#!/usr/bin/env python 

import sys, re, argparse

## From a list of GFF3 tagged vales, extract the value of
## the specified key
def get_field_value(field_list, field):
    
    val = None
    for f in field_list:
        if field in f:
            val = re.sub('.*=(.*)', r'\1', f)
            break

    return val


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help = "A GFF3 file")
    args = parser.parse_args()
    
    if args.input is '-':
        f = sys.stdin
    else:
        f = open(args.input, 'r')
    
    for line in f:
        if not line.startswith("#"):
            items = line.strip().split("\t")
            chromosome = items[0]
            entry_type = items[2]
            start      = int(items[3])
            end        = int(items[4])
            strand     = items[6]
            name       = get_field_value(items[8].strip().split(';'), 'Name')
    
            outlist = [chromosome, 
                       str(start-1),
                       str(end),
                       name,
                       '.',
                       strand]
    
            print('\t'.join(outlist))

    f.close()

