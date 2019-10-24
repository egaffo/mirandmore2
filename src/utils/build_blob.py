#!/usr/bin/env python

import sys, re, argparse
from rna import PreAnnotation
#import cPickle as pickle
import pickle
import HTSeq

def build_hairpins_annotations(target,source,env):
    #energy_re = re.compile("\(\s*(-\d+\.\d+)\)") ## this regex does not match ' ( 0.00)' limit case
    d = {}
    target = str(target[0])
    source = str(source[0])
    with open(source,"r") as f:
        while True:
            name = f.readline().strip().replace(">","")
            if not name:
                break
            seq  = f.readline().strip().replace("U","T")
            rnafoldline = f.readline().strip()
            structure = rnafoldline.split(" ")[0]
            #energy = float(energy_re.findall(rnafoldline)[0])
            energy = float(re.sub('\(|\)', '', ''.join(rnafoldline.split(' ')[1:])))
            pre = PreAnnotation(name=name,seq=seq,structure=structure,energy=energy)
            d[name] = pre
            
    with open(target,"w") as f:
        pickle.dump(d,f)

## TODO: implement using Bio.SeqIO instead of HTSeq to reduce requirements
def build_fa_blob(target,source,env):
    d = {}
    target = str(target[0])
    source = str(source[0])
    
    with open(source, "r") as f:
        stream = HTSeq.FastaReader(f)
        
        for entry in stream:
            d[entry.name] = entry.seq
            
    with open(target, "wt") as f:
        pickle.dump(d, f)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input', dest = "infile", required = True, 
                        help = '''The input file to convert into blob''')
    parser.add_argument('-t','--type', dest = "filetype", 
                        required = True, 
                        help = '''The type of file to convert. Either a FASTA or a '''\
                        '''folding structure file as output by RNAfold (ViennaRNA Package) ''', 
                        default = 'fasta', choices=['fasta', 'folding'])
    parser.add_argument('-o','--output', dest = "outfile", required = True,
                        help = '''Output blob file name''')

    args = parser.parse_args()

    env = {} # dummy environment
    if args.filetype == 'fasta':
        build_fa_blob([args.outfile], [args.infile], env)
    else:
        build_hairpins_annotations([args.outfile], [args.infile], env)

