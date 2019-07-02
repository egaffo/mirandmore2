'''
This Sconcript uses parse_mappings.pl (from mirdeep2)to parse a mirdeep2 compatible .arf file with mirdeep2 default parameters:
#parse_mappings.pl:
#   -a 0   0 mismatches allowed
#   -b 18  discard alignments < 17 nt
#   -c 25  discard alignments > 25 nt
#   -i 5   discards all alignments of uniques sequences with more than 5 alignments

When called from a SConscript it imports the following variables:
 * env
 * parse_genomic_arf_genomic_arf = genomic alignments in the mirdeep2 compatible .arf format
 * basename = basename used to name files

Returns:
 * parsed_arf = .arf file with the alignments parsed according to the mirdeep2 default parameters
'''

import os

Import('*')

try:
    env = env.Clone()
    parse_genomic_arf_genomic_arf = parse_genomic_arf_genomic_arf
    parse_genomic_arf_basename = parse_genomic_arf_basename

except NameError:
    vars = Variables('vars.py')
    vars.Add('PARSE_ARF_ALIGNMENTS', 'genomic alignments in .arf format', 
             'basename_genomic.arf')
    vars.Add('PARSE_ARF_SAMPLE', 'basename', 'basename')
    
    env = Environment(ENV=os.environ, variables=vars)
    
    Help(vars.GenerateHelpText(env))
    
    parse_genomic_arf_basename = env['PARSE_ARF_SAMPLE']
    parse_genomic_arf_genomic_arf = File(env['PARSE_ARF_ALIGNMENTS'])

source = parse_genomic_arf_genomic_arf
target = parse_genomic_arf_basename + "_parsed.arf"
parsed_arf = env.Command(target, 
                         source,
                         "parse_mappings.pl ${SOURCE} -a 0 -b 18 -c 25 -i 5 > ${TARGET}")
Return('parsed_arf')
