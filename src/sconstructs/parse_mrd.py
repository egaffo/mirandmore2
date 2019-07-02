'''
This Sconcript uses parse_mrd.py to generate a gff3 files with the predicted precursors and matures
from a output.mrd file (miRDeep2_core_algorithm.pl output) and a precursors.coords
file (output of excise_precursors_iterative_final.pl)

##parse_mrd_to_gff.py
Script which generates a gff3 files with the predicted precursors and matures
from a output.mrd file (miRDeep2_core_algorithm.pl output) and a precursors.coords
file (output of excise_precursors_iterative_final.pl)

input:
 * -i --input	       output of miRDeep2_core_algorithm.pl
 * -c --coordinates    file with the genomic coordinates of the excised precursors
 * -p --prefix         prefix to use for matures and precursors names
output:
 * -o --output     gff3 file with the predicted precursors and matures

When called from a SConscript it imports the following variables:
 * env
 * parse_mrd_MRD_FILE = output of miRDeep2_core_algorithm.pl
 * parse_mrd_PRECURSORS_COORDS = file with the genomic coordinates of the excised precursors
				                (output of excise_precursors_iterative_final.pl)
 * parse_mrd_BASENAME = basename used to name files

Returns: 
 * precursor_and_mature_gff: .gff3 file with the predicted precursors and matures
'''

import os

Import('*')

try:
    env = env.Clone()
    mrd_file = parse_mrd_MRD_FILE
    precursors_coords = parse_mrd_PRECURSORS_COORDS
    basename = parse_mrd_BASENAME

except NameError:
    vars = Variables('vars.py')
    vars.Add('MRD_FILE', 'output of miRDeep2_core_algorithm.pl', 'mm_output.mrd')
    vars.Add('PRECURSORS_COORDS', 'coordinates of the excised precursors', 'mm_excised_precursors.coords')
    vars.Add('BASENAME', 'basename', 'mm_predicted')
    
    env = Environment(ENV=os.environ,
                              variables=vars)
    
    Help(vars.GenerateHelpText(env))
    
    mrd_file = File(env['MRD_FILE'])
    precursors_coords = File(env['PRECURSORS_COORDS'])
    basename = env['BASENAME']

sources = [mrd_file, precursors_coords]
targets = basename + ".gff3"
cmdLine = "parse_mrd_to_gff.py -i ${SOURCES[0]} -c ${SOURCES[1]} -p putative_new -o ${TARGET}"

precursor_and_mature_gff = env.Command(targets, sources, cmdLine)

Return('precursor_and_mature_gff')

