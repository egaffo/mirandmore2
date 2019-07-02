'''
This Sconcript uses convert_bowtie_output.pl (from mirdeep2)to converta .bwt alignment file into a mirdeep2 compatible .arf file

When called from a SConscript it imports the following variables:
 * env
 * bwt_to_arf_ALIGNMENTS = genomic alignments in Bowtie format
 * bwt_to_arf_BASENAME = basename used to name files

Returns:
 * genomic_arf = genomic alignments in the mirdeep2 compatible .arf format
'''

import os

Import('*')

try:
    env         = env.Clone()
    alignments  = bwt_to_arf_ALIGNMENTS
    basename	= bwt_to_arf_BASENAME

except NameError:
    vars = Variables('vars.py')
    vars.Add('ALIGNMENTS', 'genomic alignments in Bowtie format', 'basename_genomic.out')
    vars.Add('BASENAME', 'basename', 'basename')
    env = Environment(ENV=os.environ, variables=vars)
    Help(vars.GenerateHelpText(env))
    basename = env['BASENAME']
    alignments = File(env['ALIGNMENTS'])

source = alignments
target = basename + "_genomic.arf"

genomic_arf = env.Command(target, source, "convert_bowtie_output.pl ${SOURCE} > ${TARGET}")

Return('genomic_arf')

