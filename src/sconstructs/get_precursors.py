'''
This Sconcript uses mirdeep2 to predict putative hairpin precursors:

When called from a SConscript it imports the following variables:
 * env
 * basename = basename used to name files
 * genome = genome in fasta format
 * alignments = genomic alignments in Bowtie format
 
Returns:
  *results = [
  	* genomic_arf = genomic alignments in the mirdeep2 compatible .arf format
 	* parsed_arf = .arf file with the alignments parsed according to the mirdeep2 default parameters
  	* putative_precursors=[
		basename+"_precursors.fa",	  	= precursors fasta sequences
		basename+"_precursors.coords",		= precursors genomic coordinates and strand
		basename+"_precursors.log",		= log with the number of precursors excised for each minimum read number threshold
		]
'''

import os

Import('*')

try:
    env = env.Clone()
    basename = get_precursors_BASENAME
    env['SAMPLE'] = get_precursors_BASENAME
    genome = get_precursors_GENOME_FASTA
    env['GENOME_FASTA'] = get_precursors_GENOME_FASTA
    alignments = get_precursors_ALIGNMENTS

except NameError as ne:
    vars = Variables('vars.py')
    #basename
    vars.Add('SAMPLE', 'basename', 'mm_excised_precursors')
    #genome sequence in fasta format
    vars.Add('GENOME_FASTA', 'genome in fasta format', 
             'MIRANDMORE_HOME/dbs/bowtie-indexes/hg38.fa')
    #genomic alignments in bowtie format
    vars.Add('ALIGNMENTS', 'genomic alignments in Bowtie format', 
             'pooled_genomic.out')

    env = Environment(variables = vars,ENV = os.environ)
    Help(vars.GenerateHelpText(env))

    basename = env['SAMPLE']
    genome = File(env['GENOME_FASTA'])
    alignments = File(env['ALIGNMENTS'])

#SRC_DIR = os.path.join(env['ENV']['MIRANDMORE_HOME'],'scons','prj')
SRC_DIR = env['SCONSCRIPT_HOME']

#mirdeep2 SConscripts
mirandmore_bwt_to_arf        = 'bwt_to_arf.py'
mirandmore_parse_genomic_arf = 'parse_genomic_arf.py'
mirandmore_excise_precursors = 'excise_precursors.py'

results = []

#mirdeep2 SConscripts calls
arf_dir = 'arfs' 
excised_precursors_dir = 'excised_precursors'

bwt_to_arf_BASENAME = basename
bwt_to_arf_ALIGNMENTS = alignments
genomic_arf = SConscript(os.path.join(arf_dir, mirandmore_bwt_to_arf),
                         variant_dir = arf_dir, src_dir = SRC_DIR,
                         duplicate = 0, 
                         exports = 'env bwt_to_arf_BASENAME bwt_to_arf_ALIGNMENTS')
results.append(genomic_arf)

parse_genomic_arf_genomic_arf = genomic_arf
parse_genomic_arf_basename = basename
parsed_arf = SConscript(os.path.join(arf_dir, mirandmore_parse_genomic_arf),
                        variant_dir = arf_dir, src_dir = SRC_DIR,
                        duplicate = 0, 
                        exports = 'env parse_genomic_arf_basename parse_genomic_arf_genomic_arf')
results.append(parsed_arf)

env_excise_precursors = env.Clone()
env_excise_precursors['ARF'] = parsed_arf
putative_precursors = SConscript(os.path.join(excised_precursors_dir, 
                                              mirandmore_excise_precursors),
                                 variant_dir = excised_precursors_dir, src_dir = SRC_DIR,
                                 duplicate = 0, 
                                 exports = 'env_excise_precursors')
results.append(putative_precursors)

Return('results') 
