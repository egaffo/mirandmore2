'''
This Sconcript uses excise_precursors_iterative_final.pl (from mirdeep2) to excise precursors from the genome with mirdeep2 default parameters: excise_precursors_iterative_final.pl
#   genome sequence in fasta format
#   parsed genomic mappings in .arf format
#   putative precursor sequences
#   putative precursor genomic coordinates
#   50000 max number of precursors to extract
#   log

When called from a SConscript it imports the following variables:
 * env
 * basename = basename used to name files
 * parsed_arf_to_excise = .arf file with the alignments parsed according to the mirdeep2 default parameters

Returns:
 *putative_precursors=[
	basename+"_precursors.fa",	  	= precursors fasta sequences
	basename+"_precursors.coords",		= precursors genomic coordinates and strand
	basename+"_precursors.log"]		= log with the number of precursors excised for each minimum read number threshold
'''

import os

Import('*')

try:
    env = env_excise_precursors.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('EXCISE_PRECURSORS_PARSED_ALIGNMENTS', 'parsed genomic alignments in .arf format', 
             'basename_parsed.arf')
    vars.Add('EXCISE_PRECURSORS_SAMPLE', 'basename', 'basename')
    vars.Add('GENOME_FASTA', 'genome in fasta format', 'MIRANDMORE_HOME/dbs/bowtie-indexes/hg38.fa')
    
    env = Environment(ENV=os.environ, variables=vars)
    
    Help(vars.GenerateHelpText(env))
    
    env['SAMPLE'] = env['EXCISE_PRECURSORS_SAMPLE']
    env['ARF'] = env['EXCISE_PRECURSORS_PARSED_ALIGNMENTS']
    genome = File(env['GENOME_FASTA'])

targets = [env['SAMPLE'] + "_precursors.fa", 
           env['SAMPLE'] + "_precursors.coords", 
           env['SAMPLE'] + "_precursors.log"]

cmdLine= '''excise_precursors_iterative_final.pl ${SOURCES[0]} '''\
         '''${SOURCES[1]} ${TARGETS[0]} ${TARGETS[1]} 50000 2> ${TARGETS[2]}'''

putative_precursors = env.Command(targets, 
				  [env['GENOME_FASTA'], env['ARF']], 
				  cmdLine)

SideEffect([env['SAMPLE'] + "_precursors.fa_all", 
           env['SAMPLE'] + "_precursors.coords_all", 
           env['SAMPLE'] + "_precursors.fa_stack"], putative_precursors)

Return('putative_precursors')
