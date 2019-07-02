'''
This Sconcript generates the alignments of unique reads vs putative precursors 

When called from a SConscript it imports the following variables:
 * env
 * mirandmore_prepare_signatures_basename = basename used to name files
 * mirandmore_prepare_signatures_collapsed_reads = fasta file with the unique reads in mirdeep2 compatible format
 * mirandmore_prepare_signatures_cpus = number of cpus for bowtie
 * mirandmore_prepare_signatures_precursors = fasta file with the sequences of the excised precursors

Returns:
 *signatures = basename + ".arf"	.arf file with unique reads aligned vs the putative precursors
                               		 sorted by chr, precursor, start and end position		
             	 
'''

import os

Import('*')

try:
    env = env.Clone()
    collapsed_reads = prepare_signatures_COLLAPSED_READS
    precursors = prepare_signatures_PRECURSORS
    cpus = prepare_signatures_CPUS
    basename = prepare_signatures_BASENAME

except NameError:
    vars = Variables('vars.py')
    vars.Add('COLLAPSED_READS',
                'parsed genomic alignments in .arf format', 'basename_parsed.arf')
    vars.Add('EXCISED_PRECURSORS_FASTA',
                'excised precursors in fasta format', 'basename_precursors.fa')
    vars.Add('CPUS', 'cpus', '1')
    vars.Add('SIGNATURES_BASENAME', 'basename', 'mm_signatures')
    
    env = Environment(ENV=os.environ,
                              variables=vars)
    
    Help(vars.GenerateHelpText(env))
    
    collapsed_reads = File(env['COLLAPSED_READS'])
    precursors = File(env['EXCISED_PRECURSORS_FASTA'])
    cpus = env['CPUS']
    basename = env['SIGNATURES_BASENAME']

#SConscripts
mirandmore_bowtie_index = 'bowtie_index.py'
mirandmore_bwt_to_arf = 'bwt_to_arf.py'
#SCONSCRIPT_HOME = os.path.join(env['ENV']['MIRANDMORE_HOME'],'scons','prj')
SCONSCRIPT_HOME = env['SCONSCRIPT_HOME']

# build the bowtie index of the excised precusrors using mirandmore_bowtie_index

bowtie_indexes = []
excised_precursor_index_dir = os.path.join('.', 'bwt_excised_precursors_index')
env_bowtie_index = env.Clone()
env_bowtie_index['FASTA_TO_INDEX'] = precursors
excised_precursor_bowtie_index = SConscript(os.path.join(excised_precursor_index_dir, 
                                                         mirandmore_bowtie_index),
                                 variant_dir = excised_precursor_index_dir,
                                 src_dir = SCONSCRIPT_HOME,
                                 duplicate = 0,
                                 exports = '''env_bowtie_index''')
bowtie_indexes.append(excised_precursor_bowtie_index)

# align the unique reads vs the excised precursors with miRDeep2 default settings for bowtie:
#  -v   1               up to 1 mismatch
#  -f                   input is in fasta format
#  --norc               dont align on the reverse complements
#  -a --best --strata   Make Bowtie guarantee that reported singleton alignments are "best" in terms of stratum
#                       (i.e. number of mismatches, or mismatches in the seed in the case of -n mode) 
#                       and in terms of the quality values at the mismatched position(s). If many valid 
#                       alignments exist and are reportable (e.g. are not disallowed via the -k option) 
#                       and they fall into more than one alignment "stratum", report only those alignments 
#                       that fall into the best stratum.
targets = [basename + '.bwt',
           basename + '_bowtie_log.txt']
sources = [collapsed_reads] + bowtie_indexes
CPUS = cpus
CmdLine = '''bowtie -f -v 1 -a --best --strata --norc $(-p ''' + CPUS + '''$) '''+ \
          '''${str(SOURCES[1].base).split(".")[0]}  ${SOURCES[0]} ${TARGETS[0]} 2>${TARGETS[1]}'''
aligned_uniques_vs_excised_precursors = env.Command(targets, sources, CmdLine)

# convert the alignments from .bwt format to .arf format using mirandmore_bwt_to_arf

bwt_to_arf_ALIGNMENTS = aligned_uniques_vs_excised_precursors[0]
bwt_to_arf_BASENAME = basename + '_unsorted'
excised_precursors_arf = SConscript(os.path.join('.', mirandmore_bwt_to_arf),
                         variant_dir = '.', src_dir = SCONSCRIPT_HOME,
                         duplicate = 0, 
                         exports = 'env bwt_to_arf_BASENAME bwt_to_arf_ALIGNMENTS')

# sort the alignments by chr, precursor, start and end position
# -v    natural sort of (version) numbers within text, used to 
#       sort the alignment by chromosome of origin of the precursor
#       and then by precursor using the precursor name.  

target = basename + '.arf'
source = excised_precursors_arf
CmdLine = 'sort -V -k6,6 -k8,8n -k9,9n ${SOURCE} > ${TARGET}' 
signatures = env.Command(target, source, CmdLine)

Return('signatures')
