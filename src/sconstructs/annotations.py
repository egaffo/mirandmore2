'''
Generates the annotation files needed by MiR&moRe:
- process annotation files
- extend annotation coordinates
- process extended annotation

Returns an array with:
* results : a dictionary with 1. processed annotation and 2. extended processed annotation
* extended_precursors : extended precursor FASTA file
'''

import os, itertools

Import('*')

try:
    env = env_annotation.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('PRECURSORS_FASTA', 'The miRNA precursor annotation', 'hairpin.fa.gz')
    vars.Add('GFF_FILE', 'The miRBase miRNA annotation in GFF3 format', 'hsa.gff3')
    vars.Add('CHRM_PREFIX_TO_REMOVE', 'Pattern to remove from chromosome names in annotation', 
             "''")
    vars.Add('EXTENDED_N', 'The number of bases to extend precursors', '30')
    vars.Add('GENOME_FASTA', 'The genome in single FASTA formatted file', 'hg38.fa')

    env = Environment(variables = vars, ENV = os.environ)
 
    Help(vars.GenerateHelpText(env))
 
    PRECURSORS_FASTA        = File(env['PRECURSORS_FASTA'])
    GFF_FILE                = File(env['GFF_FILE'])
    CHRM_PREFIX_TO_REMOVE   = env['CHRM_PREFIX_TO_REMOVE']
    EXTENDED_N                = env['EXTENDED_N']
    GENOME_FASTA            = File(env['GENOME_FASTA'])
    SPECIES                 = env['SPECIES']

#SRC_DIR = os.path.join(env['ENV']['MIRANDMORE_HOME'],'scons','prj')
SRC_DIR = env['SCONSCRIPT_HOME']

## SConscripts
mirandmore_process_annotation = 'process_annotation.py'
mirandmore_extend_hairpins = 'extend_hairpins.py'

## output directories
annotation_dir = 'annotation'
extended_precursors_dir = 'extended_annotation'

## Extend annotation coordinates
mirandmore_extend_hairpins_GFF_FILE = env['GFF_FILE']
mirandmore_extend_hairpins_EXTEND_N = env['EXTENDED_N']
mirandmore_extend_hairpins_GENOME_FASTA = env['GENOME_FASTA']

extended_precursors = SConscript(os.path.join(extended_precursors_dir, 
                                              mirandmore_extend_hairpins), 
                                 variant_dir = extended_precursors_dir, src_dir = SRC_DIR,
                                 duplicate = 0,
                                 exports = '''env '''
                                           '''mirandmore_extend_hairpins_GFF_FILE '''
                                           '''mirandmore_extend_hairpins_EXTEND_N '''
                                           '''mirandmore_extend_hairpins_GENOME_FASTA ''')

## Process both annotation and extended annotation by means of the processing SConscript
results = {}
annotation = {'PRECURSORS_FASTA': env['PRECURSORS_FASTA'],
              'CHRM_PREFIX_TO_REMOVE': env['CHRM_PREFIX_TO_REMOVE'], 
              'GFF_FILE': env['GFF_FILE']}
extended_annotation = {'PRECURSORS_FASTA':extended_precursors[1], 
                       'CHRM_PREFIX_TO_REMOVE': env['CHRM_PREFIX_TO_REMOVE'],
                       'GFF_FILE': extended_precursors[0]}

annotation_files = {annotation_dir:annotation,
                    extended_precursors_dir:extended_annotation}

for vardir, params in annotation_files.iteritems():
    
    mirandmore_process_annotation_PRECURSORS_FASTA    = params['PRECURSORS_FASTA']
    mirandmore_process_annotation_CHRM_PREFIX_TO_REMOVE = params['CHRM_PREFIX_TO_REMOVE'] 
    mirandmore_process_annotation_GFF_FILE = params['GFF_FILE']

    process_annotation_annotation = SConscript(os.path.join(vardir, 
                                                            mirandmore_process_annotation),
                                    variant_dir = vardir, src_dir = SRC_DIR, 
                                    duplicate = 0, 
                                    exports = '''env '''
                                              '''mirandmore_process_annotation_PRECURSORS_FASTA ''' 
                                              '''mirandmore_process_annotation_CHRM_PREFIX_TO_REMOVE '''
                                              '''mirandmore_process_annotation_GFF_FILE ''')
    results[vardir] = process_annotation_annotation

Clean('.', annotation_dir)
Clean('.', extended_precursors_dir)

Return('results extended_precursors')

