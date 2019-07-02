'''
'''

import os, sys, re

Import('*')

try:
    env = env.Clone()
    PRECURSORS_FASTA  = mirandmore_process_annotation_PRECURSORS_FASTA
    GFF_FILE            = mirandmore_process_annotation_GFF_FILE
    CHRM_PREFIX_TO_REMOVE = mirandmore_process_annotation_CHRM_PREFIX_TO_REMOVE

except NameError:
    vars = Variables('vars.py')
    vars.Add('PRECURSORS_FASTA', 'The miRNA precursor annotation', 'hairpin.fa.gz')
    vars.Add('GFF_FILE', 'The miRBase miRNA annotation in GFF3 format', 'hsa.gff3')
    vars.Add('CHRM_PREFIX_TO_REMOVE', 'Pattern to remove from chromosome names in annotation', 
             "''")
    env = Environment(variables = vars,ENV = os.environ)
    
    Help(vars.GenerateHelpText(env))

    PRECURSORS_FASTA  = File(env['PRECURSORS_FASTA'])
    GFF_FILE            = File(env['GFF_FILE'])
    CHRM_PREFIX_TO_REMOVE = env['CHRM_PREFIX_TO_REMOVE']

mature_table_target = 'mature-table.txt'
mature_table_sources= GFF_FILE
mature_table_cmd    = '''gff2maturetable.py -g $SOURCE -m $TARGET -c ''' + CHRM_PREFIX_TO_REMOVE
mature_table = env.Command(mature_table_target, 
                           mature_table_sources, 
                           mature_table_cmd)

foldings_target = 'hairpin.species.folded'
foldings_cmd = 'RNAfold --noPS < ${SOURCE} > ${TARGET}'
foldings = env.Command(foldings_target, 
                       [PRECURSORS_FASTA], 
                       foldings_cmd)

annotation_blob_target = 'hairpin.species.annotations.blob'
annotation_blob_cmd = 'build_blob.py -i ${SOURCES[0]} -t folding -o ${TARGETS[0]}'
annotation_blob = env.Command(annotation_blob_target, 
                              [foldings], 
                              annotation_blob_cmd)

sequences_blob_target = 'hairpin.species.blob'
sequences_blob_cmd = 'build_blob.py -i ${SOURCES[0]} -t fasta -o ${TARGETS[0]}'
sequences_blob = env.Command(sequences_blob_target, 
                             [PRECURSORS_FASTA], 
                             sequences_blob_cmd)


#env.Command(["matures.fa","matures.txt"],
#            "miRNA.dat",
#            metadata_for_matures)

Return('mature_table foldings annotation_blob sequences_blob')

