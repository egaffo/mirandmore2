import os

Import('*')

try:
    env = env_predict_novel_precursors.Clone()

except NameError:
    vars = Variables('vars.py')

    env = Environment(ENV=os.environ,
                              variables=vars)

    Help(vars.GenerateHelpText(env))

results = []

##
## Preecursor predictions
##

mirandmore_mirdeep2_dir = 'mirdeep2_precursors'


mirdeep2_predictions = []

#SConscripts
mirandmore_mirdeep2_predictions ='mirdeep2_predictions.py'

results = env['PREPROCESSING_FILES']

# pool the collapsed reads from all samples
pooling_source = [results[f]['collapsed'][0] for f in results.keys()] # all the unique read files
pooling_target = 'pooled_uniques.fa'
catenate_collapsed_cmd = 'cat ${SOURCES} > ${TARGET}' 
pooled_collapsed = env.Command(os.path.join(mirandmore_mirdeep2_dir, pooling_target), 
                               pooling_source, 
                               catenate_collapsed_cmd)
mirdeep2_predictions.append(pooled_collapsed)

# pool the genomic alignments from all samples
pooling_source = [results[f]['aligned_uniques'][0] for f in results.keys()] # all the genomic alignments files
pooling_target = 'pooled_genomic.out'
catenate_genomic_out_cmd = 'cat ${SOURCES} > ${TARGET}' 
pooled_aligned_uniques = env.Command(os.path.join(mirandmore_mirdeep2_dir, pooling_target), 
                                      pooling_source, catenate_genomic_out_cmd)
mirdeep2_predictions.append(pooled_aligned_uniques)

# miRDeep2 predictions
mirandmore_mirdeep2_predictions_POOLED_READS             = pooled_collapsed
mirandmore_mirdeep2_predictions_POOLED_ALIGNMENTS        = pooled_aligned_uniques
mirandmore_mirdeep2_predictions_GENOME_FASTA             = File(env['GENOME'])
mirandmore_mirdeep2_predictions_CPUS                     = env['CPUS']
mirandmore_mirdeep2_predictions_MRD_SCORE_THRESHOLD      = env['MRD_SCORE_THRESHOLD']
predictions = env.SConscript(os.path.join(mirandmore_mirdeep2_dir, mirandmore_mirdeep2_predictions),
                                variant_dir = mirandmore_mirdeep2_dir, 
                                src_dir = env['SCONSCRIPT_HOME'],
                                duplicate = 0,
                                exports = '''env mirandmore_mirdeep2_predictions_POOLED_READS '''
                                          '''mirandmore_mirdeep2_predictions_POOLED_ALIGNMENTS '''
                                          '''mirandmore_mirdeep2_predictions_GENOME_FASTA '''
                                          '''mirandmore_mirdeep2_predictions_CPUS '''
                                          '''mirandmore_mirdeep2_predictions_MRD_SCORE_THRESHOLD ''')

mirdeep2_predictions.append(predictions)

##
## Annotations
##

#SConscripts
mirandmore_mirdeep2_matures = 'mirdeep2_matures.py'
mirandmore_intersect_annotations = 'intersect_annotations.py'

## intersect the predicted matures with the known precursors
mirdeep2_matures_dir = 'mirdeep2_matures'
mirandmore_mirdeep2_matures_PREDICTED_GFF = predictions[4]
mirandmore_mirdeep2_matures_KNOWN_GFF = File(env['GFF_FILE'])
mirandmore_mirdeep2_matures_BASENAME = 'mirdeep2_matures'
mirdeep2_matures = env.SConscript(os.path.join(mirdeep2_matures_dir, 
                                             mirandmore_mirdeep2_matures),
                                variant_dir = mirdeep2_matures_dir, 
                                src_dir = env['SCONSCRIPT_HOME'],
                                duplicate = 0,
                                exports = '''env mirandmore_mirdeep2_matures_PREDICTED_GFF '''
                                          '''mirandmore_mirdeep2_matures_KNOWN_GFF '''
                                          '''mirandmore_mirdeep2_matures_BASENAME ''')

## intersect the predicted precursors with the known precursors
new_known_intersect_dir = 'new_known_intersect'
mirandmore_intersect_annotations_PREDICTED_GFF = predictions[4]
mirandmore_intersect_annotations_KNOWN_GFF = File(env['GFF_FILE'])
mirandmore_intersect_annotations_BASENAME = 'mm_annotations'
merged_annotations = env.SConscript(os.path.join(new_known_intersect_dir, 
                                             mirandmore_intersect_annotations),
                                variant_dir = new_known_intersect_dir, 
                                src_dir = env['SCONSCRIPT_HOME'],
                                duplicate = 0,
                                exports = '''env mirandmore_intersect_annotations_PREDICTED_GFF '''
                                          '''mirandmore_intersect_annotations_KNOWN_GFF '''
                                          '''mirandmore_intersect_annotations_BASENAME ''')

## Retrieve (novel and new) hairpin sequences
merged_precursors_fasta_cmd = '''grep miRNA_primary_transcript ${SOURCES[1]}| '''\
                              '''sed "s/miRNA_primary_transcript\(.\+\)'''\
                              '''Name=\(.\+\)/\\2\\1Name=\\2/" '''\
                              '''| bedtools getfasta -s -name -fi '''\
                              '''${SOURCES[0]} -bed - |'''\
                              '''sed -r "s@::.+:[0-9]+-[0-9]+\([-+]\)@@" - >${TARGET}'''

merged_precursors_fasta_sources = [File(env['GENOME']), merged_annotations[1]]
merged_precursors_fasta_target = 'hairpins.fa'
merged_precursors_fasta = env.Command(merged_precursors_fasta_target, 
                                      merged_precursors_fasta_sources, 
                                      merged_precursors_fasta_cmd)

## Retrieve (novel and new) mature miRNA sequences
merged_matures_fasta_cmd = '''grep "miRNA\s" ${SOURCES[1]}| '''\
                              '''sed "s/miRNA\(.\+\)'''\
                              '''Name=\([^;]\+\)/\\2\\1Name=\\2/" '''\
                              '''| bedtools getfasta -s -name -fi '''\
                              '''${SOURCES[0]} -bed - |'''\
                              '''sed -r "s@::.+:[0-9]+-[0-9]+\([-+]\)@@" - >${TARGET}'''

merged_matures_fasta = env.Command('matures.fa', 
                                      [File(env['GENOME']), merged_annotations[1]], 
                                      merged_matures_fasta_cmd)

results = {'PRECURSORS_FASTA': merged_precursors_fasta, 
          'GFF_FILE': merged_annotations[1],
          'NEW_PREMIRS_MIRS_GFF': merged_annotations[1][0]
          }

Clean('.', mirandmore_mirdeep2_dir)

Return('results')
