'''
This Sconcript uses miRDeep2 to predict putative precursors and mature miRNAs

When called from a SConscript it imports the following variables:
 * env = env
 * mirandmore_mirdeep2_predictions_POOLED_READS             = collapsed reads pooled from all samples
 * mirandmore_mirdeep2_predictions_POOLED_ALIGNMENTS        = genomic alignaments pooled from all samples
 * mirandmore_mirdeep2_predictions_GENOME_FASTA             = genome in fasta format
 * mirandmore_mirdeep2_predictions_CPUS                     = cpus
 * mirandmore_mirdeep2_predictions_MRD_SCORE_THRESHOLD      = score threshold for miRDeep2 precursors

Returns:
 * excised_precursors:
    genomic_arf                    = genomic alignments in the mirdeep2 compatible .arf format
 	parsed_arf                     = .arf file with the alignments parsed according to the mirdeep2 default parameters
    mm_excised_precursors.fa  	   = precursors fasta sequences
	mm_excised_precursors.coords   = precursors genomic coordinates and strand
	mm_excised_precursors.log	   = log with the number of precursors excised for each minimum read number threshold
 * signatures
    mm_excised.arf                 = .arf file with unique reads aligned vs the putative precursors sorted by chr, precursor, start and end position		
 * foldings
    precursors.folded              = output of RNAfold    
 * scored_precursors
    mm_output.mrd                  = output of miRDeep2_core_algorithm.pl
 * predicted_gff
    mm_predicted.gff3		   = .gff3 file with the predicted precursors and matures
'''

import os

Import('*')

try:
    env = env.Clone()
    pooled_collapsed_reads = mirandmore_mirdeep2_predictions_POOLED_READS
    pooled_genomic_alignments = mirandmore_mirdeep2_predictions_POOLED_ALIGNMENTS
    genome_fasta = mirandmore_mirdeep2_predictions_GENOME_FASTA
    cpus = mirandmore_mirdeep2_predictions_CPUS
    mrd_score_threshold = mirandmore_mirdeep2_predictions_MRD_SCORE_THRESHOLD

except NameError:
    vars = Variables('vars.py')
    vars.Add('POOLED_READS', 'collapsed reads pooled from all samples', 'pooled_uniques.fa')
    vars.Add('POOLED_ALIGNMENTS', 'genomic alignaments pooled from all samples', 'pooled_genomic.out')
    vars.Add('GENOME_FASTA', 'genome in fasta format', 'hg_19.fa')
    vars.Add('CPUS', 'cpus', '3')
    vars.Add('MRD_SCORE_THRESHOLD', 'score threshold for miRDeep2 precursors', '50')

    env = Environment(ENV=os.environ,
                              variables=vars)

    Help(vars.GenerateHelpText(env))

    pooled_collapsed_reads = File(env['POOLED_READS'])
    pooled_genomic_alignments = File(env['POOLED_ALIGNMENTS'])
    genome_fasta =  File(env['GENOME_FASTA'])
    cpus = env['CPUS']
    mrd_score_threshold =env['MRD_SCORE_THRESHOLD'] 


#SConscripts
#SCONSCRIPT_HOME = os.path.join(env['ENV']['MIRANDMORE_HOME'],'scons','prj')
SCONSCRIPT_HOME = env['SCONSCRIPT_HOME']
#mirandmore_get_precursors = 'mirandmore_get_precursors'
#mirandmore_prepare_signatures = 'mirandmore_prepare_signatures'
#mirandmore_score_precursors = 'mirandmore_score_precursors'
#mirandmore_parse_mrd = 'mirandmore_parse_mrd'

mirdeep2_predictions =[]

## excise precursors from the genome
excised_precursors_dir = 'excised_precursors'

get_precursors_BASENAME = 'mm_excised'
get_precursors_GENOME_FASTA = genome_fasta
get_precursors_ALIGNMENTS = pooled_genomic_alignments
excised_precursors = env.SConscript(os.path.join(excised_precursors_dir, 
                                             'get_precursors.py'),
                                variant_dir = excised_precursors_dir, 
                                src_dir = SCONSCRIPT_HOME,
                                duplicate = 0,
                                exports = '''env get_precursors_BASENAME '''
                                          '''get_precursors_GENOME_FASTA '''
                                          '''get_precursors_ALIGNMENTS''')
mirdeep2_predictions.append(excised_precursors)
Clean('.', excised_precursors_dir)

## prepare the read signatures of the excised precursors
signatures_dir =  'signatures'
prepare_signatures_COLLAPSED_READS = pooled_collapsed_reads
prepare_signatures_PRECURSORS = excised_precursors[2][0]
prepare_signatures_CPUS = cpus
prepare_signatures_BASENAME = 'mm_signatures'
signatures = env.SConscript(os.path.join(signatures_dir, 
                                             'prepare_signatures.py'),
                                variant_dir = signatures_dir, 
                                src_dir = SCONSCRIPT_HOME,
                                duplicate = 0,
                                exports = '''env prepare_signatures_COLLAPSED_READS '''
                                          '''prepare_signatures_PRECURSORS '''
                                          '''prepare_signatures_CPUS '''
                                          '''prepare_signatures_BASENAME ''' )
mirdeep2_predictions.append(signatures)
Clean('.', signatures_dir)

## fold the excised precursors with RNAfold
foldings_dir = 'folded_precursors'
foldings_source = excised_precursors[2][0]
foldings_target = 'precursors.folded'
foldings_cmd = 'RNAfold --noPS < ${SOURCE} > ${TARGET}'
foldings = env.Command(os.path.join(foldings_dir,foldings_target), 
                       foldings_source, 
                       foldings_cmd)
mirdeep2_predictions.append(foldings)
Clean('.', foldings_dir)

## score the excised precursors with miRDeep2 (miRDeep2_core_algorithm.pl)
predicted_precursors_dir = 'predicted_precursors'
score_precursors_SIGNATURES = signatures
score_precursors_STRUCTURES = foldings
score_precursors_THRESHOLD = mrd_score_threshold
score_precursors_BASENAME = 'mm_output.mrd'
scored_precursors = env.SConscript(os.path.join(predicted_precursors_dir, 
                                             'score_precursors.py'),
                                variant_dir = predicted_precursors_dir, 
                                src_dir = SCONSCRIPT_HOME,
                                duplicate = 0,
                                exports = '''env score_precursors_SIGNATURES '''
                                          '''score_precursors_STRUCTURES '''
					                      '''score_precursors_THRESHOLD '''
                                          '''score_precursors_BASENAME ''')
mirdeep2_predictions.append(scored_precursors)

## parse the output.mrd and the precursors coords files to generate a .gff3 file with precursors and matures 
parse_mrd_MRD_FILE = scored_precursors
parse_mrd_PRECURSORS_COORDS = excised_precursors[2][1]
parse_mrd_BASENAME = 'mm_predicted'
predicted_gff = env.SConscript(os.path.join(predicted_precursors_dir, 
                                             'parse_mrd.py'),
                                variant_dir = predicted_precursors_dir, 
                                src_dir = SCONSCRIPT_HOME,
                                duplicate = 0,
                                exports = '''env parse_mrd_MRD_FILE '''
						                  '''parse_mrd_PRECURSORS_COORDS '''
					                      '''parse_mrd_BASENAME ''')
mirdeep2_predictions.append(predicted_gff)
Clean('.', predicted_precursors_dir)

Return('mirdeep2_predictions')
