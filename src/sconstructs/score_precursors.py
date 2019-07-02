'''
This Sconcript uses miRDeep2_core_algorithm.pl to score the excised precursors using the read signatures and structures
When called from a SConscript it imports the following variables:
 * env
 * score_precursors_SIGNATURES = .arf file of read signatures (sorted by chromosome, precursor, start, end)
 * score_precursors_STRUCTURES = RNAfold structures of the excised precursors
 * score_precursors_THRESHOLD = minimum miRDeep2 score for a precursor to be included in the output
 * score_precursors_BASENAME = basename

Returns:
 * scored_precursors = miRDeep2_core_algorithm.pl output with precursors with scores > threshold
'''

import os

Import('*')

try:
    env = env.Clone()
    signatures = score_precursors_SIGNATURES
    structures = score_precursors_STRUCTURES
    threshold = score_precursors_THRESHOLD
    basename = score_precursors_BASENAME
    
except NameError:
    vars = Variables('vars.py')
    vars.Add('SCORE_PRECURSORS_SIGNATURES', '.arf file of read signatures', 
             'mm_signatures.arf')
    vars.Add('SCORE_PRECURSORS_STRUCTURES', 'RNAfold structures of the excised precursors', 
             'precursors.folded')
    vars.Add('SCORE_PRECURSORS_BASENAME', 'basename', 'mm_output.mrd')
    vars.Add('SCORE_PRECURSORS_THRESHOLD', 'miRDeep2 score threshold', '50')
        
    env = Environment(ENV=os.environ,
                              variables=vars)
    
    Help(vars.GenerateHelpText(env))
    
    signatures = File(env['SCORE_PRECURSORS_SIGNATURES'])
    structures = File(env['SCORE_PRECURSORS_STRUCTURES'])
    threshold = env['SCORE_PRECURSORS_THRESHOLD']
    basename = env['SCORE_PRECURSORS_BASENAME']

sources = [signatures,structures]
target = basename
CmdLine = 'miRDeep2_core_algorithm.pl ${SOURCES[0]} ${SOURCES[1]} -v ' + threshold + ' > ${TARGET}'
scored_precursors = env.Command(target, sources, CmdLine)
Return('scored_precursors')

