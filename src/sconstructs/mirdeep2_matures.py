'''
This Sconcript generates a .csv file with information about miRDeep2 predicted matures
1) bedtools intersect -s -wao 	to intersect all the miRDeep2 matures with the known precursors
2) mirdeep2_matures.py 		    to generate a table of miRdeep2 matures and their eventual
				                overlaps with known precursors

When called from a SConscript it imports the following variables:
 * env
 * mirandmore_mirdeep2_matures_PREDICTED_GFF    = predicted precursors in .gff3 format
 * mirandmore_mirdeep2_matures_KNOWN_GFF        = mirbase .gff annotation file
 * mirandmore_mirdeep2_matures_BASENAME         = basename used to name files

Returns:
 * basename_intersected.gff3    = output of bedtools intersect
 * basename.txt                 = table with format:
    
    name    ID  chrom   strand  start   end miRDeep2-score overlaps
    (start and end are 1 based (.gff3 format)) 
'''

import os

Import('*')

try:
    env = env.Clone()
    predicted_gff = mirandmore_mirdeep2_matures_PREDICTED_GFF
    known_gff = mirandmore_mirdeep2_matures_KNOWN_GFF
    basename = mirandmore_mirdeep2_matures_BASENAME

except NameError:
    vars = Variables('vars.py')

    vars.Add('PREDICTED_GFF', 'predicted precursors and matures in gff3 format', 
             'mm_predicted.gff3')

    vars.Add('KNOWN_GFF', 'mirbase .gff3 annotation file', 
             'hsa.gff3')

    vars.Add('BASENAME', 'basename', 'mrd_matures')

    env = Environment(ENV=os.environ,
                              variables=vars)

    Help(vars.GenerateHelpText(env))

    basename = env['BASENAME']
    known_gff = File(env['KNOWN_GFF'])
    predicted_gff = File(env['PREDICTED_GFF'])

env['SHELL'] = '/bin/bash'

results = []

sources = [predicted_gff, known_gff]
target = basename + "_intersected.gff3"
cmd_line = '''grep -P "^[^#].*miRNA[\s]" ${SOURCES[0]} | bedtools intersect -s '''\
			'''-wao -a - -b <(grep "^[^#].*miRNA_primary_transcript" ${SOURCES[1]}) > ${TARGET}'''
intersected_matures_gff = env.Command(target, sources, cmd_line)
results.append(intersected_matures_gff)

sources = intersected_matures_gff
target = basename + ".txt"
mrd_matures = env.Command(target, sources, "mirdeep2_matures.py -i ${SOURCE} -o ${TARGET}")
results.append(mrd_matures)

Return('results')

