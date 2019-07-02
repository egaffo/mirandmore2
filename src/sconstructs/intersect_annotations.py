'''
This Sconcript generates a .gff3 annotation file with the known and novel precursors:
1) bedtools intersect -s -v is used to select only the predicted precursors and matures which do not overlap any
   known precursor or mature on the same strand
2) generates a .gff3 annotation file with the known and novel precursors ( exclusion of
   predicted matures whose precursors overlapped at least one known precursor even if they didn't)
3) generates a .gff3 with all the excluded matures not overlapping any known precursor
4) generates a .gff3 file with all the excluded predicted precursors and matures
5) generates a .txt file with the numbers of novel or excluded precursors and matures

When called from a SConscript it imports the following variables:
 * env
 * mirandmore_intersect_annotations_PREDICTED_GFF = predicted precursors in .gff3 format
 * mirandmore_intersect_annotations_KNOWN_GFF = mirbase .gff annotation file
 * mirandmore_intersect_annotations_BASENAME = basename used to name files

Returns:
 * basename_intersected.gff3	 	 =  .gff3 file output of bedtools intersect -v -s
 * basename_annotations.gff3	 	 =  .gff3 file with known and putative precursors and their associated mature miRNAs
 * basename_mature_leftovers.gff3	 =  .gff3 with all the excluded matures not overlapping any known precursor = "orphan matures"
					                    (they were excluded because their precursors overlapped at least one known precursor
					                     even if they didn't)
 * basename_all_leftovers.gff3 	     =  .gff3 file with all the excluded predicted precursors and matures, including the "orphan matures"
 * basename_stats.txt                =  .txt file with statistics with format:  name\tnumber\n
                                         The file contains:
                                         	novel_precursors			number
                                        	novel_matures				number	
                                        	orphan_matures				number
                                        	discarded_overlapping_precursors	number
                                        	discarded_overlapping_matures 		number   
'''

import os

Import('*')

try:
    env = env.Clone()
    predicted_gff = mirandmore_intersect_annotations_PREDICTED_GFF
    known_gff = mirandmore_intersect_annotations_KNOWN_GFF
    basename = mirandmore_intersect_annotations_BASENAME

except NameError:
    vars = Variables('vars.py')

    vars.Add('PREDICTED_GFF', 'predicted precursors in gff3 format', 
             'mm_predicted.gff3')

    vars.Add('KNOWN_GFF', 'mirbase .gff3 annotation file', 
             'hsa.gff3')

    vars.Add('BASENAME', 'basename', 'mm_annotations')

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
dirty_gff = env.Command(target, sources, 'bedtools intersect -s -v -a ${SOURCES[0]} -b ${SOURCES[1]} > ${TARGET}')
results.append(dirty_gff)


# selects all putative precursors lines form the intersected file and passes them to sed to create REs. These REs are then used to
# extract from the intersected files the precursors and the matures associated with a precursor contained in the intersected file.
# ($|;) is used to match both the precursors lines (PRE.*;) and matures lines (PRE.*$). grep needs option -E to use REs containing |.
# The selected precursors and matures are concatenated to the gff file with the known annotations without the header.
sources = [dirty_gff, known_gff]
target = basename + ".gff3"
cmdLine = '''grep 'miRNA_primary_transcript' ${SOURCES[0]} | sed 's@.*ID=\([^;]*\).*@\\1($|;)@'| '''\
		  '''grep -E -f - ${SOURCES[0]}| cat <(grep -v '^#' ${SOURCES[1]}) - | '''\
          '''sort -k1,1 -k4,4n -k5,5rn > $TARGET'''
annotations_gff = env.Command(target, sources, cmdLine)
results.append(annotations_gff)

# selects all the excluded matures not overlapping any known precursor
# = "orphan matures" (they were excluded because their precursors overlapped
# at least one known precursor even if they didn't)
sources = [annotations_gff, dirty_gff]
target = basename + "_mature_leftovers.gff3"
# NB: grep return error code when no match is found. Wrap the command in curly braces and
# OR with true to bypass the error 
mature_leftovers_cmd = '{ grep -v -f ${SOURCES[0]} ${SOURCES[1]} > ${TARGET} || true; }'
mature_leftovers_gff = env.Command(target, sources, mature_leftovers_cmd)
results.append(mature_leftovers_gff)

# selects all the excluded precursors matures, the .gff3 file includes: 
# - orphan matures
# - precursors excluded because they overlapped known annotations
# - matures excluded because they overlapped known annotations

sources = [annotations_gff, predicted_gff]
target = basename + "_all_leftovers.gff3"
all_leftovers_gff = env.Command(target, sources, 
				'{ grep -v -f ${SOURCES[0]} ${SOURCES[1]} > ${TARGET} || true; }' )
results.append(all_leftovers_gff)

# creates a .txt file with statistics with format:  name\tnumber\n 
# The file contains:
# 	novel_precursors			number
#	novel_matures				number	
#	orphan_matures				number
#	discarded_overlapping_precursors	number
#	discarded_overlapping_matures 		number


sources = [all_leftovers_gff, predicted_gff, mature_leftovers_gff]
target = basename + "_stats.txt"
cmdLine = '''echo -n "novel_precursors\t" > ${TARGET} ; '''\
          '''grep -v -f ${SOURCES[0]} ${SOURCES[1]}| grep miRNA_primary_transcript | wc -l >> ${TARGET} ; '''\
          '''echo -n "novel_matures\t" >> ${TARGET} ; '''\
          '''grep -v -f ${SOURCES[0]} ${SOURCES[1]}| grep 'miRNA[^_]' | wc -l >> ${TARGET} ; '''\
          '''echo -n "orphan_matures\t" >> ${TARGET} ; cat ${SOURCES[2]} | wc -l >> ${TARGET} ; '''\
          '''echo -n "discarded_overlapping_precursors\t" >> ${TARGET} ; '''\
          '''grep 'miRNA_primary_transcript' ${SOURCES[0]} | wc -l >> ${TARGET} ; '''\
          '''echo -n "discarded_overlapping_matures\t" >> ${TARGET} ; '''\
          '''grep 'miRNA[^_]' ${SOURCES[0]} | grep -v -f ${SOURCES[2]} - | wc -l >> ${TARGET} '''
            
stats = env.Command(target, sources, cmdLine)
results.append(stats)

Return('results')

