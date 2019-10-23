import os
Import('*')

try:
    env = env_do_tables.Clone()
except NameError as ne:
    print 'mirandmore_do_tables: command line mode.'

    vars = Variables('vars.py')
    vars.Add('MATURE_TABLE', '', '') 
    vars.Add('HAIRPIN_ANNOTATION_BLOB', '', '')
    vars.Add('NEW_RNA_OBJECT_THRESHOLD', '', '')
    vars.Add('MIN_MORNA_LEN', '', '')
    vars.Add('SISTER_OVERHANG_LEN', '', '')
    vars.Add('SISTER_MATCHES_THRESHOLD', '', '')
    vars.Add('SAMPLE', '', '')  
    vars.Add('SPECIES', '', '')
    vars.Add('EXACT_BLOB', '', '')
    vars.Add('PRE_BLOB', '', '')

    env = Environment(variables = vars,ENV = os.environ)
    Help(vars.GenerateHelpText(env))

# discretize and quantify expression of known mirnas, new new mirnas, 
# mornas and loops using only exact alignments.
# Results are in three format:
# one R friendly (_mir_table.txt)
# one Excel friendly (_mir_table_excel.txt)
# one python friendly (_pre_processed.blob)
# Unassigned  RNA objects are listed in _assign.log
# this is probably a mess now!

OS_PATH_SEP = os.path.sep

targets = [env['SAMPLE'] + "_mir_table.txt", 
           #env['SAMPLE'] + "_mir_table_excel.txt", 
           #env['SAMPLE'] + "_pre_processed.blob", 
           env['SAMPLE'] + "_pre_summary.blob",
           env['SAMPLE'] + "_assign.log"]

#sources = [do_tables_exact_blob]#,
#           do_tables_pre_blob]

#cmdLine = '''mir_discretizer.py -p ${SOURCES[1]} -e ${SOURCES[0]}  '''\
cmdLine = '''mir_discretizer.py -e ${SOURCES[0]} '''\
          '''-m $MATURE_TABLE -a $HAIRPIN_ANNOTATION_BLOB '''\
          '''--new_rna_object_threshold $NEW_RNA_OBJECT_THRESHOLD '''\
          '''--min_morna_length $MIN_MORNA_LEN '''\
          '''--species_prefix $SPECIES --sister_overhang_len $SISTER_OVERHANG_LEN '''\
          '''--sister_matches_threshold $SISTER_MATCHES_THRESHOLD '''\
          '''${TARGETS[0].dir}%(OS_PATH_SEP)s$SAMPLE''' % locals() 

tables = env.Command(targets, 
                     env['EXACT_BLOB'], 
                     cmdLine)

Return('tables')

