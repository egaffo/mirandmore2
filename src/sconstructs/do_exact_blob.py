Import('*')

try:
    env = env_exact_blob.Clone()
except NameError as ne:
    print 'mirandmore_do_exact_blob: command line mode.'

    vars = Variables('vars.py')
    vars.Add('SAMPLE', '', '')
    vars.Add('MULTIPLE_GENOMIC_HITS_THRESHOLD', '', '')
    vars.Add('MIN_COUNT', '', '')
    vars.Add('MATURE_TABLE', '', '')
    vars.Add('HAIRPIN_ALIGNMENTS', '', '')
    vars.Add('GENOMIC_HITS_BLOB', '', '')

    env = Environment(variables = vars,ENV = os.environ)
    Help(vars.GenerateHelpText(env))

## dumps a serialized data structure with reads information for every read 
## aligning exactly with hairpin
## _exact.blob (rna.PreExactResultSet)

target = env['SAMPLE'] + "_exact.blob"

cmdLine = '''dump_exact_blob.py --mature-table $MATURE_TABLE '''\
          '''--threshold $MULTIPLE_GENOMIC_HITS_THRESHOLD '''\
          '''--genomic-hits ${SOURCES[1]} --output $TARGET --input ${SOURCES[0]} '''\
          '''-c $MIN_COUNT'''

exact_blob = env.Command(target, 
                         [env['HAIRPIN_ALIGNMENTS'], env['GENOMIC_HITS_BLOB']], 
                         cmdLine)

Return('exact_blob')
