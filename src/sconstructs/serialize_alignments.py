'''
Input:
    serialize_hairpin_alignments
    serialize_hairpin_genomic_hits_blob
Returns:
    (_pre.blob, 
    _stats.txt, 
    _matures.blob, 
    _multiple_hits_summary.txt)
'''

Import('*')

try:
    env = env_serialize_alignments.Clone()
except NameError as ne:
    print 'mirandmore_serialize_alignments: command line mode.'
    vars = Variables('vars.py')
    vars.Add('SAMPLE', '', '')
    vars.Add('MULTIPLE_GENOMIC_HITS_THRESHOLD', '', '')
    vars.Add('MIN_COUNT', '', '')
    vars.Add('ALLOWED_OVERHANG', '', '')
    vars.Add('ALLOWED_OVERHANG_MORNA', '', '')
    vars.Add('EXACT', '', '')
    vars.Add('SHORTER_OR_LONGER', '', '')
    vars.Add('MIS_1', '', '')
    vars.Add('MIS_2', '', '')
    vars.Add('FIVE_PRIME', '', '')
    vars.Add('THREE_PRIME', '', '')
    vars.Add('MIN_MORNA_LEN', '', '')
    vars.Add('MATURE_TABLE', '', '')
    vars.Add('SAM', '', '')
    vars.Add('SERIALIZE_HAIRPIN_GENOMIC_HITS_BLOB', '', '')

    env = Environment(variables = vars,ENV = os.environ)
    Help(vars.GenerateHelpText(env))

# dumps two serialized data structures with reads information by hairpin and mature. 
# Only reads aligning to matures in one of the standard categories.
#  _pre.blob (rna.PreResultSet)
#  _matures.blob (rna.MatureResultSet)

pre_blob                = env['SAMPLE'] + "_pre.blob"
stats                   = env['SAMPLE'] + "_stats.txt"
matures_blob            = env['SAMPLE'] + "_matures.blob"
multiple_hits_summary   = env['SAMPLE'] + "_multiple_hits_summary.txt"

targets = [pre_blob, 
           stats, 
           matures_blob, 
           multiple_hits_summary]

cmdLine = '''process_algn.py --mature-table $MATURE_TABLE --threshold '''\
          '''$MULTIPLE_GENOMIC_HITS_THRESHOLD --genomic-hits ${SOURCES[1]} '''\
          '''--input ${SOURCES[0]} -a $ALLOWED_OVERHANG -o $ALLOWED_OVERHANG_MORNA '''\
          '''-e $EXACT -l $SHORTER_OR_LONGER -s $MIS_1 -S $MIS_2 '''\
          '''-f $FIVE_PRIME -r $THREE_PRIME -c $MIN_COUNT '''\
          '''-n $MIN_MORNA_LEN -b ${str(TARGETS[0].abspath).split("_pre.blob")[0]}'''

blobs =  env.Command(targets, 
                     [env['SAM'], env['GENOMIC_HITS_BLOB']], 
                     cmdLine)

Return('blobs')

