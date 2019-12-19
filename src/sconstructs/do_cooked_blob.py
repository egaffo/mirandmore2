Import('*')

try:
    env = env.Clone()
    MATURE_TABLE = cooked_blob_mature_table
    base = cooked_blob_sample
    cooked_blob_exact_blob
    cooked_blob_genomic_hits_blob
except NameError as ne:
    print('mirandmore_do_cooked_blob: command line mode.')

    vars = Variables('vars.py')
    vars.Add('SAMPLE', '', '')
    vars.Add('MATURE_TABLE', '', '')
    vars.Add('EXACT_BLOB', '', '')
    vars.Add('GENOMIC_HITS_BLOB', '', '')

    env = Environment(variables = vars,ENV = os.environ)
    Help(vars.GenerateHelpText(env))

    MATURE_TABLE = env['MATURE_TABLE']
    base = env['SAMPLE']
    cooked_blob_exact_blob = env['EXACT_BLOB']
    cooked_blob_genomic_hits_blob = env['GENOMIC_HITS_BLOB']

# dumps a serialized dict with hairpin names as keys 
# and a list with number of read hits per nucleotide as values

sources = [cooked_blob_exact_blob, 
	   cooked_blob_genomic_hits_blob]

target = base + "_cooked.blob" 

cmdLine = "cook_exact_blob.py --input ${SOURCE} --output ${TARGET} --mature-table %(MATURE_TABLE)s" % locals() 

exact_blob = env.Command(target, sources, cmdLine)

Return('exact_blob')
