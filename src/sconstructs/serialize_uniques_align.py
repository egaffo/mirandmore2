Import('*')

try:
    env = env_serialize_uniques_align.Clone()
except NameError as ne:
    print('mirandmore_align_uniques: failed to import', ne)
    vars = Variables('vars.py')
    vars.Add('SAMPLE', '''The name of the sample to build the target file name. ''', '')
    vars.Add('GENOMIC_HITS', '', '')
        
    env = Environment(ENV=os.environ,
                              variables=vars)
    
    Help(vars.GenerateHelpText(env))

# serialize with pickle a dict of hits to genomic sequences with read sequence as key.
target = env['SAMPLE'] + "_genomic_hits.blob"

cmd = "genomic_hits_to_blob.py -i ${SOURCE} -o ${TARGET}"
genomic_hits_blob = env.Command(target, env['GENOMIC_HITS'], cmd)

Return('genomic_hits_blob')

