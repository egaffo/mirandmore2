Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = env_trim.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('ADAPTER', '', '')
    vars.Add('QUALITY_ENCODING', '', '')
    vars.Add('SAMPLE', 'just a name to prefix to files', '')
    
    env = Environment(ENV=os.environ,
                      variables=vars)
    
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)
    # These are the variables given from the command line when the SConscript is called
    # standalone
    base = env['SAMPLE']

# -a adapter string
# -v verbose
# -l N discard sequences shorter than N
# -c discard non-clipped sequences
# -i FILE input file
# -o FILE output file
# -z compressed output
target = env['SAMPLE'] + ".tfq.gz"
report = env['SAMPLE'] + "_trimming_report.txt"

trimmed_cmd_prefix = '''cat'''
if File(env['READS']).path.endswith('.gz'):
	trimmed_cmd_prefix = '''zcat'''

trimmed_cmd = trimmed_cmd_prefix + ''' ${SOURCE} | '''\
              '''fastx_clipper -a $ADAPTER -Q $QUALITY_ENCODING '''\
              '''-v  -l 15  -c '''\
              '''-z -i - -o ${TARGETS[0]} > ${TARGETS[1]}'''
trimmed =  env.Command([target, report], 
                        env['READS'], 
                        trimmed_cmd)

Return('trimmed')
