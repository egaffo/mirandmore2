Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = env_cutadapt.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('ADAPTER', '', '')
    vars.Add('QUALITY_ENCODING', '33 for Phread score, 64 for old Solexa', '')
    vars.Add('SAMPLE', 'just a name to prefix to files', '')
    vars.Add('MIN_LEN', 
             '', 
             '15')
    vars.Add('MAX_LEN', 
             'Max read length. Reads longer than MAX_LEN after adapter trim '\
             'are saved in a separate file', 
             '34')
    vars.Add('ADPT_OVERLAP', 'Min lenght of adapter overlap', '')
    vars.Add('CPUS', 
             '', 
             '1')
    vars.Add('CUTADAPT_EXTRA_PARAMS', '', '')
    
    env = Environment(ENV=os.environ,
                      variables=vars)
    
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print("Unknown variables:", list(unknown.keys()))
        Exit(1)

env.SetDefault(ADPT_OVERLAP = 5)
env.SetDefault(MIN_LEN = 15)

cutadapt_cmd = '''cutadapt $( -j $CPUS $) -a $ADAPTER '''\
               '''--no-indels -n 10 -O $ADPT_OVERLAP '''\
               '''--quality-base $QUALITY_ENCODING -m $MIN_LEN -M $MAX_LEN_FILTER '''\
               '''--trimmed-only '''\
               '''--too-long-output ${TARGETS[2]} '''\
               '''$CUTADAPT_EXTRA_PARAMS -o ${TARGETS[0]} '''\
               '''${SOURCES[0]} > ${TARGETS[1]}'''

cutadapt =  env.Command(['${SAMPLE}_trimmed_reads.fq.gz', 
                         '${SAMPLE}_cutadapt.log', 
                         '${SAMPLE}_long_reads.fq.gz'],
                        env['READS'], 
                        cutadapt_cmd)

Return('cutadapt')
