import os

Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = env_filter.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('MAX_LEN_FILTER', 'Maximum length of the reads. Longer reads '\
                               'will be discarded or stored in a file', '30')
    vars.Add('LONG_READS_FILE', 'The file where to store reads discarded '\
                                'because of length too large. Leave empty '\
                                'if yo do not want to store discarded reads', '')
    vars.Add('MEAN_QUAL_FILTER', 'Lowest average quality a read must have', '30')
    vars.Add('QUALITY_ENCODING', 'FASTQ encoding', 'phred')
    vars.Add('SAMPLE', 'Prefix name for results', '')
    vars.Add('READS', 'The input reads file', 'reads.fastq')
    
    env = Environment(ENV=os.environ,
                      variables=vars)
    
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

target = env['SAMPLE'] + ".fq.gz"
report = env['SAMPLE'] +"_raw_filter_report.txt"

targets = [target, report]

filtered_cmd = '''filter_trimmed.py -i ${SOURCE} '''
if not env['LONG_READS_FILE'] == '':
    filtered_cmd = filtered_cmd + '''-L ${TARGETS[2]} '''
    targets.append(env['LONG_READS_FILE'])

filtered_cmd = filtered_cmd + \
               '''-l $MAX_LEN_FILTER -q $MEAN_QUAL_FILTER '''\
               '''-e $QUALITY_ENCODING -r ${TARGETS[1]} -p 20 '''\
               '''-o ${TARGETS[0]}'''

filtered = env.Command(targets, 
                       env['READS'], 
                       filtered_cmd)

Return("filtered")
