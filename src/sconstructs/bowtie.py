import os

Import('*')

try:
    env = env_bowtie.Clone()
except NameError, ne:
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Max parallel jobs to execute', '4')
    vars.Add('BOWTIE_INDEX', 'The Bowtie index', '')
    vars.Add('READS', 'Comma separated list of input reads', '')
    vars.Add('SAMPLE', 'Name of the sample', '')
    vars.Add('BOWTIE_EXTRA_PARAMS', 'Extra parameters to be passed to Bowtie', '')

    env = Environment(ENV=os.environ,
                      variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

    env.Replace(READS = env['READS'].split(','))

env.SetDefault(BOWTIE_EXTRA_PARAMS = '')

bowtie_cmd = 'cat'
if Flatten(env['READS'])[0].endswith('.gz'):
    bowtie_cmd = 'zcat'

bowtie_cmd = bowtie_cmd + ' ${SOURCES} | bowtie -q $( -p $CPUS $) '\
                          '$BOWTIE_EXTRA_PARAMS $BOWTIE_INDEX - '\
                          '2> ${TARGETS[1]} | gzip -c > ${TARGETS[0]}'

bowtie_targets = [str(env['SAMPLE']) + '.sam.gz',
                  str(env['SAMPLE']) + '.txt']

bowtie = env.Command(bowtie_targets, 
                     env['READS'], 
                     bowtie_cmd)

Return('bowtie')
