'''

This is a SConscript script that execute the FASTQC [1] program to analyze FASTA, 
FASTQ or BAM files.

Software dependencies:
 * FASTQC

When called from a SConscript it imports the following variables:
 * env
 * fastqc_readset
  
[1] http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
'''


import os, re

Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = env_read_qual_stat.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('READS', 'FASTQ read file to process, either in plain text or gzipped (.gz)', 
             'reads.fq')
    
    env = Environment(ENV=os.environ,
                      variables=vars)
    
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print("Unknown variables:", list(unknown.keys()))
        Exit(1)

# FASTQC removes .gz, .fq, .fastq extensions (also recursively) from the input file name
# to name the output HTML and ZIP files
readset_basename = re.sub('\.fq$|\.fastq$', '', re.sub('\.gz$', '', File(env['READS']).name))

fastqc_quality_cmd = 'fastqc $SOURCE -o ${TARGETS[0].dir} > ${TARGETS[2]} 2> ${TARGETS[3]}' 
fastqc_html_target = readset_basename + '_fastqc.html'
fastqc_zip_target = readset_basename + '_fastqc.zip'

fastqc_quality = env.Command([fastqc_html_target, fastqc_zip_target,
                             '${SOURCE.filebase}_fastqc.log', '${SOURCE.filebase}_fastqc.err'], 
                             env['READS'], 
                             fastqc_quality_cmd)

Return('fastqc_quality')

