'''
This Sconcript uses fastq_to_uniques_fasta.py to collapse the reads into unique sequences in a
fasta file with miRDeep2 compatible ids
When called from a SConscript it imports the following variables:
 * env
 * collapse_reads2collapse = reads to be collapsed
 * collapse_basename = basename

Returns:
 * collapsed = unique sequences in a fasta file with miRDeep2 compatible ids
'''

import os

Import('*')

try:
    env = env_collapse.Clone()
    
except NameError:
    vars = Variables('vars.py')

    vars.Add('READS2COLLAPSE', 'parsed genomic alignments in .arf format', 
             'basename.fastq')
    vars.Add('SAMPLE', '''The name of the sample to build the target file name. '''\
                       '''Plus, the last three letters of SAMPLE will be used to '''\
                       '''identify the sample origin of collapsed readsets in FASTA headers.''', 
                       'basename')
        
    env = Environment(ENV=os.environ,
                              variables=vars)
    
    Help(vars.GenerateHelpText(env))
    
target = env['SAMPLE'] + "_uniques.fa"
collapsed = env.Command(target,
                        env['READS2COLLAPSE'],
                        """fastq_to_uniques_fasta.py -i ${SOURCE} """\
                        """-b '${SAMPLE[-3:].replace("-", "_")}' -o ${TARGET}""")
Return('collapsed')

