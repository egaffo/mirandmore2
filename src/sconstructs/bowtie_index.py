'''
Builds the Bowtie index of a fasta file with: 

bowtie-build -r -f FASTA_TO_INDEX INDEX_BASENAME
# -r    do not build the '.3.ebwt', '.4.ebwt' index files.
# -f    input is a fasta file.

When called from a SConscript it imports the following variables:
 * env
 * bowtie_index_FASTA_TO_INDEX = The fasta file to be indexed.

Returns:
 * indexes = [FASTA_BASENAME.1.ebwt,
              FASTA_BASENAME.2.ebwt,
	      FASTA_BASENAME.rev.1.ebwt,
	      FASTA_BASENAME.rev.2.ebwt,
              bwt.log, bwt.err]
'''

import os, errno

Import('*')

def SymLink(target, source, env):
    try:
        os.symlink(os.path.abspath(str(source[0])), 
                   os.path.abspath(str(target[0])))
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(os.path.abspath(str(target[0])))
            os.symlink(os.path.abspath(str(source[0])), 
                       os.path.abspath(str(target[0])))
    return None

try:
    env = env_bowtie_index.Clone()
    #fasta_to_index = bowtie_index_FASTA_TO_INDEX
except NameError as ne:
    print(ne)
    vars = Variables('vars.py')
    vars.Add('FASTA_TO_INDEX', 'The fasta file to be indexed.', 'fasta')
    
    env = Environment(variables = vars, ENV = os.environ)
 
    Help(vars.GenerateHelpText(env))
 
    #fasta_to_index = File(env['FASTA_TO_INDEX'])

current_dir = Dir('.').path
targets = ['${SOURCE.filebase}.1.ebwt', 
           '${SOURCE.filebase}.2.ebwt',
           '${SOURCE.filebase}.rev.1.ebwt', 
           '${SOURCE.filebase}.rev.2.ebwt',
           'bwt.log', 'bwt.err']

## N.B: the bowtie-build command will put index files in its execution dir,
## but if the SConscript is executed in a variant directory, we have to keep track of the 
## path by prepending the variant dir name to the index basename. Targets do not have
## to be prepended with the dir path, since it is done correctly by Scons
CmdLine = 'bowtie-build -r -f ${SOURCE} ' + \
          os.path.join(current_dir, '${SOURCE.filebase}') + ' > ${TARGETS[4]} 2> ${TARGETS[5]}'

indexes = env.Command(targets, 
                      env['FASTA_TO_INDEX'], #fasta_to_index, 
                      CmdLine)

fasta_smlnk = env.Command('${SOURCE.file}', 
                          env['FASTA_TO_INDEX'], #fasta_to_index, 
                          SymLink)

Return('indexes')

