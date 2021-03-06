from deepscons import Uri, smart_decider
import os, errno, sys

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

#vars = Variables('vars.py')
vars = Variables()
vars.Add('CPUS',
         '',
         1)

env = Environment(ENV = os.environ,
                  variables = vars)

Help(vars.GenerateHelpText(env))

env = Environment(ENV=os.environ, SHELL = '/bin/bash')

if ARGUMENTS.get('CPUS', 0):
    env['CPUS'] = ARGUMENTS.get('CPUS', 0)
    

Decider('MD5-timestamp')

#env.SetDefault(CPUS = 1)
PYTHON_VERSION = str(sys.version_info.major) + '.' + str(sys.version_info.minor)
TOOLS = os.path.join(env['ENV']['MIRANDMORE_HOME'], 'tools')
BIN = os.path.join(env['ENV']['MIRANDMORE_HOME'],'bin')
PYTHON_LIB = os.path.join(env['ENV']['MIRANDMORE_HOME'], 
                          'lib', 'python'+PYTHON_VERSION, 'site-packages')

#env.PrependENVPath('PYTHONUSERBASE', TOOLS)
env.PrependENVPath('PYTHONUSERBASE', env['ENV']['MIRANDMORE_HOME'])
env.PrependENVPath('PYTHONPATH', PYTHON_LIB)

## PIP
pip_file = 'get-pip.py'
pip_url = 'https://bootstrap.pypa.io/' + pip_file
pip_targets = [os.path.join(TOOLS, pip_file),
               os.path.join(BIN, 'pip')]
pip_cmd = ' && '.join(['wget -O ${TARGETS[0]} ' + pip_url, 
                       'python ${TARGETS[0]} --user'])
pip = env.Command(pip_targets, 
                  [], 
                  pip_cmd)

## use pip just installed and not the system one
env.PrependENVPath('PATH', os.path.join(TOOLS, 'bin'))
# BIOPYTHON
BIOPYTHON_dir = os.path.join(PYTHON_LIB, 'Bio')
BIOPYTHON_target = [os.path.join(BIOPYTHON_dir, 'SeqIO', 'FastaIO.py')]
BIOPYTHON = env.Command(BIOPYTHON_target, [pip], 
                        ['pip install --ignore-installed --user biopython'])

# HTSeq
HTSeq_dir = os.path.join(PYTHON_LIB, 'HTSeq')
HTSeq_target = [os.path.join(HTSeq_dir, '__init__.py'),
                os.path.join(BIN, 'htseq-count')]
HTSeq = env.Command(HTSeq_target, [pip], 
                    ['pip install --ignore-installed --user HTSeq'])

#env.Command(os.path.join(BIN, "${SOURCE.file}"), HTSeq[1], SymLink)

# STREAM
#stream_dir = os.path.join(PYTHON_LIB)
#stream_target = [os.path.join(stream_dir, 'stream.py')]
#stream = env.Command(stream_target, [pip], 
#                    ['pip install --ignore-installed --user stream'])


## RNAfold
RNAfold_dnwld = env.Command(os.path.join(TOOLS, "ViennaRNA-2.4.14.tar.gz"), 
            Uri("https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz"), 
            "wget -q ${SOURCE} -O ${TARGET}")
RNAfold_xtr = env.Command(os.path.join(TOOLS, "ViennaRNA-2.4.14/configure"), 
                          RNAfold_dnwld, 
                          "tar xfz ${SOURCE} -C `dirname ${SOURCE}`")
RNAfold_cfg = env.Command(os.path.join(TOOLS, "ViennaRNA-2.4.14/Makefile"), 
                          RNAfold_xtr, 
                          "cd `dirname ${SOURCE}` && ./configure")
RNAfold_make = env.Command(os.path.join(TOOLS, "ViennaRNA-2.4.14/src/bin//RNAfold"), 
                           RNAfold_cfg,
                           "cd `dirname ${SOURCE}` && make $( -j $CPUS $)")
RNAfold_link = env.Command(os.path.join(BIN, "RNAfold"),
                           RNAfold_make, 
                           SymLink)

## BOWTIE
BOWTIE_dnwld = env.Command(os.path.join(TOOLS, "bowtie-1.1.2-linux-x86_64.zip"),
            Value("http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip"),
            "wget -q ${SOURCE} -O ${TARGET}")
BOWTIE_xtr = env.Command([os.path.join(TOOLS, "bowtie-1.1.2/bowtie-build"), 
                          os.path.join(TOOLS, "bowtie-1.1.2/bowtie")],
                         BOWTIE_dnwld,
                         "unzip ${SOURCE} -d `dirname ${SOURCE}`")
BOWTIE_link = env.Command(os.path.join(BIN, "bowtie"), 
                          BOWTIE_xtr[1], 
                          SymLink)
BOWTIE_BUILD_link = env.Command(os.path.join(BIN,"bowtie-build"), 
                                BOWTIE_xtr[0], 
                                SymLink)

## FASTQC
FASTQC_zip = 'fastqc_v0.11.8.zip'
FASTQC_link = 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/' + FASTQC_zip
FASTQC_target = [os.path.join(TOOLS, FASTQC_zip),
                 os.path.join(TOOLS, 'FastQC', 'fastqc')] 
FASTQC = env.Command(FASTQC_target, [], ['wget -q -O $TARGET  ' + FASTQC_link,
                                         'unzip -u -d ${TARGETS[0].dir} ${TARGETS[0]} '\
                                         '&& chmod +x ${TARGETS[1]}'])
env.Command(os.path.join(BIN, "${SOURCE.file}"), FASTQC[1], SymLink)

## MIRDEEP2
md2_zip = "v0.1.2.zip" #'v0.0.8.zip'
md2_url = 'https://github.com/rajewsky-lab/mirdeep2/archive/' + md2_zip
md2_dir = 'mirdeep2-0.1.2' #'mirdeep2-0.0.8'
md2_dnwld = env.Command(os.path.join(TOOLS, md2_zip), 
                        Uri(md2_url), 
                        "wget -q ${SOURCE} -O ${TARGET}")
md2_xtr = env.Command([os.path.join(TOOLS, md2_dir, "src", "convert_bowtie_output.pl"),
                       os.path.join(TOOLS, md2_dir, "src", "excise_precursors_iterative_final.pl"),
                       os.path.join(TOOLS, md2_dir, "src", "parse_mappings.pl"),
                       os.path.join(TOOLS, md2_dir, "src", "miRDeep2_core_algorithm.pl")], 
                      md2_dnwld, 
                      "unzip ${SOURCE} -d `dirname ${SOURCE}`")

for l in md2_xtr:
    env.Command(os.path.join(BIN, l.name), l, SymLink)

### R PACKAGES
# getopt, plyr #data.table
env['ENV']['R_LIBS'] = os.path.join(TOOLS, "R_libs")
env['ENV']['R_INSTALL_STAGED'] = 'false' ## this is a workaround to prevent
                                    ## installation failure due to Scons 
                                    ## generating empty files and directory
                                    ## structure of targets
R_libs_targets = [os.path.join(TOOLS, 'R_libs', 'plyr', 'R', 'plyr'),
                  os.path.join(TOOLS, 'R_libs', 'getopt', 'R', 'optparse')]
R_libs = env.Command(R_libs_targets, [], 'install_R_libs.R')


# BEDTOOLS
BEDTOOLS_tar = 'bedtools-2.27.0.tar.gz'
BEDTOOLS_dir = os.path.join(TOOLS, 'bedtools2')
BEDTOOLS_link = 'https://github.com/arq5x/bedtools2/releases/download/v2.27.0/' + BEDTOOLS_tar
BEDTOOLS_target = [os.path.join(TOOLS, BEDTOOLS_tar), 
                   os.path.join(BEDTOOLS_dir, 'bin', 'bedtools')]
BEDTOOLS = env.Command(BEDTOOLS_target, [], 
                       ['wget -O $TARGET ' + BEDTOOLS_link, 
                        'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}', 
                        'cd ' + BEDTOOLS_dir + ' && make', 
                        'cd ' + Dir('.').abspath])
env.Command(os.path.join(BIN, "${SOURCE.file}"), BEDTOOLS[1], SymLink)

# SAMTOOLS
SAMTOOLS_tar = 'samtools-1.3.1.tar.bz2'
SAMTOOLS_dir = os.path.join(TOOLS, 'samtools-1.3.1')
SAMTOOLS_link = 'https://github.com/samtools/samtools/releases/download/1.3.1/' + SAMTOOLS_tar
SAMTOOLS_target = [os.path.join(TOOLS, 'samtools-1.3.1.tar.bz2'),
                   os.path.join(SAMTOOLS_dir, 'samtools')]
SAMTOOLS = env.Command(SAMTOOLS_target, [], ['wget -O $TARGET  ' + SAMTOOLS_link,
                                             'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                                             'cd ' + SAMTOOLS_dir + ' && '\
                                             'make prefix=' + SAMTOOLS_dir + ' install',
                                             'cd ' + Dir('.').abspath])
env.Command(os.path.join(BIN, "${SOURCE.file}"), SAMTOOLS[1], SymLink)


# CUTADAPT
cutadapt_dir = os.path.join(PYTHON_LIB, 'cutadapt')
cutadapt_target = [os.path.join(BIN, 'cutadapt')]
cutadapt = env.Command(cutadapt_target, [pip],
		       ['pip install --ignore-installed --user cutadapt'])
