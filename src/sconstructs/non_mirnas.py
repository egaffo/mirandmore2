'''
This script maps reads to the whole reference 
genome, eventually intersecting the alignments 
with known annotation to estimate expressed genes
'''

import os

Import('*')

try:
    env = env_non_mirna.Clone()
except NameError, ne:
    vars = Variables('vars.py')
    vars.Add('QUALITY_ENCODING', 'FASTQ encoding', 'phred')
    vars.Add('SAMPLE', 'Prefix name for results', '')
    vars.Add('READS', 'The input reads file', 'reads.fastq')
    vars.Add('ANNOTATION', 'Gene annotation in GTF/GFF/BED format', '')
    vars.Add('CPUS', 'Max parallel jobs to execute', '4')
    vars.Add('BOWTIE_INDEX', 'The Bowtie index', '')
    vars.Add('BOWTIE_EXTRA_PARAMS', 'Parameters to be passed to Bowtie.', 
             '-n 0 -l 26 -e 70 --best --strata -k 5 -m 5')

   
    env = Environment(ENV=os.environ,
                      variables=vars)
    
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

## ALIGN CLEAN READS TO REFERENCE GENOME
genome_alignments_dir = 'genome_alignments'
env_bowtie = env.Clone()

## -a and -y parameters can slow down the alignments and 
## generate large alignment files: do not use
## -k 5 : report up to 5 alignments for each read
## -m 5 : discard reads with more than 5 genomic alignments
## -n 0 : no mismatch allowed
env_bowtie.AppendUnique(BOWTIE_EXTRA_PARAMS = ['-S'])

env_bowtie['BOWTIE_INDEX'] = env_bowtie['BOWTIE_GENOME_INDEX']
if env_bowtie['BOWTIE_INDEX'].endswith('.1.ebwt'):
    env_bowtie['BOWTIE_INDEX'] = env_bowtie['BOWTIE_INDEX'][:-7]

env_bowtie.Replace(READS = [File(f).abspath for f in Flatten(env['READS'])])

genome_alignments = SConscript(os.path.join(genome_alignments_dir,
                                            'mirandmore_bowtie.scons'),
                    variant_dir = genome_alignments_dir,
                    src_dir = env['SCONSCRIPT_HOME'],
                    duplicate = 0,
                    exports = 'env_bowtie')

gene_intersect = None
if not env['ANNOTATION'] == '':
    ## INTERSECT WITH GENE ANNOTATION AND COUNT
    gene_intersect_cmd = 'zcat ${SOURCES[0]} | '\
                         'samtools view -SbuhF 4 $( -@ $CPUS $) - | '\
                         'bedtools bamtobed -i - | '\
                         'sort -k1,1 -k2,2n $( --parallel=$CPUS $) | '\
                         'bedtools coverage -sorted -a ${SOURCES[1]} -b stdin '\
                         '> $TARGET'

    gene_intersect = env.Command(env['SAMPLE'] + '.counts.gtf',
                                 [genome_alignments[0], env['ANNOTATION']],
                                 gene_intersect_cmd)

results = {'ALIGNMENTS': genome_alignments, 
           'GENE_COUNTS': gene_intersect}
Return("results")

