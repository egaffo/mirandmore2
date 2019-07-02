''' 
Aligns the collapsed reads against the genome with:

bowtie -v 0 -k 20 -p CPUS --mm -f bowtie_genome_idx SAMPLE_uniques.fa SAMPLE_genomic.out 2> genomic_alignments_stats.txt

# -v 0  no mismatches allowed  quality values are ignored mutually exclusive with -n
# -k 20 report up to 20 valid alignments 
# -f the query inputs are fasta files.
# -p number of threads to use.
# -mm Use memory-mapped I/O to load the index, the index is loaded into memory only once even if used by multiple processes

''' 

Import('*')

try:
    env = env_aligned_uniques.Clone()
except NameError as ne:
    print 'mirandmore_align_uniques: failed to import', ne
    vars = Variables('vars.py')
    vars.Add('CPUS', '', '')
    vars.Add('SAMPLE', '''The name of the sample to build the target file name. ''', '')
    vars.Add('GENOME_BWT_INDEX', '', '')
    vars.Add('UNIQUE_READS', '', '')
        
    env = Environment(ENV=os.environ,
                              variables=vars)
    
    Help(vars.GenerateHelpText(env))

targets = [env['SAMPLE'] + "_genomic.out", 
           'genomic_alignments_stats.txt']
aligned_uniqes_cmd = '''bowtie -v 0 -k 20 $(-p $CPUS --mm $) '''\
                     '''-f $GENOME_BWT_INDEX ${SOURCES[0]} '''\
                     '''${TARGETS[0]} 2> ${TARGETS[1]}'''
aligned_uniqes = env.Command(targets, 
                             [env['UNIQUE_READS']], 
                             aligned_uniqes_cmd)

Return('aligned_uniqes')
