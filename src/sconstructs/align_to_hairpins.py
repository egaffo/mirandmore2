'''
Align the reads against the extended hairpins with:

bowtie -S -n 2 -l 26 -e 70 -p CPUS --best --strata -a -y --mm BOWTIE_DB reads2align2hairpins  SAMPLE.sam 2> known_extpre_alignments_stats.txt

# -S print alignments in SAM format
# -n N Maximum number of mismatches permitted in the "seed"
# -l N seed region length ( the number of bases on the high-quality end of the read to which the -n ceiling applies. ) 
# -e N Maximum permitted total of quality values at all mismatched read positions throughout the entire alignment
# -p N number of threads
# --best Make Bowtie guarantee that reported singleton alignments are "best" in terms of stratum (i.e. number of mismatches, or mismatches in the seed in the case of -n mode) and in terms of the quality values at the mismatched position(s).
# --strata If many valid alignments exist and are reportable (e.g. are not disallowed via the -k option) and they fall into more than one alignment "stratum", report only those alignments that fall into the best stratum.
# -y Try hard slower but more accurate.
# -mm Use memory-mapped I/O to load the index, the index is loaded into memory only once even if used by multiple processes

'''
Import('*')

try:
    env = env_align_to_hairpins.Clone()
except NameError as ne:
    print 'mirandmore_align_to_hairpins:', ne, 'not set. Command line mode not implemented'
    vars = Variables('vars.py')
    vars.Add('SAMPLE', 'basename', 'mm_excised_precursors')
    vars.Add('CPUS', '', '')
    vars.Add('EXT_PRE_BWT_IDX', '''The Bowtie index of the precursors' extended sequences''', '')
    vars.Add('READS', 'The reads to be aligned', '')

    env = Environment(variables = vars,ENV = os.environ)
    Help(vars.GenerateHelpText(env))

targets = [env['SAMPLE'] + ".sam.gz", 
           'known_extpre_alignments_stats.txt']

cmd_prefix = '''cat '''
if env['READS'].name.endswith('.gz'):
	cmd_prefix = '''zcat '''

aligned_to_hairpins_cmd = cmd_prefix +\
			  ''' ${SOURCE} | bowtie -S -n 2 -l 26 -e 70 $(-p $CPUS $) '''\
			  '''--best --strata -a -y --mm $EXT_PRE_BWT_IDX - '''\
                          '''2> ${TARGETS[1]} | gzip -c > ${TARGETS[0]} '''

aligned_to_hairpins = env.Command(targets, 
                                  env['READS'], 
                                  aligned_to_hairpins_cmd)

Return('aligned_to_hairpins')

