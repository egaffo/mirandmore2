'''
'''

Import('*')

try:
    env = env.Clone()
    GFF_FILE = mirandmore_extend_hairpins_GFF_FILE
    EXTEND_N = mirandmore_extend_hairpins_EXTEND_N
    GENOME_FASTA = mirandmore_extend_hairpins_GENOME_FASTA

except NameError:
    vars = Variables('vars.py')
    vars.Add('GFF_FILE', 'The miRBase miRNA annotation in GFF3 format', 'hsa.gff3')
    vars.Add('EXTEND_N', 'The number of bases to extend precursors', '30')
    vars.Add('GENOME_FASTA', 'The genome in single FASTA formatted file', 'hg38.fa')
    env = Environment(variables = vars,ENV = os.environ)
    
    Help(vars.GenerateHelpText(env))

    GFF_FILE = File(env['GFF_FILE'])
    EXTEND_N = env['EXTEND_N']
    GENOME_FASTA = File(env['GENOME_FASTA'])


chromosome_lengths_cmd = '''fasta_len.py $SOURCE > $TARGET'''
chromosome_lengths = env.Command('chrlen.genome', GENOME_FASTA, chromosome_lengths_cmd)

extended_pre_coordinates_cmd = ''' grep miRNA_primary_transcript ${SOURCES[0]} | '''\
                               ''' bedtools slop -s -i - -g ${SOURCES[1]} -b ''' + EXTEND_N +\
                               ''' | sed 's_Name=\(.\+\);*$$_Name=\\1-ext_' > ${TARGETS[0]} '''\
                               '''&& grep -P "miRNA\\t" ${SOURCES[0]} >> ${TARGETS[0]} '''
                               ## Note on the 'sed' line: double dollar is to escape $ expansion
                               ## by Scons, as for double backslash.
                               ## Then we must append the miRNA coordinates to the annotation.
                               ## 'ID' tags are not changed, hence the file is consistent.
                               #TODO: sort the GFF by coordinates
extended_pre_coordinates = env.Command('hairpin.extended.gff3', 
                                       [GFF_FILE, chromosome_lengths],
                                       extended_pre_coordinates_cmd)

fastas_cmd = '''grep miRNA_primary_transcript ${SOURCES[1]} | '''\
             '''sed "s/miRNA_primary_transcript\(.\+\)Name=\(.\+\)/\\2\\1Name=\\2/" '''\
             ''' | bedtools getfasta -s -name -fi ${SOURCES[0]} -bed - | '''\
             '''sed "s_([-+]*)__" > ${TARGET} '''
            ## bedtools v2.27 adds -name and -name+ options. The next cmd fixes FASTA
            ## names for -name option

            # '''sed -r "s@::.+:[0-9]+-[0-9]+\([-+]\)@@" - >${TARGET}''' 
	        ## The last sed removes ::<chrom>:<start>-<end> added to the names by bedtools getfasta

fastas = env.Command('hairpin.extended.fa', 
                     [GENOME_FASTA, extended_pre_coordinates], 
                     fastas_cmd)

make_fasta_index = env.Command('${SOURCE.file}.fai', fastas, 
                               'samtools faidx $SOURCE')

Return('extended_pre_coordinates fastas make_fasta_index')
