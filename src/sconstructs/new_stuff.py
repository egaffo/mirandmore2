import os

Import("*")

try:
    env = env_new_stuff.Clone()
except NameError as ne:
    vars = Variables('vars.py')
    vars.Add('EXTENDED_N', '', '')
    vars.Add('META', '', '')
    vars.Add('NEW_PREMIRS_MIRS_GFF', '', '')
    vars.Add('EXT_HAIRPINS_FASTA', '', '')
    vars.Add('MIR_TABLES', '', '')
    vars.Add('NORMALIZED_DATA', '', '')

    env = Environment(variables = vars, ENV = os.environ)
 
    Help(vars.GenerateHelpText(env))

    env['MIR_TABLES'] = env['MIR_TABLES'].split(',')

## create new sRNA BED file, coordinates relative to extended precursors
if env['NEW_PREMIRS_MIRS_GFF'] and not env['NEW_PREMIRS_MIRS_GFF'] == '':
    cmdLine = '''get_new_srna_bed.R -m $META --extended_n $EXTENDED_N '''\
              '''--annotation_gff3 $NEW_PREMIRS_MIRS_GFF -i \"$SOURCES\" -d $TARGET.dir '''
else:
    cmdLine = '''touch ${TARGETS[0]} && touch ${TARGETS[1]}'''

get_new_srna_bed = env.Command(["new-mir-table-without-seqs.txt", "new_mir_squished_coords.bed"], 
                               env['MIR_TABLES'], 
                               cmdLine)

## retrieve FASTAs of new sRNAs
if env['NEW_PREMIRS_MIRS_GFF'] and not env['NEW_PREMIRS_MIRS_GFF'] == '':
    ## The last sed removes ::<chrom>:<start>-<end> added to the names 
    ## by bedtools getfasta -name
    #bedtoolsCmdLine = '''bedtools getfasta -fi %(EXT_HAIRPINS_FASTA)s -bed $SOURCE -s -name |'''\
    #                  '''sed -r "s@::.+:[0-9]+-[0-9]+\([-+]?\)@@" - >${TARGET}''' % locals()
    ## bedtools v2.27 adds -name and -name+ options. The next cmd fixes FASTA
    ## names for -name option
    bedtoolsCmdLine = '''bedtools getfasta -fi $EXT_HAIRPINS_FASTA '''\
                      '''-bed $SOURCE -s -name | sed "s_()__" > ${TARGET}'''
else:
    bedtoolsCmdLine = 'touch $TARGET'

new_srna_fasta = env.Command("new_srna.fa", 
                             get_new_srna_bed[1],
                             bedtoolsCmdLine)

## put FASTA into tabular format
fastaToCsvCmdLine = 'fasta_to_csv.py -i ${SOURCES[0]} -o ${TARGET}' 
fastaToCsv = env.Command("new-srna-seqs.csv", 
			 new_srna_fasta, 
			 fastaToCsvCmdLine)

## compute mean expression relative to condition
meansPerCondtionCmdLine = "means-per-condition.R -m ${SOURCES[1]} --data ${SOURCES[0]} --output ${TARGET}"
means = env.Command('means_per_condition.txt',
		    [env['NORMALIZED_DATA'], env['META']], 
		    meansPerCondtionCmdLine)

## write new stuff tables
sources_for_new_srna_table = [means[0], #means_per_condition.txt
                              fastaToCsv[0], #new-srna-seqs.csv
                              get_new_srna_bed[0]] #new-mir-table-without-seqs.txt

targets_for_new_srna_table = ["new-srna-table-with-seq.txt", 
                              "new-mir-table-with-seq.txt", 
                              "mor-table-with-seq.txt"]
newSrnaCmdLine = "create-new-srna-table.R -w ${TARGETS[0].dir}"
new_srna_tables = env.Command(targets_for_new_srna_table, 
                              sources_for_new_srna_table, 
                              newSrnaCmdLine)

Return('get_new_srna_bed new_srna_fasta fastaToCsv means new_srna_tables')
