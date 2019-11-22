import os, csv, re, itertools, collections
from collections import defaultdict

def get_matching_nodes(nodelist, rexpression):
    files = []
    for node in Flatten(nodelist):
        if re.match(rexpression, node.path):
            files.append(node)
    return files


vars = Variables('vars.py')
vars.Add('META', 
        ''''The metadata table file where you specify the project samples, conditions, etc..''', 
        'meta.csv')
vars.Add('CPUS','The number of cores to be used by multithreaded tools','1')
vars.Add('ADAPTER','Read adapter to trim. Either "classic", "truseq", "bmr", "vigneault", or '\
	'the adapter sequence itself. Only reads bearing the adapter will be kept to search '\
	'the small RNAs.', 'truseq')
vars.Add('QUALITY_ENCODING','Quality encoding of the reads.','phred')

## serialize alignemts parameters and multiple mapping
vars.Add('MULTIPLE_GENOMIC_HITS_THRESHOLD','multiple genomic hits threshold', 5)
vars.Add('MIN_COUNT','minimum number of unique reads', 10)
vars.Add('ALLOWED_OVERHANG','allowed overhang', 3)
vars.Add('ALLOWED_OVERHANG_MORNA','allowed overhang morna', 3)
vars.Add('EXACT','EXACT CONSTANT', 0)
vars.Add('SHORTER_OR_LONGER','SHORTER_OR_LONGER CONSTANT', 1)
vars.Add('MIS_1','MISMATCH_1 CONSTANT', 2)
vars.Add('MIS_2','MISMATCH 2 CONSTANT', 3)
vars.Add('FIVE_PRIME','five prime CONSTANT', 5)
vars.Add('THREE_PRIME','three prime CONSTANT', 6)
vars.Add('MIN_MORNA_LEN','min morna len', 18)
vars.Add('EXTENDED_N', 'Number of positions for extended hairpins', '30') 

## mir discretization and assignment
vars.Add('NEW_RNA_OBJECT_THRESHOLD', 'new rna object threshold',4)
vars.Add('SISTER_OVERHANG_LEN','sister overhang len',3)
vars.Add('SISTER_MATCHES_THRESHOLD','sister matches threshold',14)
vars.Add('SPECIES','Species analyzed','hsa')

## filter_trimmed.py parameters
vars.Add('MAX_LEN_FILTER','max len filter',30)
vars.Add('MEAN_QUAL_FILTER','mean qual filter',30)

## genome sequence in fasta format
vars.Add('GENOME', 'genome in fasta', 'hg38.fa')
vars.Add('PRECURSORS_FASTA', 'The miRNA precursor sequences in FASTA', 'hairpins.fa.gz')
vars.Add('GFF_FILE', 'The miRNA annotation file in GFF3', 'hsa.gff3')
vars.Add('CHRM_PREFIX_TO_REMOVE', "Pattern string to remove from annotation chromosomes."\
         "Needed if annotation and genome chromosome names differ (e.g. set 'chr' with "\
         "Ensembl genomes). Leave '' if no changes needed.", "''")

## miRDeep2 parameters
vars.Add('MRD_SCORE_THRESHOLD', 'score threshold for miRDeep2 precursors', '50')

## pre-computed Bowtie genome index
vars.Add('BOWTIE_GENOME_INDEX', 'The path and name of an already computed Bowtie genome index. '\
	 'Setting this parameter will skip the (long) generation of the index', '')

## workflow logic
vars.Add('PREDICT_NOVEL_PREMIRS', 'Wether to enable novel precursor prediction', 'True')

vars.Add('GENE_ANNOTATION', 'Gene annotation GTF file to characterize '\
         'non-miRNA reads. Leave empty to only map length discarded reads', '')
vars.Add('NM_BOWTIE_PARAMS', 'Parameters to be passed to Bowtie to align non-miRNA reads.', 
         '-n 0 -l 26 -e 70 --best --strata -k 5 -m 5 --mm')

vars.Add('NOADAPTER', 'Set "True" if adapters were already trimmed in  input reads. '\
	     'This will skip adapter trimming step.', False)

vars.Add('TRIMMER', 'The tool to use for trimming read adapter. '\
                    'Chose between cutadapt and fastxtoolkit', 
         'cutadapt')

vars.Add('CUTADAPT_EXTRA_PARAMS', '', '')
vars.Add('NON_MIRNAS', 
         'Enable characterization of reads not aligned to miRNA precursors', 
         False)


env = Environment(ENV = os.environ,
                  variables = vars)

Help(vars.GenerateHelpText(env))

SCONSCRIPT_HOME = os.path.join(env['ENV']['MIRANDMORE_HOME'],'src','sconstructs')
env['SCONSCRIPT_HOME'] = SCONSCRIPT_HOME
env['MIRANDMORE_BIN'] = os.path.join(env['ENV']['MIRANDMORE_HOME'],'bin')

env.SetDefault(NON_MIRNAS = False)

env.Replace(GENOME =  File(env['GENOME']).abspath)
if env['GENE_ANNOTATION'] != '':
    env.Replace(GENE_ANNOTATION = File(env['GENE_ANNOTATION']).abspath)
env.Replace(META = File(env['META']).abspath)

if env['NOADAPTER'] == 'True':
    env.Replace(NOADAPTER = True)
else:
    env.Replace(NOADAPTER = False)

## format read encoding and adapter parameters

if env['ADAPTER'] == "classic":
    primer = "TCGTATGCCGTCTTCTGCTTG"
elif env['ADAPTER'] == "truseq":
    primer = "TGGAATTCTCGGGTGCCAAGG"
elif env['ADAPTER'] == "bmr":
    primer = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAG"
elif env['ADAPTER'] == "vigneault":
    primer  = "AACGGGCTAATATTTATCGGTGGC"
else:
    primer = env['ADAPTER']

env.Replace(ADAPTER = primer)


# directory to store miR&moRe configuration/setting files
mirandmore_dbs_dir = 'mirandmore_dbs'
mirandmore_annotation_dir = os.path.join(mirandmore_dbs_dir, 'annotation_files')

env.Replace(GFF_FILE = File(env['GFF_FILE']).abspath)
if not env['CHRM_PREFIX_TO_REMOVE'] == "''":
    fix_annotation_chrom_names_cmd = '''sed 's#^''' + env['CHRM_PREFIX_TO_REMOVE'] +\
                                     '''##' $SOURCE > $TARGET'''
    fix_annotation_chrom_names = env.Command(os.path.join(mirandmore_annotation_dir, 
                                                'chrmfixed_${SOURCE.file}'), 
                                            File(env['GFF_FILE']), 
                                            fix_annotation_chrom_names_cmd)
    env.Replace(GFF_FILE = fix_annotation_chrom_names)

bowtie_indexes = []
if env['BOWTIE_GENOME_INDEX'] == '':
    genome_index_dir = os.path.join(mirandmore_dbs_dir, 'bwt_genome_index')
    env_bowtie_index = env.Clone()
    env_bowtie_index['FASTA_TO_INDEX'] = File(env['GENOME'])  
    genome_bowtie_index = SConscript(os.path.join(genome_index_dir, 'mirandmore_bowtie_index'),
                                     variant_dir = genome_index_dir,
                                     src_dir = SCONSCRIPT_HOME,
                                     duplicate = 0,
                                     exports = '''env_bowtie_index''')
    
    env.Replace(BOWTIE_GENOME_INDEX = genome_bowtie_index[0].abspath[:-7])
else:
    ## just use the first index file as mock
    env.Replace(BOWTIE_GENOME_INDEX = File(env['BOWTIE_GENOME_INDEX'] + ".1.ebwt").abspath[:-7])
    genome_bowtie_index = [File(env['BOWTIE_GENOME_INDEX'] + '.1.ebwt')]



bowtie_indexes.append(genome_bowtie_index)


samples_dir = 'samples'
samples = defaultdict(str)

with open(env['META']) as csvfile:
    reader = csv.DictReader(csvfile,delimiter=";")
    for row in reader:
        samples[row['sample']]=(os.path.abspath(row['file'])) 

#runs = []

## For each sample map unique reads to genome to predict novel precursors and miRNAs
mirandmore_quality                 = 'read_quality_stats.py'
if env['TRIMMER'] == 'fastxtoolkit':
    mirandmore_trim                = 'trim.py'
else:
    if env['TRIMMER'] != 'cutadapt':
        print('Trimmer', env['TRIMMER'], 'not available. Will use cutadapt.')
    mirandmore_trim                = 'cutadapt.py'
mirandmore_filter                  = 'filter.py'
mirandmore_collapse                = 'collapse.py'
mirandmore_align_uniques           = 'align_uniques.py'
mirandmore_serialize_uniques_align = 'serialize_uniques_align.py'

results          = {} ## dict (of dicts) to store results of each sample analysis 
new_matures_gffs = [] ## list of new mature annotation files (one per sample)

for sample in sorted(samples.keys()):

    readset = File(samples[sample])
    sample_dir = os.path.join(samples_dir,sample)

    results[sample] = {} ## set key to sample name and then populate with analysis step results

    ## READ QUALITY STATISTICS
    sample_quality_dir = os.path.join(sample_dir, 'read_quality_stats')    
    env_read_qual_stat = env.Clone()
    env_read_qual_stat['READS'] = readset
    quality = SConscript(os.path.join(sample_quality_dir, mirandmore_quality),
                       variant_dir = sample_quality_dir, src_dir = SCONSCRIPT_HOME, 
                       duplicate = 0, exports = 'env_read_qual_stat')

    results[sample]['quality'] = quality ## save quality results for the sample
    Clean('.', sample_quality_dir)
    
    ## TRIM READ ADAPTERS (IF NEEDED)
    trim_dir = os.path.join(sample_dir, "trimmed_reads")
    if not env['NOADAPTER']:
        env_trim = env.Clone()
        env_trim['READS'] = readset
        env_trim['SAMPLE'] = sample
        offset2enc = {'phred':33, 'solexa':64}
        env_trim.Replace(QUALITY_ENCODING = offset2enc[env['QUALITY_ENCODING']])

        if env['TRIMMER'] == 'fastxtoolkit':
            trimmed = SConscript(os.path.join(trim_dir, mirandmore_trim),
                                 variant_dir = trim_dir, 
                                 src_dir = SCONSCRIPT_HOME,
                                 duplicate = 0, 
                                 exports = 'env_trim')
        else:
            env_cutadapt = env_trim
            trimmed = SConscript(os.path.join(trim_dir, mirandmore_trim),
                                 variant_dir = trim_dir, 
                                 src_dir = SCONSCRIPT_HOME,
                                 duplicate = 0, 
                                 exports = 'env_cutadapt')

    else:
        trimmed = [readset]
    
    results[sample]['trimmed'] = trimmed
    Clean('.', trim_dir)
    
    ## APPLY QUALITY AND LENGTH FILTERS TO INPUT READS
    filter_dir = os.path.join(sample_dir, 'filtered_reads')
    env_filter = env.Clone()
    env_filter['READS'] = trimmed[0]
    env_filter['SAMPLE'] = sample
    env_filter['LONG_READS_FILE'] = ''
    if env['TRIMMER'] == 'fastxtoolkit':
        env_filter['LONG_READS_FILE'] = sample + "_length_discarded_reads.fq.gz" 
    filtered = SConscript(os.path.join(filter_dir, mirandmore_filter),
                          variant_dir = filter_dir, 
                          src_dir = SCONSCRIPT_HOME,
                          duplicate = 0, 
                          exports = 'env_filter')
    
    results[sample]['filtered'] = filtered
    Clean('.', filter_dir)
    
    if env['TRIMMER'] == 'fastxtoolkit':
        results[sample]['length_discarded'] = filtered[2]
    else:
        results[sample]['length_discarded'] = trimmed[2]
    
    ## COLLAPSE READS INTO UNIQUE SEQUENCES
    aligned_uniques_dir = os.path.join(sample_dir, "genomic_alignments")
    
    env_collapse = env.Clone()
    env_collapse['READS2COLLAPSE'] = filtered[0]
    env_collapse['SAMPLE'] = sample
    collapsed = SConscript(os.path.join(aligned_uniques_dir, mirandmore_collapse),
                           variant_dir = aligned_uniques_dir, 
			   src_dir = SCONSCRIPT_HOME,
                           duplicate = 0, 
			   exports = 'env_collapse')
    
    results[sample]['collapsed'] = collapsed
    
    ## ALIGN UNIQUE SEQUENCES TO THE WHOLE REFERENCE GENOME
    env_aligned_uniques = env.Clone()
    #env_aligned_uniques['CPUS']  = CPUS
    env_aligned_uniques['SAMPLE'] = sample
    env_aligned_uniques['GENOME_BWT_INDEX'] = str(genome_bowtie_index[0]).replace('.1.ebwt','')
    env_aligned_uniques['UNIQUE_READS'] = collapsed
    aligned_uniques = env.SConscript(os.path.join(aligned_uniques_dir, 
                                                  mirandmore_align_uniques),
                                     src_dir = SCONSCRIPT_HOME,
                                     variant_dir = aligned_uniques_dir, duplicate = 0,
                                     exports = '''env_aligned_uniques ''')
    Depends(aligned_uniques, bowtie_indexes)
    results[sample]['aligned_uniques'] = aligned_uniques

    ## PREPARE SERIALIZED OBJECT (BLOB) OF UNIQUE SEQUENCE GENOMIC ALIGNMENTS
    ## FOR LATER USE
    env_serialize_uniques_align = env.Clone()
    env_serialize_uniques_align['GENOMIC_HITS'] = aligned_uniques[0]
    env_serialize_uniques_align['SAMPLE'] = sample
    genomic_hits_blob = SConscript(os.path.join(aligned_uniques_dir, 
                                                mirandmore_serialize_uniques_align),
                                   variant_dir = aligned_uniques_dir, 
                                   src_dir = SCONSCRIPT_HOME,
                                   duplicate = 0, 
                                   exports = 'env_serialize_uniques_align')
    
    results[sample]['genomic_hits_blob'] = genomic_hits_blob
    Clean('.', aligned_uniques_dir)
    
    Clean('.', sample_dir)

Clean('.', mirandmore_dbs_dir)

results_dir = 'results'

## Precursor predictions
env['NEW_PREMIRS_MIRS_GFF'] = None
if env['PREDICT_NOVEL_PREMIRS'] == 'True':
    novel_premirs_dir = os.path.join(mirandmore_dbs_dir, 'novel_premirs')
    env_predict_novel_precursors = env.Clone()
    env_predict_novel_precursors['PREPROCESSING_FILES'] = results
    novel_premirs = SConscript(os.path.join(novel_premirs_dir, 
                                            'predict_novel_precursors.py'), 
                                   variant_dir = novel_premirs_dir, 
                                   src_dir = SCONSCRIPT_HOME,
                                   duplicate = 0,
                                   exports = '''env_predict_novel_precursors ''')

    env.Replace(PRECURSORS_FASTA = novel_premirs['PRECURSORS_FASTA'])
    env.Replace(GFF_FILE = novel_premirs['GFF_FILE'])
    env.Replace(NEW_PREMIRS_MIRS_GFF = novel_premirs['NEW_PREMIRS_MIRS_GFF'])
else:
    ## USE ONLY KNOWN MIRS
    precursors_fasta_cmd = '''grep miRNA_primary_transcript ${SOURCES[1]}| '''\
                           '''sed "s/miRNA_primary_transcript\(.\+\)'''\
                           '''Name=\(.\+\)/\\2\\1Name=\\2/" '''\
                           '''| bedtools getfasta -s -name -fi '''\
                           '''${SOURCES[0]} -bed - |'''\
                           '''sed -r "s@::.+:[0-9]+-[0-9]+\([-+]\)@@" - >${TARGET}'''
    
    precursors_fasta_sources = [File(env['GENOME']), env['GFF_FILE']]
    precursors_fasta_target = os.path.join(mirandmore_dbs_dir, 'hairpins.fa')
    precursors_fasta = env.Command(precursors_fasta_target, 
                                   precursors_fasta_sources, 
                                   precursors_fasta_cmd)
    
    env.Replace(PRECURSORS_FASTA = precursors_fasta)
    
    ## TODO?
    #matures_fasta_cmd = '''grep "miRNA\s" ${SOURCES[1]}| '''\
    #                   '''sed "s/miRNA\(.\+\)'''\
    #                   '''Name=\([^;]\+\)/\\2\\1Name=\\2/" '''\
    #                   '''| bedtools getfasta -s -name -fi '''\
    #                   '''${SOURCES[0]} -bed - |'''\
    #                   '''sed -r "s@::.+:[0-9]+-[0-9]+\([-+]\)@@" - >${TARGET}'''

    #matures_fasta = env.Command(os.path.join(mirandmore_dbs_dir, 'matures.fa'), 
    #                            [File(env['GENOME']), env['GFF_FILE'], 
    #                            matures_fasta_cmd)



## Prepare MiR&moRe annotation files and indexes
env_annotation = env.Clone()
env_annotation['GENOME_FASTA'] = File(env['GENOME'])

mirandmore_annotation = SConscript(os.path.join(mirandmore_annotation_dir, 
                                                'annotations.py'), 
                                   variant_dir = mirandmore_annotation_dir, 
                                   src_dir = SCONSCRIPT_HOME,
                                   duplicate = 0,
                                   exports = '''env_annotation ''')

env['MATURE_TABLE'] = mirandmore_annotation[0]['extended_annotation'][0][0]
env['HAIRPIN_EXTENDED_TO_FOLD'] = mirandmore_annotation[0]['extended_annotation'][1][0]
env['HAIRPIN_ANNOTATION_BLOB'] = mirandmore_annotation[0]['extended_annotation'][2][0]       
env['EXT_HAIRPINS_FASTA'] = mirandmore_annotation[1][1]

Clean('.', mirandmore_annotation_dir)

## BUILD EXTENDED HAIRPINS SEQUENCE INDEX FOR BOWTIE
ext_hairpins_index_dir = os.path.join(mirandmore_dbs_dir, 'bwt_ext_hairpins_index')

env_bowtie_index = env.Clone()
#env_bowtie_index['FASTA_TO_INDEX'] = File(mirandmore_annotation[1][1][0])
env_bowtie_index['FASTA_TO_INDEX'] = env['EXT_HAIRPINS_FASTA']
ext_hairpins_bowtie_index = SConscript(os.path.join(ext_hairpins_index_dir, 
                                                    'bowtie_index.py'),
                                 variant_dir = ext_hairpins_index_dir,
                                 src_dir = SCONSCRIPT_HOME,
                                 duplicate = 0,
                                 exports = '''env_bowtie_index ''')

env['EXT_PRE_BWT_IDX'] = str(ext_hairpins_bowtie_index[0]).replace('.1.ebwt', '')

bowtie_indexes.append(genome_bowtie_index)

Clean('.', ext_hairpins_index_dir)

## QUANTIFY MIRNA EXPRESSION FOR EACH SAMPLE
quantify_results = {}
for sample in sorted(samples.keys()):

    sample_dir = os.path.join(samples_dir,sample)
    
    env_quantify = env.Clone()
    
    env_quantify['SAMPLE'] = sample
    env_quantify['FILTERED_READS']  = results[sample]['filtered'][0]
    env_quantify['GENOMIC_HITS_BLOB'] = results[sample]['genomic_hits_blob']
    env_quantify['UNIQUE_READS'] = results[sample]['collapsed']

    quantify_sample = SConscript(os.path.join(sample_dir,'quantify.py'),
                           src_dir = SCONSCRIPT_HOME,
                           variant_dir = sample_dir, 
                           duplicate = 0,
                           exports = '''env_quantify ''')
    
    Depends(quantify_sample[0], ext_hairpins_bowtie_index)

    ## Must set explicit dependencies since the code in 'rna.py' library may hide dependencies
    Depends(quantify_sample[0], 
            [mirandmore_annotation[0][f] for f in mirandmore_annotation[0].keys()] +\
            [mirandmore_annotation[1]])
    
    quantify_results[sample] = quantify_sample[1]
    results[sample]['quantify'] = quantify_sample[0]
    
    if env['NON_MIRNAS']:
        ## CHARACTERIZE THE READS DISCARDED FOR EXCESSIVE LENGTH AND NOT
        ## ALIGNED TO MIRNA PRECURSORS:
        ## 1. align to the whole genome; 
        ## 2. intersect with the given annotation (if any)
        env_non_mirna = env.Clone()
        env_non_mirna['READS'] = [results[sample]['length_discarded'], 
                                  quantify_results[sample]['PRE_UNMAPPED']]
        env_non_mirna['SCONSCRIPT_HOME'] = SCONSCRIPT_HOME
        env_non_mirna['SAMPLE'] = sample
        env_non_mirna['ANNOTATION'] = env['GENE_ANNOTATION']
        env_non_mirna['BOWTIE_EXTRA_PARAMS'] = env['NM_BOWTIE_PARAMS'].split()
        non_mirna_dir = os.path.join(sample_dir, 'non_mirnas')
        non_mirna = SConscript(os.path.join(non_mirna_dir, 
                                            'non_mirnas.py'),
                              variant_dir = non_mirna_dir, 
                              src_dir = SCONSCRIPT_HOME,
                              duplicate = 0, exports = 'env_non_mirna')
        results[sample]['non_mirna'] = non_mirna
        Clean('.', non_mirna_dir)

if env['NON_MIRNAS']:
    ## MERGE SAMPLES' NON-MIRNA COUNTS
    sample_file_str = []
    gene_count_files = []
    for sample in samples.keys():
        count_file = results[sample]['non_mirna']['GENE_COUNTS']
        if count_file:
            count_file = count_file[0]
            sample_file_str.append("'" + sample + "' = '" + count_file.abspath + "'")
            gene_count_files.append(count_file)       
        else:
            sample_file_str.append("'" + sample + "' = 'None'")
    
    if len(gene_count_files) > 0:
        gene_count_vector = "files <- c(" + ','.join(sample_file_str) + ")"
        
        collect_gene_counts_cmd = '''Rscript -e " ''' + gene_count_vector +\
                '''; merged.gtf.counts.file <- '${TARGETS[0]}'; '''\
                '''count.matrix.file <- '${TARGETS[1]}'; '''\
                '''source(file.path('$MIRANDMORE_BIN', 'collect_gene_counts.R')) " '''
        
        collect_gene_counts_targets = [os.path.join(results_dir, f) for f in 
                                          ['non_mirna_gene_counts.gtf.csv',
                                           'non_mirna_gene_counts.csv']]
        collect_gene_counts = env.Command(collect_gene_counts_targets,
                                             gene_count_files,
                                             collect_gene_counts_cmd)


#read_quality_stats_dir = 'read_quality_stats'
### COLLECT TRIMMING STATISTICS
#env_summarize_depth = env.Clone()
#trim_results = [results[s]['trimmed'] for s in samples.keys()]
#env_summarize_depth['TRIMMING_REPORT_FILES'] = get_matching_nodes(trim_results, ".*_trimming_report\.txt")
#
#trimming_reports = SConscript(os.path.join(read_quality_stats_dir, 'mirandmore_summarize_depth'),
#                              src_dir = SCONSCRIPT_HOME,
#                              variant_dir = read_quality_stats_dir, 
#                              duplicate = 0,
#                              exports  = 'env_summarize_depth')
#
### COLLECT FILTERING STATISTICS
#filter_stats_files = []
#for sample in sorted(samples.keys()):
#	filter_stats_files.append(results[sample]['filtered'][1])
#
#collect_filter_stats_cmd = '''grep -H . ${SOURCES} > ${TARGET}'''
#collect_filter_stats = env.Command(os.path.join(read_quality_stats_dir, "filter_reports.txt"), 
#				   [filter_stats_files], 
#				   collect_filter_stats_cmd)
#
### COLLECT QUALITY READS STATISTICS RESULTS
#env_read_quality_stats = env.Clone()
#env_read_quality_stats['SAMPLES'] = samples.keys()
#
#read_quality_stats_nodes = [results[n]['quality'] for n in results.keys()]
#env_read_quality_stats['FASTQC_HTMLS'] = get_matching_nodes(read_quality_stats_nodes, '.*_fastqc\.html')
#
#read_quality_stats = SConscript(os.path.join(read_quality_stats_dir, 'mirandmore_quality_collect'),
#                                src_dir = SCONSCRIPT_HOME,
#                                variant_dir = read_quality_stats_dir, 
#                                duplicate = 0,
#                                exports = 'env_read_quality_stats')
#
#Clean('.', read_quality_stats_dir)


env['MIR_TABLES'] = [items['TABLES'][0] for items in quantify_results.values()]

## COLLECT SAMPLES' MIRNA RESULTS
env_summarize_tables = env.Clone()

env_summarize_tables['UNFILTERED_VARIANT_TABLES'] = [items['UNFILTERED_VARIANTS'] \
                                                     for items in quantify_results.values()]

summary_tables = SConscript(os.path.join(results_dir, 'summarize_tables.py'),
                            src_dir = SCONSCRIPT_HOME,
                            variant_dir = results_dir, 
                            duplicate = 0,
                            exports = '''env_summarize_tables ''')

#env['NORMALIZED_DATA'] = summary_tables[1]

## COLLECT NOVEL PREDICTIONS
env_new_stuff = env.Clone()
new_stuff = SConscript(os.path.join(results_dir,'new_stuff.py'),
                       src_dir = SCONSCRIPT_HOME,
                       variant_dir = results_dir, duplicate = 0,
                       exports = '''env_new_stuff ''')

### generate HTML with miRNA and miRNA-like expression report
#report_summary_cmd = '''Rscript -e 'results.dir <- dirname("$TARGET.abspath"); '''\
#		'''meta.file <- "${SOURCES[0].abspath}"; '''\
#		'''normalized.data.file <- "${SOURCES[1].abspath}"; '''\
#		'''normalized_variants.file <- "${SOURCES[2].abspath}"; '''\
#		'''new.srna.table.with.seq.file <- "${SOURCES[3].abspath}"; '''\
#		'''rmarkdown::render(input = "$MIRANDMORE_BIN/results_report.Rmd",'''\
#		'''output_file = "$TARGET.abspath", '''\
#		'''intermediates_dir = dirname("$TARGET.abspath") )' '''
#
#normalized_data_file = summary_tables[1] #normalized_data.txt
#normalized_variants_file = summary_tables[3] #normalized_variants.txt
#new_srna_table_with_seq_file = new_stuff[4][0] #new-srna-table-with-seq.txt
#
#report_summary = env.Command(os.path.join(results_dir, 'results_report.html'), 
#			     [env['META'],  
#			      normalized_data_file, 
#                              normalized_variants_file, 
#			      new_srna_table_with_seq_file], 
#			     report_summary_cmd)

## COLLECT PROCESSING STATISTICS
env_processing_stats = env.Clone()
env_processing_stats['RESULTS'] = results
env_processing_stats['QUANTIFY_RESULTS'] = quantify_results
processing_stats = SConscript(os.path.join(results_dir, 'processing_stats.py'),
                              src_dir = SCONSCRIPT_HOME,
                              variant_dir = results_dir, 
                              duplicate = 0,
                              exports = '''env_processing_stats get_matching_nodes''')

Clean('.', results_dir)
Clean('.', samples_dir)
