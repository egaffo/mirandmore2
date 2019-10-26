'''
'''
import os

Import('*')

try:
    env = env_quantify.Clone()
except NameError as ne:
        print 'mirandmore_quantify: command line mode. Not yet implemented', ne
	vars = Variables('vars.py')
	
	env = Environment(variables = vars,ENV = os.environ)

	Help(vars.GenerateHelpText(env))

#SRC_DIR = os.path.join(env['ENV']['MIRANDMORE_HOME'],'scons','prj')
SRC_DIR = env['SCONSCRIPT_HOME']

mirandmore_align_to_hairpins    = 'align_to_hairpins.py'
mirandmore_serialize_alignments = 'serialize_alignments.py'
mirandmore_do_exact_blob        = 'do_exact_blob.py'
mirandmore_do_cooked_blob       = 'do_cooked_blob.py'
mirandmore_do_tables            = 'do_tables.py'
mirandmore_extract_seq_of_novel_rnas = 'extract_seq_of_novel_rnas.py'
mirandmore_spit_unfiltered_variants  = 'spit_unfiltered_variants.py'

results = []
results_dict = {}

## ALIGN READS TO EXTENDED PRECURSORS
aligned_hairpins_dir = "hairpin_alignments"
env_align_to_hairpins = env.Clone()
env_align_to_hairpins['READS'] = env['FILTERED_READS']
aligned_to_hairpins = SConscript(os.path.join(aligned_hairpins_dir, mirandmore_align_to_hairpins),
                                 variant_dir = aligned_hairpins_dir, 
                                 src_dir = SRC_DIR,
                                 duplicate = 0, 
                                 exports = '''env_align_to_hairpins''')
results.append(aligned_to_hairpins)
results_dict['HAIRPIN_ALIGNMENTS'] = aligned_to_hairpins

## ALIGN TO GENOME THE READS NOT ALIGNED TO THE PRECURSORS AND
## INTERSECT ALIGNEMTS TO ANNOTATION (IF ANY) TO GET THE EXPRESSED GENES
pre_unmapped_dir = 'pre_unmapped_reads' 
sam2fastq_cmd = "zcat ${SOURCES[0]} | "\
                "samtools view $(-@ $CPUS $) -Suf 4 - | "\
                "samtools fastq - | "\
                "gzip -c > ${TARGETS[0]}"

sam2fastq = env.Command([os.path.join(pre_unmapped_dir, '${SAMPLE}_pre_unmapped.fastq.gz')],
                        aligned_to_hairpins[0],
                        sam2fastq_cmd)

results.append(sam2fastq)
results_dict['PRE_UNMAPPED'] = sam2fastq

#pre_unmapped_dir = 'pre_unmapped'
#   
#env_pre_unmapped = env.Clone()
#env_pre_unmapped['SAMPLE'] = env['SAMPLE']
#env_pre_unmapped['PRE_SAM'] = aligned_to_hairpins[0]
#env_pre_unmapped['ANNOTATION'] = env['GENE_ANNOTATION']
#env_pre_unmapped['BOWTIE_EXTRA_PARAMS'] = env['NM_BOWTIE_PARAMS'].split()
#
#sample_pre_unmapped = SConscript(os.path.join(pre_unmapped_dir, 
#                                              'mirandmore_pre_unmapped.scons'),
#                                 src_dir = SRC_DIR,
#                                 variant_dir = pre_unmapped_dir, 
#                                 duplicate = 0,
#                                 exports= 'env_pre_unmapped')
#
#results.append(sample_pre_unmapped.values())
#results_dict['PRE_UNMAPPED'] = sample_pre_unmapped
#
#Clean('.', pre_unmapped_dir)


## GENERATE PRECURSOR ALIGNMENTS' SERIALIZED OBJECT (BLOB)
env_serialize_alignments = env.Clone()
env_serialize_alignments['SAM'] = aligned_to_hairpins[0]

#hairpin_alignments_blobs = SConscript(os.path.join(aligned_hairpins_dir, 
#                                                   mirandmore_serialize_alignments),
#                                      variant_dir = aligned_hairpins_dir, 
#                                      src_dir = SRC_DIR,
#                                      duplicate = 0, 
#                                      exports = '''env_serialize_alignments ''')
#results.append(hairpin_alignments_blobs)
#results_dict['HAIRPIN_ALIGNMENTS_BLOB'] = hairpin_alignments_blobs
Clean('.', aligned_hairpins_dir)

## GENERATE SERIALIZED OBJECT (BLOB) OF EXACT ALIGNMENTS
blob_dir = "blobs"

env_exact_blob = env.Clone()
env_exact_blob['HAIRPIN_ALIGNMENTS'] = aligned_to_hairpins[0]
exact_blob = SConscript(os.path.join(blob_dir, mirandmore_do_exact_blob),
                      variant_dir = blob_dir ,src_dir = SRC_DIR,
                      duplicate = 0, 
                      exports = '''env_exact_blob ''')

results.append(exact_blob)
results_dict['HAIRPIN_ALIGNMENTS_EXACT_BLOB'] = exact_blob

#cooked_blob_genomic_hits_blob = GENOMIC_HITS_BLOB
#cooked_blob_exact_blob        = exact_blob
#cooked_blob_mature_table      = MATURE_TABLE
#cooked_blob_sample            = SAMPLE
#cooked_blob = SConscript(os.path.join(blob_dir, mirandmore_do_cooked_blob),
#                      variant_dir = blob_dir, src_dir = SRC_DIR,
#                      duplicate = 0, exports = '''env cooked_blob_genomic_hits_blob '''
#                                               '''cooked_blob_exact_blob '''
#                                               '''cooked_blob_mature_table '''
#                                               '''cooked_blob_sample ''')
#results.append(cooked_blob)

Clean('.', blob_dir)

## CREATE MIRNA EXPRESSION TABLES
table_dir = "tables"
env_do_tables = env.Clone()
env_do_tables['EXACT_BLOB'] = exact_blob
##do_tables_pre_blob = hairpin_alignments_blobs[0]
tables = SConscript(os.path.join(table_dir, mirandmore_do_tables),
                    variant_dir = table_dir, 
                    src_dir = SRC_DIR,
                    #duplicate = 0, exports = '''env do_tables_pre_blob do_tables_exact_blob '''
                    duplicate = 0, 
                    exports = '''env_do_tables ''')

results_dict['TABLES'] = tables
#[env['SAMPLE'] + "_mir_table.txt", 
#env['SAMPLE'] + "_mir_table_excel.txt", 
#env['SAMPLE'] + "_pre_summary.blob",
#env['SAMPLE'] + "_assign.log"]

results.append(tables)

Clean('.', table_dir)

## GET NUCLEOTIDE SEQUENCES OF THE NEW SMALL RNAS
sequences_dir = 'sequences'

env_extract_seq_of_novel_rnas = env.Clone()
env_extract_seq_of_novel_rnas['MIR_TABLE'] = tables[0]
novel_seqs = SConscript(os.path.join(sequences_dir, mirandmore_extract_seq_of_novel_rnas),
                        variant_dir = sequences_dir, 
                        src_dir = SRC_DIR,
                        duplicate = 0, 
                        exports = '''env_extract_seq_of_novel_rnas ''')

results.append(novel_seqs)
results_dict['NOVEL_SEQS'] = novel_seqs
Clean('.', sequences_dir)

##############################
## Include moRNAs and new sister miRNAs in unfiltered variants:
new_matures_table_source = [tables[1], env['MATURE_TABLE']]
new_matures_table_target = os.path.join(table_dir, 'mature-table-with-new.txt')
new_matures_table_cmd = 'pre_summary_to_mature_table.py -p ${SOURCES[0]} -m ${SOURCES[1]} -o ${TARGET}'
new_matures_table = env.Command(new_matures_table_target, 
                                new_matures_table_source, 
                                new_matures_table_cmd)

## process_algn.py alignments extended_mature_table.txt -> variants_pre.blob
unfiltered_variants_dir = 'unfiltered_variants'
env_serialize_alignments['MATURE_TABLE'] = new_matures_table[0]
new_hairpin_alignments_blobs = SConscript(os.path.join(unfiltered_variants_dir, 
                                                       mirandmore_serialize_alignments),
                                          variant_dir = unfiltered_variants_dir, 
                                          src_dir = SRC_DIR,
                                          duplicate = 0, 
                                          exports = '''env_serialize_alignments ''')

results.append(new_hairpin_alignments_blobs)
results_dict['HAIRPIN_ALIGNMENTS_BLOB'] = new_hairpin_alignments_blobs

Depends(new_hairpin_alignments_blobs, new_matures_table)

## unfiltered_variants_pre_blob = variants_pre.blob
env_spit_unfiltered_variants = env.Clone()
env_spit_unfiltered_variants['PRE_BLOB'] = new_hairpin_alignments_blobs[0]
unfiltered_variants = SConscript(os.path.join(unfiltered_variants_dir, 
                                              mirandmore_spit_unfiltered_variants),
                                 variant_dir = unfiltered_variants_dir, 
                                 src_dir = SRC_DIR,
                                 duplicate = 0, 
                                 exports = 'env_spit_unfiltered_variants ')

results.append(unfiltered_variants)
results_dict['UNFILTERED_VARIANTS'] = unfiltered_variants
Depends(unfiltered_variants, new_hairpin_alignments_blobs)

Clean('.', unfiltered_variants_dir)

Return('results results_dict')
