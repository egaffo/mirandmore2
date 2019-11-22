'''
This script collects and reports processing statistics
* results
'''

import os, csv, re, itertools, collections
from collections import defaultdict


Import('*')

try:
    env = env_processing_stats.Clone()
except NameError, ne:    
    vars = Variables('vars.py')

results = env['RESULTS']
quantify_results = env['QUANTIFY_RESULTS']
samples = results

read_quality_stats_dir = 'read_quality_stats'
## COLLECT TRIMMING STATISTICS
env_summarize_depth = env.Clone()
trim_results = [results[s]['trimmed'] for s in samples.keys()]
trimming_report_file = '.*_cutadapt\.log'
if env['TRIMMER'] == 'fastxtoolkit':
    trimming_report_file = ".*_trimming_report\.txt"

env_summarize_depth['TRIMMING_REPORT_FILES'] = get_matching_nodes(trim_results, 
                                                                  trimming_report_file)

trimming_reports = SConscript(os.path.join(read_quality_stats_dir,
                                           'summarize_depth.py'),
                              src_dir = env['SCONSCRIPT_HOME'],
                              variant_dir = read_quality_stats_dir, 
                              duplicate = 0,
                              exports  = 'env_summarize_depth')

## COLLECT FILTERING STATISTICS
filter_stats_files = []
for sample in sorted(samples.keys()):
	filter_stats_files.append(results[sample]['filtered'][1])

collect_filter_stats_cmd = '''grep -H . ${SOURCES} > ${TARGET}'''
collect_filter_stats = env.Command(os.path.join(read_quality_stats_dir, "filter_reports.txt"), 
				   [filter_stats_files], 
				   collect_filter_stats_cmd)

## COLLECT QUALITY READS STATISTICS RESULTS
env_read_quality_stats = env.Clone()
env_read_quality_stats['SAMPLES'] = samples.keys()

read_quality_stats_nodes = [results[n]['quality'] for n in results.keys()]
env_read_quality_stats['FASTQC_HTMLS'] = get_matching_nodes(read_quality_stats_nodes, '.*_fastqc\.html')

read_quality_stats = SConscript(os.path.join(read_quality_stats_dir,
                                             'read_quality_collect.py'),
                                src_dir = env['SCONSCRIPT_HOME'],
                                variant_dir = read_quality_stats_dir, 
                                duplicate = 0,
                                exports = 'env_read_quality_stats')

Clean('.', read_quality_stats_dir)

mirna_pre_align_list = []
mirna_seqtag_discarded_list = []
mirna_multi_gen_map_filter_list = []
mirna_pre_ass_iss_list = []
for sample,qresults in quantify_results.items():
    ## MIRNA PRECURSOR ALIGNMENTS
    mirna_pre_align_list.append(sample + '="' + qresults['HAIRPIN_ALIGNMENTS'][1].abspath + '"')
    ## MIRNA READ FILTER BY NUMBER OF READS PER SEQUENCE TAG 
    mirna_seqtag_discarded_list.append(sample + '="' + \
                                       qresults['HAIRPIN_ALIGNMENTS_BLOB'][1].abspath + '"')
    ## MIRNA FILTER BY MULTIPLE GENOMIC MAPPINGS
    mirna_multi_gen_map_filter_list.append(sample + '="' + \
                                           qresults['HAIRPIN_ALIGNMENTS_BLOB'][2].abspath + '"')
    ## PRECURSOR ASSIGNMENT ISSUES
    mirna_pre_ass_iss_list.append(sample + '="' + qresults['TABLES'][2].abspath + '"')

mirna_pre_align = 'c(' + ','.join(mirna_pre_align_list) + ')'
mirna_seqtag_discarded = 'c(' + ','.join(mirna_seqtag_discarded_list) + ')'
mirna_multi_gen_map_filter = 'c(' + ','.join(mirna_multi_gen_map_filter_list) + ')'
mirna_pre_ass_iss= 'c(' + ','.join(mirna_pre_ass_iss_list) + ')'

non_mirna_align = 'NULL'
if env['NON_MIRNAS']:
    ## NON-MIRNA ALIGNMENT STATS
    non_mirna_align_list = []
    for sample in results.keys():
        non_mirna_align_list.append(sample + '="' +\
                                    results[sample]['non_mirna']['ALIGNMENTS'][1].abspath +\
                                    '"')
    
    non_mirna_align =  'c(' + ','.join(non_mirna_align_list) + ')' 

### generate HTML result report ##
##report_summary_cmd = '''Rscript -e 'results.dir <- dirname("$TARGET.abspath"); '''\
#report_summary_cmd = '''echo 'results.dir <- dirname("${TARGETS[1].abspath}"); '''\
#		'''meta.file <- "${SOURCES[0].abspath}"; '''\
#		'''trimming.file <- "${SOURCES[1].abspath}"; '''\
#		'''filtering.file <- "${SOURCES[2].abspath}"; '''\
#        '''pre.align.files <- ''' + mirna_pre_align + '''; '''\
#        '''seqtag.files <- ''' + mirna_seqtag_discarded + '''; '''\
#        '''multi.genmap.files <- ''' + mirna_multi_gen_map_filter + '''; '''\
#        '''pre.ass.iss.files <- ''' + mirna_pre_ass_iss + '''; '''\
#        '''non.mirna.align.files <- ''' + non_mirna_align + '''; '''\
#		'''rmarkdown::render(input = "$MIRANDMORE_BIN/processing_report.Rmd",'''\
#		'''output_file = "${TARGETS[1].abspath}", '''\
#		'''intermediates_dir = dirname("${TARGETS[1].abspath}") )' > ${TARGETS[0]}'''
#
#report_summary = env.Command(['processing_report_cmd.R',
#                              'processing_report.html'], 
#            			     [env['META'], trimming_reports[0], collect_filter_stats,
#                              [quantify_results[s]['HAIRPIN_ALIGNMENTS'][1] for s in quantify_results.keys()], 
#                              [quantify_results[s]['TABLES'][3] for s in quantify_results.keys()],
#                              [results[s]['non_mirna']['ALIGNMENTS'][1] for s in results.keys()]], 
#            			     [report_summary_cmd + ' && Rscript ${TARGETS[0]}'])
#
#Return('trimming_reports collect_filter_stats read_quality_stats report_summary')
Return('trimming_reports collect_filter_stats read_quality_stats')

