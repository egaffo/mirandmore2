#!/usr/bin/env python

import sys, re, argparse

## Parse a GFF3 file and produce a table of mature miRNAs with miRNA coordinates within their
## precursor.
## * target will be the 'mature-table.txt' file (or the target name given)
## * source must be a valid GFF3 file, as provided by miRBase (e.g: hsa.gff3)
## TODO: the older GTF (a.k.a. GFF2, s.a. dme.gtf2) files only specify precursor's genomic
## coordinates. In such case you need to download and process also the mirna.dat file 
## (which is in Embl format) to extract the mature miRs' coordinates within their precursor.
def gff2maturetable(target, source, env):
    ext = str(source[0]).rpartition('.')[-1]
    pattern = 'chr'
    try:
        pattern = env['CHRM_PREFIX']
    except KeyError, e:
        sys.stderr('''WARNING: chromosome prefix not specified. '''\
                   '''Will remove default pattern '{}' from chromosome names'''.format(pattern))
    if ext.lower() == 'gff3':
        pre_id2names = {}
        table_lines  = []
        precursor_names = []
        with open(str(source[0]), 'r') as gff:
            for line in gff:
                if not line.startswith("#"):
                    items = line.split("\t")
                    chromosome = items[0].replace(pattern,"")
                    entry_type = items[2]
                    start      = int(items[3])
                    end        = int(items[4])
                    strand     = items[6]
                    Id         = re.findall(r'ID=([^;]*);',items[8])[0]
                    name       = re.findall(r'Name=([^;]*)',items[8])[0].strip()
                    

                    if entry_type == 'miRNA_primary_transcript':
                        ## check for duplicated precursor names:
                        ## for duplicated entries modify the name 
                        ## by appending the mir ID
                        if name in precursor_names:
                            sys.stderr.write("WARNING: duplicated precursor"\
                                             "name for " + name + \
                                             ". Will change into " + \
                                             name + "|" + Id + "\n")
                            name = name + "|" + Id
                        else:
                            precursor_names.append(name)

                        pre_id2names[Id] = {'name':name, 
                                            'chromosome,':chromosome, 
                                            'start':start, 
                                            'end':end}
                    elif entry_type == 'miRNA':
                        mirna_pre_id = re.findall(r'Derives_from=([^;]*)',items[8])[0].strip()
                        try:
                            ## we assume that the GFF specifies first the precursor
                            ## and then the relative miRNAs. Or, at least, all the
                            ## precursors first and then all the miRNAs
                            prename      = pre_id2names[mirna_pre_id]['name']
                            pre_start    = pre_id2names[mirna_pre_id]['start']
                            pre_end      = pre_id2names[mirna_pre_id]['end']
                            mir_in_pre_start = start - pre_start + 1
                            mir_in_pre_end   = end - pre_start + 1
                            if strand == '-':
                                mir_in_pre_start   = pre_end - end + 1
                                mir_in_pre_end     = pre_end - start + 1
                            pre_length       = pre_end - pre_start + 1
                            table_lines.append([chromosome,
                                                prename,
                                                name,
                                                strand,
                                                str(mir_in_pre_start),
                                                str(mir_in_pre_end),
                                                str(pre_length)]
                                              )
                        except KeyError, e:
                            sys.stderr.write('''WARNING: precursor ID {0} not found for miR {1}.'''\
                                       ''' Be sure that miRNA_primary_transcript feature '''\
                                       '''preceeds relative miRNA feature. '''\
                                       '''Try to sort by coordinates the GFF file; or '''\
                                       '''put all the precursor entries first.\n'''.\
                                       format(mirna_pre_id, name))
                            #continue

        with open(str(target[0]), 'w') as out:
            for line in table_lines:
                out.write('\t'.join(line) + '\n')
                
    else:
        sys.exit('''Non GFF v3 format not yet supported. Please use a valid GFF3 file.'''\
                 '''If you are sure a valid GFF3 has been passed try to rename with .gff3''')

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g','--gff', dest = "gff", required = True, 
                        help = '''A GFF v3 file with miRNA precursors and mature miRNAs '''\
                        '''annotation (features 'miRNA_primary_transcript' and 'miRNA'). '''\
                        '''The GFF must by sorted so that precursors always '''\
                        '''appear before the relative miRNAs (list all precursors '''\
                        '''first and then all miRNAs is allowed)''')
    parser.add_argument('-m','--mature_table', dest = "mature_table", 
                        required = True, help = '''The output filename''', 
                        default = 'mature-table.txt')
    parser.add_argument('-c','--chrm_prefix', dest = "chrmprefix", required = False,
                        help = '''A string pattern that will be removed from the chromosome '''\
                               '''names, e.g. 'chr' will transform 'chrX' into 'X' '''\
                               '''This is useful when you use a genome and annotation with '''\
                               '''different standards, like Ensembl genome and miRBase'''\
                               '''annotation''', 
                        default = 'chr')
    args = parser.parse_args()

    env = {'CHRM_PREFIX':args.chrmprefix}

    gff2maturetable([args.mature_table], [args.gff], env)


