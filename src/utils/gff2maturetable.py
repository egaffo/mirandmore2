#!/usr/bin/env python

import sys, re, argparse
from collections import defaultdict

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
    except KeyError as e:
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
                        except KeyError as e:
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



def get_field_value(field_list, field):
    
    val = None
    for f in field_list:
        if field in f:
            val = re.sub('.*=(.*)', r'\1', f)
            break

    return val


def maturetable2gff(mature_table, gff, prefix):
    
    precursors_pos = defaultdict(dict)

    with open(gff, 'r') as gff_features:
        for feature in gff_features:
            
            ## populate reference GFF features
            f = feature.strip().split('\t')
            chrom = f[0].strip()
            start = f[3].strip()
            end   = f[4].strip()
            strand= f[6].strip()
            #name  = re.sub('.*Name=([^;]);*', r'\1', f[8]).strip()
            name  = get_field_value(f[8].strip().split(';'), 'Name')
            ID    = get_field_value(f[8].strip().split(';'), 'ID')

            precursors_pos[name] = {'chrom': chrom,
                                    'start': start,
                                    'end':   end,
                                    'strand':strand,
                                    'ID':   ID,
                                    'gffline_list': f}

    with open(mature_table, 'r') as mt:
        for line in mt:
            l = line.strip().split('\t')
            #chrom = l[0].strip()
            start = l[4].strip()
            end   = l[5].strip()
            #strand= l[3].strip()
            name  = l[2].strip()
            pre   = l[1].strip()

            if name in precursors_pos.keys():
                outline_list = precursors_pos[name]['gffline_list']

            else:
                gen_start = int(precursors_pos[pre]['start']) + int(start) - 1
                gen_end   = int(precursors_pos[pre]['start']) + int(end) - 1
                gen_chrom = re.sub(prefix, '', precursors_pos[pre]['chrom'])
                gen_strand= precursors_pos[pre]['strand']
                gen_ID    = precursors_pos[pre]['ID']
                featuretype = 'moRNA' if 'moR' in name else 'miRNA_loop' if 'loop' in name else 'miRNA'

                outline_list = [gen_chrom, 
                                '.', 
                                featuretype,
                                str(gen_start),
                                str(gen_end),
                                '.',
                                gen_strand,
                                '.',
                                ';'.join(['ID='+name, 
                                          'Alias='+name,
                                          'Name='+name,
                                          'Derives_from='+gen_ID])]

            print('\t'.join(outline_list))

            

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
                        default = '')
    parser.add_argument('-r', '--reverse',
                        dest = 'reverse', 
                        required = False,
                        action = 'store_true',
                        help = '''Reverse the conversion: from a mature-table.txt '''\
                               '''formatted file return a GFF with the genomic '''\
                               '''coordinates. The original GFF is required '''\
                               '''since it bears the genomic coordinates of '''\
                               '''each precursor. This is useful to get genomic '''\
                               '''coordinates of newly predicte RNAs, such as moRNAs. '''\
                               '''Output will be to stdout''')

    args = parser.parse_args()

    env = {'CHRM_PREFIX':args.chrmprefix}

    if args.reverse:
        maturetable2gff(args.mature_table, args.gff, args.chrmprefix)
    else:
        gff2maturetable([args.mature_table], [args.gff], env)


