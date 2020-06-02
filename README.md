MiR&moRe2
=========

A bioinformatics software pipeline to characterize microRNA and microRNA-offset RNAs from RNA-seq data

Installation
============

MiR&moRe2 runs on Linux platforms (tested on Ubuntu 16.04 LTS). It requires Python3, R, Java and Perl to be installed in your workstation. Other tools will be installed during miR&moRe2 installation.  
N.B: to install some dependency program you need the required tools in your system, f.i. 'make' and a C compiler.  
In a fresh install Ubuntu Server 16.04 LTS you might need to install Java and the following packages:  

    apt install python3-dev build-essential zlib1g-dev python3-distutil libbz2-dev liblzma-dev libncurses-dev r-base-core libssl-dev libcurl4-openssl-dev

Clone the miR&moRe2 GIT repository from [GitHub][mnm_link] or download and extract the package: 

    git clone git@github.com:egaffo/mirandmore2.git

Enter the miR&moRe2 directory and run the installation script:

    cd mirandmore  
    ./install_mnm2  

How to use
==========

Prepare your project directory and files by downloading:

  - [miRBase][mirbase_link] annotation in GFF3 format.
  - reference genome in FASTA

Specify the samples to analyze
------------------------------

A semicolon separated values file must be created to declare the sample names and the relative FASTQ read files. It must contain 2 colons named 'sample' and 'file'. An optional column 'condition' can be included. Write one row for each sample.  
Example 'meta.csv':

    sample;file;condition  
    MM_115;data/MM_115_small.fastq;A  
    MM_115_replica1;data/MM_115_smallx2.fastq;B  
    MM_115_replica2;data/MM_115_smallx2.fastq;B  

Specify and save project parameters
-----------------------------------

Write a 'vars.py' file, which must be located in the project directory (i.e. the directory in which you will run miR&moRe). In it you have to declare parameters of non default values, as Python variables.
Example 'vars.py': 

    CPUS = '8'  
    GFF_FILE = 'path/to/mirbase_hsa.gff3'  
    GENOME = 'path/to/hg38.fa'  
    ADAPTER = 'TCGTATGCCGTCTTCTGCTTGT'  
    META = 'meta.csv'  
    
Run the analysis
----------------

In you project directory then you'll have two files:
    
    vars.py
    meta.csv

Enter your project directory and execute miR&moRe

    cd myprj
    /path/to/mirnadmore/mirandmore 


Speed up your analysis
======================

Use a precompiled genome index
------------------------------

MiR&moRe2 automatically builds the genome index for Bowtie when only the FASTA genome is given. Building the index will take hours, but you can skip this step by providing a previously built index. 
For instance, you can download pre-built Bowtie indexes from the [Bowtie website][bowtie_index_link]. Then you'll have to set the BOWTIE_GENOME_INDEX parameter in the vars.py file. For instance:

    BOWTIE_GENOME_INDEX = '/home/enrico/annotation/indexes/GCA_000001405.15_GRCh38_no_alt_analysis_set'

N.B: set the reference genome parameter with the FASTA sequence used to build the genome index.      

Already adapter trimmed input reads
-----------------------------------

You can skip the adapter trimming step by setting the 

    NO_ADAPTER = 'True' 
    
in the vars.py file.  

Parallel tasks execution
------------------------

To speed up your analysis, you can tell miR&moRe2 to run multiple tasks in parallel by setting the '-j N\_JOBS' parameter when launching mirandmore. _N\_JOBS_ is the number of parallel jobs to execute.  Example:

    /path/to/mirnadmore/mirandmore "-j 2"
    
__Caveat:__ many tasks of miR&moRe2 use multiprocess computing. When setting the parallel task option you may have these tasks running at the same time and you'll have to balance the CPUS parameter and the -j option in order not to overload your machine. For instance, if you want each parallel task to use 8 cores, then you'll set CPU = '8' in the vars.py file, and if you want 4 tasks to run in parallel you'll set '-j 4' when _calling mirandmore_. In this case you'll need 32 cores in your machine as it may happen that all 4 tasks using 8 cores each will run in parallel.

Output
==================

The 'results' directory will be generated after running miR&moRe2. In it you'll find:

  - _raw\_data.txt_ table with mature small RNA read count (rows) in each sample (columns)  
  - _raw_variants.txt_  table with small RNA isoform read count (rows) in each sample (columns), with sequence, position in the extended precursor and type of isoform
  - _mnm2\_sRNA.gff3_  genome annotation of the detected small RNAs
  - _mnm2\_sRNA.sequences_  sequences of the detected small RNAs in tabulat format
  - _read\_quality_stats/depth.csv_  table with the input raw read amount and number of reads passing adapter and min length filter
  - _read\_quality_stats/filter\_reports.txt_  report of other read filters
  - _read\_quality_stats/samplename\_fastqc.html_  HTML file to inspect each sample input read quality as obtained with [FASTQC][fastqc_link]  

Parameters
================

The following variables can be set in the vars.py file  

META: 'The metadata table file where you specify the project samples, conditions, etc..
    default: meta.csv

CPUS: The number of cores to be used by multithreaded tools
    default: 1

ADAPTER: Read adapter to trim. The adapter sequence itself. Only reads bearing the adapter will be kept to search the small RNAs.
    default: TGGAATTCTCGGGTGCCAAGG

GENOME: genome in FASTA
    default: hg38.fa

GFF_FILE: The miRNA annotation file in GFF3
    default: hsa.gff3

BOWTIE_GENOME_INDEX: The path and name of an already computed Bowtie genome index. Setting this parameter will skip the (long) generation of the index
    default: 

PREDICT_NOVEL_PREMIRS: Wether to enable novel precursor prediction with miRDeep2
    default: False

MRD_SCORE_THRESHOLD: score threshold for miRDeep2 precursors
    default: 50
    
MAX_LEN_FILTER: maximum length allowed for reads
    default: 30

MEAN_QUAL_FILTER: mean Phred read quality. Mean across all base qualities of the read.
    default: 30

MIN_COUNT: minimum number of unique reads
    default: 10

NOADAPTER: Set "True" if adapters were already trimmed in  input reads. This will skip adapter trimming step.
    default: False

CHRM_PREFIX_TO_REMOVE: Pattern string to remove from annotation chromosomes. Needed if annotation and genome chromosome names differ (e.g. set 'chr' with Ensembl genomes). Leave '' if no changes needed.
    default: ''

MULTIPLE_GENOMIC_HITS_THRESHOLD: multiple genomic hits threshold
    default: 5

TRIMMER: The tool to use for trimming read adapter. Only cutadapt is currently supported
    default: cutadapt

CUTADAPT_EXTRA_PARAMS: additional parameters to be passed to Cutadapt. For instance use '-u 3' for SMARTER small RNA-seq library.
    default: 

MIN_MORNA_LEN: minimun lenght for moRNAs
    default: 18

EXTENDED_N: Number of positions for extended hairpins at each end
    default: 30

MORFILTER: Whether to account (permissive) or not (conservative) sequences that align also to miRNAs to estimate moRNA expression. Default is to account sequences aligned only to moRNAs to prevent miRNA sequences to increase moRNA read count.   
    default: conservative

ALLOWED_OVERHANG: The number of bases allowed to overhang the small RNA coordinates, either at 3' and 5'. This parameter controls the isoform length assigned to each small RNA. For instance, if miR-xx-3p starts at base 10 and ends at base 30 on its precursor, reads aligned from base 7 and/or aligned up to base 33 will be accounted, by default, as miR-xx-3p isomiRs. Conversely, reads aligned before base 7 and/or after base 33 will not increase miR-xx-3p read count.  
    default: 3

NON_MIRNAS: (Experimental) Enable characterization of reads not aligned to miRNA precursors
    default: False

GENE_ANNOTATION: (Experimental) Gene annotation GTF file to characterize non-miRNA reads. Leave empty to only map length discarded reads
    default: 

NM_BOWTIE_PARAMS: (Experimental) Parameters to be passed to Bowtie to align non-miRNA reads.
    default: -n 0 -l 26 -e 70 --best --strata -k 5 -m 5 --mm

[mnm_link]: https://github.com/egaffo/mirandmore2
[mirbase_link]: ftp://mirbase.org/pub/mirbase/CURRENT/genomes/
[bowtie_index_link]: ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/
[fastqc_link]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
