#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

#option_list <- list(
#  make_option(c("-o", "--output"), action="store", type="character",
#              default="ld_counts_", help="Output files' prefix"),
#  make_option(c("-g", "--gtf_list"), action="store", type="character", default="",
#              help=paste("The file with the set of miRNAs to output. Must be a text ",
#                         "file with a line for each miRNA. Default prints all ",
#                         "miRNAs.", sep=""))#,
#)
#parser <- OptionParser(usage="%prog [options] ",
#                       option_list=option_list,
#                       description="")
#arguments <- parse_args(parser, positional_arguments=TRUE)
#opt <- arguments$options
#
#if(length(arguments$args) != 1) {
#  cat("Incorrect number of required positional arguments\n\n")
#  print_help(parser)
#  stop()
#} else {
#  
#  unfvar      <- arguments$args[1]
#  gtf.list     <- opt$gtf_list
#  output      <- opt$output

## input files and output names
#files <- c('CR_bio3' = 'CR_bio3.counts.gtf', 'CR_bio2' = 'CR_bio2.counts.gtf', ...)
#merged.gtf.counts.file <- "merged.counts.gtf.csv"
#count.matrix.file <- "counts.csv"

counts <- rbindlist(sapply(files, FUN = fread, 
			   USE.NAMES=T, simplify=F),use.names=T, idcol="sample")
counts[, sample_id := gsub(".*samples/([^/]*)/.*", "\\1", sample)]
counts[, gene_id := gsub('.*gene_id "([^"]*)".*', "\\1", V9)]
counts[, gene_name := gsub('.*gene_name "([^"]*)".*', "\\1", V9)]

write.csv(counts[V10 > 0], merged.gtf.counts.file, row.names = F)

samples <- unique(counts$sample_id)

counts.data.matrix <- dcast(counts[, .(gene_id, gene_name, sample_id, V10,
                                       chr = V1, start = V4, end = V5, strand = V7)], 
                            formula = gene_id + gene_name + chr + start + end+strand ~ sample_id, 
                            value.var="V10", 
                            fill=0)
counts.data.matrix <- counts.data.matrix[gene_id %in% counts[V10 > 0, as.character(unique(gene_id))]]

write.csv(counts.data.matrix[order(rowMeans(counts.data.matrix[, samples, with=F]), 
                                   decreasing=T),], 
	  count.matrix.file, row.names = F)

#}


