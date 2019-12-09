#!/usr/bin/env Rscript

library(getopt)
# library(plyr)
library(data.table)

opt_spec <- matrix(ncol=4, byrow=TRUE, c(
    'meta',  'm', 1, 'character',
    'verbose', 'v' , 1, 'logical',
    'input','i',1,'character',
    'destination','d',1,'character',
    'unfiltered_variants', 'u', '1', 'character'))

opt = getopt(opt_spec)


if (is.null(opt$meta)) {
    stop("please specify a meta file\n")
}

get.rna.type <- function(x){
    rnatype <- ""
    # rnatype <- "known."
    # if(grepl("putative_new", x)) rnatype <- "new.pre."
    if(grepl("moR", x)) rnatype <- paste0(rnatype, "moR")
    if(grepl("loop", x)) rnatype <- paste0(rnatype, "loop")
    if(!grepl("moR|loop", rnatype)) rnatype <- paste0(rnatype, "miR")
    rnatype <- paste0(rnatype,
                      ifelse(grepl("putative_new", x),
                             ".new.pre",
                             ".known"))
    rnatype
}

variants.files <- unlist(strsplit(x = opt$unfiltered_variants, split = " "))
names(variants.files) <-
    sapply(variants.files,
           function(x)sub("_unfiltered_variants.txt", "", basename(x)))

raw.vars.m <- rbindlist(lapply(variants.files, fread),
                        use.names = T,
                        idcol = "sample_id")

raw.vars.m[, rna.type := sapply(X = mature,
                                FUN = get.rna.type,
                                USE.NAMES = F),
           by = mature]

## select/filter only pure moRNAs: discard sequences common to any miR variant
all.mor.seqs <- raw.vars.m[grepl("moR", mature), ]

keep.mor.seqs <-
    raw.vars.m[,.(srna.types = paste0(unique(rna.type),
                                      collapse = "|")),
               by = seq][srna.types %in% c("moR.known", "moR.new.pre"),
                         unique(seq)]

real.mors <- raw.vars.m[seq %in% keep.mor.seqs, ]

## purge moRs from data and then reinsert only the pure moRs
raw.vars.m <- raw.vars.m[!grepl("moR", rna.type), ]
raw.vars.m <- rbindlist(list(raw.vars.m, real.mors), use.names = T)

raw.vars <- dcast(raw.vars.m,
                  pre + mature + seq + start + end + category + type ~ sample_id,
                  value.var = "count",
                  fill = 0)

## save as table
fwrite(x = raw.vars,
       file = file.path(opt$destination, "raw_variants.txt"),
       sep = "\t")

## compute sRNA expression melting the isoforms
raw.pre.m <- raw.vars.m[, .(reads = sum(count)),
                        by = .(sample_id, pre, mature)][reads > 0]
raw.m <- raw.pre.m[, .(reads = as.integer(ceiling(mean(reads)))),
                   by = .(sample_id, name = mature)]

raw.mat.dt <- dcast(raw.m,
                    name ~ sample_id,
                    value.var = "reads",
                    fill = 0)

## set mock column
## TODO: set proper values
# raw.mat.dt$new <- FALSE

## save as table
fwrite(x = raw.mat.dt,
       file = file.path(opt$destination, "raw_data.txt"),
       sep = "\t")



# source(file.path(Sys.getenv("MIRANDMORE_HOME"),"lib","collect_lib.R"))
# meta <- read.table(opt$meta,stringsAsFactors=F,sep=";",header=T)
# sample.names <- unlist(strsplit(opt$input," "))
# new_flag = prepare.flag.vector(sample.names)
# count.data <- prepare.count.data(sample.names)
# count.data <- count.data[,c("name",meta$sample)]
# #count.data.cds.wrapper <- count.data.to.cds(count.data,meta$condition,"name")
# #normalized.data <- normalize.count.data(count.data.cds.wrapper)
# #normalized.data <- cbind(name=normalized.data[,1],new=new_flag,normalized.data[,-1, drop = FALSE])
# count.data <- cbind(name=count.data[,1],new=new_flag,count.data[,-1, drop = FALSE])
# #count.data.cds  <- count.data.cds.wrapper$cds
# #variants.data <-  prepare.variants(sample.names)
# variants.data <-  prepare.variants(opt$unfiltered_variants)
# variants.meta <- c("pre", "mature", "seq", "start", "end", "category", "type")
# colnames(variants.data) <-  gsub("\\.","-",colnames(variants.data))
# variants.data <- variants.data[,c(variants.meta,meta$sample)]
# #variants.cds.wrapper <- count.data.to.cds(variants.data,meta$condition,"type")
# #variants.normalized  <-  normalize.count.data(variants.cds.wrapper)
# #variants.cds <- variants.cds.wrapper$cds
# write.table(count.data,file.path(opt$destination,"raw_data.txt"),sep=";",quote=F,row.names=F)
# #write.table(normalized.data,file.path(opt$destination,"normalized_data.txt"),sep=";",quote=F,row.names=F)
# write.table(variants.data,file.path(opt$destination,"raw_variants.txt"),sep=";",row.names=F,quote=F)
# #write.table(variants.normalized,file.path(opt$destination,"normalized_variants.txt"),sep=";",quote=F,row.names=F)
# #save(list=c("meta","variants.normalized","variants.cds","variants.data","count.data.cds","count.data","normalized.data"),
# #     file="summary.rda")
