#!/usr/bin/env Rscript

library(getopt)
library(plyr)

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

source(file.path(Sys.getenv("MIRANDMORE_HOME"),"lib","collect_lib.R"))
meta <- read.table(opt$meta,stringsAsFactors=F,sep=";",header=T)
sample.names <- unlist(strsplit(opt$input," ")) 
new_flag = prepare.flag.vector(sample.names)
count.data <- prepare.count.data(sample.names)
count.data <- count.data[,c("name",meta$sample)]
count.data.cds.wrapper <- count.data.to.cds(count.data,meta$condition,"name")
normalized.data <- normalize.count.data(count.data.cds.wrapper)
normalized.data <- cbind(name=normalized.data[,1],new=new_flag,normalized.data[,-1, drop = FALSE])
count.data <- cbind(name=count.data[,1],new=new_flag,count.data[,-1, drop = FALSE])
count.data.cds  <- count.data.cds.wrapper$cds
#variants.data <-  prepare.variants(sample.names)
variants.data <-  prepare.variants(opt$unfiltered_variants)
variants.meta <- c("pre", "mature", "seq", "start", "end", "category", "type")
colnames(variants.data) <-  gsub("\\.","-",colnames(variants.data))
variants.data <- variants.data[,c(variants.meta,meta$sample)]
variants.cds.wrapper <- count.data.to.cds(variants.data,meta$condition,"type")
variants.normalized  <-  normalize.count.data(variants.cds.wrapper)
variants.cds <- variants.cds.wrapper$cds
write.table(count.data,file.path(opt$destination,"raw_data.txt"),sep=";",quote=F,row.names=F)
write.table(normalized.data,file.path(opt$destination,"normalized_data.txt"),sep=";",quote=F,row.names=F)
write.table(variants.data,file.path(opt$destination,"raw_variants.txt"),sep=";",row.names=F,quote=F)
write.table(variants.normalized,file.path(opt$destination,"normalized_variants.txt"),sep=";",quote=F,row.names=F)
#save(list=c("meta","variants.normalized","variants.cds","variants.data","count.data.cds","count.data","normalized.data"), 
#     file="summary.rda")
