#!/usr/bin/env Rscript

library(getopt)
library(plyr)

opt_spec <- matrix(ncol=4, byrow=TRUE, 
                   c('meta',  'm', 1, 'character',
                     'data','d',1,'character',
                     'output','o',1,'character'))

opt = getopt(opt_spec)

meta <- read.table(opt$meta,sep=";",header=T,stringsAsFactors=F)
data <- read.table(opt$data,sep=";",header=T,stringsAsFactors=F,quote="",check.names=F)

meta.names <- c("name","new")
sample.names <- meta[,"sample"]
mat <- as.matrix(data[,sample.names])
means <- apply(mat,1,function(x) tapply(x,meta$condition, mean))

if(length(unique(meta$condition)) > 1){
    means <- t(means)
}

means.df <- cbind(data[,meta.names],means)
write.table(means.df,file=opt$output,sep=";",quote=F,row.names=F)
