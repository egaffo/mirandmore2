#!/usr/bin/env Rscript

library(getopt)
library(plyr)

opt_spec <- matrix(ncol=4, byrow=TRUE, c(
    'meta',  'm', 1, 'character',
    #'species','s',1,'character',
    'annotation_gff3', 'a', 1, 'character',
    'extended_n','n',1,'integer',
    'input','i',1,'character',
    'destination','d',1,'character'))

pre.process.positions <- function(df){
    df[,c("pre","name","new","start","end")]
}

merge.positions <- function(dfl) {
    it <- int.generator()
    data <- Reduce(function(x,y) {
        merge(x,y,all=T,by=c("pre","name"),suffixes=c(it(),it()))
    }, dfl)
    data
}

transform.coords  <- function(df) {
    abs.start <-ifelse(df$strand=="+",
                       df$start + df$start.pre - EXTENDED_N,
                       df$end.pre + EXTENDED_N  - df$end +1)
    abs.end <- ifelse(df$strand=="+",
                      df$end + df$start.pre - EXTENDED_N-1,
                      df$end.pre + EXTENDED_N  - df$start)
    df$position = paste(df$chromosome,":",abs.start,"-",abs.end,sep="")
    cbind(df,abs.start,abs.end)
}


opt = getopt(opt_spec)


if (is.null(opt$meta)) {
    stop("please specify a meta file\n")
}

source(file.path(Sys.getenv("MIRANDMORE_HOME"),"lib","collect_lib.R"))
meta <- read.table(opt$meta,stringsAsFactors=F,sep=";",header=T)
#SPECIES <- opt$species
annotation_gff3 <- opt$annotation_gff3
EXTENDED_N <- opt$extended_n
#source(file.path(Sys.getenv("MIRANDMORE_HOME"),"etc","config.R"))
#meta = read.csv("meta.csv",sep=";")
#sample.names <- meta$sample
sample.names <- unlist(strsplit(opt$input," ")) 
#filenames <- file.path(sample.names,  paste( sample.names, "_mir_table.txt", sep=""))
filenames <- sample.names
l_ply(filenames,load_file)
ll <- mget(paste(meta$sample,"_mir_table",sep="") ,envir=.GlobalEnv)
ll <- lapply(ll, pre.process.positions)
data <- merge.positions(ll)
starts <- as.matrix(data[,  grepl("start",colnames(data))])
actual.starts <-  apply(starts,1, min, na.rm=T)
ends <-   as.matrix(data[,  grepl("end",colnames(data))])
new_ <- as.matrix(data[,grep("new",colnames(data))])
actual.ends <-  apply(ends,1, max, na.rm=T)
new_ <- apply(new_, 1, all,na.rm=T)
new.data <-  data.frame(data[,c("pre","name")], new=new_, start=actual.starts-1, end=actual.ends-1)


#gff <- file.path(Sys.getenv("MIRANDMORE_HOME"),"annotations","mirbase", paste(SPECIES,".gff3",sep=""))
gff <- file.path(annotation_gff3)
the.gff <- read.table(gff,sep="\t",stringsAsFactors=F) 
the.gff <- the.gff[,c("V1","V4","V5","V7","V9")] 
colnames(the.gff) <- c("chromosome","start","end","strand","pre")
the.gff$pre <- sub(".*?Name=([^;]*).*","\\1", the.gff$pre,perl=T)

new.data.bis <- subset(new.data, new==TRUE)[, c("pre", "start", "end", "name")]
if(nrow(new.data.bis) > 0){
    new.data.bis <- merge(new.data.bis, the.gff, by="pre", suffixes=c("",".pre"))
    new.data.bis <- transform.coords(new.data.bis)
}else{
    new.data.bis <- data.frame(pre = character(),
                               name = character(),
                               strand = character(),
                               position = character())
}
new.data.bis <- new.data.bis[,c("pre","name","strand","position")]

write.table(new.data.bis,file.path(opt$destination,"new-mir-table-without-seqs.txt"),quote=F,sep="\t",row.names=F,col.names=T)

new.data <- merge(new.data, the.gff[, c("pre", "strand")], by="pre", suffixes=c("",".pre"))
new.data$pre <- sapply(new.data$pre, function(x){paste0(x, "-ext")})
new.data <- subset(new.data, new==TRUE)[, c("pre", "start", "end", "name", "strand")]
write.table(new.data, file.path(opt$destination,"new_mir_squished_coords.bed"), quote=F, sep="\t", row.names=F, col.names=F)
