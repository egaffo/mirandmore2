#require(DESeq2)
require(plyr)

a <- function(name,value) {
  assign(name,value,envir=.GlobalEnv)
}

load_file <- function(name,query="\\.txt",replace="") {
    df_name <- basename(sub(query,replace,name))
    if (inherits( try( get(df_name,envir=.GlobalEnv),T ), "try-error"))
        a(df_name,read.table(name,sep=";",quote="",header=TRUE,stringsAsFactors=F))
}

pre.process <- function(df) {
  df <- arrange(df,name)
  df <- ddply(df,"name",summarize,expression=round(sum(expression)))
  df
}

pre.process.flag <- function(df) {
    df <- arrange(df,name)
    df <- ddply(df,"name",summarize,new = new[1])
    df 
}


int.generator <- function() {
  start = 0
  return(function(){
    start <<- start + 1
    start
  })}

colgenerator <- function(cols) {
  #cols <- c("red","blue","green")
  i <- 0
  return(function(){
    i <<- i + 1
    return (cols[i])
  })}


na.to.zero <- function(v) {
  v[is.na(v)] <- 0
  v
}


round.vec <- function(vec) {
  if (class(vec) == "numeric") {
    round(vec)
  } else {
    vec 
  }    
}    

merge.variants <- function(dfl) {
  it <- int.generator()
  keys <- c("pre","mature","seq","start","end","category","type")
  data <- Reduce(function(x,y) {
    merge(x,y,all=T,by=keys,suffixes=c(it(),it()))
  },dfl)
  to.keep <- c(keys,grep("count",names(data),value=T))
  data<-subset(data,T,to.keep) 
  names(data) <- c(keys,names(dfl))
  data <- lapply(data,na.to.zero)
  data <- lapply(data,round.vec)
  names(data) <- sub("_unfiltered_variants","",names(data))
    
  as.data.frame(data,stringsAsFactors=F)
}

merge.counts <- function(dfl,mangler=NULL) {
  it <- int.generator()
  data <- Reduce(function(x,y) {
    merge(x,y,all=T,by=c("name"),suffixes=c(it(),it()))
  },dfl)
  n <- length(colnames(data))
  if (is.null(mangler))
    mangler <- function(x)x
  colnames(data)[2:n] <- mangler(names(dfl))
  data.no.na <- lapply(data[-1],na.to.zero)
  data <- cbind(data[1],data.no.na)  
  data
}

#normalize <- function(cds) {
  #t( t(counts(cds) / sizeFactors(cds) ))
#}

mir.table.name.mangler <- function(vec) {
  sub("_mir.*","",vec)
}

arrange.and.sort.adts <- function(names) {
  l_ply(names,function(x){
    df = get(x,envir=.GlobalEnv)
    df = transform(df,n=nchar(seq))
    df = arrange(df,desc(n),desc(freq))
    assign(x,df,envir=.GlobalEnv)
  })
}

aggregate.samples.by.f <- function(sample.names,aggregator) {
  tablenames <- basename(sub("\\.txt","",sample.names))
  l_ply(sample.names,load_file)
  mir.tables <- mget(tablenames,envir=.GlobalEnv)
  mir.tables.def <- lapply(mir.tables,aggregator)
  count.data <- merge.counts(mir.tables.def,mir.table.name.mangler)
  count.data
}

prepare.count.data <- function(sample.names) {
    aggregate.samples.by.f(sample.names,pre.process)
}


prepare.count.data.flag <- function(sample.names) {
    aggregate.samples.by.f(sample.names,pre.process.flag)
}


prepare.flag.vector <- function(sample.names) {
    data <- prepare.count.data.flag(sample.names)
    mat <- as.matrix(data[,-1])
    apply(mat,1,sum) > 0
}


prepare.variants <- function(unfiltered_variants_filepaths) {
    #sample.names <- sub("_mir_table.txt","",sample.names) 
    #filenames <-  paste(sample.names, "_unfiltered_variants.txt", sep="")
    filenames <- unlist(strsplit(unfiltered_variants_filepaths, " ")) 
    tablenames <- basename(sub("\\.txt","", filenames))
    l_ply(filenames,load_file)
    variants.dfl = mget(tablenames, envir=.GlobalEnv)
    variants.data <- merge.variants(variants.dfl)
    variants.data
}  


#to.cds <- function(data,conditions) {
#  count.mat <- as.matrix(data)
#  cds <- newCountDataSet(count.mat,conditions)
#  cds <- estimateSizeFactors(cds)
#  cds
#}
#
#count.data.to.cds  <- function(count.data, conditions, last) {
#    i <- match(last, colnames(count.data))
#    meta.data.names <- colnames(count.data)[1:i]
#    meta.data <- as.data.frame(count.data[, 1:i], stringsAsFactors = F)
#    colnames(meta.data) <- meta.data.names
#    data <- count.data[, -c(1:i), drop = FALSE]
#    if (last == "name")
#        rownames(data) <- count.data[, "name"]
#    cds <- to.cds(data, conditions)
#    list(cds = cds, meta.data = meta.data)
#}
#
#normalize.count.data <- function(obj) {
#    normalized.data <- counts(obj$cds,normalized=T)
#    cbind(obj$meta.data,as.data.frame(normalized.data,stringsAsFactors=F))
#}

filter.out.mir.tablef <- function(df) {
  return (df[,c("name","expression")])
}

