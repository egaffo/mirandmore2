#!/usr/bin/env Rscript

library(getopt)
library(plyr)

opt_spec <- matrix(ncol = 4, byrow = TRUE, 
                   c('workingdir', 'w', 1, 'character'))

opt = getopt(opt_spec)

data <- tryCatch(read.table(file.path(opt$workingdir, "means_per_condition.txt"), 
                            sep = ";", header = T, stringsAsFactors = F, quote = ""), 
                 error = function(e) setNames(data.frame(matrix(ncol = 2, nrow = 0)), 
                                            c("name", "new")))

seqs <- tryCatch(read.table(file.path(opt$workingdir, "new-srna-seqs.csv"),
                            sep = ",", header = T, stringsAsFactors = F, quote = ""), 
                 error = function(e) setNames(data.frame(matrix(ncol = 2, nrow = 0)), 
                                            c("name", "seq")))

annotations <- tryCatch(read.table(file.path(opt$workingdir, "new-mir-table-without-seqs.txt"),
                                   sep = "\t", header = T, stringsAsFactors = F, quote = ""), 
                        error = function(e) setNames(data.frame(matrix(ncol = 4, nrow = 0)), 
                                                   c("pre", "name", "strand", "position")))

tmp = merge(annotations, seqs, by = "name")
tmp = merge(tmp, data, by = "name")

write.table(tmp, file = file.path(opt$workingdir, "new-srna-table-with-seq.txt"),
            sep = ";", row.names = F, quote = F)

mirnas <- subset(tmp, !grepl("moR|loop", name, perl = T))
write.table(tmp, file = file.path(opt$workingdir, "new-mir-table-with-seq.txt"),
            sep = ";", row.names = F, quote = F)

mornas <- subset(tmp, grepl("moR", name, perl = T))
write.table(mornas, file = file.path(opt$workingdir, "mor-table-with-seq.txt"),
            sep = ";", row.names = F, quote = F)
