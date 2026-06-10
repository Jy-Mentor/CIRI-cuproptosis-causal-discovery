#!/usr/bin/env Rscript
gpl_file <- "D:/反向网络药理学/L1 数据集/bulk/GSE97537(24H)/GPL1355-10794 (1).txt"
gpl <- read.table(gpl_file, header=TRUE, sep="\t", quote="", 
                  comment.char="", stringsAsFactors=FALSE, fill=TRUE, nrows=5)
cat("总列数:", ncol(gpl), "\n")
for (i in seq_len(ncol(gpl))) {
  cat("  [", i, "] ", colnames(gpl)[i], "\n", sep="")
}
print(gpl[1:3, 1:15])
