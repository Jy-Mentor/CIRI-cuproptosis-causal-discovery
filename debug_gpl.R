gpl <- read.table("D:/反向网络药理学/L1 数据集/bulk/GSE97537(24H)/GPL1355-10794 (1).txt", header=TRUE, sep="\t", quote="", skip=17, comment.char="", stringsAsFactors=FALSE, fill=TRUE, nrows=3)
print(colnames(gpl))
cat("Gene Symbol" %in% colnames(gpl))
