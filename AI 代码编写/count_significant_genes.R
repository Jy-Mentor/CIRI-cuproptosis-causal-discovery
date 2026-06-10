# 读取GSE61616.top.table (1).tsv文件
data <- read.table("GSE61616.top.table (1).tsv", header = TRUE, sep = "\t")

# 统计显著差异基因数量（FC≥0.2 及 adj.P.Val < 0.05）
# 注意：logFC的绝对值≥0.2对应的FC≈1.15（因为log2(1.15)≈0.2）
significant_genes <- data[abs(data$logFC) >= 0.2 & data$adj.P.Val < 0.05, ]
count <- nrow(significant_genes)

# 输出结果
cat("显著差异基因数量（FC≥0.2 及 校正P<0.05）：", count, "\n")

# 输出上调和下调基因数量
up_regulated <- significant_genes[significant_genes$logFC > 0, ]
down_regulated <- significant_genes[significant_genes$logFC < 0, ]
cat("上调基因数量：", nrow(up_regulated), "\n")
cat("下调基因数量：", nrow(down_regulated), "\n")

# 输出前10个显著差异基因
if (count > 0) {
  cat("\n前10个显著差异基因：\n")
  print(head(significant_genes[, c("Gene.symbol", "adj.P.Val", "logFC")], 10))
}
