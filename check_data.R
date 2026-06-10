# 检查数据分布
library(dplyr)

data <- read.table("C:/Users/Jy-Mentor-7/Downloads/GSE97537.top.table (1).tsv", 
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("数据维度:", nrow(data), "行\n\n")

cat("=== logFC 分布 ===\n")
print(summary(data$logFC))

cat("\n=== logFC 范围 ===\n")
cat("Min:", min(data$logFC, na.rm = TRUE), "\n")
cat("Max:", max(data$logFC, na.rm = TRUE), "\n")

cat("\n=== 正负分布 ===\n")
cat("logFC > 0:", sum(data$logFC > 0, na.rm = TRUE), "\n")
cat("logFC < 0:", sum(data$logFC < 0, na.rm = TRUE), "\n")
cat("logFC = 0:", sum(data$logFC == 0, na.rm = TRUE), "\n")

cat("\n=== 不同阈值下的显著基因数 ===\n")
for(th in c(0.5, 1, 1.5, 2)) {
  up <- sum(data$logFC > th & data$adj.P.Val < 0.05, na.rm = TRUE)
  down <- sum(data$logFC < -th & data$adj.P.Val < 0.05, na.rm = TRUE)
  cat(sprintf("|logFC| > %.1f: 上调=%d, 下调=%d\n", th, up, down))
}

cat("\n=== Padj 分布 ===\n")
cat("Padj < 0.05:", sum(data$adj.P.Val < 0.05, na.rm = TRUE), "\n")
cat("Padj < 0.01:", sum(data$adj.P.Val < 0.01, na.rm = TRUE), "\n")

cat("\n=== 前10个最大logFC ===\n")
print(head(data[order(-data$logFC), c("Gene.symbol", "logFC", "adj.P.Val")], 10))

cat("\n=== 前10个最小logFC ===\n")
print(head(data[order(data$logFC), c("Gene.symbol", "logFC", "adj.P.Val")], 10))
