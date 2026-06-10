#!/usr/bin/env Rscript

# 测试读取10x数据文件
file_path <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE/GSM5319987_sham1_barcodes.tsv"

cat("尝试读取文件: ", file_path, "\n")

# 尝试直接读取
tryCatch({
  content <- readLines(file_path, n = 10)
  cat("成功读取前10行:\n")
  print(content)
}, error = function(e) {
  cat("读取失败: ", e$message, "\n")
})

# 尝试使用Read10X
cat("\n尝试使用Read10X读取整个样本:\n")
tryCatch({
  library(Seurat)
  sample_path <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE/GSM5319987_sham1"
  counts <- Read10X(sample_path)
  cat("成功读取，基因数: ", nrow(counts), " 细胞数: ", ncol(counts), "\n")
}, error = function(e) {
  cat("Read10X失败: ", e$message, "\n")
})
