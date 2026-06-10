# GSE23160 对照组问题诊断
library(limma)
library(dplyr)

setwd("D:/反向网络药理学/L1 数据集/bulk/GSE23160(主验证集时序差异分析，2h,8h,24h)")

# 读取数据
expr <- read.delim("GSE23160_series_matrix.txt.gz", 
                   comment.char = "!", 
                   header = TRUE, 
                   check.names = FALSE, 
                   row.names = 1, 
                   stringsAsFactors = FALSE)
expr <- as.matrix(expr)

# 读取样本信息
con <- gzfile("GSE23160_series_matrix.txt.gz", "rt")
lines <- readLines(con, warn = FALSE)
close(con)

# 提取所有Sample相关行
sample_accession <- grep("^!Sample_geo_accession", lines, value = TRUE)
sample_title <- grep("^!Sample_title", lines, value = TRUE)
sample_characteristics <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)

cat("=== 样本完整信息 ===\n")
for(i in 1:min(8, length(sample_title))) {
  title <- strsplit(sample_title[i], "\t")[[1]][-1]
  cat("\n样本", i, ":\n")
  cat("  Accession:", strsplit(sample_accession[i], "\t")[[1]][-1], "\n")
  cat("  Title:", title, "\n")
}

# 关键：检查Series_matrix.txt.gz中的数据处理信息
cat("\n=== 数据预处理信息 ===\n")
data_processing <- grep("!Series_data_processing", lines, value = TRUE)
for(dp in data_processing) {
  cat("  ", dp, "\n")
}

# 检查数据范围（确认是否已log2）
cat("\n=== 数据范围检查 ===\n")
cat("原始数据范围:", range(expr, na.rm = TRUE), "\n")
cat("如果范围>20，说明是线性尺度；如果<20，说明已是log2尺度\n")

# 检查Sham样本在每个时间点的情况
cat("\n=== 关键问题：Sham是否有时间点？ ===\n")
cat("GEO的series_matrix通常将Sham作为共同对照\n")
cat("但有些研究会在不同时间点设置不同的Sham对照\n")

# 尝试用GEO2R的标准方法
# GEO2R读取的是处理后的值，直接使用

# 测试不同数据转换方式下的DEG数量
test_deg_count <- function(expr_matrix, region, time_point, sham_indices, ir_indices) {
  expr_sub <- expr_matrix[, c(sham_indices, ir_indices)]
  
  group <- factor(c(rep("Sham", length(sham_indices)), rep("IR", length(ir_indices))))
  design <- model.matrix(~0 + group)
  colnames(design) <- c("Sham", "IR")
  
  fit <- lmFit(expr_sub, design)
  contrast.matrix <- makeContrasts(IRvsSham = IR - Sham, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  results <- topTable(fit2, coef = "IRvsSham", number = Inf, sort.by = "P")
  
  # 计数
  c(
    log2_P0.05_FC0.585 = sum(results$P.Value < 0.05 & abs(results$logFC) >= 0.585),
    log2_P0.05_FC0 = sum(results$P.Value < 0.05 & abs(results$logFC) >= 0),
    log2_P0.05_FC1.0 = sum(results$P.Value < 0.05 & abs(results$logFC) >= 1.0)
  )
}

# Cortex样本索引（1-16）
cortex_sham <- 1:4
cortex_ir_2h <- 5:8
cortex_ir_8h <- 9:12
cortex_ir_24h <- 13:16

# Striatum样本索引（17-32）
striatum_sham <- 17:20
striatum_ir_2h <- 21:24
striatum_ir_8h <- 25:28
striatum_ir_24h <- 29:32

# 测试1：原始数据
cat("\n=== 测试1：原始数据 ===\n")
cat("Cortex 2h vs Sham:\n")
print(test_deg_count(expr, "Cortex", "2h", cortex_sham, cortex_ir_2h))

cat("Cortex 8h vs Sham:\n")
print(test_deg_count(expr, "Cortex", "8h", cortex_sham, cortex_ir_8h))

cat("Cortex 24h vs Sham:\n")
print(test_deg_count(expr, "Cortex", "24h", cortex_sham, cortex_ir_24h))

# 测试2：log2转换后
cat("\n=== 测试2：log2(expr+1) ===\n")
expr_log2 <- log2(expr + 1)
cat("Cortex 2h vs Sham:\n")
print(test_deg_count(expr_log2, "Cortex", "2h", cortex_sham, cortex_ir_2h))

cat("Cortex 8h vs Sham:\n")
print(test_deg_count(expr_log2, "Cortex", "8h", cortex_sham, cortex_ir_8h))

cat("Cortex 24h vs Sham:\n")
print(test_deg_count(expr_log2, "Cortex", "24h", cortex_sham, cortex_ir_24h))

# 测试3：log2转换 + 分位数标准化
cat("\n=== 测试3：log2(expr+1) + 分位数标准化 ===\n")
expr_norm <- normalizeBetweenArrays(expr_log2, method = "quantile")
cat("Cortex 2h vs Sham:\n")
print(test_deg_count(expr_norm, "Cortex", "2h", cortex_sham, cortex_ir_2h))

cat("Cortex 8h vs Sham:\n")
print(test_deg_count(expr_norm, "Cortex", "8h", cortex_sham, cortex_ir_8h))

cat("Cortex 24h vs Sham:\n")
print(test_deg_count(expr_norm, "Cortex", "24h", cortex_sham, cortex_ir_24h))

# 测试4：分位数标准化（线性尺度）
cat("\n=== 测试4：分位数标准化（线性尺度）===\n")
expr_norm_linear <- normalizeBetweenArrays(expr, method = "quantile")
cat("Cortex 2h vs Sham:\n")
print(test_deg_count(expr_norm_linear, "Cortex", "2h", cortex_sham, cortex_ir_2h))

cat("Cortex 8h vs Sham:\n")
print(test_deg_count(expr_norm_linear, "Cortex", "8h", cortex_sham, cortex_ir_8h))

cat("Cortex 24h vs Sham:\n")
print(test_deg_count(expr_norm_linear, "Cortex", "24h", cortex_sham, cortex_ir_24h))

# 测试5：如果Sham有各自时间点（假设Sham样本按顺序分配）
cat("\n=== 测试5：假设Sham按时间点分配（4个Sham分给3个时间点）===\n")
cat("这个假设不太可能，因为4个Sham不能整除3个时间点\n")

cat("\n=== 诊断完成 ===\n")
