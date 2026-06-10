# 诊断脚本：检查数据问题和对照组设置
library(limma)
library(ggplot2)
library(dplyr)

setwd("D:/反向网络药理学/L1 数据集/bulk/GSE23160(主验证集时序差异分析，2h,8h,24h)")

# ==================== 1. 检查数据基本信息 ====================
expr <- read.delim("GSE23160_series_matrix.txt.gz", 
                   comment.char = "!", 
                   header = TRUE, 
                   check.names = FALSE, 
                   row.names = 1, 
                   stringsAsFactors = FALSE)
expr <- as.matrix(expr)

cat("=== 数据基本信息 ===\n")
cat("数据范围:", range(expr, na.rm=TRUE), "\n")
cat("中位数:", median(expr, na.rm=TRUE), "\n")
cat("均值:", mean(expr, na.rm=TRUE), "\n")
cat("样本数:", ncol(expr), "\n")
cat("探针数:", nrow(expr), "\n")

cat("\n=== 检查数据是否在log2尺度 ===\n")
cat("如果数据>100，可能是线性尺度；如果<20，可能是log2尺度\n")
cat("前5个探针的表达值范围:\n")
for(i in 1:5) {
  cat("  探针", rownames(expr)[i], ": range =", range(expr[i,]), "\n")
}

# ==================== 2. 解析样本信息 ====================
con <- gzfile("GSE23160_series_matrix.txt.gz", "rt")
lines <- readLines(con, warn = FALSE)
close(con)

# 找Sample_title
title_lines <- grep("^!Sample_title", lines, value = TRUE)
sample_titles <- unlist(strsplit(title_lines, "\t"))[-1]
sample_titles <- gsub('"', '', sample_titles)

# 找Sample_characteristics_ch1
char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
cat("\n=== 样本特征信息（前5个样本）===\n")
for(i in 1:min(5, length(char_lines))) {
  cat("  样本", i, ":\n")
  cat("    Title:", sample_titles[i], "\n")
  # 找对应的characteristics
  idx <- grep(paste0("!Sample_geo_accession.*", strsplit(lines[grep("!Sample_geo_accession", lines)[i]], "\t")[2]), lines)
  # 直接打印前20个characteristics行
  if(i <= 5) {
    char_idx <- grep("^!Sample_characteristics_ch1", lines)
    if(length(char_idx) >= i) {
      cat("    Characteristics:", lines[char_idx[i]], "\n")
    }
  }
}

cat("\n=== 所有样本标题 ===\n")
for(i in seq_along(sample_titles)) {
  cat("  ", i, ": ", sample_titles[i], "\n", sep="")
}

# ==================== 3. 检查Sham样本是否有时间点信息 ====================
cat("\n=== Sham样本详细信息 ===\n")
sham_titles <- sample_titles[grepl("Sham", sample_titles, ignore.case = TRUE)]
for(t in sham_titles) {
  cat("  ", t, "\n")
}

# 查找所有Sample_characteristics
cat("\n=== 所有Sample_characteristics_ch1 ===\n")
char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
for(i in seq_along(char_lines)) {
  cat("  样本", i, ": ", char_lines[i], "\n", sep="")
}

# ==================== 4. 尝试不同的分组方式 ====================
cat("\n=== 重新解析样本分组 ===\n")

sample_info <- data.frame(
  sample_id = colnames(expr),
  title = sample_titles,
  stringsAsFactors = FALSE
)

# 方法1：按原始标题解析
sample_info$brain_region <- ifelse(grepl("^Cortex-", sample_titles, ignore.case = TRUE), "Cortex",
                            ifelse(grepl("^Striatum-", sample_titles, ignore.case = TRUE), "Striatum", "Unknown"))

sample_info$treatment <- ifelse(grepl("-Sham-", sample_titles, ignore.case = TRUE), "Sham", "IR")

sample_info$time_group <- ifelse(sample_info$treatment == "Sham", "Ctrl",
                          ifelse(grepl("-2h-", sample_titles, ignore.case = TRUE), "2h",
                          ifelse(grepl("-8h-", sample_titles, ignore.case = TRUE), "8h",
                          ifelse(grepl("-24h-", sample_titles, ignore.case = TRUE), "24h", "Unknown"))))

cat("当前分组方式1（Sham=Ctrl）:\n")
print(table(sample_info$brain_region, sample_info$treatment, sample_info$time_group))

# 方法2：假设Sham也有时间点（根据文件顺序或其他信息）
# GSE23160的Sham通常是所有时间点的共同对照
# 但如果Sham是特定时间点的对照，则需要重新分组

cat("\n=== 尝试方法2：Sham作为每个时间点的对照（当前方法）===\n")
cat("每个时间点: IR(4) vs Sham(4)\n")
cat("总样本: Cortex(16), Striatum(16)\n")

# ==================== 5. 对比Zhao文献的方法 ====================
cat("\n=== 尝试Zhao文献的方法：所有IR vs Sham（不区分时间点）===\n")

for(region in c("Cortex", "Striatum")) {
  cat("\n---", region, "---\n")
  
  region_idx <- sample_info$brain_region == region
  expr_region <- expr[, region_idx]
  sample_region <- sample_info[region_idx, ]
  
  # Zhao方法：所有IR vs Sham
  design <- model.matrix(~0 + treatment, data = sample_region)
  colnames(design) <- c("IR", "Sham")
  
  fit <- lmFit(expr_region, design)
  contrast.matrix <- makeContrasts(IRvsSham = IR - Sham, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  results <- topTable(fit2, coef = "IRvsSham", number = Inf, sort.by = "P")
  
  # 测试不同阈值
  for(fc_thresh in c(0, 0.263, 0.585, 1.0)) {
    sig <- results %>% filter(P.Value < 0.05 & abs(logFC) >= fc_thresh)
    cat(sprintf("  P<0.05, |log2FC|>=%.3f: %d DEGs (上调: %d, 下调: %d)\n", 
                fc_thresh, nrow(sig), sum(sig$logFC > 0), sum(sig$logFC < 0)))
  }
}

cat("\n=== 诊断完成 ===\n")
