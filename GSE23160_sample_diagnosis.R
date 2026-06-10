# GSE23160 样本信息诊断脚本
series_matrix_file <- "D:/反向网络药理学/L1 数据集/bulk/GSE23160(主验证集时序差异分析，2h,8h,24h)/GSE23160_series_matrix.txt.gz"

con <- gzfile(series_matrix_file, "rt")
lines <- readLines(con, warn = FALSE)
close(con)

# 提取所有样本相关元数据行
meta_lines <- grep("^!Sample_", lines, value = TRUE)

# 查看有哪些样本字段
fields <- unique(sub("^!Sample_([^ ]+).*", "\\1", meta_lines))
cat("样本字段列表:\n")
print(fields)

cat("\n\n========================================\n")

# 提取前3个样本的所有信息（用于调试）
sample_title <- grep("^!Sample_title", lines, value = TRUE)
sample_char <- grep("^!Sample_characteristics", lines, value = TRUE)
sample_source <- grep("^!Sample_source_name", lines, value = TRUE)

cat("!Sample_title:\n")
cat(sample_title, sep = "\n")

cat("\n\n!Sample_characteristics:\n")
if(length(sample_char) > 0) {
  cat(sample_char, sep = "\n")
} else {
  cat("未找到Sample_characteristics字段\n")
}

cat("\n\n!Sample_source_name:\n")
if(length(sample_source) > 0) {
  cat(sample_source, sep = "\n")
} else {
  cat("未找到Sample_source_name字段\n")
}

# 检查是否有其他可能包含分组信息的字段
other_fields <- grep("^!Sample_", lines, value = TRUE)
cat("\n\n所有!Sample_开头的行（前3个样本）:\n")
sample_count <- length(grep("^!Sample_title", lines, value = TRUE))
cat(paste0("总样本数: ", sample_count, "\n\n"))

# 提取每个样本的关键信息
for(i in 1:min(6, sample_count)) {
  cat(paste0("===== 样本 ", i, " =====\n"))
  for(field in c("!Sample_title", "!Sample_characteristics_ch1", "!Sample_characteristics_ch2", 
                 "!Sample_characteristics_ch3", "!Sample_source_name", "!Sample_organism")) {
    idx <- grep(paste0("^", field), lines)
    if(length(idx) >= i) {
      cat(paste0(field, ": ", lines[idx[i]], "\n"))
    }
  }
  cat("\n")
}
