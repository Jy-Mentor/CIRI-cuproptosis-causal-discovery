#!/usr/bin/env Rscript
# 检查VCF文件的FORMAT格式

library(data.table)

# 检查VCF格式
check_vcf_format <- function() {
  vcf_file <- "D:/EQTL/ieu-a-83.vcf"
  
  # 读取前150行（包含header和一些数据行）
  lines <- readLines(vcf_file, n = 150)
  
  # 找到header行
  header_line <- lines[grep("^#CHROM", lines)]
  cat("Header line:\n")
  cat(header_line, "\n\n")
  
  # 找到FORMAT列的索引
  header_parts <- strsplit(header_line, "\t")[[1]]
  format_idx <- which(header_parts == "FORMAT")
  cat(sprintf("FORMAT column index: %d\n\n", format_idx))
  
  # 找到第一个数据行
  data_lines <- lines[!grepl("^#", lines)]
  if (length(data_lines) > 0) {
    first_data_line <- data_lines[1]
    cat("First data line:\n")
    cat(first_data_line, "\n\n")
    
    # 提取FORMAT字段
    data_parts <- strsplit(first_data_line, "\t")[[1]]
    if (length(data_parts) >= format_idx) {
      format_field <- data_parts[format_idx]
      cat(sprintf("FORMAT field: %s\n", format_field))
      
      # 提取FORMAT标签
      format_tags <- strsplit(format_field, ":")[[1]]
      cat("FORMAT tags:\n")
      print(format_tags)
      
      # 提取OUTCOME字段（FORMAT后面的列）
      if (length(data_parts) > format_idx) {
        outcome_field <- data_parts[format_idx + 1]
        cat(sprintf("\nOUTCOME field: %s\n", outcome_field))
        
        # 分割OUTCOME字段
        outcome_tags <- strsplit(outcome_field, ":")[[1]]
        cat("OUTCOME values:\n")
        print(outcome_tags)
        
        # 显示标签和值的对应关系
        cat("\nTag-value mapping:\n")
        for (i in seq_along(format_tags)) {
          if (i <= length(outcome_tags)) {
            cat(sprintf("%s: %s\n", format_tags[i], outcome_tags[i]))
          }
        }
      }
    }
  }
}

# 运行检查
if (sys.nframe() == 0) {
  check_vcf_format()
}
