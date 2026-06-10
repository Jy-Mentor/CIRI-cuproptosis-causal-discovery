#!/usr/bin/env Rscript
# 测试VCF解析功能

library(data.table)
library(futile.logger)

# 设置日志
flog.appender(appender.console())
flog.threshold(INFO)

# 快速解析GWAS VCF文件（向量化版本）
parse_gwas_vcf_fast <- function(vcf_file, nrows = Inf) {
  flog.info(sprintf("快速解析VCF文件: %s", vcf_file))
  
  # 找到header行数
  preview <- readLines(vcf_file, n = 1000)
  header_line <- which(grepl("^#CHROM", preview))[1]
  if (is.na(header_line)) stop("无法找到VCF header")
  
  flog.info(sprintf("Header位于第%d行，开始流式读取...", header_line))
  
  # 使用fread从header开始读取
  vcf_data <- fread(
    vcf_file,
    skip = header_line - 1,
    sep = "\t",
    header = TRUE,
    nrows = nrows,  # 测试时可用1000限制行数
    select = c(3, 4, 5, 9, 10),  # ID, REF, ALT, FORMAT, OUTCOME列（按位置选更快）
    col.names = c("ID", "REF", "ALT", "FORMAT", "OUTCOME"),
    showProgress = TRUE
  )
  
  flog.info(sprintf("读取完成: %d行，开始解析FORMAT字段...", nrow(vcf_data)))
  
  # 向量化解析（关键优化）
  # 1. 提取FORMAT标签（仅第一行）
  format_tags <- strsplit(vcf_data$FORMAT[1], ":")[[1]]
  
  # 2. 找到标签索引，处理标签不存在的情况
  get_tag_index <- function(tag) {
    idx <- which(format_tags == tag)
    if (length(idx) == 0) return(NA) else return(idx)
  }
  
  tag_idx <- list(
    ES = get_tag_index("ES"),
    SE = get_tag_index("SE"),
    LP = get_tag_index("LP"),
    AF = get_tag_index("AF")
  )
  
  # 3. 一次性分割所有OUTCOME字段（data.table内置tstrsplit是C级速度）
  outcome_parts <- tstrsplit(vcf_data$OUTCOME, ":", fixed = TRUE)
  
  # 4. 提取数值（向量化，无循环）
  vcf_data[, beta.outcome := as.numeric(outcome_parts[[tag_idx$ES]])]
  vcf_data[, se.outcome := as.numeric(outcome_parts[[tag_idx$SE]])]
  vcf_data[, LP := as.numeric(outcome_parts[[tag_idx$LP]])]
  
  # 处理AF标签不存在的情况
  if (!is.na(tag_idx$AF)) {
    vcf_data[, AF := as.numeric(outcome_parts[[tag_idx$AF]])]
  } else {
    vcf_data[, AF := NA]  # 如果没有AF标签，设为NA
  }
  
  vcf_data[, pval.outcome := 10^(-LP)]
  
  # 4. 清理和重命名
  vcf_data <- vcf_data[!is.na(beta.outcome) & !is.na(se.outcome) & !is.na(pval.outcome),
                       .(SNP = ID, REF, ALT, beta.outcome, se.outcome, pval.outcome, AF)]
  
  flog.info(sprintf("VCF解析完成: %d行有效数据", nrow(vcf_data)))
  return(vcf_data)
}

# 测试函数
test_vcf_parse <- function() {
  vcf_file <- "D:/EQTL/ieu-a-83.vcf"
  
  # 测试1：文件是否存在
  if (!file.exists(vcf_file)) {
    flog.error(sprintf("VCF文件不存在: %s", vcf_file))
    return(FALSE)
  }
  
  # 测试2：读取前100行
  flog.info("测试1：读取前100行...")
  tryCatch({
    preview <- readLines(vcf_file, n = 100)
    flog.info(sprintf("前100行读取成功，共%d行", length(preview)))
  }, error = function(e) {
    flog.error(sprintf("读取前100行失败: %s", e$message))
    return(FALSE)
  })
  
  # 测试3：解析前1000行
  flog.info("测试2：解析前1000行...")
  tryCatch({
    vcf_data <- parse_gwas_vcf_fast(vcf_file, nrows = 1000)
    flog.info(sprintf("解析成功，返回%d行数据", nrow(vcf_data)))
    flog.info("前5行数据:")
    print(head(vcf_data, 5))
  }, error = function(e) {
    flog.error(sprintf("解析失败: %s", e$message))
    return(FALSE)
  })
  
  # 测试4：文件大小
  file_info <- file.info(vcf_file)
  flog.info(sprintf("文件大小: %.2f GB", file_info$size / 1e9))
  
  return(TRUE)
}

# 运行测试
if (sys.nframe() == 0) {
  flog.info("开始VCF解析测试...")
  success <- test_vcf_parse()
  if (success) {
    flog.info("VCF解析测试成功！")
  } else {
    flog.error("VCF解析测试失败！")
  }
}
