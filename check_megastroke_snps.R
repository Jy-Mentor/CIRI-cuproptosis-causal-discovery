#!/usr/bin/env Rscript
# 检查 MEGASTROKE 数据中 PRKCQ 的 SNP 是否存在

suppressPackageStartupMessages({
  library(data.table)
})

cat("检查 MEGASTROKE 数据中的 PRKCQ SNP\n")
cat("=========================================\n\n")

# 加载 MEGASTROKE 数据
megastroke_file <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"
cat("加载 MEGASTROKE 数据...\n")
outcome_data <- fread(megastroke_file, sep = "\t", stringsAsFactors = FALSE)

cat(sprintf("总 SNP 数：%d\n\n", nrow(outcome_data)))

# PRKCQ 的 SNP
prkcq_snps <- c("rs2007252", "rs2255088", "rs658230")
cat("PRKCQ 的 SNP:\n")
print(prkcq_snps)
cat("\n")

# 检查这些 SNP 是否在 MEGASTROKE 中
cat("检查匹配...\n")
matched <- outcome_data[outcome_data$hm_rsid %in% prkcq_snps, ]
cat(sprintf("匹配的 SNP 数：%d\n", nrow(matched)))

if (nrow(matched) > 0) {
  cat("\n匹配的数据:\n")
  print(matched[, c("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "hm_effect_allele_frequency")])
} else {
  cat("\n未找到匹配的 SNP！\n")
  
  # 检查所有可能的 rsID 格式
  cat("\n检查 MEGASTROKE 的 rsID 格式...\n")
  cat("前 10 个 hm_rsid:\n")
  print(head(outcome_data$hm_rsid, 10))
  
  cat("\n检查是否有 rs 前缀的 SNP:\n")
  rs_snps <- outcome_data[grepl("^rs", outcome_data$hm_rsid), ]
  cat(sprintf("有 rs 前缀的 SNP 数：%d\n", nrow(rs_snps)))
  
  if (nrow(rs_snps) > 0) {
    cat("\n前 10 个 rsSNP:\n")
    print(head(rs_snps$hm_rsid, 10))
  }
}
