#!/usr/bin/env Rscript
# ================================================================================
# 检查 MEGASTROKE 数据列名
# ================================================================================

suppressPackageStartupMessages({
  library(data.table)
})

cat("======================================================================\n")
cat("检查 MEGASTROKE 数据列名\n")
cat("======================================================================\n\n")

MEGASTROKE_FILE <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"

cat("加载 MEGASTROKE 数据...\n")
outcome_data <- fread(MEGASTROKE_FILE, sep = "\t", stringsAsFactors = FALSE, nrows = 10)

cat("MEGASTROKE 数据列名:\n")
print(names(outcome_data))
cat("\n")

cat("前 3 行数据:\n")
print(head(outcome_data, 3))
cat("\n")

cat("检查关键列:\n")
cat("hm_rsid:", "hm_rsid" %in% names(outcome_data), "\n")
cat("hm_beta:", "hm_beta" %in% names(outcome_data), "\n")
cat("hm_se:", "hm_se" %in% names(outcome_data), "\n")
cat("hm_effect_allele:", "hm_effect_allele" %in% names(outcome_data), "\n")
cat("hm_other_allele:", "hm_other_allele" %in% names(outcome_data), "\n")
cat("hm_effect_allele_frequency:", "hm_effect_allele_frequency" %in% names(outcome_data), "\n")
cat("hm_pval:", "hm_pval" %in% names(outcome_data), "\n")
