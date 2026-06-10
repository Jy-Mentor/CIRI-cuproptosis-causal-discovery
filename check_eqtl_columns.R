#!/usr/bin/env Rscript
# 检查 eQTLGen 数据的列名

suppressPackageStartupMessages({
  library(readxl)
})

CLUMPED_FILE <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.001.xlsx"
eqtlgen_data <- read_excel(CLUMPED_FILE)

cat("eQTLGen 数据列名:\n")
print(names(eqtlgen_data))
cat("\n")

cat("前 3 行数据:\n")
print(head(eqtlgen_data, 3))
cat("\n")

cat("检查 SNP 列:\n")
cat("SNP 列是否存在:", "SNP" %in% names(eqtlgen_data), "\n")
cat("SNP 列类型:", class(eqtlgen_data$SNP), "\n")
cat("前 5 个 SNP:", head(eqtlgen_data$SNP, 5), "\n")
