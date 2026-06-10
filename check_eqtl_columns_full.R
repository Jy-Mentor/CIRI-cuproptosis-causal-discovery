#!/usr/bin/env Rscript
# ================================================================================
# 诊断 eQTLGen 数据列名
# ================================================================================

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

cat("检查关键列:\n")
cat("SNP:", "SNP" %in% names(eqtlgen_data), "\n")
cat("beta.exposure:", "beta.exposure" %in% names(eqtlgen_data), "\n")
cat("se.exposure:", "se.exposure" %in% names(eqtlgen_data), "\n")
cat("effect_allele.exposure:", "effect_allele.exposure" %in% names(eqtlgen_data), "\n")
cat("other_allele.exposure:", "other_allele.exposure" %in% names(eqtlgen_data), "\n")
cat("eaf.exposure:", "eaf.exposure" %in% names(eqtlgen_data), "\n")
cat("pval.exposure:", "pval.exposure" %in% names(eqtlgen_data), "\n")
