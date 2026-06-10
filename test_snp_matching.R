#!/usr/bin/env Rscript
# 测试单个基因的 SNP 匹配和等位基因情况

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(readxl)
  library(data.table)
})

cat("测试 PRKCQ 基因的 SNP 匹配情况\n")
cat("=========================================\n\n")

# 加载 eQTLGen 数据
eqtlgen_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-8_kb_10000_r2_0.001.xlsx"
eqtlgen_data <- read_excel(eqtlgen_file)

# 获取 PRKCQ 的 SNP
prkcq_snps <- eqtlgen_data[eqtlgen_data$gene == "PRKCQ" | eqtlgen_data$gene == "ENSG00000065675", ]
cat("PRKCQ 暴露数据:\n")
print(prkcq_snps[, c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "beta.exposure", "se.exposure")])
cat("\n")

# 格式化暴露数据
exposure <- data.frame(
  SNP = prkcq_snps$SNP,
  BETA = prkcq_snps$beta.exposure,
  SE = prkcq_snps$se.exposure,
  EFFECT_ALLELE = prkcq_snps$effect_allele.exposure,
  OTHER_ALLELE = prkcq_snps$other_allele.exposure,
  EAF = prkcq_snps$eaf.exposure,
  PVAL = prkcq_snps$pval.exposure,
  stringsAsFactors = FALSE
)

exposure_fmt <- format_data(exposure, type = "exposure")
cat("格式化后的暴露数据:\n")
print(exposure_fmt[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")])
cat("\n")

# 加载 MEGASTROKE 数据
megastroke_file <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"
outcome_data <- fread(megastroke_file, sep = "\t", stringsAsFactors = FALSE)

# 提取匹配的 SNP
outcome_matched <- outcome_data[outcome_data$hm_rsid %in% prkcq_snps$SNP, ]
cat("MEGASTROKE 匹配的结局数据:\n")
print(outcome_matched[, c("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "hm_effect_allele_frequency")])
cat("\n")

# 格式化结局数据
outcome_data_df <- as.data.frame(outcome_data)
outcome_fmt <- format_data(
  outcome_data_df,
  type = "outcome",
  snp_col = "hm_rsid",
  beta_col = "hm_beta",
  se_col = "hm_beta",
  effect_allele_col = "hm_effect_allele",
  other_allele_col = "hm_other_allele",
  eaf_col = "hm_effect_allele_frequency",
  pval_col = "hm_beta"
)

# 提取匹配的结局 SNP
outcome_matched_fmt <- outcome_fmt[outcome_fmt$SNP %in% prkcq_snps$SNP, ]
cat("格式化后的结局数据:\n")
print(outcome_matched_fmt[, c("SNP", "beta.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")])
cat("\n")

# 尝试 harmonise
cat("尝试 harmonise_data...\n")
dat <- harmonise_data(exposure_fmt, outcome_fmt)
cat("Harmonised 数据:\n")
if (nrow(dat) > 0) {
  print(dat[, c("SNP", "beta.exposure", "beta.outcome", "effect_allele.exposure", "effect_allele.outcome", "other_allele.exposure", "other_allele.outcome", "samplesize.outcome")])
} else {
  cat("  无数据返回！\n")
}
cat("\n")

cat("问题诊断:\n")
cat("1. 暴露和结局的 SNP ID 是否匹配？", length(intersect(exposure_fmt$SNP, outcome_fmt$SNP)), "/", nrow(exposure_fmt), "\n")
cat("2. harmonise_data 返回的数据行数：", nrow(dat), "\n")
if (nrow(dat) == 0) {
  cat("3. 可能原因：等位基因不匹配或所有 SNP 都被过滤\n")
}
