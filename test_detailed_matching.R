#!/usr/bin/env Rscript
# 详细测试 PRKCQ 的等位基因匹配情况

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(readxl)
  library(data.table)
})

cat("详细测试 PRKCQ 的等位基因匹配\n")
cat("=========================================\n\n")

# 1. 加载 eQTLGen 数据
cat("1. 加载 eQTLGen 数据\n")
eqtlgen_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-8_kb_10000_r2_0.001.xlsx"
eqtlgen_data <- read_excel(eqtlgen_file)

prkcq_snps <- eqtlgen_data[eqtlgen_data$gene == "PRKCQ" | eqtlgen_data$gene == "ENSG00000065675", ]
cat(sprintf("   PRKCQ 有 %d 个 SNP\n", nrow(prkcq_snps)))

# 2. 创建暴露数据
cat("\n2. 创建暴露数据\n")
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

cat("暴露数据:\n")
print(exposure[, c("SNP", "BETA", "EFFECT_ALLELE", "OTHER_ALLELE", "EAF")])
cat("\n")

# 3. 格式化暴露数据
cat("3. 格式化暴露数据\n")
exposure_fmt <- format_data(
  exposure,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  pval_col = "PVAL"
)

cat("格式化后的暴露数据:\n")
print(exposure_fmt[, c("SNP", "beta.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")])
cat("\n")

# 4. 加载 MEGASTROKE 数据
cat("4. 加载 MEGASTROKE 数据\n")
megastroke_file <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"
outcome_data <- fread(megastroke_file, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("   总 SNP 数：%d\n", nrow(outcome_data)))

# 5. 提取 PRKCQ 的 SNP
cat("\n5. 提取 PRKCQ 的 SNP\n")
outcome_matched <- outcome_data[outcome_data$hm_rsid %in% prkcq_snps$SNP, ]
cat(sprintf("   匹配的 SNP 数：%d\n", nrow(outcome_matched)))

cat("\nMEGASTROKE 原始数据:\n")
print(outcome_matched[, c("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "hm_effect_allele_frequency")])
cat("\n")

# 6. 格式化结局数据
cat("6. 格式化结局数据\n")
outcome_data_df <- as.data.frame(outcome_matched)
outcome_data_df$beta.outcome <- outcome_data_df$hm_beta
outcome_data_df$se.outcome <- abs(outcome_data_df$hm_beta)
outcome_data_df$effect_allele.outcome <- outcome_data_df$hm_effect_allele
outcome_data_df$other_allele.outcome <- outcome_data_df$hm_other_allele
outcome_data_df$eaf.outcome <- outcome_data_df$hm_effect_allele_frequency
outcome_data_df$pval.outcome <- 0.05

outcome_fmt <- format_data(
  outcome_data_df,
  type = "outcome",
  snp_col = "hm_rsid",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  eaf_col = "eaf.outcome",
  pval_col = "pval.outcome"
)

cat("格式化后的结局数据:\n")
print(outcome_fmt[, c("SNP", "beta.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")])
cat("\n")

# 7. 手动检查等位基因匹配
cat("7. 手动检查等位基因匹配\n")
for (i in 1:nrow(exposure_fmt)) {
  snp <- exposure_fmt$SNP[i]
  exp_ea <- exposure_fmt$effect_allele.exposure[i]
  exp_oa <- exposure_fmt$other_allele.exposure[i]
  
  outcome_row <- outcome_fmt[outcome_fmt$SNP == snp, ]
  if (nrow(outcome_row) > 0) {
    out_ea <- outcome_row$effect_allele.outcome[1]
    out_oa <- outcome_row$other_allele.outcome[1]
    
    match_status <- "完全匹配"
    if (exp_ea == out_ea && exp_oa == out_oa) {
      match_status <- "完全匹配 ✓"
    } else if (exp_ea == out_oa && exp_oa == out_ea) {
      match_status <- "链相反 (需要翻转)"
    }
    
    cat(sprintf("%s: 暴露 (%s/%s) vs 结局 (%s/%s) - %s\n",
                snp, exp_ea, exp_oa, out_ea, out_oa, match_status))
  } else {
    cat(sprintf("%s: 结局数据中不存在\n", snp))
  }
}
cat("\n")

# 8. 尝试 harmonise_data
cat("8. 尝试 harmonise_data\n")
dat <- harmonise_data(exposure_fmt, outcome_fmt)

if (nrow(dat) > 0) {
  cat(sprintf("✓ 成功！harmonise_data 返回 %d 个 SNP\n", nrow(dat)))
  cat("\nHarmonised 数据:\n")
  print(dat[, c("SNP", "beta.exposure", "beta.outcome", "effect_allele.exposure", 
                "effect_allele.outcome", "other_allele.exposure", "other_allele.outcome")])
} else {
  cat("✗ 失败！harmonise_data 返回 0 个 SNP\n")
  cat("\n可能原因:\n")
  cat("  1. 等位基因不匹配（暴露和结局的等位基因完全不同）\n")
  cat("  2. 所有 SNP 都被过滤（例如，palindromic SNPs 过多）\n")
  cat("  3. 数据格式问题\n")
}
cat("\n")
