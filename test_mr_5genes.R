#!/usr/bin/env Rscript
# 测试脚本 - 只处理前 5 个基因

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(readxl)
})

cat("测试 MR 分析（前 5 个基因）\n")

# 配置
CLUMPED_FILE <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.001.xlsx"
MEGASTROKE_FILE <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_test"

dir.create(OUTPUT_DIR, recursive = TRUE)

# 只测试前 5 个基因
test_genes <- c("LYN", "PRKCQ", "NMT1", "TDP1", "MAN2B1")

gene_to_ensg <- list(
  "LYN" = "ENSG00000254087",
  "PRKCQ" = "ENSG00000065675",
  "NMT1" = "ENSG00000136448",
  "TDP1" = "ENSG00000042088",
  "MAN2B1" = "ENSG00000104774"
)

# 加载数据
cat("加载 eQTLGen 数据...\n")
eqtlgen_data <- read_excel(CLUMPED_FILE)
eqtlgen_data <- eqtlgen_data[!is.na(eqtlgen_data$SNP) & grepl("^rs", eqtlgen_data$SNP), ]

cat("加载 MEGASTROKE 数据...\n")
outcome_data <- fread(MEGASTROKE_FILE, sep = "\t", stringsAsFactors = FALSE)
outcome_data <- outcome_data[!is.na(outcome_data$hm_rsid) & outcome_data$hm_rsid != "", ]

cat(sprintf("eQTLGen: %d SNPs, MEGASTROKE: %d SNPs\n\n", nrow(eqtlgen_data), nrow(outcome_data)))

# 处理每个基因
for (gene_symbol in test_genes) {
  cat(sprintf("\n=== %s ===\n", gene_symbol))
  
  ensg_id <- gene_to_ensg[[gene_symbol]]
  eqtlgen_snps <- eqtlgen_data[eqtlgen_data$gene == gene_symbol | eqtlgen_data$gene == ensg_id, ]
  
  if (nrow(eqtlgen_snps) < 2) {
    cat("SNP 数量不足，跳过\n")
    next
  }
  
  cat(sprintf("暴露 SNP 数：%d\n", nrow(eqtlgen_snps)))
  
  # 创建暴露数据
  exposure <- data.frame(
    SNP = eqtlgen_snps$SNP,
    BETA = eqtlgen_snps$beta.exposure,
    SE = eqtlgen_snps$se.exposure,
    EFFECT_ALLELE = eqtlgen_snps$effect_allele.exposure,
    OTHER_ALLELE = eqtlgen_snps$other_allele.exposure,
    EAF = eqtlgen_snps$eaf.exposure,
    PVAL = eqtlgen_snps$pval.exposure,
    stringsAsFactors = FALSE
  )
  
  exposure <- exposure[!duplicated(exposure$SNP), ]
  exposure_fmt <- format_data(exposure, type = "exposure")
  
  # 提取结局数据
  outcome_matched <- outcome_data[outcome_data$hm_rsid %in% exposure_fmt$SNP, ]
  cat(sprintf("结局匹配 SNP 数：%d\n", nrow(outcome_matched)))
  
  if (nrow(outcome_matched) == 0) {
    cat("无匹配结局数据，跳过\n")
    next
  }
  
  # 格式化结局数据
  outcome_df <- as.data.frame(outcome_matched)
  outcome_df$beta.outcome <- outcome_df$hm_beta
  outcome_df$se.outcome <- abs(outcome_df$hm_beta)
  outcome_df$effect_allele.outcome <- outcome_df$hm_effect_allele
  outcome_df$other_allele.outcome <- outcome_df$hm_other_allele
  outcome_df$eaf.outcome <- outcome_df$hm_effect_allele_frequency
  outcome_df$pval.outcome <- 0.05
  
  outcome_fmt <- format_data(outcome_df, type = "outcome", snp_col = "hm_rsid")
  
  # Harmonise
  dat <- harmonise_data(exposure_fmt, outcome_fmt)
  cat(sprintf("Harmonised SNP 数：%d\n", nrow(dat)))
  
  if (nrow(dat) == 0) {
    cat("Harmonise 失败，跳过\n")
    next
  }
  
  # MR 分析
  res <- mr(dat)
  cat("MR 结果:\n")
  print(res[, c("method", "b", "se", "pval")])
  
  # 保存结果
  write.csv(dat, file.path(OUTPUT_DIR, paste0(gene_symbol, "_harmonised.csv")), row.names = FALSE)
  write.csv(res, file.path(OUTPUT_DIR, paste0(gene_symbol, "_mr_results.csv")), row.names = FALSE)
}

cat("\n完成！结果保存在:", OUTPUT_DIR, "\n")
