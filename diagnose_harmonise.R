#!/usr/bin/env Rscript
# ================================================================================
# 诊断 Harmonise 失败问题
# ================================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(readxl)
})

cat("======================================================================\n")
cat("诊断 Harmonise 失败问题\n")
cat("======================================================================\n\n")

# 配置
CLUMPED_FILE <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.001.xlsx"
MEGASTROKE_FILE <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"

# 加载数据
cat("加载 eQTLGen 数据...\n")
eqtlgen_data <- read_excel(CLUMPED_FILE)
eqtlgen_data <- eqtlgen_data[!is.na(eqtlgen_data$SNP) & grepl("^rs", eqtlgen_data$SNP), ]

cat("加载 MEGASTROKE 数据...\n")
outcome_data <- fread(MEGASTROKE_FILE, sep = "\t", stringsAsFactors = FALSE)
outcome_data <- outcome_data[!is.na(outcome_data$hm_rsid) & outcome_data$hm_rsid != "", ]

cat(sprintf("eQTLGen: %d SNPs, MEGASTROKE: %d SNPs\n\n", nrow(eqtlgen_data), nrow(outcome_data)))

# 测试基因：LYN
gene_symbol <- "LYN"
ensg_id <- "ENSG00000254087"

cat(sprintf("测试基因：%s (ENSG: %s)\n\n", gene_symbol, ensg_id))

# 获取 eQTL 数据
eqtlgen_snps <- eqtlgen_data[eqtlgen_data$gene == gene_symbol | eqtlgen_data$gene == ensg_id, ]
cat(sprintf("找到 %d 个 eQTL SNPs\n", nrow(eqtlgen_snps)))

if (nrow(eqtlgen_snps) > 0) {
  cat("\n暴露数据前 3 行:\n")
  print(head(eqtlgen_snps[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", 
                               "other_allele.exposure", "eaf.exposure", "pval.exposure")], 3))
  
  # 准备暴露数据
  exposure <- eqtlgen_snps[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", 
                                "other_allele.exposure", "eaf.exposure", "pval.exposure")]
  names(exposure) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", 
                       "other_allele.exposure", "eaf.exposure", "pval.exposure")
  
  exposure <- exposure[!duplicated(exposure$SNP) & !is.na(exposure$SNP) & exposure$SNP != "", ]
  cat(sprintf("\n去重后暴露 SNP 数：%d\n", nrow(exposure)))
  
  # 格式化暴露数据
  cat("\n格式化暴露数据...\n")
  exposure_fmt <- format_data(exposure, type = "exposure")
  cat(sprintf("格式化后暴露数据维度：%d x %d\n", nrow(exposure_fmt), ncol(exposure_fmt)))
  cat("暴露数据列名:\n")
  print(names(exposure_fmt))
  cat("\n暴露数据前 3 行:\n")
  print(head(exposure_fmt[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", 
                               "other_allele.exposure", "eaf.exposure")], 3))
  
  # 提取结局数据
  cat("\n提取匹配的结局数据...\n")
  outcome_matched <- outcome_data[outcome_data$hm_rsid %in% exposure_fmt$SNP, ]
  cat(sprintf("匹配的结局 SNP 数：%d\n", nrow(outcome_matched)))
  
  if (nrow(outcome_matched) > 0) {
    cat("\n结局数据前 3 行:\n")
    print(head(outcome_matched[, c("hm_rsid", "hm_beta", "hm_se", "hm_effect_allele", 
                                    "hm_other_allele", "hm_effect_allele_frequency")], 3))
    
    # 准备结局数据
    outcome_df <- as.data.frame(outcome_matched)
    outcome_df$beta.outcome <- as.numeric(outcome_df$hm_beta)
    outcome_df$se.outcome <- abs(as.numeric(outcome_df$hm_beta))
    outcome_df$effect_allele.outcome <- as.character(outcome_df$hm_effect_allele)
    outcome_df$other_allele.outcome <- as.character(outcome_df$hm_other_allele)
    outcome_df$eaf.outcome <- as.numeric(outcome_df$hm_effect_allele_frequency)
    outcome_df$pval.outcome <- 0.05
    
    cat("\n格式化结局数据...\n")
    outcome_fmt <- format_data(outcome_df, type = "outcome", snp_col = "hm_rsid")
    cat(sprintf("格式化后结局数据维度：%d x %d\n", nrow(outcome_fmt), ncol(outcome_fmt)))
    cat("结局数据列名:\n")
    print(names(outcome_fmt))
    cat("\n结局数据前 3 行:\n")
    print(head(outcome_fmt[, c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", 
                                "other_allele.outcome", "eaf.outcome")], 3))
    
    # 检查两个数据集的 SNP 重叠
    cat("\n检查 SNP 重叠:\n")
    cat(sprintf("暴露 SNP 数：%d\n", nrow(exposure_fmt)))
    cat(sprintf("结局 SNP 数：%d\n", nrow(outcome_fmt)))
    cat(sprintf("共同 SNP 数：%d\n", length(intersect(exposure_fmt$SNP, outcome_fmt$SNP))))
    
    # Harmonise
    cat("\n尝试 Harmonise...\n")
    dat <- tryCatch({
      harmonise_data(exposure_fmt, outcome_fmt)
    }, error = function(e) {
      cat(sprintf("Harmonise 失败：%s\n", e$message))
      return(NULL)
    })
    
    if (!is.null(dat)) {
      cat(sprintf("Harmonise 成功！维度：%d x %d\n", nrow(dat), ncol(dat)))
      cat("\nHarmonised 数据前 3 行:\n")
      print(head(dat[, c("SNP", "beta.exposure", "beta.outcome", "effect_allele.exposure", 
                          "effect_allele.outcome")], 3))
    }
  } else {
    cat("无匹配的结局数据\n")
  }
} else {
  cat("无 eQTL 数据\n")
}

cat("\n完成！\n")
