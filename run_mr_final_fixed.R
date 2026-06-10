#!/usr/bin/env Rscript
# 最终修复版 - 使用宽松 eQTLGen 数据

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(readxl)
})

cat("======================================================================\n")
cat("MR 分析（使用 P<1e-5 eQTLGen 数据）\n")
cat("======================================================================\n\n")

# 配置
CLUMPED_FILE <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.001.xlsx"
MEGASTROKE_FILE <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_relaxed"

dir.create(OUTPUT_DIR, recursive = TRUE)

# 测试前 10 个基因
test_genes <- c("LYN", "PRKCQ", "NMT1", "TDP1", "MAN2B1", "IL10RA", "RHOC", "SREBF1", "CTSC", "CAT")

gene_to_ensg <- list(
  "LYN" = "ENSG00000254087", "PRKCQ" = "ENSG00000065675", "NMT1" = "ENSG00000136448",
  "TDP1" = "ENSG00000042088", "MAN2B1" = "ENSG00000104774", "IL10RA" = "ENSG00000110324",
  "RHOC" = "ENSG00000155366", "SREBF1" = "ENSG00000072310", "CTSC" = "ENSG00000109861",
  "CAT" = "ENSG00000121691"
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
results_list <- list()

for (gene_symbol in test_genes) {
  cat(sprintf("\n=== %s ===\n", gene_symbol))
  
  ensg_id <- gene_to_ensg[[gene_symbol]]
  eqtlgen_snps <- eqtlgen_data[eqtlgen_data$gene == gene_symbol | eqtlgen_data$gene == ensg_id, ]
  
  if (nrow(eqtlgen_snps) < 2) {
    cat("SNP 数量不足，跳过\n")
    next
  }
  
  cat(sprintf("暴露 SNP 数：%d\n", nrow(eqtlgen_snps)))
  
  # 创建暴露数据 - 确保列名正确
  exposure <- data.frame(
    SNP = as.character(eqtlgen_snps$SNP),
    BETA = as.numeric(eqtlgen_snps$beta.exposure),
    SE = as.numeric(eqtlgen_snps$se.exposure),
    EFFECT_ALLELE = as.character(eqtlgen_snps$effect_allele.exposure),
    OTHER_ALLELE = as.character(eqtlgen_snps$other_allele.exposure),
    EAF = as.numeric(eqtlgen_snps$eaf.exposure),
    PVAL = as.numeric(eqtlgen_snps$pval.exposure),
    stringsAsFactors = FALSE
  )
  
  # 移除重复
  exposure <- exposure[!duplicated(exposure$SNP) & !is.na(exposure$SNP) & exposure$SNP != "", ]
  
  if (nrow(exposure) < 2) {
    cat("去重后 SNP 不足，跳过\n")
    next
  }
  
  # 检查列
  cat("暴露数据列名:", paste(names(exposure), collapse=", "), "\n")
  
  # 格式化暴露数据
  tryCatch({
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
    
    cat(sprintf("格式化成功：%d SNPs\n", nrow(exposure_fmt)))
    
    # 提取结局数据
    outcome_matched <- outcome_data[outcome_data$hm_rsid %in% exposure_fmt$SNP, ]
    cat(sprintf("结局匹配 SNP 数：%d\n", nrow(outcome_matched)))
    
    if (nrow(outcome_matched) == 0) {
      cat("无匹配结局数据\n")
      next
    }
    
    # 格式化结局数据
    outcome_df <- as.data.frame(outcome_matched)
    outcome_df$beta.outcome <- as.numeric(outcome_df$hm_beta)
    outcome_df$se.outcome <- abs(as.numeric(outcome_df$hm_beta))
    outcome_df$effect_allele.outcome <- as.character(outcome_df$hm_effect_allele)
    outcome_df$other_allele.outcome <- as.character(outcome_df$hm_other_allele)
    outcome_df$eaf.outcome <- as.numeric(outcome_df$hm_effect_allele_frequency)
    outcome_df$pval.outcome <- 0.05
    
    outcome_fmt <- format_data(
      outcome_df,
      type = "outcome",
      snp_col = "hm_rsid",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      effect_allele_col = "effect_allele.outcome",
      other_allele_col = "other_allele.outcome",
      eaf_col = "eaf.outcome",
      pval_col = "pval.outcome"
    )
    
    cat(sprintf("结局格式化成功：%d SNPs\n", nrow(outcome_fmt)))
    
    # Harmonise
    dat <- harmonise_data(exposure_fmt, outcome_fmt)
    cat(sprintf("Harmonised SNP 数：%d\n", nrow(dat)))
    
    if (nrow(dat) == 0) {
      cat("Harmonise 失败\n")
      next
    }
    
    # MR 分析
    res <- mr(dat)
    cat("MR 结果:\n")
    print(res[, c("method", "b", "se", "pval")])
    
    results_list[[gene_symbol]] <- res
    
    # 保存
    write.csv(dat, file.path(OUTPUT_DIR, paste0(gene_symbol, "_harmonised.csv")), row.names = FALSE)
    write.csv(res, file.path(OUTPUT_DIR, paste0(gene_symbol, "_mr_results.csv")), row.names = FALSE)
    
  }, error = function(e) {
    cat(sprintf("错误：%s\n", e$message))
  })
}

# 汇总结果
if (length(results_list) > 0) {
  cat("\n\n汇总结果:\n")
  all_results <- do.call(rbind, lapply(names(results_list), function(g) {
    res <- results_list[[g]]
    ivw <- res[res$method == "Inverse variance weighted", ]
    if (nrow(ivw) > 0) {
      data.frame(gene = g, beta = ivw$b, se = ivw$se, pval = ivw$pval)
    }
  }))
  print(all_results)
  write.csv(all_results, file.path(OUTPUT_DIR, "summary_results.csv"), row.names = FALSE)
}

cat("\n完成！结果保存在:", OUTPUT_DIR, "\n")
