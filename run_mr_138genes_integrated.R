#!/usr/bin/env Rscript
# ================================================================================
# 138 基因 MR 分析（使用整合的 eQTL 数据：eQTLGen + GTEx 全血 + GTEx 脑组织）
# ================================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

cat("======================================================================\n")
cat("138 基因 MR 分析（使用整合的 eQTL 数据）\n")
cat("======================================================================\n\n")

# 配置
CLUMPED_FILE <- "D:/EQTL/clump/merged_eqtl_all_sources_p_1e-05_simple.tsv"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_138genes_integrated"
MEGASTROKE_FILE <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 加载数据
cat("加载整合的 eQTL 数据...\n")
eqtl_data <- fread(CLUMPED_FILE, sep = "\t")
cat(sprintf("  总 SNP 数：%d\n", nrow(eqtl_data)))

cat("加载 MEGASTROKE 数据...\n")
outcome_data <- fread(MEGASTROKE_FILE, sep = "\t", stringsAsFactors = FALSE)
outcome_data <- outcome_data[!is.na(outcome_data$hm_rsid) & outcome_data$hm_rsid != "", ]
cat(sprintf("  总 SNP 数：%d\n\n", nrow(outcome_data)))

# 138 个目标基因
target_genes <- c(
  "LYN", "PRKCQ", "NMT1", "TDP1", "MAN2B1", "IL10RA", "RHOC", "SREBF1",
  "KCNA5", "HIF1A", "CTSC", "CAT", "FABP4", "STAT5A", "FABP2", "B2M",
  "RBM39", "HBS1L", "CHFR", "NUDCD2", "TCN2", "SCN9A", "JAK1", "GPX1",
  "CTSB", "CASP8", "FABP5", "XDH", "MB", "POLR2D", "HSD17B10", "MAPKAPK2",
  "SEC13", "PCTP", "ZEB1", "RELA", "IRF1", "GFAP", "CPT2", "BRD3",
  "NR3C1", "F3", "C3", "ITGA1", "CITED2", "HIBADH", "SAT2", "TSPO",
  "PTGS1", "IMPDH2", "FLT4", "CPT1A", "AKT1", "CCR5", "PTPRF", "HPGDS",
  "PTPRJ", "CASK", "MGAT1", "IGFBP2", "TOP2A", "PPARG", "IL6", "EPHX1",
  "CP", "AIF1", "PLA2G4A", "ALDH9A1", "S100A6", "DDC", "CUL4B", "BST1",
  "CNDP2", "TNF", "PARP1", "IKBKB", "EGFR", "COL1A1", "ADRB1", "SPHK1",
  "GCH1", "ACADVL", "STARD13", "CTSD", "PDCD6IP", "PTPRC", "TGFB1", "PABPC1",
  "HTR2C", "CTSS", "CNR2", "ACTA2", "FNTA", "RENBP", "CCNA2", "PTGR1",
  "LEF1", "SAT1", "XRCC6", "TBXAS1", "NR1H3", "HTR2B", "CTSL", "CDK4",
  "CXCR3", "TIMP1", "OAZ1", "STK4", "ZHX2", "MKNK2", "SERPINB10", "ACADM",
  "STAT3", "NFKB1", "HSPA5", "CTSK", "CCND1", "PTPN2", "PTPN6", "PA2G4",
  "HSD17B4", "ACAD11", "PDCD6", "PARP12", "SERPINB1", "STAT1", "NFE2L2",
  "HMOX1", "CTSF", "CCL2", "MAOB", "ICAM1", "FDX1", "LIAS", "LIPT1",
  "DLAT", "PDHB", "PDHX", "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
  "HIF1A", "MTOR", "NFKB1", "GPX4"
)

target_genes <- unique(target_genes)
cat(sprintf("目标基因数：%d\n\n", length(target_genes)))

# 统计每个基因的 SNP 数量
cat("统计每个基因的 SNP 数量...\n")
gene_snp_counts <- list()
for (gene in target_genes) {
  snps <- eqtl_data[eqtl_data$gene == gene, ]
  gene_snp_counts[[gene]] <- nrow(snps)
}

snp_counts_df <- data.frame(
  gene = names(gene_snp_counts),
  snp_count = unname(unlist(gene_snp_counts)),
  stringsAsFactors = FALSE
)

cat(sprintf("有 eQTL 的基因数：%d\n", sum(snp_counts_df$snp_count > 0)))
cat(sprintf("平均每个基因的 SNP 数：%.2f\n", mean(snp_counts_df$snp_count[snp_counts_df$snp_count > 0])))
cat(sprintf("中位数：%.0f\n", median(snp_counts_df$snp_count[snp_counts_df$snp_count > 0])))
cat(sprintf("范围：%d - %d\n\n", min(snp_counts_df$snp_count[snp_counts_df$snp_count > 0]), max(snp_counts_df$snp_count)))

# 处理每个基因
results_list <- list()
skipped <- data.frame(gene = character(), reason = character(), stringsAsFactors = FALSE)
failed <- data.frame(gene = character(), reason = character(), n_snps = integer(), stringsAsFactors = FALSE)

cat("开始处理基因...\n\n")

for (i in seq_along(target_genes)) {
  gene_symbol <- target_genes[i]
  cat(sprintf("[%d/%d] %s: ", i, length(target_genes), gene_symbol))
  
  # 获取 eQTL 数据
  gene_eqtl <- eqtl_data[eqtl_data$gene == gene_symbol, ]
  
  if (nrow(gene_eqtl) < 2) {
    cat(sprintf("无 eQTL 数据，跳过\n"))
    skipped <- rbind(skipped, data.frame(gene = gene_symbol, reason = "无 eQTL 数据", stringsAsFactors = FALSE))
    next
  }
  
  # 直接使用 eQTL 数据（已经是 TwoSampleMR 格式）
  exposure <- gene_eqtl[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", 
                            "other_allele.exposure", "eaf.exposure", "pval.exposure")]
  
  exposure <- exposure[!duplicated(exposure$SNP) & !is.na(exposure$SNP) & exposure$SNP != "", ]
  
  if (nrow(exposure) < 2) {
    cat("去重后 SNP 不足，跳过\n")
    skipped <- rbind(skipped, data.frame(gene = gene_symbol, reason = "去重后 SNP 不足", stringsAsFactors = FALSE))
    next
  }
  
  # 暴露数据已经是正确格式，直接使用
  exposure_fmt <- as.data.frame(exposure)
  exposure_fmt$type <- "exposure"
  exposure_fmt$mr_keep.exposure <- TRUE
  exposure_fmt$id.exposure <- "Integrated_eQTL"
  exposure_fmt$exposure <- gene_symbol
  
  # 提取结局数据（使用 rsID 匹配）
  outcome_matched <- outcome_data[outcome_data$hm_rsid %in% exposure_fmt$SNP, ]
  
  if (nrow(outcome_matched) == 0) {
    cat("无匹配结局数据，跳过\n")
    skipped <- rbind(skipped, data.frame(gene = gene_symbol, reason = "无匹配结局数据", stringsAsFactors = FALSE))
    next
  }
  
  # 格式化结局数据
  outcome_df <- as.data.frame(outcome_matched)
  outcome_df$beta.outcome <- as.numeric(outcome_df$hm_beta)
  outcome_df$se.outcome <- as.numeric(outcome_df$standard_error)
  outcome_df$effect_allele.outcome <- as.character(outcome_df$hm_effect_allele)
  outcome_df$other_allele.outcome <- as.character(outcome_df$hm_other_allele)
  outcome_df$eaf.outcome <- as.numeric(outcome_df$hm_effect_allele_frequency)
  outcome_df$pval.outcome <- as.numeric(outcome_df$p_value)
  
  # 结局数据格式化
  outcome_fmt <- as.data.frame(outcome_df[, c("hm_rsid", "beta.outcome", "se.outcome", 
                                               "effect_allele.outcome", "other_allele.outcome", 
                                               "eaf.outcome", "pval.outcome")])
  names(outcome_fmt)[1] <- "SNP"
  outcome_fmt$type <- "outcome"
  outcome_fmt$mr_keep.outcome <- TRUE
  outcome_fmt$id.outcome <- "MEGASTROKE"
  outcome_fmt$outcome <- "Stroke"
  
  # Harmonise
  dat <- tryCatch({
    harmonise_data(exposure_fmt, outcome_fmt)
  }, error = function(e) {
    cat(sprintf("Harmonise 失败：%s\n", e$message))
    failed <- rbind(failed, data.frame(gene = gene_symbol, reason = paste("Harmonise:", e$message), n_snps = nrow(exposure_fmt), stringsAsFactors = FALSE))
    return(NULL)
  })
  
  if (is.null(dat) || nrow(dat) == 0) {
    cat("Harmonise 失败\n")
    next
  }
  
  # MR 分析
  res <- tryCatch({
    mr(dat, method_list = c("mr_ivw", "mr_weighted_median"))
  }, error = function(e) {
    cat(sprintf("MR 分析失败：%s\n", e$message))
    failed <- rbind(failed, data.frame(gene = gene_symbol, reason = paste("MR:", e$message), n_snps = nrow(dat), stringsAsFactors = FALSE))
    return(NULL)
  })
  
  if (is.null(res) || nrow(res) == 0) {
    cat("MR 分析失败\n")
    next
  }
  
  # 提取 IVW 结果
  ivw_res <- res[res$method == "Inverse variance weighted", ]
  
  if (nrow(ivw_res) == 0) {
    cat("无 IVW 结果\n")
    next
  }
  
  # 计算 F 统计量
  f_stat <- mean((exposure_fmt$beta.exposure / exposure_fmt$se.exposure)^2, na.rm = TRUE)
  
  # 创建结果
  result <- data.frame(
    gene = gene_symbol,
    beta = ivw_res$b[1],
    se = ivw_res$se[1],
    or = exp(ivw_res$b[1]),
    pval = ivw_res$pval[1],
    f_stat = f_stat,
    nsnp = nrow(dat),
    stringsAsFactors = FALSE
  )
  
  results_list[[gene_symbol]] <- result
  
  # 保存
  write.csv(dat, file.path(OUTPUT_DIR, paste0(gene_symbol, "_harmonised.csv")), row.names = FALSE)
  write.csv(res, file.path(OUTPUT_DIR, paste0(gene_symbol, "_mr_results.csv")), row.names = FALSE)
  
  cat(sprintf("✓ OR=%.3f, P=%.3f, F=%.2f\n", result$or, result$pval, f_stat))
}

# 汇总结果
cat("\n\n汇总结果:\n")
if (length(results_list) > 0) {
  all_results <- do.call(rbind, results_list)
  all_results <- all_results[order(all_results$pval), ]
  
  # FDR 校正
  all_results$fdr <- p.adjust(all_results$pval, method = "fdr")
  
  write.csv(all_results, file.path(OUTPUT_DIR, "summary_results.csv"), row.names = FALSE)
  write.csv(all_results, file.path(OUTPUT_DIR, "mr_results_fdr_corrected.csv"), row.names = FALSE)
  
  if (nrow(skipped) > 0) {
    write.csv(skipped, file.path(OUTPUT_DIR, "skipped_genes.csv"), row.names = FALSE)
  }
  
  if (nrow(failed) > 0) {
    write.csv(failed, file.path(OUTPUT_DIR, "failed_genes.csv"), row.names = FALSE)
  }
  
  cat(sprintf("成功分析：%d 个基因\n", nrow(all_results)))
  cat(sprintf("跳过：%d 个基因\n", nrow(skipped)))
  cat(sprintf("失败：%d 个基因\n", nrow(failed)))
  
  cat("\nTop 10 最显著的基因:\n")
  print(head(all_results[, c("gene", "beta", "or", "pval", "fdr", "f_stat", "nsnp")], 10))
  
  if (any(all_results$pval < 0.05)) {
    cat("\n显著基因 (P < 0.05):\n")
    print(all_results[all_results$pval < 0.05, c("gene", "beta", "or", "pval", "fdr", "f_stat", "nsnp")])
  }
  
  if (any(all_results$fdr < 0.05)) {
    cat("\nFDR 显著基因 (FDR < 0.05):\n")
    print(all_results[all_results$fdr < 0.05, c("gene", "beta", "or", "pval", "fdr", "f_stat", "nsnp")])
  }
} else {
  cat("无成功分析的基因\n")
}

cat("\n完成！结果保存在:", OUTPUT_DIR, "\n")
