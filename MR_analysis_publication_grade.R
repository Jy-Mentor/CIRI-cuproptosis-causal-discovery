#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# 孟德尔随机化分析 - 9个Hub基因 vs 脑卒中风险
# 版本: 2.0 (Publication Grade)
# 日期: 2025-01-24
# 作者: 标准化MR分析流程
# 
# 核心科学假设:
# - NFKB1 (炎症核心): 高表达 → 增加卒中风险 (OR > 1)
# - FDX1 (铜死亡驱动): 高表达 → 增加卒中风险 (OR > 1) [铜死亡促进细胞死亡]
# - STAT3 (转录因子): 高表达 → 增加卒中风险 (OR > 1)
# =============================================================================

setwd("D:/EQTL")
set.seed(42)  # 确保结果可重复

# 版本控制 - 记录关键包版本
 cat("=" ,rep("=", 70), "\n", sep="")
cat("MR分析 - 9个Hub基因 vs 脑卒中 (Publication Grade v2.0)\n")
cat("=" ,rep("=", 70), "\n\n", sep="")
cat("分析日期:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("随机种子: 42\n\n")

# 加载包并记录版本
suppressPackageStartupMessages({
  if (!require("TwoSampleMR", quietly = TRUE)) install.packages("TwoSampleMR")
  if (!require("ieugwasr", quietly = TRUE)) install.packages("ieugwasr")
  if (!require("data.table", quietly = TRUE)) install.packages("data.table")
  if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!require("readr", quietly = TRUE)) install.packages("readr")
  if (!require("MRPRESSO", quietly = TRUE)) install.packages("MRPRESSO")
  
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(MRPRESSO)
})

# 记录包版本
cat("包版本信息:\n")
cat("- TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
cat("- ieugwasr:", as.character(packageVersion("ieugwasr")), "\n")
cat("- MRPRESSO:", as.character(packageVersion("MRPRESSO")), "\n\n")

# =============================================================================
# 配置
# =============================================================================
# 基因列表 - 修正为9个基因
p0_genes <- c("NFKB1", "FDX1", "STAT3")  # 核心层
p1_genes <- c("HIF1A", "HMOX1", "GPX4", "TNF", "IL6", "AGER")  # 补充层
all_genes <- c(p0_genes, p1_genes)  # 共9个基因

# 数据路径
eqtl_dir <- "D:/EQTL/rawdata"
outcome_file <- "D:/EQTL/finngen_R12_I9_STR"
ld_ref_file <- "D:/EQTL/clump/EUR"  # 1000G欧洲人群LD参考
output_dir <- paste0("D:/EQTL/MR_9genes_Results_", format(Sys.time(), "%Y%m%d"))

# 创建带日期的输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 检查LD参考文件
if (!file.exists(paste0(ld_ref_file, ".bed"))) {
  cat("警告: 本地LD参考文件未找到:", ld_ref_file, "\n")
  cat("将使用在线LD参考 (需要IEU GWAS token)\n\n")
  use_local_ld <- FALSE
} else {
  cat("使用本地LD参考:", ld_ref_file, "\n\n")
  use_local_ld <- TRUE
}

cat("=" ,rep("=", 70), "\n", sep="")
cat("分析配置:\n")
cat("- P0层基因 (核心):", paste(p0_genes, collapse = ", "), "\n")
cat("- P1层基因 (补充):", paste(p1_genes, collapse = ", "), "\n")
cat("- 总计:", length(all_genes), "个基因\n\n")

cat("科学假设:\n")
cat("- NFKB1: 高表达 → 炎症 → 卒中风险↑ (OR > 1)\n")
cat("- FDX1: 高表达 → 铜死亡 → 卒中风险↑ (OR > 1) [注意:铜死亡驱动细胞死亡]\n")
cat("- STAT3: 高表达 → 促炎 → 卒中风险↑ (OR > 1)\n\n")

cat("数据路径:\n")
cat("- eQTL目录:", eqtl_dir, "\n")
cat("- 结局数据:", outcome_file, "\n")
cat("- 输出目录:", output_dir, "\n")
cat("=" ,rep("=", 70), "\n\n", sep="")

# =============================================================================
# 步骤1: 读取结局数据 (FinnGen R12 I9_STR)
# =============================================================================
cat(rep("-", 70), "\n", sep="")
cat("步骤1: 读取结局数据\n")
cat(rep("-", 70), "\n", sep="")

if (!file.exists(outcome_file)) {
  stop("错误: 结局文件不存在: ", outcome_file)
}

# 读取FinnGen数据
# FinnGen格式: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls
outcome_data <- fread(outcome_file, 
                      select = c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt_cases", "af_alt_controls"),
                      col.names = c("SNP", "effect_allele", "other_allele", "beta.outcome", "se.outcome", 
                                   "pval.outcome", "eaf.outcome", "eaf.outcome.controls"))

cat("原始结局数据:", nrow(outcome_data), "行\n")

# 移除缺失SNP
outcome_data <- outcome_data[!is.na(SNP) & SNP != "", ]
cat("过滤后数据:", nrow(outcome_data), "行 (移除缺失SNP)\n")

# 添加必要的列
outcome_data$phenotype <- "Ischemic_Stroke_FinnGen_R12"
outcome_data$units <- "log odds"
outcome_data$samplesize.outcome <- 377277  # FinnGen R12总样本量
outcome_data$ncase.outcome <- 19862  # I9_STR病例数
outcome_data$ncontrol.outcome <- 357415  # 对照数
outcome_data$prevalence.outcome <- 0.0526  # 卒中患病率

# 确认FinnGen beta是log(OR) - FinnGen使用logistic回归
# 对于二分类结局(I9_STR)，beta已经是log(OR)
cat("数据格式确认: FinnGen I9_STR为二分类结局，beta = log(OR)\n")
cat("样本量:", format(unique(outcome_data$samplesize.outcome)[1], big.mark=","),
    "病例:", format(unique(outcome_data$ncase.outcome)[1], big.mark=","),
    "对照:", format(unique(outcome_data$ncontrol.outcome)[1], big.mark=","), "\n\n")

# =============================================================================
# MR分析主函数 (Publication Grade)
# =============================================================================
run_mr_publication_grade <- function(gene_symbol) {
  
  result_summary <- list(
    gene = gene_symbol,
    success = FALSE,
    exposure_nsnps_5e8 = 0,
    exposure_nsnps_1e5 = 0,
    clumped_nsnps = 0,
    harmonized_nsnps = 0,
    mr_ivw = NULL,
    mr_egger = NULL,
    mr_wald = NULL,
    heterogeneity = NULL,
    pleiotropy = NULL,
    steiger = NULL,
    mrpresso = NULL,
    leaveoneout = NULL,
    error = NULL
  )
  
  cat("\n" ,rep("=", 70), "\n", sep="")
  cat("分析基因:", gene_symbol, "\n")
  cat(rep("=", 70), "\n", sep="")
  
  # ---------------------------------------------------------------------------
  # 1. 读取eQTL暴露数据
  # ---------------------------------------------------------------------------
  cat("\n[1/10] 读取eQTL暴露数据...\n")
  
  # 使用Ensembl ID精确匹配
  ensembl_id <- switch(gene_symbol,
    "NFKB1" = "ENSG00000109320",
    "FDX1" = "ENSG00000137731", 
    "STAT3" = "ENSG00000168610",
    "HIF1A" = "ENSG00000100644",
    "HMOX1" = "ENSG00000100292",
    "GPX4" = "ENSG00000167468",
    "TNF" = "ENSG00000232810",
    "IL6" = "ENSG00000136244",
    "AGER" = "ENSG00000204305",
    NULL
  )
  
  if (is.null(ensembl_id)) {
    result_summary$error <- "Unknown gene - no Ensembl ID mapping"
    return(result_summary)
  }
  
  # 精确匹配eQTL文件
  eqtl_file <- file.path(eqtl_dir, paste0("eqtl-a-", ensembl_id, "_full_outcome.rds"))
  
  if (!file.exists(eqtl_file)) {
    # 尝试模糊匹配
    alt_files <- list.files(eqtl_dir, pattern = paste0(".*", gene_symbol, ".*\\.rds$"), full.names = TRUE)
    if (length(alt_files) > 0) {
      eqtl_file <- alt_files[1]
      cat("  精确匹配失败，使用模糊匹配:", basename(eqtl_file), "\n")
    } else {
      result_summary$error <- "eQTL file not found"
      return(result_summary)
    }
  }
  
  cat("  eQTL文件:", basename(eqtl_file), "\n")
  
  exposure_dat <- tryCatch(readRDS(eqtl_file), error = function(e) {
    result_summary$error <- paste("Read error:", conditionMessage(e))
    return(NULL)
  })
  
  if (is.null(exposure_dat)) return(result_summary)
  
  cat("  原始SNP数:", nrow(exposure_dat), "\n")
  
  # ---------------------------------------------------------------------------
  # 2. 筛选显著SNP (严格阈值 p < 5e-8)
  # ---------------------------------------------------------------------------
  cat("\n[2/10] 筛选显著SNP (p < 5e-8)...\n")
  
  exposure_5e8 <- exposure_dat[exposure_dat$pval < 5e-8, ]
  result_summary$exposure_nsnps_5e8 <- nrow(exposure_5e8)
  cat("  p<5e-8 SNP数:", nrow(exposure_5e8), "\n")
  
  # 如果不足，使用宽松阈值 p < 1e-5 (但有明确记录)
  if (nrow(exposure_5e8) < 1) {
    cat("  警告: 全基因组显著SNP不足，使用阈值 p < 1e-5\n")
    exposure_sig <- exposure_dat[exposure_dat$pval < 1e-5, ]
    threshold_used <- "1e-5"
  } else {
    exposure_sig <- exposure_5e8
    threshold_used <- "5e-8"
  }
  
  result_summary$exposure_nsnps_1e5 <- nrow(exposure_sig)
  cat("  使用阈值:", threshold_used, "SNP数:", nrow(exposure_sig), "\n")
  
  if (nrow(exposure_sig) == 0) {
    result_summary$error <- "No significant SNPs found"
    return(result_summary)
  }
  
  # ---------------------------------------------------------------------------
  # 3. LD Clumping (使用本地1000G参考)
  # ---------------------------------------------------------------------------
  cat("\n[3/10] LD Clumping (r2<0.001, kb=10000)...\n")
  
  exposure_clumped <- tryCatch({
    if (use_local_ld) {
      # 使用本地LD参考
      clump_data(exposure_sig, 
                 clump_kb = 10000, 
                 clump_r2 = 0.001,
                 clump_p1 = 1,
                 clump_p2 = 1,
                 bfile = ld_ref_file,
                 plink_bin = "plink")  # 需要PLINK在PATH中
    } else {
      # 使用在线LD参考
      cat("  使用在线LD参考 (需要IEU token)\n")
      clump_data(exposure_sig, 
                 clump_kb = 10000, 
                 clump_r2 = 0.001)
    }
  }, error = function(e) {
    cat("  Clumping错误:", conditionMessage(e), "\n")
    cat("  尝试不clumping继续 (结果可能受LD影响)\n")
    return(exposure_sig)
  })
  
  result_summary$clumped_nsnps <- nrow(exposure_clumped)
  cat("  Clumping后SNP数:", nrow(exposure_clumped), "\n")
  
  # ---------------------------------------------------------------------------
  # 4. 提取结局数据
  # ---------------------------------------------------------------------------
  cat("\n[4/10] 提取结局数据...\n")
  
  outcome_dat <- outcome_data[outcome_data$SNP %in% exposure_clumped$SNP, ]
  cat("  匹配结局SNP数:", nrow(outcome_dat), "\n")
  
  if (nrow(outcome_dat) == 0) {
    result_summary$error <- "No overlapping SNPs"
    return(result_summary)
  }
  
  # ---------------------------------------------------------------------------
  # 5. 数据格式化
  # ---------------------------------------------------------------------------
  cat("\n[5/10] 数据格式化...\n")
  
  # 暴露数据格式化
  exposure_formatted <- format_data(
    exposure_clumped,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    phenotype_col = "Phenotype",
    samplesize_col = "samplesize"
  )
  
  # 结局数据格式化
  outcome_formatted <- format_data(
    outcome_dat,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    phenotype_col = "phenotype",
    samplesize_col = "samplesize.outcome",
    ncase_col = "ncase.outcome",
    ncontrol_col = "ncontrol.outcome"
  )
  
  # ---------------------------------------------------------------------------
  # 6. Harmonization
  # ---------------------------------------------------------------------------
  cat("\n[6/10] 数据Harmonization...\n")
  
  dat <- harmonise_data(
    exposure_dat = exposure_formatted,
    outcome_dat = outcome_formatted,
    action = 2  # 尝试推断回文SNP
  )
  
  dat <- dat[dat$mr_keep == TRUE, ]
  result_summary$harmonized_nsnps <- nrow(dat)
  cat("  Harmonization后有效SNP:", nrow(dat), "\n")
  
  if (nrow(dat) == 0) {
    result_summary$error <- "No SNPs after harmonization"
    return(result_summary)
  }
  
  # ---------------------------------------------------------------------------
  # 7. Steiger Filtering (检验方向性)
  # ---------------------------------------------------------------------------
  cat("\n[7/10] Steiger Filtering (方向性检验)...\n")
  
  steiger_result <- tryCatch({
    dat_steiger <- steiger_filtering(dat)
    n_steiger_pass <- sum(dat_steiger$steiger_dir == TRUE, na.rm = TRUE)
    cat("  通过Steiger检验的SNP:", n_steiger_pass, "/", nrow(dat_steiger), "\n")
    result_summary$steiger <- dat_steiger
    dat_steiger
  }, error = function(e) {
    cat("  Steiger检验错误:", conditionMessage(e), "\n")
    return(dat)
  })
  
  # 使用通过Steiger检验的SNP继续分析
  if (!is.null(result_summary$steiger)) {
    dat <- dat[dat$steiger_dir == TRUE | is.na(dat$steiger_dir), ]
    cat("  继续分析SNP数:", nrow(dat), "\n")
  }
  
  # ---------------------------------------------------------------------------
  # 8. MR分析 (优先IVW)
  # ---------------------------------------------------------------------------
  cat("\n[8/10] MR分析...\n")
  
  # 方法选择: SNP>=3用IVW，否则用Wald ratio
  if (nrow(dat) >= 3) {
    methods_list <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
    cat("  使用多SNP方法: IVW + Egger + Weighted Median\n")
  } else {
    methods_list <- c("mr_wald_ratio")
    cat("  使用单SNP方法: Wald ratio\n")
  }
  
  res <- mr(dat, method_list = methods_list)
  
  # 提取主结果 (优先IVW)
  if ("Inverse variance weighted" %in% res$method) {
    main_res <- res[res$method == "Inverse variance weighted", ]
    result_summary$mr_ivw <- main_res
  } else {
    main_res <- res[res$method == "Wald ratio", ]
    result_summary$mr_wald <- main_res
  }
  
  cat("\n  MR结果:\n")
  print(res[, c("method", "nsnp", "b", "se", "pval")])
  
  # ---------------------------------------------------------------------------
  # 9. 敏感性分析
  # ---------------------------------------------------------------------------
  cat("\n[9/10] 敏感性分析...\n")
  
  # 9a. 异质性检验
  if (nrow(dat) >= 3) {
    het <- mr_heterogeneity(dat)
    result_summary$heterogeneity <- het
    cat("  异质性检验:\n")
    print(het[, c("method", "Q", "Q_df", "Q_pval")])
  }
  
  # 9b. 多效性检验 (MR-Egger intercept)
  if (nrow(dat) >= 4) {
    pleio <- mr_pleiotropy_test(dat)
    result_summary$pleiotropy <- pleio
    cat("\n  MR-Egger Intercept:", round(pleio$egger_intercept, 4), 
        "(p =", format.pval(pleio$pval, eps = 0.001), ")\n")
  }
  
  # 9c. MR-PRESSO
  if (nrow(dat) >= 4) {
    cat("\n  MR-PRESSO分析...\n")
    mrpresso <- tryCatch({
      mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                data = dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE)
    }, error = function(e) {
      cat("    MR-PRESSO错误:", conditionMessage(e), "\n")
      return(NULL)
    })
    result_summary$mrpresso <- mrpresso
    if (!is.null(mrpresso)) {
      cat("    异常值检测完成\n")
    }
  }
  
  # 9d. 留一法分析
  if (nrow(dat) > 3) {
    cat("\n  留一法敏感性分析...\n")
    loo <- mr_leaveoneout(dat)
    result_summary$leaveoneout <- loo
    cat("    完成 (", nrow(loo), " iterations)\n")
  }
  
  # ---------------------------------------------------------------------------
  # 10. 保存结果
  # ---------------------------------------------------------------------------
  cat("\n[10/10] 保存结果...\n")
  
  write.csv(res, file.path(output_dir, paste0("MR_result_", gene_symbol, ".csv")), row.names = FALSE)
  
  # 保存敏感性分析结果
  if (!is.null(result_summary$heterogeneity)) {
    write.csv(result_summary$heterogeneity, 
              file.path(output_dir, paste0("Heterogeneity_", gene_symbol, ".csv")), row.names = FALSE)
  }
  
  result_summary$success <- TRUE
  cat("  完成!\n")
  
  return(result_summary)
}

# =============================================================================
# 运行所有基因分析
# =============================================================================
cat("\n" ,rep("=", 70), "\n", sep="")
cat("开始MR分析 (共", length(all_genes), "个基因)\n")
cat(rep("=", 70), "\n\n", sep="")

all_results <- list()

for (gene in all_genes) {
  result <- run_mr_publication_grade(gene)
  all_results[[gene]] <- result
  cat("\n")
}

# =============================================================================
# 汇总结果
# =============================================================================
cat("\n" ,rep("=", 70), "\n", sep="")
cat("MR分析汇总\n")
cat(rep("=", 70), "\n\n", sep="")

summary_list <- list()

for (gene in names(all_results)) {
  res <- all_results[[gene]]
  
  if (res$success && (!is.null(res$mr_ivw) || !is.null(res$mr_wald))) {
    
    # 优先使用IVW结果
    if (!is.null(res$mr_ivw)) {
      main_res <- res$mr_ivw
      method_used <- "IVW"
    } else {
      main_res <- res$mr_wald
      method_used <- "Wald ratio"
    }
    
    or <- exp(main_res$b)
    or_lower <- exp(main_res$b - 1.96 * main_res$se)
    or_upper <- exp(main_res$b + 1.96 * main_res$se)
    
    # 科学预期判断 (修正版)
    # NFKB1, FDX1, STAT3 都是促炎/促死亡基因，预期 OR > 1
    expectation <- ""
    if (gene %in% c("NFKB1", "FDX1", "STAT3", "TNF", "IL6")) {
      # 这些是促炎/促死亡基因，高表达应增加风险
      expectation <- ifelse(or > 1, "✓ 符合预期 (风险性)", "✗ 与预期相反")
    } else if (gene %in% c("GPX4", "HMOX1")) {
      # 这些是抗氧化基因，高表达应降低风险
      expectation <- ifelse(or < 1, "✓ 符合预期 (保护性)", "✗ 与预期相反")
    }
    
    summary_list[[gene]] <- data.frame(
      Gene = gene,
      Layer = ifelse(gene %in% p0_genes, "P0 (Core)", "P1 (Supplementary)"),
      Threshold = ifelse(res$exposure_nsnps_5e8 >= 1, "5e-8", "1e-5"),
      NSNP_Exp = res$exposure_nsnps_1e5,
      NSNP_Clump = res$clumped_nsnps,
      NSNP_Final = res$harmonized_nsnps,
      Beta = main_res$b,
      SE = main_res$se,
      OR = or,
      OR_95L = or_lower,
      OR_95U = or_upper,
      P_value = main_res$pval,
      Method = method_used,
      Expectation = expectation,
      stringsAsFactors = FALSE
    )
  } else {
    summary_list[[gene]] <- data.frame(
      Gene = gene,
      Layer = ifelse(gene %in% p0_genes, "P0 (Core)", "P1 (Supplementary)"),
      Threshold = NA,
      NSNP_Exp = res$exposure_nsnps_1e5,
      NSNP_Clump = res$clumped_nsnps,
      NSNP_Final = 0,
      Beta = NA, SE = NA, OR = NA, OR_95L = NA, OR_95U = NA,
      P_value = NA,
      Method = "Failed",
      Expectation = res$error,
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

cat("汇总结果:\n")
print(summary_df[, c("Gene", "OR", "P_value", "NSNP_Final", "Method", "Expectation")])

write.csv(summary_df, file.path(output_dir, "MR_summary_9genes.csv"), row.names = FALSE)
cat("\n汇总结果已保存:\n")
cat("-", file.path(output_dir, "MR_summary_9genes.csv"), "\n")

# 显著结果
cat("\n")
cat(rep("-", 70), "\n", sep="")
cat("显著结果 (P < 0.05):\n")
cat(rep("-", 70), "\n", sep="")

sig_results <- summary_df[!is.na(summary_df$P_value) & summary_df$P_value < 0.05, ]
if (nrow(sig_results) > 0) {
  for (i in 1:nrow(sig_results)) {
    direction <- ifelse(sig_results$OR[i] > 1, "风险性", "保护性")
    cat(sig_results$Gene[i], 
        ": OR =", round(sig_results$OR[i], 3),
        "[", round(sig_results$OR_95L[i], 3), "-", round(sig_results$OR_95U[i], 3), "]",
        "(", direction, "), P =", format.pval(sig_results$P_value[i], eps = 0.001),
        sig_results$Expectation[i], "\n")
  }
} else {
  cat("未发现显著关联 (P < 0.05)\n")
}

# P0层验证
cat("\n")
cat(rep("-", 70), "\n", sep="")
cat("P0层核心基因验证:\n")
cat(rep("-", 70), "\n", sep="")

p0_results <- summary_df[summary_df$Layer == "P0 (Core)", ]
for (i in 1:nrow(p0_results)) {
  gene <- p0_results$Gene[i]
  or <- p0_results$OR[i]
  pval <- p0_results$P_value[i]
  
  if (!is.na(or)) {
    cat(gene, ": OR =", round(or, 3), 
        "[", round(p0_results$OR_95L[i], 3), "-", round(p0_results$OR_95U[i], 3), "]",
        ", P =", format.pval(pval, eps = 0.001),
        p0_results$Expectation[i], "\n")
  } else {
    cat(gene, ": 分析失败 (", p0_results$Expectation[i], ")\n")
  }
}

# =============================================================================
# 标准MR森林图 (Publication Grade)
# =============================================================================
cat("\n")
cat(rep("-", 70), "\n", sep="")
cat("生成标准MR森林图\n")
cat(rep("-", 70), "\n", sep="")

if (nrow(summary_df[!is.na(summary_df$OR), ]) > 0) {
  
  plot_data <- summary_df[!is.na(summary_df$OR), ]
  plot_data$Significance <- ifelse(plot_data$P_value < 0.05, "Significant", 
                                   ifelse(plot_data$P_value < 0.1, "Trend", "NS"))
  plot_data$Label <- paste0(
    plot_data$Gene, "\n",
    "(n=", plot_data$NSNP_Final, ", ", plot_data$Method, ")"
  )
  
  # 标准MR森林图
  p <- ggplot(plot_data, aes(y = reorder(Gene, OR))) +
    # 效应点
    geom_point(aes(x = OR, color = Significance, size = -log10(P_value + 1e-10)), 
               shape = 15, stroke = 1) +
    # 误差线
    geom_errorbarh(aes(xmin = OR_95L, xmax = OR_95U, color = Significance), 
                   height = 0.2, linewidth = 0.8) +
    # 参考线 OR=1
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    # 颜色设置
    scale_color_manual(values = c("Significant" = "#D55E00", "Trend" = "#E69F00", "NS" = "#999999"),
                       name = "Significance") +
    scale_size_continuous(name = "-log10(P)", range = c(2, 6)) +
    # 坐标轴
    scale_x_log10(breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2),
                  labels = c("0.50", "0.75", "1.00", "1.25", "1.50", "2.00")) +
    # 标签
    labs(
      title = "Mendelian Randomization: 9 Hub Genes and Ischemic Stroke Risk",
      subtitle = paste0("Data source: GTEx eQTL + FinnGen R12 (n=", 
                        format(unique(outcome_data$ncase.outcome)[1], big.mark=","), " cases)"),
      x = "Odds Ratio (95% CI)",
      y = "Gene (n SNP, Method)",
      caption = paste0("Analysis date: ", format(Sys.time(), "%Y-%m-%d"), 
                       " | Seed: 42 | Threshold: 5e-8/1e-5")
    ) +
    # 主题
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      plot.caption = element_text(size = 9, color = "gray50")
    ) +
    # 添加OR数值标注
    geom_text(aes(x = OR_95U + 0.05, 
                  label = sprintf("%.2f [%.2f-%.2f]", OR, OR_95L, OR_95U)),
              hjust = 0, size = 3, color = "black")
  
  # 保存
  ggsave(file.path(output_dir, "MR_forest_plot_publication.png"), 
         p, width = 12, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "MR_forest_plot_publication.pdf"), 
         p, width = 12, height = 10)
  
  cat("森林图已保存:\n")
  cat("- PNG:", file.path(output_dir, "MR_forest_plot_publication.png"), "\n")
  cat("- PDF:", file.path(output_dir, "MR_forest_plot_publication.pdf"), "\n")
}

# =============================================================================
# 完成
# =============================================================================
cat("\n")
cat("=" ,rep("=", 70), "\n", sep="")
cat("分析完成!\n")
cat(rep("=", 70), "\n", sep="")
cat("输出目录:", output_dir, "\n")
cat("所有结果文件已保存 (带日期标记，可追溯)\n")
cat(rep("=", 70), "\n", sep="")

# 创建分析日志
log_file <- file.path(output_dir, "analysis_log.txt")
sink(log_file)
cat("MR Analysis Log\n")
cat("===============\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Seed: 42\n")
cat("R Version:", R.version.string, "\n")
cat("TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
cat("\nGenes analyzed:", paste(all_genes, collapse = ", "), "\n")
cat("\nSummary:\n")
print(summary_df[, c("Gene", "OR", "P_value", "Method")])
sink()
