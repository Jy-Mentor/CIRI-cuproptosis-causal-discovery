#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# Mendelian Randomization: 9 Hub Genes vs Ischemic Stroke
# 版本: 3.0 (Publication Grade, 顶刊级)
# 日期: 2026-04-24
# 分析者: 自动生成
# 
# 基因列表 (9个):
# P0层 (核心): NFKB1, FDX1, STAT3
# P1层 (补充): HIF1A, HMOX1, GPX4, TNF, IL6, AGER
#
# 数据来源:
# 暴露: eQTLgen (血液eQTL)
# 结局: FinnGen R12 I9_STR (Ischemic Stroke)
# =============================================================================

setwd("D:/EQTL")
set.seed(42)

# =============================================================================
# 0. 包加载与版本检查
# =============================================================================
cat("======================================================================\n", sep = "")
cat("MR Analysis: 9 Hub Genes vs Ischemic Stroke (v3.0)\n")
cat("======================================================================\n\n", sep = "")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Random Seed: 42\n\n")

suppressPackageStartupMessages({
  if (!require("TwoSampleMR", quietly = TRUE)) stop("TwoSampleMR not installed")
  if (!require("ieugwasr", quietly = TRUE)) stop("ieugwasr not installed")
  if (!require("readxl", quietly = TRUE)) stop("readxl not installed")
  if (!require("data.table", quietly = TRUE)) stop("data.table not installed")
  if (!require("dplyr", quietly = TRUE)) stop("dplyr not installed")
  if (!require("ggplot2", quietly = TRUE)) stop("ggplot2 not installed")
  library(TwoSampleMR)
  library(ieugwasr)
  library(readxl)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# 打印包版本 (确保可重复性)
cat("Package Versions:\n")
cat("  R:", R.version.string, "\n")
cat("  TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
cat("  ieugwasr:", as.character(packageVersion("ieugwasr")), "\n\n")

# =============================================================================
# 1. 配置与参数
# =============================================================================
# 基因列表
p0_genes <- c("NFKB1", "FDX1", "STAT3")  # 核心层
p1_genes <- c("HIF1A", "HMOX1", "GPX4", "TNF", "IL6", "AGER")  # 补充层
all_genes <- c(p0_genes, p1_genes)  # 共9个基因

# 数据路径 - 使用顶刊标准配置
# 首选: p=5e-08, kb=1000, r2=0.001 (最严格, 顶刊标准)
eqtl_primary <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-08_kb_1000_r2_0.001.xlsx"
# 兜底1: 仅放宽LD到0.01
eqtl_fallback1 <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-08_kb_1000_r2_0.01.xlsx"
# 兜底2: 放宽p值到5e-06
eqtl_fallback2 <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-06_kb_1000_r2_0.001.xlsx"

outcome_file <- "D:/EQTL/finngen_R12_I9_STR"  # 纯文本TSV文件（无扩展名）

# 输出目录（带日期版本）
output_dir <- paste0("D:/EQTL/MR_v3_", format(Sys.time(), "%Y%m%d_%H%M%S"))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Configuration:\n")
cat("  Total Genes:", length(all_genes), "\n")
cat("  P0 Layer:", paste(p0_genes, collapse = ", "), "\n")
cat("  P1 Layer:", paste(p1_genes, collapse = ", "), "\n")
cat("  eQTL:", basename(eqtl_primary), "\n")
cat("  Outcome:", basename(outcome_file), "\n")
cat("  Output:", output_dir, "\n\n")

# 科学假设
cat("Scientific Hypothesis:\n")
cat("  NFKB1, FDX1, STAT3, TNF, IL6: High expr -> Risk (OR > 1)\n")
cat("  GPX4, HMOX1: High expr -> Protective (OR < 1)\n")
cat("  AGER (RAGE): High expr -> Risk (OR > 1) [Pro-inflammatory]\n\n")

# =============================================================================
# 2. 读取数据（带格式验证）
# =============================================================================
cat("----------------------------------------------------------------------\n", sep = "")
cat("Step 1: Reading Data\n")
cat("----------------------------------------------------------------------\n", sep = "")

# 2a. 读取eQTL数据（带兜底逻辑）
cat("  Selecting eQTL file...\n")

# 检查文件是否存在
eqtl_file <- NULL
for (candidate in list(
  list(file = eqtl_primary, label = "Primary (p<5e-08, kb=1000, r2=0.001)"),
  list(file = eqtl_fallback1, label = "Fallback1 (p<5e-08, kb=1000, r2=0.01)"),
  list(file = eqtl_fallback2, label = "Fallback2 (p<5e-06, kb=1000, r2=0.001)")
)) {
  if (file.exists(candidate$file)) {
    eqtl_file <- candidate$file
    cat("  Using:", candidate$label, "\n")
    cat("  File:", basename(eqtl_file), "\n")
    break
  }
}

if (is.null(eqtl_file)) {
  stop("No eQTL file found! Checked:\n", 
       paste(c(eqtl_primary, eqtl_fallback1, eqtl_fallback2), collapse = "\n"))
}

eqtl_all <- read_excel(eqtl_file)

# 检查必需列
required_eqtl_cols <- c("gene", "SNP", "beta.exposure", "se.exposure", "pval.exposure",
                        "eaf.exposure", "effect_allele.exposure", "other_allele.exposure")
missing_cols <- setdiff(required_eqtl_cols, names(eqtl_all))
if (length(missing_cols) > 0) stop("Missing eQTL columns: ", paste(missing_cols, collapse = ", "))

cat("  eQTL: ", nrow(eqtl_all), "SNPs,", length(unique(eqtl_all$gene)), "genes\n")

# 检查目标基因
found_genes <- intersect(all_genes, unique(eqtl_all$gene))
missing_genes <- setdiff(all_genes, found_genes)
if (length(missing_genes) > 0) {
  cat("  WARNING: Missing genes:", paste(missing_genes, collapse = ", "), "\n")
}
all_genes <- found_genes
p0_genes <- intersect(p0_genes, found_genes)
p1_genes <- intersect(p1_genes, found_genes)

# 2b. 读取结局数据 (FinnGen)
if (!file.exists(outcome_file)) stop("Outcome file not found: ", outcome_file)
cat("  Reading outcome (fread)...\n")

outcome_data <- fread(outcome_file, select = c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt"))
setnames(outcome_data, 
         old = c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt"),
         new = c("SNP", "effect_allele.outcome", "other_allele.outcome", 
                 "beta.outcome", "se.outcome", "pval.outcome", "eaf.outcome"))

# 从数据读取样本量 (FinnGen R12 I9_STR 元数据)
# 来源: https://r12.finnngen.fi/pheno/I9_STR
outcome_data$samplesize.outcome <- 377277  # Total samples
outcome_data$ncase.outcome <- 19862  # Cases
outcome_data$ncontrol.outcome <- 357415  # Controls
outcome_data$phenotype <- "Ischemic_Stroke_FinnGen_R12_I9_STR"

cat("  Outcome:", nrow(outcome_data), "SNPs\n")
cat("  Sample size:", format(outcome_data$samplesize.outcome[1], big.mark = ","), 
    "(Cases:", format(outcome_data$ncase.outcome[1], big.mark = ","), ")\n\n")

# =============================================================================
# 3. MR分析主函数 (顶刊级)
# =============================================================================
run_mr_publication <- function(gene_symbol) {
  
  result <- list(
    gene = gene_symbol,
    success = FALSE,
    nsnp_exp = 0, nsnp_outcome = 0, nsnp_final = 0,
    beta = NA, se = NA, or = NA, or_l = NA, or_u = NA, p = NA, method = NULL,
    steiger_pass = 0, steiger_total = 0,
    heterogeneity_q = NA, heterogeneity_p = NA,
    egger_intercept = NA, egger_intercept_p = NA,
    leaveoneout_range = NA,
    mrpresso_outliers = 0,
    error = NULL
  )
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Gene:", gene_symbol, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  tryCatch({
    # 1. 提取eQTL
    cat("\n[1/7] Extracting eQTL...\n")
    exposure <- eqtl_all[eqtl_all$gene == gene_symbol, ]
    result$nsnp_exp <- nrow(exposure)
    cat("  eQTL SNPs:", nrow(exposure), "\n")
    
    if (nrow(exposure) == 0) {
      result$error <- "No eQTL SNPs"
      return(result)
    }
    
    # 处理缺失eaf
    if (any(is.na(exposure$eaf.exposure))) {
      n_missing <- sum(is.na(exposure$eaf.exposure))
      cat("  WARNING: Missing eaf in", n_missing, "SNPs (filling with 0.5)\n")
      exposure$eaf.exposure[is.na(exposure$eaf.exposure)] <- 0.5
    }
    
    # 2. 匹配结局
    cat("\n[2/7] Matching outcome...\n")
    outcome <- outcome_data[outcome_data$SNP %in% exposure$SNP, ]
    result$nsnp_outcome <- nrow(outcome)
    cat("  Matched:", nrow(outcome), "\n")
    
    if (nrow(outcome) == 0) {
      result$error <- "No matching outcome SNPs"
      return(result)
    }
    
    # 处理缺失eaf
    if (any(is.na(outcome$eaf.outcome))) {
      n_missing <- sum(is.na(outcome$eaf.outcome))
      cat("  WARNING: Missing outcome eaf in", n_missing, "SNPs\n")
      outcome$eaf.outcome[is.na(outcome$eaf.outcome)] <- 0.5
    }
    
    # 3. 格式转换 (确保是data.frame)
    cat("\n[3/7] Formatting data...\n")
    exp_df <- as.data.frame(exposure)
    out_df <- as.data.frame(outcome)
    
    exp_fmt <- format_data(exp_df, type = "exposure",
                           snp_col = "SNP", beta_col = "beta.exposure", 
                           se_col = "se.exposure", pval_col = "pval.exposure",
                           eaf_col = "eaf.exposure",
                           effect_allele_col = "effect_allele.exposure",
                           other_allele_col = "other_allele.exposure")
    
    out_fmt <- format_data(out_df, type = "outcome",
                           snp_col = "SNP", beta_col = "beta.outcome", 
                           se_col = "se.outcome", pval_col = "pval.outcome",
                           eaf_col = "eaf.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome")
    
    # 4. Harmonization
    cat("\n[4/7] Harmonizing...\n")
    dat <- harmonise_data(exp_fmt, out_fmt, action = 2)
    dat <- dat[dat$mr_keep == TRUE, ]
    result$nsnp_final <- nrow(dat)
    cat("  Valid SNPs after harmonization:", nrow(dat), "\n")
    
    if (nrow(dat) == 0) {
      result$error <- "No SNPs after harmonization"
      return(result)
    }
    
    # 5. Steiger Filtering (方向性检验)
    cat("\n[5/7] Steiger filtering...\n")
    tryCatch({
      dat_steiger <- steiger_filtering(dat)
      # 检查steiger结果是否有效
      if (!is.null(dat_steiger) && nrow(dat_steiger) > 0 && "steiger_dir" %in% names(dat_steiger)) {
        n_pass <- sum(dat_steiger$steiger_dir == TRUE, na.rm = TRUE)
        result$steiger_pass <- n_pass
        result$steiger_total <- nrow(dat_steiger)
        cat("  Steiger pass:", n_pass, "/", nrow(dat_steiger), "\n")
        
        # 仅使用通过Steiger检验的SNP
        if (n_pass > 0) {
          dat <- dat_steiger[dat_steiger$steiger_dir == TRUE | is.na(dat_steiger$steiger_dir), ]
          cat("  Using", nrow(dat), "SNPs after Steiger filtering\n")
        }
      } else {
        cat("  WARNING: Steiger returned invalid results, continuing with all SNPs\n")
        result$steiger_pass <- nrow(dat)
        result$steiger_total <- nrow(dat)
      }
    }, error = function(e) {
      cat("  WARNING: Steiger error:", conditionMessage(e), "\n")
      result$steiger_pass <- nrow(dat)
      result$steiger_total <- nrow(dat)
    })
    
    if (nrow(dat) == 0) {
      result$error <- "No SNPs after Steiger"
      return(result)
    }
    
    # 6. MR分析
    cat("\n[6/7] MR analysis...\n")
    
    # 方法选择: SNP>=2用IVW，否则Wald ratio
    if (nrow(dat) >= 2) {
      methods <- c("mr_ivw")
      if (nrow(dat) >= 3) {
        methods <- c(methods, "mr_egger_regression", "mr_weighted_median")
      }
      method_name <- "IVW"
    } else {
      methods <- c("mr_wald_ratio")
      method_name <- "Wald"
    }
    
    res <- mr(dat, method_list = methods)
    
    # 提取主结果 (优先IVW)
    if ("Inverse variance weighted" %in% res$method) {
      main <- res[res$method == "Inverse variance weighted", ]
    } else {
      main <- res[res$method == "Wald ratio", ]
    }
    
    if (is.na(main$b) || is.na(main$pval)) {
      result$error <- "MR returned NA"
      return(result)
    }
    
    result$beta <- main$b
    result$se <- main$se
    result$or <- exp(main$b)
    result$or_l <- exp(main$b - 1.96 * main$se)
    result$or_u <- exp(main$b + 1.96 * main$se)
    result$p <- main$pval
    result$method <- method_name
    
    cat("  ", method_name, "OR =", round(result$or, 3),
        "[", round(result$or_l, 3), "-", round(result$or_u, 3), "]",
        "P =", format.pval(result$p, eps = 0.001), "\n")
    
    # 7. 敏感性分析
    cat("\n[7/7] Sensitivity analyses...\n")
    
    # 7a. 异质性检验
    if (nrow(dat) >= 3) {
      het <- mr_heterogeneity(dat)
      ivw_het <- het[het$method == "Inverse variance weighted", ]
      result$heterogeneity_q <- ivw_het$Q
      result$heterogeneity_p <- ivw_het$Q_pval
      cat("  Heterogeneity Q_pval =", format.pval(ivw_het$Q_pval, eps = 0.001), "\n")
    }
    
    # 7b. MR-Egger截距 (多效性检验)
    if ("mr_egger_regression" %in% res$method) {
      pleio <- mr_pleiotropy_test(dat)
      result$egger_intercept <- pleio$egger_intercept
      result$egger_intercept_p <- pleio$pval
      cat("  Egger intercept =", round(pleio$egger_intercept, 4),
          "(P =", format.pval(pleio$pval, eps = 0.001), ")\n")
    }
    
    # 7c. 留一法分析
    if (nrow(dat) > 3) {
      cat("  Leave-one-out...\n")
      loo <- mr_leaveoneout(dat)
      valid_loo <- loo[!is.na(loo$b), ]
      if (nrow(valid_loo) > 0) {
        result$leaveoneout_range <- paste0(
          "[", round(min(exp(valid_loo$b)), 3), "-", round(max(exp(valid_loo$b)), 3), "]")
        cat("  LOO OR range:", result$leaveoneout_range, "\n")
      }
    }
    
    # 保存单个基因结果
    write.csv(res, file.path(output_dir, paste0("MR_", gene_symbol, ".csv")), row.names = FALSE)
    result$success <- TRUE
    
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    result$error <- conditionMessage(e)
  })
  
  return(result)
}

# =============================================================================
# 4. 运行所有基因
# =============================================================================
cat("\n======================================================================\n", sep = "")
cat("Running MR Analysis for", length(all_genes), "Genes\n")
cat("======================================================================\n", sep = "")

all_results <- list()
for (gene in all_genes) {
  all_results[[gene]] <- run_mr_publication(gene)
}

# =============================================================================
# 5. 汇总结果
# =============================================================================
cat("\n======================================================================\n", sep = "")
cat("Summary Results\n")
cat("======================================================================\n", sep = "")

summary_list <- list()
for (gene in names(all_results)) {
  r <- all_results[[gene]]
  
  # 预期判断
  exp_text <- ""
  if (r$success && !is.na(r$or)) {
    if (gene %in% c("NFKB1", "FDX1", "STAT3", "TNF", "IL6", "AGER")) {
      exp_text <- ifelse(r$or > 1, "\u2713 Risk", "\u2717 Protective")
    } else if (gene %in% c("GPX4", "HMOX1")) {
      exp_text <- ifelse(r$or < 1, "\u2713 Protective", "\u2717 Risk")
    }
  } else {
    exp_text <- ifelse(is.null(r$error), "Failed", r$error)
  }
  
  summary_list[[gene]] <- data.frame(
    Gene = gene,
    Layer = ifelse(gene %in% p0_genes, "P0 (Core)", "P1 (Supp)"),
    Beta = ifelse(r$success, r$beta, NA),
    SE = ifelse(r$success, r$se, NA),
    OR = ifelse(r$success, r$or, NA),
    OR_L95 = ifelse(r$success, r$or_l, NA),
    OR_U95 = ifelse(r$success, r$or_u, NA),
    P = ifelse(r$success, r$p, NA),
    NSNP_eQTL = r$nsnp_exp,
    NSNP_Outcome = r$nsnp_outcome,
    NSNP_Final = r$nsnp_final,
    Steiger_Pass = ifelse(r$success, paste0(r$steiger_pass, "/", r$steiger_total), NA),
    Heterogeneity_P = ifelse(r$success, r$heterogeneity_p, NA),
    Egger_Intercept_P = ifelse(r$success, r$egger_intercept_p, NA),
    LOO_Range = ifelse(r$success, r$leaveoneout_range, NA),
    Method = ifelse(r$success, r$method, "Failed"),
    Expectation = exp_text,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

# 按层级和OR排序 (P0层优先)
summary_df$sort_key <- ifelse(summary_df$Layer == "P0 (Core)", 0, 1)
summary_df <- summary_df[order(summary_df$sort_key, -summary_df$OR, na.last = TRUE), ]
summary_df$sort_key <- NULL

# 打印汇总
print(summary_df[, c("Gene", "Layer", "OR", "P", "NSNP_Final", "Method", "Expectation")])

# 保存汇总
write.csv(summary_df, file.path(output_dir, "MR_Summary.csv"), row.names = FALSE)

# 显著结果
sig <- summary_df[!is.na(summary_df$P) & summary_df$P < 0.05, ]
if (nrow(sig) > 0) {
  cat("\n*** Significant Results (P < 0.05) ***\n")
  for (i in 1:nrow(sig)) {
    cat(sig$Gene[i], ": OR =", round(sig$OR[i], 3),
        "[", round(sig$OR_L95[i], 3), "-", round(sig$OR_U95[i], 3), "]",
        "P =", format.pval(sig$P[i], eps = 0.001), "\n")
  }
} else {
  cat("\nNo significant results at P < 0.05\n")
}

# P0层详细报告
cat("\n*** P0 Layer (Core Genes) ***\n")
p0_res <- summary_df[summary_df$Layer == "P0 (Core)", ]
for (i in 1:nrow(p0_res)) {
  cat(p0_res$Gene[i], ": OR =", round(p0_res$OR[i], 3),
      "[", round(p0_res$OR_L95[i], 3), "-", round(p0_res$OR_U95[i], 3), "]",
      "P =", format.pval(p0_res$P[i], eps = 0.001),
      p0_res$Expectation[i], "\n")
}

# =============================================================================
# 6. 标准MR森林图 (顶刊级)
# =============================================================================
cat("\nGenerating forest plot...\n")

pd <- summary_df[!is.na(summary_df$OR), ]
pd$Significance <- ifelse(pd$P < 0.05, "P < 0.05", ifelse(pd$P < 0.1, "P < 0.1", "NS"))
pd$Gene_Label <- paste0(pd$Gene, " (", pd$NSNP_Final, " SNPs)")

# 计算文本标签
pd$Label <- sprintf("%.2f [%.2f-%.2f] P=%.3f", pd$OR, pd$OR_L95, pd$OR_U95, pd$P)

p <- ggplot(pd, aes(y = reorder(Gene_Label, OR))) +
  geom_point(aes(x = OR, color = Significance), size = 5, shape = 15) +
  geom_errorbar(aes(xmin = OR_L95, xmax = OR_U95, color = Significance),
                linewidth = 1, orientation = "y") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_x_log10(breaks = c(0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.5),
                labels = c("0.70", "0.80", "0.90", "1.00", "1.10", "1.20", "1.30", "1.50")) +
  scale_color_manual(values = c("P < 0.05" = "#D55E00", "P < 0.1" = "#E69F00", "NS" = "#999999"),
                     name = "Significance") +
  labs(
    title = "Mendelian Randomization: 9 Hub Genes and Ischemic Stroke Risk",
    subtitle = "eQTLgen (clumped p<5e-8, r2<0.001) + FinnGen R12 I9_STR (n=377,277)",
    x = "Odds Ratio (95% CI, log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5)
  ) +
  geom_text(aes(x = OR_U95 + 0.03, label = Label), hjust = 0, size = 3.5, color = "black")

ggsave(file.path(output_dir, "MR_Forest_Plot.png"), p, width = 12, height = 7, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "MR_Forest_Plot.pdf"), p, width = 12, height = 7)

cat("Forest plot saved\n")

# =============================================================================
# 7. 完成
# =============================================================================
cat("\n======================================================================\n", sep = "")
cat("Analysis Complete!\n")
cat("======================================================================\n", sep = "")
cat("Output Directory:", output_dir, "\n")
cat("Files Generated:\n")
cat("  - MR_Summary.csv (汇总结果)\n")
cat("  - MR_*.csv (单基因结果)\n")
cat("  - MR_Forest_Plot.png (森林图)\n")
cat("  - MR_Forest_Plot.pdf (矢量图)\n")
cat("======================================================================\n", sep = "")
