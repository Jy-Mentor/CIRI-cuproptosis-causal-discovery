#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# MR分析 - 使用本地eQTL和GWAS数据
# 8个Hub基因: NFKB1, FDX1, STAT3, HIF1A, HMOX1, GPX4, TNF, IL6, AGER
# 结局: 脑卒中 (使用FinnGen R12 I9_STR数据)

setwd("D:/EQTL")

# 加载包
suppressPackageStartupMessages({
  if (!require("TwoSampleMR", quietly = TRUE)) install.packages("TwoSampleMR")
  if (!require("ieugwasr", quietly = TRUE)) install.packages("ieugwasr")
  if (!require("data.table", quietly = TRUE)) install.packages("data.table")
  if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!require("readr", quietly = TRUE)) install.packages("readr")
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

cat(rep("=", 70), "\n", sep = "")
cat("本地数据MR分析: 8个Hub基因 vs 脑卒中\n")
cat(rep("=", 70), "\n\n", sep = "")

# ============================================
# 配置
# ============================================
# 基因列表
p0_genes <- c("NFKB1", "FDX1", "STAT3")
p1_genes <- c("HIF1A", "HMOX1", "GPX4", "TNF", "IL6", "AGER")
all_genes <- c(p0_genes, p1_genes)

# 数据路径
eqtl_dir <- "D:/EQTL/rawdata"
outcome_file <- "D:/EQTL/finngen_R12_I9_STR/finngen_R12_I9_STR.tsv"
clump_dir <- "D:/EQTL/clump"
output_dir <- "D:/EQTL/MR_8genes_Results"

# 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("分析配置:\n")
cat("- P0层基因:", paste(p0_genes, collapse = ", "), "\n")
cat("- P1层基因:", paste(p1_genes, collapse = ", "), "\n")
cat("- eQTL目录:", eqtl_dir, "\n")
cat("- 结局数据:", outcome_file, "\n")
cat("- 输出目录:", output_dir, "\n\n")

# ============================================
# 读取结局数据 (脑卒中 GWAS)
# ============================================
cat(rep("-", 70), "\n", sep = "")
cat("步骤1: 读取结局数据 (FinnGen R12 I9_STR)\n")
cat(rep("-", 70), "\n", sep = "")

if (!file.exists(outcome_file)) {
  cat("错误: 结局文件不存在:", outcome_file, "\n")
  cat("请确保FinnGen脑卒中GWAS数据已下载\n")
  quit(status = 1)
}

# 读取结局数据
outcome_data <- fread(outcome_file, select = c("rsid", "chrom", "pos", "ref", "alt", 
                                                 "beta", "se", "pval", "n_cases", "n_controls"))
cat("结局数据行数:", nrow(outcome_data), "\n")
cat("列名:", paste(names(outcome_data), collapse = ", "), "\n\n")

# 标准化列名以匹配TwoSampleMR格式
setnames(outcome_data, old = c("rsid", "ref", "alt", "beta", "se", "pval"),
         new = c("SNP", "effect_allele", "other_allele", "beta.outcome", "se.outcome", "pval.outcome"))

# ============================================
# MR分析主函数
# ============================================
run_mr_local <- function(gene_symbol) {
  cat(rep("=", 60), "\n", sep = "")
  cat("分析基因:", gene_symbol, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  result <- list(
    gene = gene_symbol,
    success = FALSE,
    exposure_nsnps = 0,
    harmonized_nsnps = 0,
    mr_result = NULL,
    heterogeneity = NULL,
    pleiotropy = NULL,
    error = NULL
  )
  
  # 1. 读取eQTL暴露数据
  cat("\n1. 读取eQTL暴露数据...\n")
  
  # 查找eQTL文件
  eqtl_pattern <- paste0("eqtl-a-.*", gene_symbol, "_full_outcome.rds")
  eqtl_files <- list.files(eqtl_dir, pattern = eqtl_pattern, full.names = TRUE)
  
  # 如果找不到精确匹配，尝试模糊匹配
  if (length(eqtl_files) == 0) {
    eqtl_files <- list.files(eqtl_dir, pattern = paste0(".*", gene_symbol, ".*"), full.names = TRUE)
  }
  
  if (length(eqtl_files) == 0) {
    cat("   未找到eQTL文件，尝试通过Ensembl ID查找...\n")
    result$error <- "No eQTL file found"
    return(result)
  }
  
  cat("   找到eQTL文件:", basename(eqtl_files[1]), "\n")
  
  # 读取eQTL数据
  tryCatch({
    exposure_dat <- readRDS(eqtl_files[1])
    cat("   eQTL数据行数:", nrow(exposure_dat), "\n")
    
    # 筛选显著SNP (p < 5e-8)
    exposure_sig <- exposure_dat[exposure_dat$pval < 5e-8, ]
    cat("   显著SNP数 (p<5e-8):", nrow(exposure_sig), "\n")
    
    if (nrow(exposure_sig) == 0) {
      cat("   警告: 无全基因组显著SNP，尝试p<1e-5...\n")
      exposure_sig <- exposure_dat[exposure_dat$pval < 1e-5, ]
      cat("   p<1e-5 SNP数:", nrow(exposure_sig), "\n")
    }
    
    result$exposure_nsnps <- nrow(exposure_sig)
    
    # 2. LD Clumping (去除连锁不平衡)
    cat("\n2. LD Clumping...\n")
    
    if (nrow(exposure_sig) > 0) {
      # 读取预计算的clumping结果或使用ieugwasr
      tryCatch({
        exposure_clumped <- clump_data(exposure_sig, 
                                        clump_kb = 10000, 
                                        clump_r2 = 0.001,
                                        clump_p1 = 1,
                                        clump_p2 = 1)
        cat("   Clumping后SNP数:", nrow(exposure_clumped), "\n")
      }, error = function(e) {
        cat("   Clumping失败，使用原始数据:", conditionMessage(e), "\n")
        exposure_clumped <- exposure_sig
      })
    } else {
      result$error <- "No significant SNPs"
      return(result)
    }
    
    # 3. 提取结局数据中的对应SNP
    cat("\n3. 提取结局数据...\n")
    
    outcome_dat <- outcome_data[outcome_data$SNP %in% exposure_clumped$SNP, ]
    cat("   匹配的结局SNP数:", nrow(outcome_dat), "\n")
    
    if (nrow(outcome_dat) == 0) {
      result$error <- "No overlapping SNPs"
      return(result)
    }
    
    # 准备暴露数据格式
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
      phenotype_col = "Phenotype"
    )
    
    # 准备结局数据格式
    outcome_formatted <- format_data(
      outcome_dat,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      pval_col = "pval.outcome",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele"
    )
    
    # 4. Harmonization
    cat("\n4. 数据Harmonization...\n")
    dat <- harmonise_data(
      exposure_dat = exposure_formatted,
      outcome_dat = outcome_formatted
    )
    
    # 过滤有效SNP
    dat <- dat[dat$mr_keep == TRUE, ]
    result$harmonized_nsnps <- nrow(dat)
    cat("   Harmonization后有效SNP:", nrow(dat), "\n")
    
    if (nrow(dat) == 0) {
      result$error <- "No SNPs after harmonization"
      return(result)
    }
    
    # 5. 执行MR分析
    cat("\n5. 执行MR分析...\n")
    
    methods_list <- c("mr_wald_ratio")
    if (nrow(dat) >= 3) {
      methods_list <- c(methods_list, "mr_ivw")
    }
    if (nrow(dat) >= 4) {
      methods_list <- c(methods_list, "mr_egger_regression")
    }
    
    res <- mr(dat, method_list = methods_list)
    result$mr_result <- res
    
    cat("\n   MR结果:\n")
    print(res[, c("exposure", "method", "nsnp", "b", "se", "pval")])
    
    # 6. 敏感性分析
    if (nrow(dat) >= 3) {
      cat("\n6. 异质性检验...\n")
      het <- mr_heterogeneity(dat)
      result$heterogeneity <- het
      print(het[, c("method", "Q", "Q_df", "Q_pval")])
    }
    
    if ("mr_egger_regression" %in% res$method) {
      cat("\n7. MR-Egger多效性检验...\n")
      pleio <- mr_pleiotropy_test(dat)
      result$pleiotropy <- pleio
      cat("   Intercept:", round(pleio$egger_intercept, 4), 
          "(p =", format.pval(pleio$pval, eps = 0.001), ")\n")
    }
    
    # 8. 保存结果
    write.csv(res, file.path(output_dir, paste0("MR_result_", gene_symbol, ".csv")), row.names = FALSE)
    
    result$success <- TRUE
    
  }, error = function(e) {
    cat("   错误:", conditionMessage(e), "\n")
    result$error <- conditionMessage(e)
  })
  
  return(result)
}

# ============================================
# 运行所有基因分析
# ============================================
cat(rep("=", 70), "\n", sep = "")
cat("开始MR分析\n")
cat(rep("=", 70), "\n\n", sep = "")

all_results <- list()

for (gene in all_genes) {
  result <- run_mr_local(gene)
  all_results[[gene]] <- result
  cat("\n")
}

# ============================================
# 汇总结果
# ============================================
cat(rep("=", 70), "\n", sep = "")
cat("MR分析汇总\n")
cat(rep("=", 70), "\n\n", sep = "")

summary_list <- list()

for (gene in names(all_results)) {
  res <- all_results[[gene]]
  
  if (res$success && !is.null(res$mr_result)) {
    main_res <- res$mr_result[res$mr_result$method %in% c("Inverse variance weighted", "Wald ratio"), ]
    if (nrow(main_res) > 0) {
      main_res <- main_res[1, ]
      
      or <- exp(main_res$b)
      or_lower <- exp(main_res$b - 1.96 * main_res$se)
      or_upper <- exp(main_res$b + 1.96 * main_res$se)
      
      # 判断预期
      expectation <- ""
      if (gene == "NFKB1") {
        expectation <- ifelse(or > 1, "✓ 符合预期 (OR>1)", "✗ 与预期相反")
      } else if (gene == "FDX1") {
        expectation <- ifelse(or < 1, "✓ 符合预期 (OR<1)", "✗ 与预期相反")
      }
      
      summary_list[[gene]] <- data.frame(
        Gene = gene,
        Layer = ifelse(gene %in% p0_genes, "P0 (Core)", "P1 (Supplementary)"),
        Beta = main_res$b,
        SE = main_res$se,
        OR = or,
        OR_lower = or_lower,
        OR_upper = or_upper,
        P_value = main_res$pval,
        NSNP = main_res$nsnp,
        Method = main_res$method,
        Expectation = expectation,
        stringsAsFactors = FALSE
      )
    }
  } else {
    summary_list[[gene]] <- data.frame(
      Gene = gene,
      Layer = ifelse(gene %in% p0_genes, "P0 (Core)", "P1 (Supplementary)"),
      Beta = NA,
      SE = NA,
      OR = NA,
      OR_lower = NA,
      OR_upper = NA,
      P_value = NA,
      NSNP = 0,
      Method = "Failed",
      Expectation = res$error,
      stringsAsFactors = FALSE
    )
  }
}

# 合并结果
summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

cat("汇总结果:\n")
print(summary_df[, c("Gene", "Layer", "OR", "P_value", "NSNP", "Expectation")])

# 保存汇总
write.csv(summary_df, file.path(output_dir, "MR_summary_8genes.csv"), row.names = FALSE)
cat("\n汇总结果已保存至:", file.path(output_dir, "MR_summary_8genes.csv"), "\n")

# 显著结果
cat("\n")
cat(rep("-", 70), "\n", sep = "")
cat("显著结果 (P < 0.05):\n")
cat(rep("-", 70), "\n", sep = "")

sig_results <- summary_df[!is.na(summary_df$P_value) & summary_df$P_value < 0.05, ]
if (nrow(sig_results) > 0) {
  for (i in 1:nrow(sig_results)) {
    direction <- ifelse(sig_results$OR[i] > 1, "风险性", "保护性")
    cat(sig_results$Gene[i], ": OR =", round(sig_results$OR[i], 3), 
        "(", direction, "), P =", format.pval(sig_results$P_value[i], eps = 0.001), 
        sig_results$Expectation[i], "\n")
  }
} else {
  cat("未发现显著关联 (P < 0.05)\n")
}

# P0层验证
cat("\n")
cat(rep("-", 70), "\n", sep = "")
cat("P0层核心基因验证:\n")
cat(rep("-", 70), "\n", sep = "")

p0_results <- summary_df[summary_df$Layer == "P0 (Core)", ]
for (i in 1:nrow(p0_results)) {
  gene <- p0_results$Gene[i]
  or <- p0_results$OR[i]
  pval <- p0_results$P_value[i]
  
  if (!is.na(or)) {
    cat(gene, ": OR =", round(or, 3), ", P =", format.pval(pval, eps = 0.001), 
        p0_results$Expectation[i], "\n")
  } else {
    cat(gene, ": 分析失败 (", p0_results$Expectation[i], ")\n")
  }
}

# ============================================
# 可视化
# ============================================
cat("\n")
cat(rep("-", 70), "\n", sep = "")
cat("生成可视化\n")
cat(rep("-", 70), "\n", sep = "")

if (nrow(summary_df[!is.na(summary_df$OR), ]) > 0) {
  
  plot_data <- summary_df[!is.na(summary_df$OR), ]
  plot_data$Significance <- ifelse(plot_data$P_value < 0.05, "Significant", 
                                   ifelse(plot_data$P_value < 0.1, "Trend", "NS"))
  
  p <- ggplot(plot_data, aes(x = reorder(Gene, OR), y = OR, color = Significance)) +
    geom_point(aes(size = -log10(P_value)), position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    coord_flip() +
    scale_color_manual(values = c("Significant" = "red", "Trend" = "orange", "NS" = "gray")) +
    labs(
      title = "MR分析: 8个Hub基因与脑卒中风险",
      subtitle = "使用FinnGen R12 I9_STR GWAS数据",
      x = "基因",
      y = "OR (95% CI)",
      color = "显著性",
      size = "-log10(P)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_dir, "MR_forest_plot.png"), p, width = 10, height = 8, dpi = 300)
  cat("森林图已保存至:", file.path(output_dir, "MR_forest_plot.png"), "\n")
}

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("分析完成!\n")
cat(rep("=", 70), "\n", sep = "")
cat("所有结果保存在:", output_dir, "\n")
cat(rep("=", 70), "\n", sep = "")
