#!/usr/bin/env Rscript
# MR分析 - 9 Hub Genes vs 脑卒中

setwd("D:/EQTL")
set.seed(42)

cat("=" ,rep("=", 70), "\n", sep="")
cat("MR分析: 9 Hub Genes vs Ischemic Stroke\n")
cat("=" ,rep("=", 70), "\n\n", sep="")

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

# 基因
p0_genes <- c("NFKB1", "FDX1", "STAT3")
p1_genes <- c("HIF1A", "HMOX1", "GPX4", "TNF", "IL6", "AGER")
all_genes <- c(p0_genes, p1_genes)

# 路径
eqtl_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-8_kb_10000_r2_0.001.xlsx"
outcome_file <- "D:/EQTL/finngen_R12_I9_STR"
output_dir <- paste0("D:/EQTL/MR_Results_", format(Sys.time(), "%Y%m%d"))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("输出目录:", output_dir, "\n\n")

# 读取数据
cat("[1/2] 读取eQTL...\n")
eqtl_all <- read_excel(eqtl_file)
cat("  eQTL SNPs:", nrow(eqtl_all), "\n")

cat("[2/2] 读取FinnGen...\n")
outcome_data <- fread(outcome_file, select = c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt"))
setnames(outcome_data, c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt"),
         c("SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", 
           "se.outcome", "pval.outcome", "eaf.outcome"))
outcome_data$phenotype <- "Stroke"
outcome_data$samplesize.outcome <- 377277
outcome_data$ncase.outcome <- 19862
cat("  Outcome SNPs:", nrow(outcome_data), "\n\n")

# MR函数
run_mr <- function(gene) {
  res <- list(gene = gene, success = FALSE, nsnp = 0, beta = NA, se = NA, or = NA, p = NA, error = NULL)
  
  cat("\n", gene, "\n", rep("-", 40), "\n", sep = "")
  
  tryCatch({
    # eQTL
    exposure <- eqtl_all[eqtl_all$gene == gene, ]
    res$nsnp <- nrow(exposure)
    cat("  eQTL SNPs:", nrow(exposure), "\n")
    
    if (nrow(exposure) == 0) {
      res$error <- "No eQTL"
      return(res)
    }
    
    # 匹配结局
    outcome <- outcome_data[outcome_data$SNP %in% exposure$SNP, ]
    cat("  Matched:", nrow(outcome), "\n")
    
    if (nrow(outcome) == 0) {
      res$error <- "No match"
      return(res)
    }
    
    # 转换为data.frame
    exp_df <- as.data.frame(exposure)
    out_df <- as.data.frame(outcome)
    
    # 格式化
    exp_fmt <- format_data(exp_df, type = "exposure",
                           snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure",
                           pval_col = "pval.exposure", eaf_col = "eaf.exposure",
                           effect_allele_col = "effect_allele.exposure",
                           other_allele_col = "other_allele.exposure")
    
    out_fmt <- format_data(out_df, type = "outcome",
                           snp_col = "SNP", beta_col = "beta.outcome", se_col = "se.outcome",
                           pval_col = "pval.outcome", eaf_col = "eaf.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome")
    
    # Harmonization
    dat <- harmonise_data(exp_fmt, out_fmt)
    dat <- dat[dat$mr_keep == TRUE, ]
    cat("  After harmonization:", nrow(dat), "\n")
    
    if (nrow(dat) == 0) {
      res$error <- "No valid SNPs"
      return(res)
    }
    
    # MR分析 - 根据SNP数量选择方法
    if (nrow(dat) == 1) {
      res_mr <- mr(dat, method_list = c("mr_wald_ratio"))
      main <- res_mr[res_mr$method == "Wald ratio", ]
      method_name <- "Wald"
    } else {
      res_mr <- mr(dat, method_list = c("mr_ivw"))
      main <- res_mr[res_mr$method == "Inverse variance weighted", ]
      method_name <- "IVW"
    }
    
    # 检查结果是否有效
    if (is.na(main$b) || is.na(main$pval)) {
      res$error <- "MR returned NA"
      return(res)
    }
    
    res$beta <- main$b
    res$se <- main$se
    res$or <- exp(main$b)
    res$or_l <- exp(main$b - 1.96 * main$se)
    res$or_u <- exp(main$b + 1.96 * main$se)
    res$p <- main$pval
    res$method <- method_name
    res$success <- TRUE
    
    cat("  ", method_name, "OR =", round(res$or, 3),
        "[", round(res$or_l, 3), "-", round(res$or_u, 3), "]",
        "P =", format.pval(res$p, eps = 0.001), "\n")
    
    write.csv(res_mr, file.path(output_dir, paste0(gene, ".csv")), row.names = FALSE)
    
  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
    res$error <- conditionMessage(e)
  })
  
  return(res)
}

# 运行分析
cat("=" ,rep("=", 70), "\n", sep="")
cat("Running MR analysis...\n")
cat(rep("=", 70), "\n", sep = "")

results <- list()
for (gene in all_genes) {
  results[[gene]] <- run_mr(gene)
}

# 汇总
cat("\n")
cat("=" ,rep("=", 70), "\n", sep="")
cat("Summary Results\n")
cat(rep("=", 70), "\n\n", sep = "")

summary_list <- list()
for (gene in names(results)) {
  r <- results[[gene]]
  
  exp_text <- ""
  if (r$success && !is.na(r$or)) {
    if (gene %in% c("NFKB1", "FDX1", "STAT3", "TNF", "IL6")) {
      exp_text <- ifelse(r$or > 1, "✓ Risk", "✗ Protective")
    } else if (gene %in% c("GPX4", "HMOX1")) {
      exp_text <- ifelse(r$or < 1, "✓ Protective", "✗ Risk")
    } else {
      exp_text <- ""
    }
  } else {
    exp_text <- ifelse(is.null(r$error), "Failed", r$error)
  }
  
  summary_list[[gene]] <- data.frame(
    Gene = gene,
    Layer = ifelse(gene %in% p0_genes, "P0", "P1"),
    Beta = ifelse(r$success, r$beta, NA),
    SE = ifelse(r$success, r$se, NA),
    OR = ifelse(r$success, r$or, NA),
    OR_L95 = ifelse(r$success, r$or_l, NA),
    OR_U95 = ifelse(r$success, r$or_u, NA),
    P = ifelse(r$success, r$p, NA),
    NSNP = r$nsnp,
    Method = ifelse(r$success, ifelse(is.null(r$method), "Unknown", r$method), "Failed"),
    Expectation = exp_text,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

print(summary_df[, c("Gene", "OR", "P", "NSNP", "Expectation")])
write.csv(summary_df, file.path(output_dir, "Summary.csv"), row.names = FALSE)

# 显著结果
sig <- summary_df[!is.na(summary_df$P) & summary_df$P < 0.05, ]
if (nrow(sig) > 0) {
  cat("\nSignificant results (P < 0.05):\n")
  for (i in 1:nrow(sig)) {
    cat(sig$Gene[i], ": OR =", round(sig$OR[i], 3),
        "[", round(sig$OR_L95[i], 3), "-", round(sig$OR_U95[i], 3), "]", 
        "P =", format.pval(sig$P[i]), "\n")
  }
}

# P0层
cat("\nP0 Layer:\n")
p0 <- summary_df[summary_df$Layer == "P0", ]
for (i in 1:nrow(p0)) {
  cat(p0$Gene[i], ": OR =", round(p0$OR[i], 3), "P =", format.pval(p0$P[i]), p0$Expectation[i], "\n")
}

# 森林图
if (sum(!is.na(summary_df$OR)) > 0) {
  pd <- summary_df[!is.na(summary_df$OR), ]
  pd$Sig <- ifelse(pd$P < 0.05, "Sig", ifelse(pd$P < 0.1, "Trend", "NS"))
  
  p <- ggplot(pd, aes(y = reorder(Gene, OR))) +
    geom_point(aes(x = OR, color = Sig), size = 4) +
    geom_errorbarh(aes(xmin = OR_L95, xmax = OR_U95, color = Sig), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_x_log10() +
    scale_color_manual(values = c("Sig" = "red", "Trend" = "orange", "NS" = "gray")) +
    labs(title = "MR: 9 Genes vs Stroke", x = "OR (95% CI)", y = NULL) +
    theme_minimal()
  
  ggsave(file.path(output_dir, "Forest.png"), p, width = 10, height = 6, dpi = 300)
  cat("\nForest plot saved\n")
}

cat("\nDone! Results in:", output_dir, "\n")
