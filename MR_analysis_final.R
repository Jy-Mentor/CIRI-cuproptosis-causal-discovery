#!/usr/bin/env Rscript
# MR分析 - 使用eQTLgen clumping数据 + FinnGen结局

setwd("D:/EQTL")
set.seed(42)

cat("=" ,rep("=", 70), "\n", sep="")
cat("MR分析: 9 Hub Genes vs 脑卒中 (eQTLgen + FinnGen)\n")
cat("=" ,rep("=", 70), "\n\n", sep="")

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(readxl)
  library(dplyr)
  library(ggplot2)
})

# 基因配置
p0_genes <- c("NFKB1", "FDX1", "STAT3")
p1_genes <- c("HIF1A", "HMOX1", "GPX4", "TNF", "IL6", "AGER")
all_genes <- c(p0_genes, p1_genes)

# 文件路径
eqtl_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-8_kb_10000_r2_0.001.xlsx"
outcome_file <- "D:/EQTL/finngen_R12_I9_STR"
output_dir <- paste0("D:/EQTL/MR_Final_Results_", format(Sys.time(), "%Y%m%d"))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("配置:\n")
cat("- 基因:", length(all_genes), "个 (P0:", length(p0_genes), ", P1:", length(p1_genes), ")\n")
cat("- eQTL:", basename(eqtl_file), "\n")
cat("- 结局: FinnGen R12 I9_STR\n\n")

# 读取eQTL
cat("[1/3] 读取eQTL数据...\n")
eqtl_all <- read_excel(eqtl_file)
cat("  eQTL SNP数:", nrow(eqtl_all), "\n")
cat("  基因数:", length(unique(eqtl_all$gene)), "\n\n")

# 检查基因
found <- intersect(all_genes, unique(eqtl_all$gene))
missing <- setdiff(all_genes, found)
cat("  找到基因:", paste(found, collapse = ", "), "\n")
if (length(missing) > 0) cat("  缺失:", paste(missing, collapse = ", "), "\n")
cat("\n")

all_genes <- found
p0_genes <- intersect(p0_genes, found)
p1_genes <- intersect(p1_genes, found)

# 读取FinnGen (使用fread更快)
cat("[2/3] 读取FinnGen结局数据...\n")
library(data.table)
outcome_data <- fread(outcome_file, select = c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt"))

# 重命名列为TwoSampleMR格式
setnames(outcome_data, c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt"),
         c("SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome", 
           "se.outcome", "pval.outcome", "eaf.outcome"))

# 添加必需列
outcome_data$phenotype <- "Ischemic_Stroke"
outcome_data$samplesize.outcome <- 377277
outcome_data$ncase.outcome <- 19862
outcome_data$ncontrol.outcome <- 357415

cat("  结局SNP数:", nrow(outcome_data), "\n\n")

# MR函数
run_mr <- function(gene) {
  result <- list(gene = gene, success = FALSE, error = NULL,
                 nsnp_exp = 0, nsnp_match = 0, nsnp_final = 0,
                 beta = NA, se = NA, or = NA, or_l = NA, or_u = NA, p = NA)
  
  cat("\n分析:", gene, "\n")
  cat(rep("-", 50), "\n", sep = "")
  
  tryCatch({
    # 1. 提取eQTL
    exposure <- eqtl_all[eqtl_all$gene == gene, ]
    result$nsnp_exp <- nrow(exposure)
    cat("  eQTL SNP:", nrow(exposure), "\n")
    
    if (nrow(exposure) == 0) {
      result$error <- "No eQTL SNPs"
      return(result)
    }
    
    # 2. 匹配结局SNP
    outcome_match <- outcome_data[outcome_data$SNP %in% exposure$SNP, ]
    result$nsnp_match <- nrow(outcome_match)
    cat("  匹配FinnGen:", nrow(outcome_match), "\n")
    
    if (nrow(outcome_match) == 0) {
      result$error <- "No matching SNPs in outcome"
      return(result)
    }
    
    # 3. 格式化暴露数据
    # 转换为data.frame (TwoSampleMR要求)
    exposure_df <- as.data.frame(exposure)
    outcome_df <- as.data.frame(outcome_match)
    
    exposure_fmt <- format_data(
      exposure_df,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "beta.exposure",
      se_col = "se.exposure",
      pval_col = "pval.exposure",
      eaf_col = "eaf.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure",
      phenotype_col = "exposure"
    )
    
    # 4. 格式化结局数据
    outcome_fmt <- format_data(
      outcome_df,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      pval_col = "pval.outcome",
      eaf_col = "eaf.outcome",
      effect_allele_col = "effect_allele.outcome",
      other_allele_col = "other_allele.outcome"
    )
    
    # 5. Harmonization
    dat <- harmonise_data(exposure_fmt, outcome_fmt)
    dat <- dat[dat$mr_keep == TRUE, ]
    result$nsnp_final <- nrow(dat)
    cat("  Harmonization后:", nrow(dat), "\n")
    
    if (nrow(dat) == 0) {
      result$error <- "No SNPs after harmonization"
      return(result)
    }
    
    # 6. MR分析
    if (nrow(dat) >= 3) {
      methods <- c("mr_ivw", "mr_egger_regression")
      method_name <- "IVW"
    } else {
      methods <- c("mr_wald_ratio")
      method_name <- "Wald"
    }
    
    res <- mr(dat, method_list = methods)
    
    if (method_name == "IVW") {
      main <- res[res$method == "Inverse variance weighted", ]
    } else {
      main <- res[res$method == "Wald ratio", ]
    }
    
    result$beta <- main$b
    result$se <- main$se
    result$or <- exp(main$b)
    result$or_l <- exp(main$b - 1.96 * main$se)
    result$or_u <- exp(main$b + 1.96 * main$se)
    result$p <- main$pval
    result$method <- method_name
    
    cat("  结果: OR =", round(result$or, 3),
        "[", round(result$or_l, 3), "-", round(result$or_u, 3), "]",
        "P =", format.pval(result$p, eps = 0.001), "\n")
    
    # 敏感性
    if (nrow(dat) >= 3) {
      het <- mr_heterogeneity(dat)
      cat("  异质性P:", format.pval(het$Q_pval[1], eps = 0.001), "\n")
    }
    
    if (nrow(dat) >= 4) {
      pleio <- mr_pleiotropy_test(dat)
      cat("  Egger截距P:", format.pval(pleio$pval, eps = 0.001), "\n")
    }
    
    # 保存
    write.csv(res, file.path(output_dir, paste0("MR_", gene, ".csv")), row.names = FALSE)
    result$success <- TRUE
    
  }, error = function(e) {
    cat("  错误:", conditionMessage(e), "\n")
    result$error <- conditionMessage(e)
  })
  
  return(result)
}

# 运行分析
cat("=" ,rep("=", 70), "\n", sep="")
cat("[3/3] 开始MR分析\n")
cat(rep("=", 70), "\n", sep="")

results <- list()
for (gene in all_genes) {
  results[[gene]] <- run_mr(gene)
}

# 汇总
cat("\n")
cat("=" ,rep("=", 70), "\n", sep="")
cat("汇总结果\n")
cat(rep("=", 70), "\n\n", sep="")

summary_df <- data.frame()
for (gene in names(results)) {
  r <- results[[gene]]
  
  exp_text <- ""
  if (r$success) {
    if (gene %in% c("NFKB1", "FDX1", "STAT3", "TNF", "IL6")) {
      exp_text <- ifelse(r$or > 1, "✓ 风险性", "✗ 相反(保护性)")
    } else if (gene %in% c("GPX4", "HMOX1")) {
      exp_text <- ifelse(r$or < 1, "✓ 保护性", "✗ 相反(风险性)")
    }
  } else {
    exp_text <- r$error
  }
  
  summary_df <- rbind(summary_df, data.frame(
    Gene = gene,
    Layer = ifelse(gene %in% p0_genes, "P0", "P1"),
    Beta = r$beta,
    SE = r$se,
    OR = r$or,
    OR_L95 = r$or_l,
    OR_U95 = r$or_u,
    P = r$p,
    NSNP_eQTL = r$nsnp_exp,
    NSNP_Match = r$nsnp_match,
    NSNP_Final = r$nsnp_final,
    Method = ifelse(r$success, r$method, "Failed"),
    Expectation = exp_text,
    stringsAsFactors = FALSE
  ))
}

print(summary_df[, c("Gene", "OR", "P", "NSNP_Final", "Expectation")])
write.csv(summary_df, file.path(output_dir, "MR_Summary.csv"), row.names = FALSE)

# 显著结果
sig <- summary_df[!is.na(summary_df$P) & summary_df$P < 0.05, ]
if (nrow(sig) > 0) {
  cat("\n显著结果 (P < 0.05):\n")
  cat(rep("-", 50), "\n", sep = "")
  for (i in 1:nrow(sig)) {
    direction <- ifelse(sig$OR[i] > 1, "风险性", "保护性")
    cat(sig$Gene[i], ": OR =", round(sig$OR[i], 3),
        "[", round(sig$OR_L95[i], 3), "-", round(sig$OR_U95[i], 3), "]",
        "P =", format.pval(sig$P[i], eps = 0.001),
        "(", direction, ")\n")
  }
} else {
  cat("\n无显著结果 (P < 0.05)\n")
}

# P0层详细结果
cat("\n")
cat("P0层核心基因:\n")
cat(rep("-", 50), "\n", sep = "")
p0_res <- summary_df[summary_df$Layer == "P0", ]
for (i in 1:nrow(p0_res)) {
  cat(p0_res$Gene[i], ":")
  if (!is.na(p0_res$OR[i])) {
    cat(" OR =", round(p0_res$OR[i], 3),
        "P =", format.pval(p0_res$P[i], eps = 0.001))
  }
  cat(" ", p0_res$Expectation[i], "\n")
}

# 森林图
if (nrow(summary_df[!is.na(summary_df$OR), ]) > 0) {
  cat("\n生成森林图...\n")
  
  pd <- summary_df[!is.na(summary_df$OR), ]
  pd$Signif <- ifelse(pd$P < 0.05, "Significant",
                      ifelse(pd$P < 0.1, "Trend", "NS"))
  
  p <- ggplot(pd, aes(y = reorder(Gene, OR))) +
    geom_point(aes(x = OR, color = Signif), size = 5, shape = 15) +
    geom_errorbarh(aes(xmin = OR_L95, xmax = OR_U95, color = Signif), 
                   height = 0.3, linewidth = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_x_log10(breaks = c(0.5, 0.8, 1, 1.25, 1.5, 2)) +
    scale_color_manual(values = c("Significant" = "#D55E00", 
                                  "Trend" = "#E69F00", 
                                  "NS" = "#999999")) +
    labs(
      title = "Mendelian Randomization: 9 Hub Genes and Ischemic Stroke",
      subtitle = "eQTLgen (clumped) + FinnGen R12 I9_STR",
      x = "Odds Ratio (95% CI)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  ggsave(file.path(output_dir, "Forest_Plot.png"), p, width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "Forest_Plot.pdf"), p, width = 10, height = 8)
  cat("  森林图已保存\n")
}

cat("\n")
cat("=" ,rep("=", 70), "\n", sep="")
cat("完成! 结果目录:", output_dir, "\n")
cat("=" ,rep("=", 70), "\n", sep="")
