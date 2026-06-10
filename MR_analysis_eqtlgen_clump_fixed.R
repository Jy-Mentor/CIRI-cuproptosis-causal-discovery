#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# MR分析 - 使用eQTLgen clumping数据
# 9个Hub基因 vs 脑卒中风险

setwd("D:/EQTL")
set.seed(42)

cat(rep("=", 70), "\n", sep="")
cat("MR分析 - 9个Hub基因 vs 脑卒中 (eQTLgen)\n")
cat(rep("=", 70), "\n\n", sep="")

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(readxl)
  library(dplyr)
  library(ggplot2)
})

# 配置
p0_genes <- c("NFKB1", "FDX1", "STAT3")
p1_genes <- c("HIF1A", "HMOX1", "GPX4", "TNF", "IL6", "AGER")
all_genes <- c(p0_genes, p1_genes)

eqtl_clump_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-8_kb_10000_r2_0.001.xlsx"
outcome_file <- "D:/EQTL/finngen_R12_I9_STR"
output_dir <- paste0("D:/EQTL/MR_9genes_Results_", format(Sys.time(), "%Y%m%d"))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("基因:", length(all_genes), "个\n")
cat("eQTL:", eqtl_clump_file, "\n\n")

# 读取eQTL
cat("读取eQTL数据...\n")
eqtl_all <- read_excel(eqtl_clump_file)
cat("总SNP:", nrow(eqtl_all), "，基因:", length(unique(eqtl_all$gene)), "\n\n")

# 检查基因
found_genes <- intersect(all_genes, unique(eqtl_all$gene))
missing <- setdiff(all_genes, found_genes)
cat("找到:", paste(found_genes, collapse = ", "), "\n")
if (length(missing) > 0) cat("缺失:", paste(missing, collapse = ", "), "\n")
cat("\n")

all_genes <- found_genes
p0_genes <- intersect(p0_genes, found_genes)
p1_genes <- intersect(p1_genes, found_genes)

# 读取FinnGen
cat("读取FinnGen...\n")
outcome_data <- read.table(outcome_file, header = TRUE, sep = "\t",
                           colClasses = c("rsids" = "character"))
outcome_data$SNP <- outcome_data$rsids
outcome_data$beta.outcome <- outcome_data$beta
outcome_data$se.outcome <- outcome_data$sebeta
outcome_data$pval.outcome <- outcome_data$pval
outcome_data$effect_allele <- outcome_data$ref
outcome_data$other_allele <- outcome_data$alt
outcome_data$eaf.outcome <- outcome_data$af_alt
outcome_data$phenotype <- "Stroke"
outcome_data$samplesize.outcome <- 377277
cat("结局SNP:", nrow(outcome_data), "\n\n")

# MR函数
run_mr <- function(gene) {
  res_list <- list(gene = gene, success = FALSE, nsnp = 0, error = NULL)
  
  cat("\n分析:", gene, "\n")
  
  tryCatch({
    # 提取eQTL
    exposure <- eqtl_all[eqtl_all$gene == gene, ]
    res_list$nsnp <- nrow(exposure)
    cat("  eQTL SNP:", nrow(exposure), "\n")
    
    if (nrow(exposure) == 0) {
      res_list$error <- "No SNPs"
      return(res_list)
    }
    
    # 提取结局
    outcome <- outcome_data[outcome_data$SNP %in% exposure$SNP, ]
    cat("  匹配结局:", nrow(outcome), "\n")
    
    if (nrow(outcome) == 0) {
      res_list$error <- "No overlap"
      return(res_list)
    }
    
    # Harmonization
    dat <- harmonise_data(exposure, outcome)
    dat <- dat[dat$mr_keep == TRUE, ]
    cat("  Harmonization后:", nrow(dat), "\n")
    
    if (nrow(dat) == 0) {
      res_list$error <- "No valid SNPs"
      return(res_list)
    }
    
    # MR
    if (nrow(dat) >= 3) {
      methods <- c("mr_ivw", "mr_egger_regression")
    } else {
      methods <- c("mr_wald_ratio")
    }
    
    res <- mr(dat, method_list = methods)
    
    if ("Inverse variance weighted" %in% res$method) {
      main <- res[res$method == "Inverse variance weighted", ]
      method_used <- "IVW"
    } else {
      main <- res[res$method == "Wald ratio", ]
      method_used <- "Wald"
    }
    
    or <- exp(main$b)
    or_l <- exp(main$b - 1.96 * main$se)
    or_u <- exp(main$b + 1.96 * main$se)
    
    cat("  OR:", round(or, 3), "[", round(or_l, 3), "-", round(or_u, 3), "]",
        "P:", format.pval(main$pval, eps = 0.001), "\n")
    
    # 敏感性
    if (nrow(dat) >= 3) {
      het <- mr_heterogeneity(dat)
      cat("  异质性Q_pval:", format.pval(het$Q_pval[het$method == "Inverse variance weighted"], eps = 0.001), "\n")
    }
    
    if (nrow(dat) >= 4) {
      pleio <- mr_pleiotropy_test(dat)
      cat("  Egger截距P:", format.pval(pleio$pval, eps = 0.001), "\n")
    }
    
    write.csv(res, file.path(output_dir, paste0("MR_", gene, ".csv")), row.names = FALSE)
    
    res_list$success <- TRUE
    res_list$beta <- main$b
    res_list$se <- main$se
    res_list$or <- or
    res_list$or_l <- or_l
    res_list$or_u <- or_u
    res_list$p <- main$pval
    res_list$method <- method_used
    
  }, error = function(e) {
    cat("  错误:", conditionMessage(e), "\n")
    res_list$error <- conditionMessage(e)
  })
  
  return(res_list)
}

# 运行
cat(rep("=", 70), "\n", sep="")
cat("开始分析\n")
cat(rep("=", 70), "\n\n", sep="")

results <- list()
for (gene in all_genes) {
  results[[gene]] <- run_mr(gene)
}

# 汇总
cat("\n")
cat(rep("=", 70), "\n", sep="")
cat("汇总结果\n")
cat(rep("=", 70), "\n\n", sep="")

summary_df <- data.frame()
for (gene in names(results)) {
  r <- results[[gene]]
  if (r$success) {
    exp_text <- ""
    if (gene %in% c("NFKB1", "FDX1", "STAT3", "TNF", "IL6")) {
      exp_text <- ifelse(r$or > 1, "✓ 风险性", "✗ 相反")
    } else if (gene %in% c("GPX4", "HMOX1")) {
      exp_text <- ifelse(r$or < 1, "✓ 保护性", "✗ 相反")
    }
    
    summary_df <- rbind(summary_df, data.frame(
      Gene = gene,
      Layer = ifelse(gene %in% p0_genes, "P0", "P1"),
      OR = r$or, OR_L = r$or_l, OR_U = r$or_u,
      P = r$p, NSNP = r$nsnp, Method = r$method,
      Expectation = exp_text
    ))
  } else {
    summary_df <- rbind(summary_df, data.frame(
      Gene = gene, Layer = ifelse(gene %in% p0_genes, "P0", "P1"),
      OR = NA, OR_L = NA, OR_U = NA, P = NA,
      NSNP = r$nsnp, Method = "Failed", Expectation = r$error
    ))
  }
}

print(summary_df[, c("Gene", "OR", "P", "NSNP", "Expectation")])
write.csv(summary_df, file.path(output_dir, "MR_summary.csv"), row.names = FALSE)

# 显著结果
sig <- summary_df[!is.na(summary_df$P) & summary_df$P < 0.05, ]
if (nrow(sig) > 0) {
  cat("\n显著 (P<0.05):\n")
  for (i in 1:nrow(sig)) {
    cat(sig$Gene[i], ": OR=", round(sig$OR[i], 3),
        "[", round(sig$OR_L[i], 3), "-", round(sig$OR_U[i], 3), "], P=", format.pval(sig$P[i]), "\n")
  }
}

# 森林图
if (nrow(summary_df[!is.na(summary_df$OR), ]) > 0) {
  pd <- summary_df[!is.na(summary_df$OR), ]
  pd$Sig <- ifelse(pd$P < 0.05, "Sig", ifelse(pd$P < 0.1, "Trend", "NS"))
  
  p <- ggplot(pd, aes(y = reorder(Gene, OR))) +
    geom_point(aes(x = OR, color = Sig), size = 4) +
    geom_errorbarh(aes(xmin = OR_L, xmax = OR_U, color = Sig), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_x_log10() +
    scale_color_manual(values = c("Sig" = "red", "Trend" = "orange", "NS" = "gray")) +
    labs(title = "MR: Hub Genes vs Stroke", x = "OR (95% CI)", y = "Gene") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "forest.png"), p, width = 10, height = 6, dpi = 300)
  cat("\n森林图已保存\n")
}

cat("\n完成! 结果在:", output_dir, "\n")
