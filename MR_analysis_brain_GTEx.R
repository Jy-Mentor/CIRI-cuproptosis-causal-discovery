#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# Mendelian Randomization: 9 Hub Genes vs Ischemic Stroke
# 版本: 4.0 (Brain GTEx Edition)
# 日期: 2026-04-24
# 
# 暴露: GTEx Brain eQTL (在线获取)
# 结局: FinnGen R12 I9_STR
# =============================================================================

setwd("D:/EQTL")
set.seed(42)

# =============================================================================
# 0. 包加载与版本检查
# =============================================================================
cat("======================================================================\n", sep = "")
cat("MR Analysis: 9 Hub Genes vs Ischemic Stroke (Brain GTEx v4.0)\n")
cat("======================================================================\n\n", sep = "")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Random Seed: 42\n\n")

suppressPackageStartupMessages({
  if (!require("TwoSampleMR", quietly = TRUE)) stop("TwoSampleMR not installed")
  if (!require("ieugwasr", quietly = TRUE)) stop("ieugwasr not installed")
  if (!require("data.table", quietly = TRUE)) stop("data.table not installed")
  if (!require("dplyr", quietly = TRUE)) stop("dplyr not installed")
  if (!require("ggplot2", quietly = TRUE)) stop("ggplot2 not installed")
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# 打印包版本
cat("Package Versions:\n")
cat("  R:", R.version.string, "\n")
cat("  TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
cat("  ieugwasr:", as.character(packageVersion("ieugwasr")), "\n\n")

# =============================================================================
# 0.5 设置OpenGWAS API Token
# =============================================================================
cat("Setting up OpenGWAS API authentication...\n")
# 设置JWT token (从 https://api.opengwas.io/ 获取)
# Token有效期至: 2026-05-08 06:34:56 UTC
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3NzAxMjQ5NiwiZXhwIjoxNzc4MjIyMDk2fQ.sz-zieAniOZLfPGDiQRo3W6z5PQ6HfEvsfSjvKCVLOlUJumnpM-jE9Se6QAqXezfH-ffEEbM0heC3pQu6MjBF1ttqaRpIdOv8S7_pq-s-HLVQUJZ3Ge5Ort66N89riaj9nRS-FdGidJAft58CEmjk-Q1lpJpJbOSGJiKCRQ_Man5vbCBiTH3p_e490wbapqqpxaXWs4ki-xYUDlQXuc6u2lZIh7_6TSei6Chmew0EWy8wcVgSTJGAmAiT8jrFTq5ydxsPvzZMnMZgEezAfrHIi9i3zJqbTzggs8Okqy2owfzLzcE7MKt_EcUSHYvMzRsL3YcbILipkodoy-sbzNmpw")

# 验证token
tryCatch({
  user_info <- ieugwasr::user()
  cat("  Token verified successfully!\n")
  cat("  User:", user_info$user, "\n")
  cat("  Tier:", user_info$tier, "\n")
  cat("  Remaining queries:", user_info$allowance - user_info$used, "\n\n")
}, error = function(e) {
  cat("  WARNING: Token verification failed:", conditionMessage(e), "\n")
  cat("  Continuing with local data fallback...\n\n")
})

# =============================================================================
# 1. 配置与参数
# =============================================================================
# 基因列表
p0_genes <- c("NFKB1", "FDX1", "STAT3")
p1_genes <- c("HIF1A", "HMOX1", "GPX4", "TNF", "IL6", "AGER")
all_genes <- c(p0_genes, p1_genes)

# GTEx脑组织选择
# 可选: "Brain_Frontal_Cortex", "Brain_Putamen", "Brain_Hippocampus", 
#       "Brain_Cerebellum", "Brain_Cortex", "Brain_Amygdala"
brain_tissue <- "Brain_Frontal_Cortex"  # 默认使用额叶皮层
alt_tissues <- c("Brain_Putamen", "Brain_Cortex", "Brain_Hippocampus")

outcome_file <- "D:/EQTL/finngen_R12_I9_STR"

# 输出目录
output_dir <- paste0("D:/EQTL/MR_BrainGTEx_", brain_tissue, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Configuration:\n")
cat("  Total Genes:", length(all_genes), "\n")
cat("  P0 Layer:", paste(p0_genes, collapse = ", "), "\n")
cat("  P1 Layer:", paste(p1_genes, collapse = ", "), "\n")
cat("  Brain Tissue:", brain_tissue, "\n")
cat("  Alternative tissues:", paste(alt_tissues, collapse = ", "), "\n")
cat("  Outcome:", basename(outcome_file), "\n")
cat("  Output:", output_dir, "\n\n")

cat("Scientific Hypothesis:\n")
cat("  NFKB1, FDX1, STAT3, TNF, IL6, AGER: High expr -> Risk (OR > 1)\n")
cat("  GPX4, HMOX1: High expr -> Protective (OR < 1)\n\n")

# =============================================================================
# 2. 获取GTEx脑组织eQTL数据
# =============================================================================
cat("----------------------------------------------------------------------\n", sep = "")
cat("Step 1: Fetching GTEx Brain eQTL Data from IEU GWAS Database\n")
cat("----------------------------------------------------------------------\n", sep = "")

# 将基因名转换为GTEx ID格式 (ENSG)
# 使用gwasglue或手动映射
cat("  Querying gene information...\n")

# 从IEU获取GTEx脑组织基因eQTL数据
get_gene_eqtl <- function(gene_symbol, tissue = brain_tissue) {
  cat("  Fetching", gene_symbol, "from", tissue, "...\n")
  
  tryCatch({
    # 方法1: 使用gwasinfo搜索GTEx脑组织数据
    cat("    Querying IEU GWAS database...\n")
    
    # 获取所有可用数据集信息
    gwas_info <- ieugwasr::gwasinfo()
    
    # 搜索GTEx脑组织中该基因的数据
    # GTEx格式示例: "eqtl-a-ENSG00000109332" (Gene: NFKB1)
    gtex_pattern <- paste0("gtex.*", tissue, "|GTEx.*", tissue, "|", tissue)
    gene_pattern <- gene_symbol
    
    gtex_matches <- gwas_info[grepl(tissue, gwas_info$trait, ignore.case = TRUE), ]
    
    if (nrow(gtex_matches) > 0) {
      cat("    Found", nrow(gtex_matches), "GTEx datasets for", tissue, "\n")
      
      # 搜索特定基因
      gene_matches <- gtex_matches[grepl(gene_symbol, gtex_matches$trait, ignore.case = TRUE) |
                                     grepl(gene_symbol, gtex_matches$id, ignore.case = TRUE), ]
      
      if (nrow(gene_matches) > 0) {
        cat("    Found gene", gene_symbol, "in", nrow(gene_matches), "datasets\n")
        cat("    Using dataset:", gene_matches$id[1], "\n")
        
        # 提取instruments (token已通过环境变量OPENGWAS_JWT自动使用)
         exp_dat <- extract_instruments(
           outcomes = gene_matches$id[1],
           p1 = 5e-08,
           clump = TRUE,
           r2 = 0.001,
           kb = 1000
         )
        
        if (!is.null(exp_dat) && nrow(exp_dat) > 0) {
          cat("    SUCCESS: Found", nrow(exp_dat), "instruments\n")
          return(exp_dat)
        }
      }
    }
    
    # 方法2: 尝试直接查询eqtl-a格式
    cat("    Trying direct eQTL query...\n")
    
    # 首先获取基因ID映射
    # 使用genes()函数或手动映射常用基因
    gene_to_ensg <- list(
      "NFKB1" = "ENSG00000109332",
      "FDX1" = "ENSG00000137730",
      "STAT3" = "ENSG00000168610",
      "HIF1A" = "ENSG00000100644",
      "HMOX1" = "ENSG00000100292",
      "GPX4" = "ENSG00000167468",
      "TNF" = "ENSG00000232810",
      "IL6" = "ENSG00000136244",
      "AGER" = "ENSG00000204305"
    )
    
    if (gene_symbol %in% names(gene_to_ensg)) {
      ensg_id <- gene_to_ensg[[gene_symbol]]
      eqtl_id <- paste0("eqtl-a-", ensg_id)
      
      cat("    Trying ID:", eqtl_id, "\n")
      
      exp_dat <- extract_instruments(
        outcomes = eqtl_id,
        p1 = 5e-08,
        clump = TRUE,
        r2 = 0.001,
        kb = 1000,
        access_token = Sys.getenv("OPENGWAS_JWT")
      )
      
      if (!is.null(exp_dat) && nrow(exp_dat) > 0) {
        cat("    SUCCESS: Found", nrow(exp_dat), "instruments\n")
        return(exp_dat)
      }
    }
    
    cat("    No GTEx data found for", gene_symbol, "\n")
    return(NULL)
    
  }, error = function(e) {
    cat("    ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# 测试获取一个基因
cat("\nTesting data retrieval for NFKB1...\n")
test_dat <- get_gene_eqtl("NFKB1", brain_tissue)

if (is.null(test_dat)) {
  cat("\nWARNING: Could not retrieve GTEx data directly.\n")
  cat("Trying alternative approach using local eQTLgen with brain-specific annotations...\n")
  cat("\nNOTE: For true brain-specific MR, you need:\n")
  cat("  1. Download GTEx V8 eQTL data from https://gtexportal.org/\n")
  cat("  2. Or use pre-processed GTEx brain eQTL files\n")
  cat("  3. Current script will fall back to eQTLgen with brain gene annotation\n\n")
  
  use_local <- TRUE
} else {
  cat("\nSUCCESS: GTEx brain data accessible!\n")
  use_local <- FALSE
}

# =============================================================================
# 3. 读取结局数据 (FinnGen)
# =============================================================================
cat("\n----------------------------------------------------------------------\n", sep = "")
cat("Step 2: Reading Outcome Data (FinnGen)\n")
cat("----------------------------------------------------------------------\n", sep = "")

if (!file.exists(outcome_file)) stop("Outcome file not found: ", outcome_file)

cat("  Reading outcome...\n")
outcome_data <- fread(outcome_file, select = c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt"))
setnames(outcome_data, 
         old = c("rsids", "ref", "alt", "beta", "sebeta", "pval", "af_alt"),
         new = c("SNP", "effect_allele.outcome", "other_allele.outcome", 
                 "beta.outcome", "se.outcome", "pval.outcome", "eaf.outcome"))

outcome_data$samplesize.outcome <- 377277
outcome_data$ncase.outcome <- 19862
outcome_data$ncontrol.outcome <- 357415
outcome_data$phenotype <- "Ischemic_Stroke_FinnGen_R12_I9_STR"

cat("  Outcome:", nrow(outcome_data), "SNPs\n")
cat("  Sample size:", format(outcome_data$samplesize.outcome[1], big.mark = ","), "\n\n")

# =============================================================================
# 4. MR分析主函数
# =============================================================================
run_mr_brain <- function(gene_symbol) {
  
  result <- list(
    gene = gene_symbol,
    success = FALSE,
    nsnp_exp = 0, nsnp_outcome = 0, nsnp_final = 0,
    beta = NA, se = NA, or = NA, or_l = NA, or_u = NA, p = NA, method = NULL,
    tissue = brain_tissue,
    error = NULL
  )
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Gene:", gene_symbol, "| Tissue:", brain_tissue, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  tryCatch({
    # 1. 获取暴露数据
    cat("\n[1/5] Fetching GTEx eQTL...\n")
    
    if (use_local) {
      # 回退到本地eQTLgen数据
      cat("  Using local eQTLgen data as proxy...\n")
      eqtl_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-08_kb_1000_r2_0.001.xlsx"
      eqtl_all <- readxl::read_excel(eqtl_file)
      exposure <- eqtl_all[eqtl_all$gene == gene_symbol, ]
      
      if (nrow(exposure) == 0) {
        result$error <- "No eQTL SNPs in local data"
        return(result)
      }
      
      # 格式转换
      exp_dat <- format_data(
        as.data.frame(exposure),
        type = "exposure",
        snp_col = "SNP",
        beta_col = "beta.exposure",
        se_col = "se.exposure",
        pval_col = "pval.exposure",
        eaf_col = "eaf.exposure",
        effect_allele_col = "effect_allele.exposure",
        other_allele_col = "other_allele.exposure"
      )
    } else {
      # 使用在线GTEx数据
      exp_dat <- get_gene_eqtl(gene_symbol, brain_tissue)
      if (is.null(exp_dat) || nrow(exp_dat) == 0) {
        result$error <- "No GTEx instruments found"
        return(result)
      }
    }
    
    result$nsnp_exp <- nrow(exp_dat)
    cat("  eQTL SNPs:", nrow(exp_dat), "\n")
    
    # 2. 提取结局数据
    cat("\n[2/5] Extracting outcome data...\n")
    out_dat <- outcome_data[outcome_data$SNP %in% exp_dat$SNP, ]
    
    if (nrow(out_dat) == 0) {
      result$error <- "No matching outcome SNPs"
      return(result)
    }
    
    out_fmt <- format_data(
      as.data.frame(out_dat),
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      pval_col = "pval.outcome",
      eaf_col = "eaf.outcome",
      effect_allele_col = "effect_allele.outcome",
      other_allele_col = "other_allele.outcome"
    )
    
    result$nsnp_outcome <- nrow(out_fmt)
    cat("  Matched:", nrow(out_fmt), "SNPs\n")
    
    # 3. 数据协调
    cat("\n[3/5] Harmonizing data...\n")
    dat <- harmonise_data(exp_dat, out_fmt, action = 2)
    dat <- dat[dat$mr_keep == TRUE, ]
    result$nsnp_final <- nrow(dat)
    cat("  Valid SNPs:", nrow(dat), "\n")
    
    if (nrow(dat) == 0) {
      result$error <- "No SNPs after harmonization"
      return(result)
    }
    
    # 4. Steiger filtering
    cat("\n[4/5] Steiger filtering...\n")
    tryCatch({
      dat_steiger <- steiger_filtering(dat)
      if (!is.null(dat_steiger) && "steiger_dir" %in% names(dat_steiger)) {
        n_pass <- sum(dat_steiger$steiger_dir == TRUE, na.rm = TRUE)
        cat("  Steiger pass:", n_pass, "/", nrow(dat_steiger), "\n")
        if (n_pass > 0) {
          dat <- dat_steiger[dat_steiger$steiger_dir == TRUE | is.na(dat_steiger$steiger_dir), ]
        }
      }
    }, error = function(e) {
      cat("  WARNING: Steiger error, continuing with all SNPs\n")
    })
    
    # 5. MR分析
    cat("\n[5/5] MR analysis...\n")
    
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
    
    if ("Inverse variance weighted" %in% res$method) {
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
    
    cat("  ", method_name, "OR =", round(result$or, 3),
        "[", round(result$or_l, 3), "-", round(result$or_u, 3), "]",
        "P =", format.pval(result$p, eps = 0.001), "\n")
    
    write.csv(res, file.path(output_dir, paste0("MR_", gene_symbol, "_", brain_tissue, ".csv")), row.names = FALSE)
    result$success <- TRUE
    
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    result$error <- conditionMessage(e)
  })
  
  return(result)
}

# =============================================================================
# 5. 运行所有基因
# =============================================================================
cat("\n======================================================================\n", sep = "")
cat("Running MR Analysis for", length(all_genes), "Genes\n")
cat("Tissue:", brain_tissue, "\n")
cat("======================================================================\n", sep = "")

all_results <- list()
for (gene in all_genes) {
  all_results[[gene]] <- run_mr_brain(gene)
}

# =============================================================================
# 6. 汇总结果
# =============================================================================
cat("\n======================================================================\n", sep = "")
cat("Summary Results -", brain_tissue, "\n")
cat("======================================================================\n", sep = "")

summary_list <- list()
for (gene in names(all_results)) {
  r <- all_results[[gene]]
  
  # 预期判断
  exp_text <- ""
  if (r$success && !is.na(r$or)) {
    if (gene %in% c("NFKB1", "FDX1", "STAT3", "TNF", "IL6", "AGER")) {
      exp_text <- ifelse(r$or > 1, "✓ Risk", "✗ Protective")
    } else if (gene %in% c("GPX4", "HMOX1")) {
      exp_text <- ifelse(r$or < 1, "✓ Protective", "✗ Risk")
    }
  } else {
    exp_text <- ifelse(is.null(r$error), "Failed", r$error)
  }
  
  summary_list[[gene]] <- data.frame(
    Gene = gene,
    Layer = ifelse(gene %in% p0_genes, "P0 (Core)", "P1 (Supp)"),
    Tissue = r$tissue,
    Beta = ifelse(r$success, r$beta, NA),
    SE = ifelse(r$success, r$se, NA),
    OR = ifelse(r$success, r$or, NA),
    OR_L95 = ifelse(r$success, r$or_l, NA),
    OR_U95 = ifelse(r$success, r$or_u, NA),
    P = ifelse(r$success, r$p, NA),
    NSNP_Final = r$nsnp_final,
    Method = ifelse(r$success, r$method, "Failed"),
    Expectation = exp_text,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

# 按层级排序
summary_df$sort_key <- ifelse(summary_df$Layer == "P0 (Core)", 0, 1)
summary_df <- summary_df[order(summary_df$sort_key, -summary_df$OR, na.last = TRUE), ]
summary_df$sort_key <- NULL

print(summary_df[, c("Gene", "Layer", "OR", "P", "NSNP_Final", "Method", "Expectation")])

write.csv(summary_df, file.path(output_dir, "MR_Summary_BrainGTEx.csv"), row.names = FALSE)

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

# 与血液结果对比
if (file.exists("D:/EQTL/MR_v3_20260424_141437/MR_Summary.csv")) {
  cat("\n*** Comparison with Blood eQTLgen ***\n")
  blood_res <- read.csv("D:/EQTL/MR_v3_20260424_141437/MR_Summary.csv")
  
  comparison <- merge(
    summary_df[, c("Gene", "OR", "P")],
    blood_res[, c("Gene", "OR", "P")],
    by = "Gene",
    suffixes = c("_Brain", "_Blood")
  )
  
  cat("Gene\t\tBrain OR\tBlood OR\tDirection Agreement\n")
  for (i in 1:nrow(comparison)) {
    agree <- ifelse(
      (comparison$OR_Brain[i] > 1 & comparison$OR_Blood[i] > 1) |
      (comparison$OR_Brain[i] < 1 & comparison$OR_Blood[i] < 1),
      "YES", "NO"
    )
    cat(sprintf("%-10s\t%.3f\t\t%.3f\t\t%s\n", 
                comparison$Gene[i], 
                comparison$OR_Brain[i], 
                comparison$OR_Blood[i],
                agree))
  }
}

# =============================================================================
# 7. 可视化
# =============================================================================
cat("\nGenerating forest plot...\n")

pd <- summary_df[!is.na(summary_df$OR), ]
if (nrow(pd) > 0) {
  pd$Significance <- ifelse(pd$P < 0.05, "P < 0.05", ifelse(pd$P < 0.1, "P < 0.1", "NS"))
  pd$Gene_Label <- paste0(pd$Gene, " (", pd$NSNP_Final, " SNPs)")
  pd$Label <- sprintf("%.2f [%.2f-%.2f] P=%.3f", pd$OR, pd$OR_L95, pd$OR_U95, pd$P)
  
  p <- ggplot(pd, aes(y = reorder(Gene_Label, OR))) +
    geom_point(aes(x = OR, color = Significance), size = 5, shape = 15) +
    geom_errorbar(aes(xmin = OR_L95, xmax = OR_U95, color = Significance),
                  linewidth = 1, orientation = "y") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_x_log10(breaks = c(0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.5),
                  labels = c("0.70", "0.80", "0.90", "1.00", "1.10", "1.20", "1.30", "1.50")) +
    scale_color_manual(values = c("P < 0.05" = "#D55E00", "P < 0.1" = "#E69F00", "NS" = "#999999")) +
    labs(
      title = paste("MR:", brain_tissue, "eQTL vs Ischemic Stroke"),
      subtitle = paste("GTEx + FinnGen R12 (n=377,277)"),
      x = "Odds Ratio (95% CI, log scale)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    ) +
    geom_text(aes(x = OR_U95 + 0.03, label = Label), hjust = 0, size = 3.5)
  
  ggsave(file.path(output_dir, "MR_Forest_BrainGTEx.png"), p, width = 12, height = 7, dpi = 300, bg = "white")
  cat("Forest plot saved\n")
}

# =============================================================================
# 8. 完成
# =============================================================================
cat("\n======================================================================\n", sep = "")
cat("Analysis Complete!\n")
cat("======================================================================\n", sep = "")
cat("Output Directory:", output_dir, "\n")
cat("Files Generated:\n")
cat("  - MR_Summary_BrainGTEx.csv\n")
cat("  - MR_*_", brain_tissue, ".csv (单基因结果)\n", sep="")
cat("  - MR_Forest_BrainGTEx.png\n")
cat("======================================================================\n", sep = "")
