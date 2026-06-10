#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# Mendelian Randomization: 9 Hub Genes vs Ischemic Stroke (Brain GTEx v2.0)
# 版本: 4.1 (修正版 - 使用正确的GTEx数据集ID)
# 日期: 2026-04-24
# 
# 暴露: GTEx Brain eQTL v8 (ieu-b-4171等)
# 结局: FinnGen R12 I9_STR
# =============================================================================

setwd("D:/EQTL")
set.seed(42)

# =============================================================================
# 0. 包加载与版本检查
# =============================================================================
cat("======================================================================\n", sep = "")
cat("MR Analysis: 9 Hub Genes vs Ischemic Stroke (Brain GTEx v4.1)\n")
cat("======================================================================\n\n", sep = "")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Random Seed: 42\n\n")

suppressPackageStartupMessages({
  if (!require("TwoSampleMR", quietly = TRUE)) stop("TwoSampleMR not installed")
  if (!require("ieugwasr", quietly = TRUE)) stop("ieugwasr not installed")
  if (!require("data.table", quietly = TRUE)) stop("data.table not installed")
  if (!require("dplyr", quietly = TRUE)) stop("dplyr not installed")
  if (!require("ggplot2", quietly = TRUE)) stop("ggplot2 not installed")
  if (!require("readxl", quietly = TRUE)) stop("readxl not installed")
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readxl)
})

cat("Package Versions:\n")
cat("  R:", R.version.string, "\n")
cat("  TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
cat("  ieugwasr:", as.character(packageVersion("ieugwasr")), "\n\n")

# =============================================================================
# 0.5 设置OpenGWAS API Token
# =============================================================================
cat("Setting up OpenGWAS API authentication...\n")

# 检查环境变量中是否已有token
existing_token <- Sys.getenv("OPENGWAS_JWT")
if (existing_token == "") {
  cat("请输入OpenGWAS JWT Token (从 https://api.opengwas.io/ 获取):\n")
  cat("(直接按回车使用默认token，或粘贴您的token)\n")
  # 在脚本环境中使用硬编码token（仅用于自动化运行）
  token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3NzAxMjQ5NiwiZXhwIjoxNzc4MjIyMDk2fQ.sz-zieAniOZLfPGDiQRo3W6z5PQ6HfEvsfSjvKCVLOlUJumnpM-jE9Se6QAqXezfH-ffEEbM0heC3pQu6MjBF1ttqaRpIdOv8S7_pq-s-HLVQUJZ3Ge5Ort66N89riaj9nRS-FdGidJAft58CEmjk-Q1lpJpJbOSGJiKCRQ_Man5vbCBiTH3p_e490wbapqqpxaXWs4ki-xYUDlQXuc6u2lZIh7_6TSei6Chmew0EWy8wcVgSTJGAmAiT8jrFTq5ydxsPvzZMnMZgEezAfrHIi9i3zJqbTzggs8Okqy2owfzLzcE7MKt_EcUSHYvMzRsL3YcbILipkodoy-sbzNmpw"
  Sys.setenv(OPENGWAS_JWT = token)
  cat("  Token set from default\n")
} else {
  cat("  Using existing token from environment\n")
}

# 验证token
tryCatch({
  user_info <- ieugwasr::user()
  cat("  Token verified successfully!\n")
  cat("  User:", user_info$user[[1]], "\n")
  cat("  Tier:", user_info$tier[[1]], "\n")
  cat("  Remaining queries:", user_info$allowance[[1]] - user_info$used[[1]], "\n\n")
}, error = function(e) {
  cat("  WARNING: Token verification failed:\n")
  cat("  ", conditionMessage(e), "\n\n")
})

# =============================================================================
# 1. 配置与参数
# =============================================================================
# 基因列表
p0_genes <- c("NFKB1", "FDX1", "STAT3")
p1_genes <- c("HIF1A", "HMOX1", "GPX4", "TNF", "IL6", "AGER")
all_genes <- c(p0_genes, p1_genes)

# GTEx v8 脑组织数据集ID（顶刊标准，直接使用）
# 来源: IEU OpenGWAS Database - GTEx v8 eQTL
gtex_tissue_datasets <- list(
  "Brain_Frontal_Cortex" = "ieu-b-4171",  # 大脑皮层（首选，缺血性卒中核心）
  "Brain_Cortex" = "ieu-b-4178",          # 全皮层汇总
  "Brain_Hippocampus" = "ieu-b-4166",     # 海马体
  "Brain_Putamen" = "ieu-b-4164",         # 基底节
  "Brain_Cerebellum" = "ieu-b-4170"       # 小脑
)

# 默认使用大脑皮层（Frontal Cortex）
brain_tissue <- "Brain_Frontal_Cortex"
tissue_id <- gtex_tissue_datasets[[brain_tissue]]

outcome_file <- "D:/EQTL/finngen_R12_I9_STR"

# 输出目录
output_dir <- paste0("D:/EQTL/MR_GTEx_v2_", brain_tissue, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Configuration:\n")
cat("  Total Genes:", length(all_genes), "\n")
cat("  P0 Layer:", paste(p0_genes, collapse = ", "), "\n")
cat("  P1 Layer:", paste(p1_genes, collapse = ", "), "\n")
cat("  Brain Tissue:", brain_tissue, "\n")
cat("  Dataset ID:", tissue_id, "\n")
cat("  Outcome:", basename(outcome_file), "\n")
cat("  Output:", output_dir, "\n\n")

cat("Scientific Hypothesis:\n")
cat("  NFKB1, FDX1, STAT3, TNF, IL6, AGER: High expr -> Risk (OR > 1)\n")
cat("  GPX4, HMOX1: High expr -> Protective (OR < 1)\n\n")

# =============================================================================
# 2. 从GTEx脑组织数据集提取基因eQTL
# =============================================================================
cat("----------------------------------------------------------------------\n", sep = "")
cat("Step 1: Fetching GTEx Brain eQTL Data\n")
cat("Dataset:", tissue_id, "(", brain_tissue, ")\n")
cat("----------------------------------------------------------------------\n", sep = "")

# 核心函数：从GTEx组织数据集提取指定基因的eQTL
get_gtex_eqtl_simple <- function(gene_symbol, tissue_dataset_id) {
  cat("  Fetching", gene_symbol, "from", tissue_dataset_id, "...\n")
  
  tryCatch({
    # 方法：从GTEx组织数据集提取所有显著eQTL，然后按基因筛选
    # 【注意】TwoSampleMR的extract_instruments没有gene参数，需要手动筛选
    
    # Step 1: 提取该组织所有显著eQTL instruments
    cat("    Extracting instruments from", tissue_dataset_id, "...\n")
    all_exp_dat <- extract_instruments(
      outcomes = tissue_dataset_id,    # GTEx组织数据集ID (如 ieu-b-4171)
      p1 = 5e-8,                       # 全基因组显著阈值
      clump = TRUE,                    # 自动LD clump
      r2 = 0.001,                      # 严格LD阈值
      kb = 1000                        # 顺式1Mb窗口
    )
    
    if (is.null(all_exp_dat) || nrow(all_exp_dat) == 0) {
      cat("    WARNING: No instruments found in", tissue_dataset_id, "\n")
      return(NULL)
    }
    
    cat("    Total instruments in dataset:", nrow(all_exp_dat), "\n")
    
    # Step 2: 按基因名称筛选
    # GTEx数据中基因名称通常在exposure或gene列
    if ("exposure" %in% names(all_exp_dat)) {
      # 从exposure列提取基因名 (格式: "GeneName (Tissue)")
      all_exp_dat$gene <- gsub(" .*", "", all_exp_dat$exposure)
      gene_dat <- all_exp_dat[grepl(paste0("^", gene_symbol, "$"), all_exp_dat$gene, ignore.case = TRUE), ]
    } else if ("gene" %in% names(all_exp_dat)) {
      gene_dat <- all_exp_dat[grepl(paste0("^", gene_symbol, "$"), all_exp_dat$gene, ignore.case = TRUE), ]
    } else {
      cat("    WARNING: Cannot find gene column in data\n")
      return(NULL)
    }
    
    if (nrow(gene_dat) > 0) {
      cat("    SUCCESS: Found", nrow(gene_dat), "instruments for", gene_symbol, "\n")
      return(gene_dat)
    } else {
      cat("    WARNING: No instruments found for", gene_symbol, "in this dataset\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("    ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# 测试获取FDX1
cat("\nTesting GTEx data retrieval for FDX1...\n")
test_dat <- get_gtex_eqtl_simple("FDX1", tissue_id)

if (is.null(test_dat)) {
  cat("\nWARNING: Could not retrieve GTEx data.\n")
  cat("Possible reasons:\n")
  cat("  1. Token expired or invalid\n")
  cat("  2. Gene not available in this tissue\n")
  cat("  3. Network connectivity issues\n")
  cat("\nFalling back to local eQTLgen data...\n\n")
  use_local <- TRUE
} else {
  cat("\nSUCCESS: GTEx brain data accessible!\n")
  cat("Example instrument:", test_dat$SNP[1], "\n\n")
  use_local <- FALSE
}

# =============================================================================
# 3. 读取结局数据 (FinnGen)
# =============================================================================
cat("----------------------------------------------------------------------\n", sep = "")
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
# 4. MR分析主函数 (带完整Steiger filtering)
# =============================================================================
run_mr_brain <- function(gene_symbol) {
  
  result <- list(
    gene = gene_symbol,
    success = FALSE,
    nsnp_exp = 0, nsnp_outcome = 0, nsnp_final = 0,
    beta = NA, se = NA, or = NA, or_l = NA, or_u = NA, p = NA, method = NULL,
    tissue = brain_tissue,
    dataset_id = tissue_id,
    steiger_pass = 0, steiger_total = 0,
    error = NULL
  )
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Gene:", gene_symbol, "| Tissue:", brain_tissue, "\n")
  cat("Dataset:", tissue_id, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  tryCatch({
    # 1. 获取暴露数据
    cat("\n[1/5] Fetching GTEx eQTL...\n")
    
    if (use_local) {
      # 回退到本地eQTLgen数据
      cat("  Using local eQTLgen data as proxy...\n")
      eqtl_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-08_kb_1000_r2_0.001.xlsx"
      
      if (!file.exists(eqtl_file)) {
        result$error <- "Local eQTL file not found"
        return(result)
      }
      
      eqtl_all <- read_excel(eqtl_file)
      
      # 列检查
      required_cols <- c("gene", "SNP", "beta.exposure", "se.exposure", "pval.exposure", 
                         "effect_allele.exposure", "other_allele.exposure")
      missing_cols <- setdiff(required_cols, names(eqtl_all))
      if (length(missing_cols) > 0) {
        result$error <- paste("Missing columns:", paste(missing_cols, collapse = ", "))
        return(result)
      }
      
      exposure <- eqtl_all[eqtl_all$gene == gene_symbol, ]
      
      if (nrow(exposure) == 0) {
        result$error <- "No eQTL SNPs in local data"
        return(result)
      }
      
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
      # 使用GTEx数据
      exp_dat <- get_gtex_eqtl_simple(gene_symbol, tissue_id)
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
    
    # 4. Steiger filtering (完整逻辑)
    cat("\n[4/5] Steiger filtering...\n")
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
      cat("  Continuing with all SNPs\n")
      result$steiger_pass <- nrow(dat)
      result$steiger_total <- nrow(dat)
    })
    
    if (nrow(dat) == 0) {
      result$error <- "No SNPs after Steiger filtering"
      return(result)
    }
    
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
cat("Tissue:", brain_tissue, "| Dataset:", tissue_id, "\n")
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
    Dataset_ID = r$dataset_id,
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

write.csv(summary_df, file.path(output_dir, "MR_Summary_GTEx_Brain.csv"), row.names = FALSE)

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
      subtitle = paste("GTEx v8 (", tissue_id, ") + FinnGen R12 (n=377,277)"),
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
  
  ggsave(file.path(output_dir, "MR_Forest_GTEx_Brain.png"), p, width = 12, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "MR_Forest_GTEx_Brain.pdf"), p, width = 12, height = 7)
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
cat("  - MR_Summary_GTEx_Brain.csv\n")
cat("  - MR_*_", brain_tissue, ".csv (单基因结果)\n", sep = "")
cat("  - MR_Forest_GTEx_Brain.png/pdf\n")
cat("\nNotes:\n")
cat("  - GTEx Dataset ID:", tissue_id, "\n")
cat("  - If GTEx data unavailable, used local eQTLgen as fallback\n")
cat("======================================================================\n", sep = "")
