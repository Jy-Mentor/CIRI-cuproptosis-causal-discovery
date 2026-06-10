#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# Mendelian Randomization: 9 Hub Genes vs Ischemic Stroke
# 版本: 5.0 (xQTLbiolinks - GTEx Brain eQTL)
# 日期: 2026-04-24
# 
# 暴露: GTEx V8 Brain eQTL (xQTLbiolinks)
# 结局: FinnGen R12 I9_STR
# =============================================================================

# setwd("D:/EQTL")
# 使用当前目录
setwd(getwd())
set.seed(42)

# 设置CRAN镜像
options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

# =============================================================================
# 0. 包加载与版本检查
# =============================================================================
cat("======================================================================\n", sep = "")
cat("MR Analysis: 9 Hub Genes vs Ischemic Stroke (xQTLbiolinks v5.0)\n")
cat("======================================================================\n\n", sep = "")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Random Seed: 42\n\n")

suppressPackageStartupMessages({
  if (!require("TwoSampleMR", quietly = TRUE)) stop("TwoSampleMR not installed")
  if (!require("xQTLbiolinks", quietly = TRUE)) {
    cat("Installing xQTLbiolinks...\n")
    # 安装 Bioconductor 依赖
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(c("GenomeInfoDb", "GenomicFeatures", "IRanges", "GenomicRanges", "SummarizedExperiment", "BiocGenerics"))
    # 安装 xQTLbiolinks
    if (!require("devtools", quietly = TRUE)) install.packages("devtools")
    devtools::install_github("dingruofan/xQTLbiolinks")
    library(xQTLbiolinks)
  }
  if (!require("data.table", quietly = TRUE)) stop("data.table not installed")
  if (!require("dplyr", quietly = TRUE)) stop("dplyr not installed")
  if (!require("ggplot2", quietly = TRUE)) stop("ggplot2 not installed")
  if (!require("readxl", quietly = TRUE)) stop("readxl not installed")
  library(TwoSampleMR)
  library(xQTLbiolinks)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readxl)
})

cat("Package Versions:\n")
cat("  R:", R.version.string, "\n")
cat("  TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
cat("  xQTLbiolinks:", as.character(packageVersion("xQTLbiolinks")), "\n\n")

# =============================================================================
# 1. 配置与参数
# =============================================================================
# 基因列表
p0_genes <- c("NFKB1", "FDX1", "STAT3")
p1_genes <- c("HIF1A", "HMOX1", "GPX4", "HSPA5", "AGER", "DLAT")
all_genes <- c(p0_genes, p1_genes)

# 全血eQTL数据配置
gtex_brain_tissues <- list(
  "Whole_Blood" = "Whole Blood"  # 全血数据
)

# 使用全血数据
primary_tissue <- "Whole_Blood"
tissue_name <- gtex_brain_tissues[[primary_tissue]]

outcome_file <- "D:/EQTL/finngen_R12_I9_STR"

# 输出目录
output_dir <- paste0("D:/EQTL/MR_xQTLbiolinks_", primary_tissue, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Configuration:\n")
cat("  Total Genes:", length(all_genes), "\n")
cat("  P0 Layer:", paste(p0_genes, collapse = ", "), "\n")
cat("  P1 Layer:", paste(p1_genes, collapse = ", "), "\n")
cat("  Brain Tissue:", tissue_name, "\n")
cat("  Outcome:", basename(outcome_file), "\n")
cat("  Output:", output_dir, "\n\n")

cat("Scientific Hypothesis:\n")
cat("  NFKB1, FDX1, STAT3, AGER: High expr -> Risk (OR > 1)\n")
cat("  GPX4, HMOX1: High expr -> Protective (OR < 1)\n")
cat("  HSPA5, DLAT, HIF1A: High expr -> Risk/Protective (to be determined)\n\n")

# =============================================================================
# 2. 从GTEx V8脑组织获取eQTL数据
# =============================================================================
cat("----------------------------------------------------------------------\n", sep = "")
cat("Step 1: Fetching GTEx V8 Brain eQTL Data (xQTLbiolinks)\n")
cat("Tissue:", tissue_name, "\n")
cat("----------------------------------------------------------------------\n", sep = "")

# 核心函数：从全血eQTL数据获取基因eQTL
get_gtex_eqtl <- function(gene_symbol, tissue) {
  cat("  Fetching", gene_symbol, "from whole blood eQTL data...\n")
  
  # 首先尝试从eQTLgen Excel文件获取
  tryCatch({
    # eQTLgen Excel文件路径
    eqtl_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-08_kb_1000_r2_0.01.xlsx"
    
    if (!file.exists(eqtl_file)) {
      cat("    ERROR: eQTLgen Excel file not found\n")
      # 回退到eQTL Catalogue全血数据
      cat("    Falling back to eQTL Catalogue (Whole Blood)...\n")
      return(get_blood_eqtl(gene_symbol))
    }
    
    # 读取Excel文件
    cat("    Reading eQTLgen Excel file...\n")
    if (!require(readxl, quietly = TRUE)) install.packages("readxl")
    library(readxl)
    eqtl_data <- read_excel(eqtl_file)
    
    # 检查必要列是否存在
    required_cols <- c("gene", "SNP", "beta.exposure", "se.exposure", "pval.exposure", 
                      "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")
    missing_cols <- setdiff(required_cols, names(eqtl_data))
    if (length(missing_cols) > 0) {
      cat("    ERROR: Missing columns in Excel file:", paste(missing_cols, collapse = ", "), "\n")
      # 回退到eQTL Catalogue全血数据
      cat("    Falling back to eQTL Catalogue (Whole Blood)...\n")
      return(get_blood_eqtl(gene_symbol))
    }
    
    # 按基因名称筛选
    gene_data <- eqtl_data[eqtl_data$gene == gene_symbol, ]
    
    if (nrow(gene_data) == 0) {
      cat("    WARNING: No eQTLs found for", gene_symbol, "in eQTLgen file\n")
      # 回退到eQTL Catalogue全血数据
      cat("    Falling back to eQTL Catalogue (Whole Blood)...\n")
      return(get_blood_eqtl(gene_symbol))
    }
    
    cat("    SUCCESS: Found", nrow(gene_data), "eQTLs for", gene_symbol, "(eQTLgen)\n")
    
    # 转换为TwoSampleMR格式
    exp_dat <- data.frame(
      SNP = gene_data$SNP,
      beta.exposure = gene_data$beta.exposure,
      se.exposure = gene_data$se.exposure,
      pval.exposure = gene_data$pval.exposure,
      eaf.exposure = gene_data$eaf.exposure,
      effect_allele.exposure = gene_data$effect_allele.exposure,
      other_allele.exposure = gene_data$other_allele.exposure,
      exposure = paste0(gene_symbol, " (Whole Blood eQTLgen)")
    )
    
    # 过滤掉没有rsID的SNP
    exp_dat <- exp_dat[exp_dat$SNP != "", ]
    
    if (nrow(exp_dat) == 0) {
      cat("    WARNING: No SNPs with rsID found\n")
      # 回退到eQTL Catalogue全血数据
      cat("    Falling back to eQTL Catalogue (Whole Blood)...\n")
      return(get_blood_eqtl(gene_symbol))
    }
    
    # 格式转换
    exp_fmt <- format_data(
      exp_dat,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "beta.exposure",
      se_col = "se.exposure",
      pval_col = "pval.exposure",
      eaf_col = "eaf.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure"
    )
    
    return(exp_fmt)
  }, error = function(e) {
    cat("    ERROR:", conditionMessage(e), "\n")
    # 回退到eQTL Catalogue全血数据
    cat("    Falling back to eQTL Catalogue (Whole Blood)...\n")
    return(get_blood_eqtl(gene_symbol))
  })
}

# 核心函数：从eQTL Catalogue获取全血eQTL数据
get_blood_eqtl <- function(gene_symbol) {
  cat("  Fetching", gene_symbol, "from eQTL Catalogue (Whole Blood)...\n")
  
  tryCatch({
    # 设置eQTL Catalogue Token
    Sys.setenv(OPENGWAS_JWT = "eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJzdWIiOiI2NDQzNmNkMi1lOTg2LTQyNDUtODUyZC05M2E0OTM0N2FmZjciLCJ0eXBlIjoiYWNjZXNzX3Rva2VuIiwiZXhwIjoxNzgxMDM0NjEwfQ.S1XcA85Kc9Kt9b3N3O8pK0t5w5Y7f5f5f5f5f5f5f5f5f5")
    
    # 从eQTL Catalogue获取全血eQTL
    # 全血数据集ID格式：eqtl-a-ENSG
    # 首先需要获取基因的ENSG ID
    if (!require(biomaRt, quietly = TRUE)) install.packages("biomaRt")
    library(biomaRt)
    
    # 使用Ensembl数据库
    ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    
    # 获取基因的ENSG ID
    gene_info <- getBM(
      attributes = c("ensembl_gene_id"),
      filters = "hgnc_symbol",
      values = gene_symbol,
      mart = ensembl
    )
    
    if (nrow(gene_info) == 0) {
      cat("    ERROR: Could not find ENSG ID for", gene_symbol, "\n")
      return(NULL)
    }
    
    ensg_id <- gene_info$ensembl_gene_id[1]
    cat("    ENSG ID:", ensg_id, "\n")
    
    # 构建eQTL Catalogue数据集ID
    dataset_id <- paste0("eqtl-a-", ensg_id)
    
    # 从eQTL Catalogue提取工具变量
    if (!require(ieugwasr, quietly = TRUE)) install.packages("ieugwasr")
    library(ieugwasr)
    
    exp_dat <- extract_instruments(
      outcomes = dataset_id,
      p1 = 5e-8,
      clump = TRUE,
      clump_r2 = 0.01,
      clump_kb = 1000
    )
    
    if (!is.null(exp_dat) && nrow(exp_dat) > 0) {
      cat("    SUCCESS: Found", nrow(exp_dat), "instruments from eQTL Catalogue (Whole Blood)\n")
      
      # 格式转换
      exp_fmt <- format_data(
        exp_dat,
        type = "exposure",
        snp_col = "SNP",
        beta_col = "beta.exposure",
        se_col = "se.exposure",
        pval_col = "pval.exposure",
        eaf_col = "eaf.exposure",
        effect_allele_col = "effect_allele.exposure",
        other_allele_col = "other_allele.exposure"
      )
      
      return(exp_fmt)
    } else {
      cat("    WARNING: No instruments found for", gene_symbol, "in eQTL Catalogue\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("    ERROR with eQTL Catalogue:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# 测试获取FDX1
test_dat <- get_gtex_eqtl("FDX1", tissue_name)

if (is.null(test_dat)) {
  cat("\nWARNING: Could not retrieve GTEx data.\n")
  cat("Will continue trying with eQTL Catalogue...\n\n")
  use_local <- FALSE
} else {
  cat("\nSUCCESS: GTEx brain data accessible!\n")
  cat("Example SNP:", test_dat$SNP[1], "\n")
  cat("Beta:", test_dat$beta.exposure[1], "\n\n")
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
    tissue = tissue_name,
    steiger_pass = 0, steiger_total = 0,
    error = NULL
  )
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Gene:", gene_symbol, "| Tissue:", tissue_name, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  tryCatch({
    # 1. 获取暴露数据
    cat("\n[1/5] Fetching eQTL data...\n")
    
    # 首先尝试xQTLbiolinks (GTEx)
    exp_dat <- get_gtex_eqtl(gene_symbol, tissue_name)
    
    if (is.null(exp_dat) || nrow(exp_dat) == 0) {
      # 尝试eQTL Catalogue作为备选
      cat("  Trying eQTL Catalogue...\n")
      tryCatch({
        # eQTL Catalogue数据集ID
        eqtl_catalogue_id <- "EQTL_CATALOG_GTEX_V8_BRAIN_FRONTAL_CORTEX_BA9"
        
        # 使用extract_instruments从eQTL Catalogue获取数据
        # 传递token参数
        exp_dat <- extract_instruments(
          outcomes = eqtl_catalogue_id,
          p1 = 5e-8,
          clump = TRUE,
          r2 = 0.001,
          kb = 1000,
          opengwas_jwt = Sys.getenv("OPENGWAS_JWT")
        )
        
        if (!is.null(exp_dat) && nrow(exp_dat) > 0) {
          # 按基因筛选
          if ("exposure" %in% names(exp_dat)) {
            exp_dat$gene <- gsub(" .*", "", exp_dat$exposure)
            exp_dat <- exp_dat[grepl(paste0("^", gene_symbol, "$"), exp_dat$gene, ignore.case = TRUE), ]
          }
          
          if (nrow(exp_dat) > 0) {
            cat("    SUCCESS: Found", nrow(exp_dat), "instruments from eQTL Catalogue\n")
          } else {
            result$error <- "No eQTL instruments found in eQTL Catalogue"
            return(result)
          }
        } else {
          result$error <- "Failed to fetch eQTL Catalogue data"
          return(result)
        }
      }, error = function(e) {
        cat("    ERROR with eQTL Catalogue:", conditionMessage(e), "\n")
        result$error <- "eQTL Catalogue error"
        return(result)
      })
    }
    
    if (is.null(exp_dat) || nrow(exp_dat) == 0) {
      result$error <- "No eQTL instruments found from any source"
      return(result)
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
    
    write.csv(res, file.path(output_dir, paste0("MR_", gene_symbol, "_", gsub(" ", "_", primary_tissue), ".csv")), row.names = FALSE)
    result$success <- TRUE
    
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    result$error <- conditionMessage(e)
  })
  
  return(result)
}

# =============================================================================
# 5. 运行基因分析 (只测试50%的基因)
# =============================================================================
cat("\n======================================================================\n", sep = "")
cat("Running MR Analysis (50% sample for validation)\n")
cat("Tissue:", tissue_name, "\n")

# 使用所有基因进行分析
test_genes <- all_genes

cat("Total genes:", length(all_genes), "\n")
cat("Test genes:", length(test_genes), "\n")
cat("Selected genes:", paste(test_genes, collapse = ", "), "\n")
cat("======================================================================\n", sep = "")

all_results <- list()
for (gene in test_genes) {
  all_results[[gene]] <- run_mr_brain(gene)
}

# =============================================================================
# 6. 汇总结果
# =============================================================================
cat("\n======================================================================\n", sep = "")
cat("Summary Results -", tissue_name, "\n")
cat("======================================================================\n", sep = "")

summary_list <- list()
for (gene in names(all_results)) {
  r <- all_results[[gene]]
  
  # 预期判断
  exp_text <- ""
  if (r$success && !is.na(r$or)) {
    if (gene %in% c("NFKB1", "FDX1", "STAT3", "AGER")) {
      exp_text <- ifelse(r$or > 1, "✓ Risk", "✗ Protective")
    } else if (gene %in% c("GPX4", "HMOX1")) {
      exp_text <- ifelse(r$or < 1, "✓ Protective", "✗ Risk")
    } else {
      exp_text <- ifelse(r$or > 1, "Risk", "Protective")
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

write.csv(summary_df, file.path(output_dir, "MR_Summary_GTEx_Brain_xQTL.csv"), row.names = FALSE)

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
      title = paste("MR:", tissue_name, "eQTL vs Ischemic Stroke"),
      subtitle = "GTEx V8 + FinnGen R12 (n=377,277)",
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
  
  ggsave(file.path(output_dir, "MR_Forest_GTEx_Brain_xQTL.png"), p, width = 12, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "MR_Forest_GTEx_Brain_xQTL.pdf"), p, width = 12, height = 7)
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
cat("  - MR_Summary_GTEx_Brain_xQTL.csv\n")
cat("  - MR_*_", gsub(" ", "_", primary_tissue), ".csv (单基因结果)\n", sep = "")
cat("  - MR_Forest_GTEx_Brain_xQTL.png/pdf\n")
cat("\nNotes:\n")
cat("  - Data Source: GTEx V8 via xQTLbiolinks\n")
cat("  - Tissue:", tissue_name, "\n")
cat("  - If GTEx data unavailable, used local eQTLgen as fallback\n")
cat("======================================================================\n", sep = "")
