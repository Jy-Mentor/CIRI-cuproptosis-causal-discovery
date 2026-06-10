#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# MR分析 - 8个Hub基因 vs 脑卒中风险
# 使用真实TwoSampleMR包和IEU GWAS数据
# 注意: 从2024年5月起，IEU GWAS API需要JWT token认证

# 设置工作目录
setwd("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙")

# ============================================
# 获取IEU GWAS Token (必需)
# ============================================
# 1. 访问 https://api.opengwas.io/
# 2. 使用Google账号登录
# 3. 获取JWT token
# 4. 设置环境变量: Sys.setenv(OPENGWAS_JWT = "your_token_here")
# 或使用: ieugwasr::get_access_token()

# 检查token
cat(rep("=", 70), "\n", sep = "")
cat("MR分析准备\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("IEU GWAS API认证状态检查...\n")
cat("从2024年5月起，访问OpenGWAS API需要JWT token。\n")
cat("请按以下步骤获取token:\n")
cat("1. 访问 https://api.opengwas.io/\n")
cat("2. 使用Google账号登录\n")
cat("3. 复制JWT token\n")
cat("4. 在R中设置: Sys.setenv(OPENGWAS_JWT = 'your_token')\n\n")

# 尝试加载包
suppressPackageStartupMessages({
  if (!require("TwoSampleMR", quietly = TRUE)) {
    install.packages("TwoSampleMR", repos = "https://cloud.r-project.org/")
  }
  if (!require("ieugwasr", quietly = TRUE)) {
    install.packages("ieugwasr", repos = "https://cloud.r-project.org/")
  }
  if (!require("dplyr", quietly = TRUE)) {
    install.packages("dplyr", repos = "https://cloud.r-project.org/")
  }
  if (!require("readr", quietly = TRUE)) {
    install.packages("readr", repos = "https://cloud.r-project.org/")
  }
  library(TwoSampleMR)
  library(ieugwasr)
  library(dplyr)
  library(readr)
})

# 尝试获取token
jwt_token <- Sys.getenv("OPENGWAS_JWT")
if (jwt_token == "") {
  cat("警告: 未设置OPENGWAS_JWT环境变量\n")
  cat("尝试交互式获取token...\n")
  tryCatch({
    # 尝试获取token (需要交互式环境)
    jwt_token <- get_access_token()
    cat("成功获取token\n")
  }, error = function(e) {
    cat("无法获取token:", conditionMessage(e), "\n")
    cat("\n请先获取token后再运行此脚本。\n")
    cat("或者使用离线数据模式 (见脚本底部注释)\n")
    quit(status = 1)
  })
} else {
  cat("使用已设置的JWT token\n")
}

# 创建输出目录
output_dir <- "MR_analysis_real_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat(rep("=", 70), "\n", sep = "")
cat("孟德尔随机化分析: 8个Hub基因 vs 脑卒中风险\n")
cat(rep("=", 70), "\n\n", sep = "")

# ============================================
# 基因配置
# ============================================
p0_genes <- c("NFKB1", "FDX1", "STAT3")
p1_genes <- c("HIF1A", "HMOX1", "GPX4", "TNF", "IL6", "AGER")
all_genes <- c(p0_genes, p1_genes)

cat("分析基因列表:\n")
cat("P0层 (必须):", paste(p0_genes, collapse = ", "), "\n")
cat("P1层 (补充):", paste(p1_genes, collapse = ", "), "\n\n")

# ============================================
# 搜索可用的eQTL数据源
# ============================================
cat(rep("-", 70), "\n", sep = "")
cat("步骤1: 搜索可用的eQTL数据源\n")
cat(rep("-", 70), "\n", sep = "")

tryCatch({
  ao <- available_outcomes()
  cat("IEU GWAS数据库中可用研究数:", nrow(ao), "\n")
  
  # 查找eQTL研究
  eqtl_studies <- ao[grep("eqtl", ao$trait, ignore.case = TRUE), ]
  cat("eQTL相关研究:", nrow(eqtl_studies), "个\n")
  
  # 查找特定基因eQTL
  for (gene in head(all_genes, 3)) {
    gene_eqtl <- eqtl_studies[grep(gene, eqtl_studies$trait, ignore.case = TRUE), ]
    if (nrow(gene_eqtl) > 0) {
      cat("\n", gene, "相关eQTL研究:\n")
      print(gene_eqtl[1:min(3, nrow(gene_eqtl)), c("id", "trait", "sample_size")])
    }
  }
}, error = function(e) {
  cat("无法获取eQTL列表:", conditionMessage(e), "\n")
  cat("请检查token是否有效\n")
})

# ============================================
# 定义eQTL和结局数据源
# ============================================
eqtl_source <- "eqtl-a-"  # GTEx v8 eQTL前缀
stroke_outcome <- "ebi-a-GCST006906"  # MEGASTROKE all stroke

cat("\n使用的数据源:\n")
cat("- eQTL: GTEx v8 (", eqtl_source, "*)\n", sep = "")
cat("- 结局: MEGASTROKE全脑卒中 (", stroke_outcome, ")\n\n", sep = "")

# ============================================
# MR分析函数
# ============================================
run_mr_for_gene <- function(gene_symbol) {
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
  
  tryCatch({
    # 1. 提取暴露数据 (eQTL)
    exposure_id <- paste0(eqtl_source, gene_symbol)
    cat("\n1. 提取暴露数据:", exposure_id, "\n")
    
    exposure_dat <- extract_instruments(
      outcomes = exposure_id,
      p1 = 5e-08,
      clump = TRUE,
      r2 = 0.001,
      kb = 10000
    )
    
    if (is.null(exposure_dat) || nrow(exposure_dat) == 0) {
      cat("   未找到eQTL数据\n")
      result$error <- "No eQTL instruments found"
      return(result)
    }
    
    result$exposure_nsnps <- nrow(exposure_dat)
    cat("   找到", nrow(exposure_dat), "个SNP工具变量\n")
    
    # 2. 提取结局数据
    cat("\n2. 提取结局数据 (MEGASTROKE)...\n")
    outcome_dat <- extract_outcome_data(
      snps = exposure_dat$SNP,
      outcomes = stroke_outcome
    )
    
    if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
      cat("   未找到结局数据\n")
      result$error <- "No outcome data found"
      return(result)
    }
    
    cat("   找到", nrow(outcome_dat), "个SNP的结局数据\n")
    
    # 3. 数据Harmonization
    cat("\n3. 数据Harmonization...\n")
    dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat
    )
    
    # 过滤有效SNP
    dat <- dat[dat$mr_keep == TRUE, ]
    
    result$harmonized_nsnps <- nrow(dat)
    cat("   Harmonization后:", nrow(dat), "个有效SNP\n")
    
    if (nrow(dat) == 0) {
      result$error <- "No SNPs after harmonization"
      return(result)
    }
    
    # 4. 执行MR分析
    cat("\n4. 执行MR分析...\n")
    
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
    
    # 5. 异质性检验
    if (nrow(dat) >= 3) {
      cat("\n5. 异质性检验...\n")
      het <- mr_heterogeneity(dat)
      result$heterogeneity <- het
      print(het[, c("method", "Q", "Q_df", "Q_pval")])
    }
    
    # 6. 多效性检验
    if ("mr_egger_regression" %in% res$method) {
      cat("\n6. MR-Egger多效性检验...\n")
      pleio <- mr_pleiotropy_test(dat)
      result$pleiotropy <- pleio
      cat("   Intercept:", round(pleio$egger_intercept, 4), 
          "(p =", format.pval(pleio$pval, eps = 0.001), ")\n")
    }
    
    # 7. 单个SNP分析
    cat("\n7. 单个SNP分析:\n")
    single_snp <- mr_singlesnp(dat)
    print(single_snp[, c("SNP", "b", "se", "p")])
    
    # 8. 留一法分析
    if (nrow(dat) > 3) {
      cat("\n8. 留一法分析...\n")
      loo <- mr_leaveoneout(dat)
      print(loo[, c("SNP", "b", "se", "p")])
    }
    
    result$success <- TRUE
    
    # 保存结果
    write.csv(res, file.path(output_dir, paste0("MR_result_", gene_symbol, ".csv")), row.names = FALSE)
    
  }, error = function(e) {
    cat("   错误:", conditionMessage(e), "\n")
    result$error <- conditionMessage(e)
  })
  
  return(result)
}

# ============================================
# 运行MR分析
# ============================================
cat(rep("=", 70), "\n", sep = "")
cat("开始MR分析\n")
cat(rep("=", 70), "\n\n", sep = "")

all_results <- list()

for (gene in all_genes) {
  result <- run_mr_for_gene(gene)
  all_results[[gene]] <- result
  cat("\n")
}

# ============================================
# 汇总结果
# ============================================
cat(rep("=", 70), "\n", sep = "")
cat("MR分析汇总结果\n")
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
      Error = ifelse(is.null(res$error), "Unknown", res$error),
      stringsAsFactors = FALSE
    )
  }
}

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  cat("汇总结果:\n")
  print(summary_df[, c("Gene", "Layer", "OR", "P_value", "NSNP", "Method")])
  
  write.csv(summary_df, file.path(output_dir, "MR_summary_all_genes.csv"), row.names = FALSE)
  cat("\n汇总结果已保存至:", file.path(output_dir, "MR_summary_all_genes.csv"), "\n")
  
  # 显著结果
  sig_results <- summary_df[!is.na(summary_df$P_value) & summary_df$P_value < 0.05, ]
  if (nrow(sig_results) > 0) {
    cat("\n显著结果 (P < 0.05):\n")
    for (i in 1:nrow(sig_results)) {
      direction <- ifelse(sig_results$OR[i] > 1, "风险性", "保护性")
      cat(sig_results$Gene[i], ": OR =", round(sig_results$OR[i], 3), 
          "(", direction, "), P =", format.pval(sig_results$P_value[i], eps = 0.001), "\n")
    }
  }
  
  # P0层验证
  cat("\nP0层核心基因:\n")
  p0_results <- summary_df[summary_df$Layer == "P0 (Core)", ]
  print(p0_results[, c("Gene", "OR", "P_value")])
}

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("分析完成\n")
cat(rep("=", 70), "\n", sep = "")
cat("所有结果保存在:", output_dir, "\n")
cat(rep("=", 70), "\n", sep = "")

# ============================================
# 离线数据模式说明
# ============================================
cat("\n")
cat(rep("-", 70), "\n", sep = "")
cat("离线数据模式说明\n")
cat(rep("-", 70), "\n", sep = "")
cat("如果无法获取IEU GWAS token，可以使用本地数据模式:\n\n")
cat("1. 从以下数据库下载eQTL和GWAS summary statistics:\n")
cat("   - GTEx Portal: https://gtexportal.org/\n")
cat("   - MEGASTROKE: https://www.megastroke.org/\n")
cat("   - GWAS Catalog: https://www.ebi.ac.uk/gwas/\n\n")
cat("2. 使用本地文件格式进行MR分析 (见TwoSampleMR文档)\n")
cat("   - 暴露数据: 包含SNP, effect_allele, other_allele, beta, se, pval, eaf\n")
cat("   - 结局数据: 同上\n\n")
cat("3. 使用read_exposure_data()和read_outcome_data()读取本地文件\n")
cat(rep("-", 70), "\n", sep = "")
