#!/usr/bin/env Rscript
# QC检查脚本 - 5个成功分析的基因

library(data.table)
library(TwoSampleMR)
library(dplyr)

output_dir <- "D:/EQTL/mr_results_p1e-05_fixed"
genes <- c("NFKB1", "STAT3", "HIF1A", "HSPA5", "HMOX1")

cat("========== QC检查 - 5个基因 ==========\n\n")

# 1. 检查异质性和多效性
cat("【1. 异质性与多效性检查】\n\n")

for (gene in genes) {
  cat(paste0("========== ", gene, " ==========\n"))
  
  # 读取harmonised数据
  harm_file <- file.path(output_dir, paste0(gene, "_harmonised_data.csv"))
  if (!file.exists(harm_file)) {
    cat("  文件不存在\n\n")
    next
  }
  
  dat <- fread(harm_file)
  cat("  SNP数:", nrow(dat), "\n")
  
  if (nrow(dat) < 2) {
    cat("  SNP不足，跳过QC\n\n")
    next
  }
  
  # 异质性检验
  hetero <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  if (!is.null(hetero)) {
    ivw_q <- hetero$Q_pval[hetero$method == "Inverse variance weighted"]
    if (length(ivw_q) > 0) {
      cat("  Cochran's Q P值:", format(ivw_q, digits = 3, scientific = TRUE), "\n")
      if (ivw_q < 0.05) {
        cat("  ⚠️ 存在显著异质性\n")
      } else {
        cat("  ✓ 无异质性\n")
      }
    }
  }
  
  # 多效性检验
  pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  if (!is.null(pleio) && !is.na(pleio$pval)) {
    cat("  Egger intercept P值:", format(pleio$pval, digits = 3), "\n")
    if (pleio$pval < 0.05) {
      cat("  ⚠️ 存在定向多效性\n")
    } else {
      cat("  ✓ 无多效性\n")
    }
  } else {
    cat("  Egger intercept: 无法计算 (SNP不足)\n")
  }
  
  cat("\n")
}

# 2. 留一法分析
cat("\n【2. 留一法分析】\n\n")

for (gene in genes) {
  cat(paste0("========== ", gene, " ==========\n"))
  
  harm_file <- file.path(output_dir, paste0(gene, "_harmonised_data.csv"))
  if (!file.exists(harm_file)) {
    cat("  文件不存在\n\n")
    next
  }
  
  dat <- fread(harm_file)
  
  if (nrow(dat) < 2) {
    cat("  SNP不足，无法做留一法\n\n")
    next
  }
  
  # 留一法
  loo <- tryCatch(mr_leaveoneout(dat), error = function(e) NULL)
  
  if (!is.null(loo) && nrow(loo) > 0) {
    # 保存结果
    write.csv(loo, file.path(output_dir, paste0(gene, "_loo.csv")), row.names = FALSE)
    
    # 检查稳定性
    ivw_loo <- loo[loo$method == "Inverse variance weighted", ]
    if (nrow(ivw_loo) > 1) {
      all_beta <- ivw_loo$b[ivw_loo$SNP == "All"]
      if (length(all_beta) > 0) {
        single_loo <- ivw_loo[ivw_loo$SNP != "All", ]
        sign_changes <- sum(sign(single_loo$b) != sign(all_beta), na.rm = TRUE)
        
        if (sign_changes > 0) {
          cat("  ⚠️ 存在", sign_changes, "个离群SNP!\n")
        } else {
          cat("  ✓ 结果稳定\n")
        }
      }
    }
  }
  
  cat("\n")
}

# 3. 汇总结果
cat("\n【3. 汇总结果】\n\n")

summary_df <- data.frame(
  Gene = character(),
  nSNP = integer(),
  Beta = numeric(),
  OR = numeric(),
  P_value = numeric(),
  Q_pval = numeric(),
  Egger_pval = numeric(),
  Reliability = character(),
  stringsAsFactors = FALSE
)

for (gene in genes) {
  harm_file <- file.path(output_dir, paste0(gene, "_harmonised_data.csv"))
  if (!file.exists(harm_file)) next
  
  dat <- fread(harm_file)
  n_snp <- nrow(dat)
  
  # MR结果
  res <- tryCatch(mr(dat), error = function(e) NULL)
  if (!is.null(res)) {
    ivw <- res[res$method == "Inverse variance weighted", ]
    if (nrow(ivw) > 0) {
      beta <- ivw$b[1]
      pval <- ivw$pval[1]
      
      # QC结果
      hetero <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
      q_pval <- NA
      if (!is.null(hetero)) {
        q_pval <- hetero$Q_pval[hetero$method == "Inverse variance weighted"][1]
      }
      
      pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
      egger_pval <- NA
      if (!is.null(pleio)) {
        egger_pval <- pleio$pval[1]
      }
      
      # 可靠性评级
      reliability <- "⭐⭐⭐"
      if (!is.na(q_pval) && q_pval < 0.05) reliability <- "⭐⭐"
      if (!is.na(q_pval) && q_pval < 0.01) reliability <- "⭐"
      
      summary_df <- rbind(summary_df, data.frame(
        Gene = gene,
        nSNP = n_snp,
        Beta = round(beta, 4),
        OR = round(exp(beta), 3),
        P_value = format(pval, digits = 3, scientific = TRUE),
        Q_pval = ifelse(is.na(q_pval), "NA", format(q_pval, digits = 3, scientific = TRUE)),
        Egger_pval = ifelse(is.na(egger_pval), "NA", format(egger_pval, digits = 3)),
        Reliability = reliability
      ))
    }
  }
}

cat("汇总表:\n")
print(summary_df, row.names = FALSE)

cat("\n========== 分析完成 ==========\n")
