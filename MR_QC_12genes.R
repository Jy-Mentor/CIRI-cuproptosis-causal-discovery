#!/usr/bin/env Rscript
# QC检查脚本 - 12个成功分析的基因

library(data.table)
library(TwoSampleMR)
library(dplyr)

output_dir <- "D:/EQTL/mr_results_p1e-05_fixed"
genes <- c("NFKB1", "STAT3", "HIF1A", "HSPA5", "HMOX1", 
           "NFE2L2", "LIAS", "IKBKB", "PARP1", "CASP8", 
           "MTOR", "PTPRC")

cat("========== QC检查 - 12个基因 ==========\n\n")

# 存储汇总结果
summary_list <- list()

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
  n_snp <- nrow(dat)
  cat("  SNP数:", n_snp, "\n")
  
  if (n_snp < 2) {
    cat("  SNP不足，跳过QC\n\n")
    next
  }
  
  # 初始化结果
  q_pval <- NA
  egger_pval <- NA
  hetero_status <- "NA"
  pleio_status <- "NA"
  
  # 异质性检验
  hetero <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  if (!is.null(hetero)) {
    ivw_q <- hetero$Q_pval[hetero$method == "Inverse variance weighted"]
    if (length(ivw_q) > 0) {
      q_pval <- ivw_q[1]
      cat("  Cochran's Q P值:", format(q_pval, digits = 3, scientific = TRUE), "\n")
      if (!is.na(q_pval)) {
        if (q_pval < 0.05) {
          cat("  ⚠️ 存在显著异质性\n")
          hetero_status <- "Significant"
        } else if (q_pval < 0.1) {
          cat("  ⚡ 边缘异质性\n")
          hetero_status <- "Borderline"
        } else {
          cat("  ✓ 无异质性\n")
          hetero_status <- "None"
        }
      }
    }
  }
  
  # 多效性检验
  pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  if (!is.null(pleio) && !is.na(pleio$pval)) {
    egger_pval <- pleio$pval[1]
    cat("  Egger intercept P值:", format(egger_pval, digits = 3), "\n")
    if (egger_pval < 0.05) {
      cat("  ⚠️ 存在定向多效性\n")
      pleio_status <- "Significant"
    } else {
      cat("  ✓ 无多效性\n")
      pleio_status <- "None"
    }
  } else {
    cat("  Egger intercept: 无法计算 (SNP不足)\n")
  }
  
  # MR结果
  res <- tryCatch(mr(dat, method_list = c("mr_ivw")), error = function(e) NULL)
  beta <- NA
  or <- NA
  pval <- NA
  if (!is.null(res)) {
    ivw <- res[res$method == "Inverse variance weighted", ]
    if (nrow(ivw) > 0) {
      beta <- ivw$b[1]
      or <- exp(beta)
      pval <- ivw$pval[1]
    }
  }
  
  # 保存结果
  summary_list[[gene]] <- data.frame(
    Gene = gene,
    nSNP = n_snp,
    Beta = ifelse(is.na(beta), NA, round(beta, 4)),
    OR = ifelse(is.na(or), NA, round(or, 3)),
    P_value = ifelse(is.na(pval), NA, format(pval, digits = 3, scientific = TRUE)),
    Q_pval = ifelse(is.na(q_pval), "NA", format(q_pval, digits = 3, scientific = TRUE)),
    Hetero_status = hetero_status,
    Egger_pval = ifelse(is.na(egger_pval), "NA", format(egger_pval, digits = 3)),
    Pleio_status = pleio_status,
    stringsAsFactors = FALSE
  )
  
  cat("\n")
}

# 2. 留一法分析
cat("\n【2. 留一法分析】\n\n")

loo_summary <- list()

for (gene in genes) {
  cat(paste0("========== ", gene, " ==========\n"))
  
  harm_file <- file.path(output_dir, paste0(gene, "_harmonised_data.csv"))
  if (!file.exists(harm_file)) {
    cat("  文件不存在\n\n")
    loo_summary[[gene]] <- data.frame(Gene = gene, Stability = "File missing", Outlier_SNPs = "", stringsAsFactors = FALSE)
    next
  }
  
  dat <- fread(harm_file)
  
  if (nrow(dat) < 2) {
    cat("  SNP不足，无法做留一法\n\n")
    loo_summary[[gene]] <- data.frame(Gene = gene, Stability = "NA (SNP<2)", Outlier_SNPs = "", stringsAsFactors = FALSE)
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
          outlier_snps <- paste(single_loo$SNP[sign(single_loo$b) != sign(all_beta)], collapse = ", ")
          cat("  离群SNP:", outlier_snps, "\n")
          loo_summary[[gene]] <- data.frame(Gene = gene, Stability = "UNSTABLE", Outlier_SNPs = outlier_snps, stringsAsFactors = FALSE)
        } else {
          cat("  ✓ 结果稳定，无离群SNP\n")
          # 显示最大影响
          beta_diff <- abs(single_loo$b - all_beta)
          max_diff_idx <- which.max(beta_diff)
          cat("  最大影响SNP:", single_loo$SNP[max_diff_idx], 
              "(Beta变化:", round(beta_diff[max_diff_idx], 4), ")\n")
          loo_summary[[gene]] <- data.frame(Gene = gene, Stability = "Stable", Outlier_SNPs = "", stringsAsFactors = FALSE)
        }
      } else {
        loo_summary[[gene]] <- data.frame(Gene = gene, Stability = "No All row", Outlier_SNPs = "", stringsAsFactors = FALSE)
      }
    } else {
      loo_summary[[gene]] <- data.frame(Gene = gene, Stability = "Insufficient data", Outlier_SNPs = "", stringsAsFactors = FALSE)
    }
  } else {
    cat("  留一法失败\n")
    loo_summary[[gene]] <- data.frame(Gene = gene, Stability = "Failed", Outlier_SNPs = "", stringsAsFactors = FALSE)
  }
  
  cat("\n")
}

# 3. 汇总结果
cat("\n【3. 汇总结果】\n\n")

# 合并汇总表
summary_df <- do.call(rbind, summary_list)
loo_df <- do.call(rbind, loo_summary)

# 添加可靠性评级
summary_df$Reliability <- "⭐⭐⭐"
summary_df$Reliability[summary_df$Hetero_status == "Significant"] <- "⭐"
summary_df$Reliability[summary_df$Hetero_status == "Borderline"] <- "⭐⭐"
summary_df$Reliability[summary_df$Pleio_status == "Significant"] <- "⭐"

# 合并留一法结果
summary_df <- merge(summary_df, loo_df, by = "Gene", all.x = TRUE)

cat("完整QC汇总表:\n")
cat(strrep("=", 120), "\n", sep = "")
print(summary_df, row.names = FALSE)
cat(strrep("=", 120), "\n", sep = "")

# 保存汇总表
write.csv(summary_df, file.path(output_dir, "MR_QC_summary_12genes.csv"), row.names = FALSE)

cat("\n")

# 4. 分层总结
cat("\n【4. 分层总结】\n\n")

cat("⭐⭐⭐ 高可靠性基因 (无异质性，无多效性):\n")
high_rel <- summary_df[summary_df$Reliability == "⭐⭐⭐", ]
if (nrow(high_rel) > 0) {
  cat("  ", paste(high_rel$Gene, collapse = ", "), "\n")
} else {
  cat("  无\n")
}

cat("\n⭐⭐ 中等可靠性基因 (边缘异质性):\n")
med_rel <- summary_df[summary_df$Reliability == "⭐⭐", ]
if (nrow(med_rel) > 0) {
  cat("  ", paste(med_rel$Gene, collapse = ", "), "\n")
} else {
  cat("  无\n")
}

cat("\n⭐ 低可靠性基因 (显著异质性或多效性):\n")
low_rel <- summary_df[summary_df$Reliability == "⭐", ]
if (nrow(low_rel) > 0) {
  cat("  ", paste(low_rel$Gene, collapse = ", "), "\n")
} else {
  cat("  无\n")
}

cat("\n结果不稳定基因 (留一法检测):\n")
unstable <- summary_df[summary_df$Stability == "UNSTABLE", ]
if (nrow(unstable) > 0) {
  for (i in 1:nrow(unstable)) {
    cat("  ", unstable$Gene[i], ":", unstable$Outlier_SNPs[i], "\n")
  }
} else {
  cat("  无\n")
}

cat("\n========== 分析完成 ==========\n")
cat("QC汇总表已保存:", file.path(output_dir, "MR_QC_summary_12genes.csv"), "\n")
