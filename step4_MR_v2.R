# ============================================================================
# 第4步：孟德尔随机化（MR）人群验证
# 验证BCP核心基因与缺血性脑卒中的因果关系
# ============================================================================

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

cat("=== 第4步：孟德尔随机化（MR）人群验证 ===\n\n")

# 检查TwoSampleMR包
if(!require(TwoSampleMR, quietly=TRUE)) {
  cat("安装TwoSampleMR包...\n")
  if(!require(remotes, quietly=TRUE)) install.packages("remotes")
  remotes::install_github("MRCIEU/TwoSampleMR")
  library(TwoSampleMR)
}

if(!require(ieugwasr, quietly=TRUE)) {
  remotes::install_github("MRCIEU/ieugwasr")
  library(ieugwasr)
}

library(dplyr)
library(ggplot2)

result_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/causal_analysis_results"
mr_dir <- file.path(result_dir, "MR_results")
dir.create(mr_dir, showWarnings=FALSE, recursive=TRUE)

# ============================================================================
# 1. 定义分析基因
# ============================================================================
cat("【4.1】定义BCP核心基因...\n")

# 8个BCP核心基因
bcp_genes <- c("IL6", "STAT3", "NFKB1", "TGFB1", "AGER", "PTGS2", "TLR4", "FDX1")

# 检查D:/EQTL中已分析的结果
existing_csv <- "D:/EQTL/MR_10000kb_Results/mr_main_results.csv"
if(file.exists(existing_csv)) {
  existing_res <- read.csv(existing_csv, stringsAsFactors=FALSE)
  existing_genes <- unique(existing_res$gene)
  cat(sprintf("  D:/EQTL中已分析的基因: %s\n", paste(existing_genes, collapse=", ")))
  
  # 合并结果
  new_genes <- setdiff(bcp_genes, existing_genes)
} else {
  existing_genes <- c()
  new_genes <- bcp_genes
}

cat(sprintf("  总共 %d 个BCP基因\n", length(bcp_genes)))
cat(sprintf("  新分析: %s\n", paste(new_genes, collapse=", ")))

# ============================================================================
# 2. MR分析 - 使用IEU GWAS在线数据
# ============================================================================
cat("\n【4.2】执行MR分析（使用IEU GWAS在线数据）...\n")

mr_results_all <- list()

for(gene in new_genes) {
  cat(sprintf("\n  分析基因: %s...\n", gene))
  
  tryCatch({
    # 1. 获取eQTL instruments（血液表达）
    # 尝试不同的eQTL数据源
    eqtl_id <- paste0("eqtl-a-", gene)
    
    cat("    提取eQTL instruments...\n")
    exposure_dat <- extract_instruments(
      outcomes = eqtl_id,
      p1 = 5e-08,
      clump = TRUE,
      r2 = 0.001,
      kb = 10000
    )
    
    if(is.null(exposure_dat) || nrow(exposure_dat) == 0) {
      # 尝试其他eQTL数据集
      eqtl_id2 <- paste0("eqtlgen-", gene)
      exposure_dat <- extract_instruments(
        outcomes = eqtl_id2,
        p1 = 5e-08,
        clump = TRUE,
        r2 = 0.001,
        kb = 10000
      )
      
      if(is.null(exposure_dat) || nrow(exposure_dat) == 0) {
        cat(sprintf("    [跳过] %s 没有找到eQTL instruments\n", gene))
        next
      }
    }
    
    cat(sprintf("    找到 %d 个eQTL SNPs\n", nrow(exposure_dat)))
    
    # 2. 获取结局数据（缺血性脑卒中）
    # 使用FinnGen R12 缺血性脑卒中
    cat("    提取结局数据（缺血性脑卒中）...\n")
    outcome_dat <- extract_outcome_data(
      snps = exposure_dat$SNP,
      outcomes = "finn-b-I9_STR"  # FinnGen缺血性脑卒中
    )
    
    if(is.null(outcome_dat) || nrow(outcome_dat) == 0) {
      cat(sprintf("    [跳过] 没有找到结局数据\n"))
      next
    }
    
    cat(sprintf("    找到 %d 个结局SNPs\n", nrow(outcome_dat)))
    
    # 3. 协调数据
    cat("    协调数据...\n")
    dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat,
      action = 2
    )
    
    if(nrow(dat) < 3) {
      cat(sprintf("    [跳过] SNP数量不足 (%d < 3)\n", nrow(dat)))
      next
    }
    
    cat(sprintf("    协调后: %d 个SNPs\n", nrow(dat)))
    
    # 4. 执行MR分析
    cat("    执行MR分析...\n")
    res <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw", "mr_weighted_median"))
    
    # 5. 敏感性分析
    cat("    敏感性分析...\n")
    het <- mr_heterogeneity(dat)
    pleio <- mr_pleiotropy_test(dat)
    
    # 保存结果
    mr_results_all[[gene]] <- list(
      gene = gene,
      mr_results = res,
      heterogeneity = het,
      pleiotropy = pleio,
      n_snps = nrow(dat)
    )
    
    # 输出结果
    cat("    MR结果:\n")
    for(i in 1:nrow(res)) {
      method <- res$method[i]
      b <- res$b[i]
      p <- res$pval[i]
      sig <- ifelse(p < 0.05, "***", ifelse(p < 0.1, "*", ""))
      cat(sprintf("      %s: b=%.3f, p=%.4f %s\n", method, b, p, sig))
    }
    
    # 保存到文件
    write.csv(res, file.path(mr_dir, paste0("MR_results_", gene, ".csv")), row.names=FALSE)
    cat(sprintf("    ✅ %s 完成\n", gene))
    
  }, error = function(e) {
    cat(sprintf("    [错误] %s: %s\n", gene, conditionMessage(e)))
  })
}

# ============================================================================
# 3. 汇总结果
# ============================================================================
cat("\n【4.3】MR分析汇总...\n")

if(length(mr_results_all) > 0) {
  cat(sprintf("\n  成功分析 %d 个基因:\n\n", length(mr_results_all)))
  
  # 创建汇总表
  summary_list <- list()
  for(gene in names(mr_results_all)) {
    res <- mr_results_all[[gene]]$mr_results
    
    # IVW结果
    ivw_row <- res[res$method == "Inverse variance weighted", ]
    if(nrow(ivw_row) == 0) ivw_row <- res[1, ]  # 如果没有IVW，使用第一行
    
    # Weighted median结果
    wm_row <- res[res$method == "Weighted median", ]
    if(nrow(wm_row) == 0) wm_row <- data.frame(b=NA, pval=NA)
    
    summary_list[[gene]] <- data.frame(
      gene = gene,
      n_snps = mr_results_all[[gene]]$n_snps,
      ivw_beta = ivw_row$b,
      ivw_se = ivw_row$se,
      ivw_p = ivw_row$pval,
      wm_beta = ifelse(nrow(wm_row) > 0, wm_row$b, NA),
      wm_p = ifelse(nrow(wm_row) > 0, wm_row$pval, NA),
      stringsAsFactors = FALSE
    )
  }
  
  summary_df <- do.call(rbind, summary_list)
  
  cat("  MR分析结果汇总:\n")
  print(summary_df, row.names=FALSE)
  
  write.csv(summary_df, file.path(mr_dir, "MR_summary_results.csv"), row.names=FALSE)
  
  # 识别显著基因
  sig_ivw <- summary_df$gene[summary_df$ivw_p < 0.05]
  sig_wm <- summary_df$gene[summary_df$wm_p < 0.05]
  suggestive <- summary_df$gene[summary_df$ivw_p < 0.1 | summary_df$wm_p < 0.1]
  
  if(length(sig_ivw) > 0) {
    cat(sprintf("\n  ✅ IVW显著 (p<0.05): %s\n", paste(sig_ivw, collapse=", ")))
  }
  if(length(sig_wm) > 0) {
    cat(sprintf("  ✅ Weighted median显著 (p<0.05): %s\n", paste(sig_wm, collapse=", ")))
  }
  if(length(suggestive) > 0 && length(sig_ivw) == 0 && length(sig_wm) == 0) {
    cat(sprintf("  ⚠️ 提示性关联 (p<0.1): %s\n", paste(suggestive, collapse=", ")))
  }
  if(length(sig_ivw) == 0 && length(sig_wm) == 0 && length(suggestive) == 0) {
    cat("\n  ⚠️ 没有发现显著因果关联 (p < 0.05)\n")
  }
  
  # 保存所有结果
  saveRDS(mr_results_all, file.path(mr_dir, "all_MR_results.rds"))
  
  # 可视化
  cat("\n【4.4】生成可视化...\n")
  pdf(file.path(mr_dir, "MR_forest_plot.pdf"), width=12, height=8)
  
  res_all <- bind_rows(lapply(mr_results_all, function(x) x$mr_results))
  if(nrow(res_all) > 0) {
    res_all$significant <- ifelse(res_all$pval < 0.05, "p<0.05", "n.s.")
    
    p <- ggplot(res_all, aes(x=exposure, y=b, ymin=b-1.96*se, ymax=b+1.96*se, color=method, shape=significant)) +
      geom_pointrange(position=position_dodge(width=0.5), size=0.8) +
      geom_hline(yintercept=0, linetype="dashed", color="red") +
      coord_flip() +
      scale_shape_manual(values=c("p<0.05"=16, "n.s."=1)) +
      labs(title="MR Analysis: BCP Genes vs Ischemic Stroke",
           subtitle="FinnGen R12 Ischemic Stroke (I9_STR)",
           x="Gene", y="Effect Size (Beta)") +
      theme_bw() +
      theme(legend.position="bottom")
    print(p)
  }
  
  dev.off()
  cat("  ✅ 森林图已保存\n")
  
} else {
  cat("  [警告] 没有成功完成任何MR分析\n")
}

cat("\n"); cat(rep("=", 60), sep=""); cat("\n")
cat("第4步完成！\n")
cat("结果目录:", mr_dir, "\n")
cat("\n关键发现:\n")
if(exists("sig_ivw") && length(sig_ivw) > 0) {
  cat(sprintf("  - IVW显著基因: %s\n", paste(sig_ivw, collapse=", ")))
}
if(exists("sig_wm") && length(sig_wm) > 0) {
  cat(sprintf("  - Weighted median显著基因: %s\n", paste(sig_wm, collapse=", ")))
}
cat("  - 在人群水平验证了BCP基因与缺血性脑卒中的关系\n")
cat(rep("=", 60), sep=""); cat("\n")
