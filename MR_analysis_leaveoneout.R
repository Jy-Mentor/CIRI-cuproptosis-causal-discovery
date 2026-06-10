#!/usr/bin/env Rscript
# 留一法分析脚本

library(data.table)
library(TwoSampleMR)

output_dir <- "D:/EQTL/mr_results_p1e-05_fixed"

cat("========== 留一法（Leave-one-out）分析 ==========\n\n")

# 获取所有基因的结果文件
genes <- c("NFKB1", "FDX1", "STAT3", "HIF1A", "HMOX1", "GPX4", "HSPA5", "AGER")

for (gene in genes) {
  cat("==========", gene, "==========\n")
  
  # 读取harmonised数据
  harm_file <- file.path(output_dir, paste0(gene, "_harmonised_data.csv"))
  if (!file.exists(harm_file)) {
    cat("  文件不存在，跳过\n\n")
    next
  }
  
  dat <- fread(harm_file)
  cat("  SNP数:", nrow(dat), "\n")
  
  if (nrow(dat) < 2) {
    cat("  SNP不足2个，无法做留一法\n\n")
    next
  }
  
  # 留一法分析
  loo <- tryCatch({
    mr_leaveoneout(dat)
  }, error = function(e) {
    cat("  留一法失败:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(loo) || nrow(loo) == 0) {
    cat("  留一法无结果\n\n")
    next
  }
  
  # 保存留一法结果
  write.csv(loo, file.path(output_dir, paste0(gene, "_loo.csv")), row.names = FALSE)
  
  # 检查是否有某个SNP去掉后结果反转
  ivw_loo <- loo[loo$method == "Inverse variance weighted", ]
  if (nrow(ivw_loo) > 1) {
    # 获取完整分析的beta（SNP列为"All"的行）
    all_beta <- ivw_loo$b[ivw_loo$SNP == "All"]
    if (length(all_beta) > 0) {
      # 检查单个SNP去掉后的结果
      single_loo <- ivw_loo[ivw_loo$SNP != "All", ]
      sign_changes <- sum(sign(single_loo$b) != sign(all_beta), na.rm = TRUE)
      
      if (sign_changes > 0) {
        cat("  ⚠️ 警告: 存在", sign_changes, "个离群SNP主导结果!\n")
        # 显示具体是哪个SNP
        outlier_snps <- single_loo$SNP[sign(single_loo$b) != sign(all_beta)]
        cat("  离群SNP:", paste(outlier_snps, collapse = ", "), "\n")
      } else {
        cat("  ✓ 结果稳定，无明显离群SNP\n")
      }
      
      # 显示最大影响
      beta_diff <- abs(single_loo$b - all_beta)
      max_diff_idx <- which.max(beta_diff)
      cat("  最大影响SNP:", single_loo$SNP[max_diff_idx], 
          "(Beta变化:", round(beta_diff[max_diff_idx], 4), ")\n")
    }
  }
  
  cat("\n")
}

cat("========== 分析完成 ==========\n")
cat("留一法结果已保存到各基因的_loo.csv文件\n")
