#!/usr/bin/env Rscript
# MR 分析结果导出为 Excel
library(readxl)
library(writexl)
library(dplyr)

cat("\n=== MR 分析结果导出为 Excel ===\n\n")

# 设置路径
results_dir <- "D:/下载/MR_batch_results"
output_file <- "D:/下载/MR_batch_results/MR_results_summary.xlsx"

# 读取主结果文件
cat("读取 MR 分析结果...\n")
mr_file <- file.path(results_dir, "20260508/MR_results_main.csv")

if (!file.exists(mr_file)) {
  cat("错误：找不到主结果文件:", mr_file, "\n")
  quit(status = 1)
}

mr_results <- read.csv(mr_file, stringsAsFactors = FALSE)

cat("成功读取", nrow(mr_results), "个基因的结果\n\n")

# 创建 Excel 工作簿
cat("创建 Excel 汇总报告...\n")

# 按 P 值排序
mr_results_sorted <- mr_results %>%
  arrange(pval)

# 添加显著性标记
mr_results_sorted$significance <- ifelse(
  mr_results_sorted$pval < 0.05,
  ifelse(mr_results_sorted$fdr_qval < 0.05, "FDR 显著", "P 显著"),
  "不显著"
)

# 从 OR_95CI 中提取 OR 值（格式："1.048 (0.999-1.100)"）
extract_or <- function(or_ci_str) {
  or_val <- as.numeric(gsub(" \\(.*\\)", "", or_ci_str))
  return(or_val)
}

mr_results_sorted$OR <- extract_or(mr_results_sorted$OR_95CI)

# 添加 OR 方向标记
mr_results_sorted$effect_direction <- ifelse(
  mr_results_sorted$OR > 1,
  "风险因素 (OR>1)",
  ifelse(mr_results_sorted$OR < 1, "保护因素 (OR<1)", "无效应 (OR=1)")
)

# 导出到 Excel
sheets_list <- list(
  "所有结果_按 P 值排序" = mr_results_sorted,
  "显著基因_P<0.05" = mr_results_sorted %>% filter(pval < 0.05),
  "FDR 显著基因" = mr_results_sorted %>% filter(fdr_sig == TRUE),
  "保护因素" = mr_results_sorted %>% filter(OR < 1 & pval < 0.05),
  "风险因素" = mr_results_sorted %>% filter(OR > 1 & pval < 0.05)
)

# 写入 Excel
write_xlsx(sheets_list, path = output_file)

cat("\n✓ Excel 报告已创建:", output_file, "\n\n")

# 打印摘要
cat("=== 结果摘要 ===\n")
cat("总基因数:", nrow(mr_results_sorted), "\n")
cat("P 值显著基因数:", sum(mr_results_sorted$pval < 0.05), "\n")
cat("FDR 校正后显著基因数:", sum(mr_results_sorted$fdr_sig), "\n")
cat("保护因素基因数:", sum(mr_results_sorted$OR < 1 & mr_results_sorted$pval < 0.05), "\n")
cat("风险因素基因数:", sum(mr_results_sorted$OR > 1 & mr_results_sorted$pval < 0.05), "\n\n")

# 显示 Top 10 显著基因
cat("=== Top 10 最显著基因 ===\n")
top10 <- head(mr_results_sorted %>% filter(pval < 0.05), 10)
if (nrow(top10) > 0) {
  print(top10[, c("gene", "pval", "OR_95CI", "effect_direction", "significance")])
} else {
  cat("无 P<0.05 的显著基因\n")
}

cat("\n✓ 导出完成！\n\n")
