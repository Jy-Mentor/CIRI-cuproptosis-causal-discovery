#!/usr/bin/env Rscript
# ================================================================================
# 138 基因 MR 分析（使用本地 MEGASTROKE 文件）
# 数据源：MEGASTROKE GWAS summary statistics (GCST006906)
# ================================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

cat("======================================================================\n")
cat("138 基因 MR 分析（使用本地 MEGASTROKE 文件）\n")
cat("======================================================================\n\n")

# 配置
EXPOSURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_matched/matched_data"
OUTPUT_DIR