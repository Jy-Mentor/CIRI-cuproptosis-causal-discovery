#!/usr/bin/env Rscript
# ================================================================================
# 138 基因 MR 分析（使用本地 MEGASTROKE 数据）
# 参考 GitHub 权威项目最佳实践
# - 使用 GWAS Catalog 的 MEGASTROKE 数据（GCST006906）
# - 双源 eQTL 数据（eQTLGen + GTEx）
# - 无需 API 调用，直接使用本地文件
# ================================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

cat("======================================================================\n")
cat("138 基因 MR 分析（使用本地 MEGASTROKE 数据）\n")
cat("======================================================================\n\n")

# 配置
EXPOSURE_DIR