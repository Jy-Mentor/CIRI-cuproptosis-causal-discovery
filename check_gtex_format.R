#!/usr/bin/env Rscript
# ================================================================================
# 检查 GTEx 数据格式
# ================================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(arrow)
  library(data.table)
})

cat("======================================================================\n")
cat("检查 GTEx 数据格式\n")
cat("======================================================================\n\n")

# 文件路径
GTEx_DIR <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL"

# 检查全血 eQTL
cat("1. 检查全血 eQTL 数据...\n")
blood_genes_file <- file.path(GTEx_DIR, "Whole_Blood.v11.eGenes.txt")
blood_eqtl_file <- file.path(GTEx_DIR, "Whole_Blood.v11.eQTLs.signif_pairs.parquet")

if (file.exists(blood_genes_file)) {
  cat("  读取 eGenes 文件...\n")
  blood_genes <- fread(blood_genes_file, sep = "\t", header = TRUE)
  cat(sprintf("  eGenes 文件维度：%d x %d\n", nrow(blood_genes), ncol(blood_genes)))
  cat("  列名:\n")
  print(names(blood_genes))
  cat("\n")
}

if (file.exists(blood_eqtl_file)) {
  cat("  读取 eQTLs parquet 文件...\n")
  blood_eqtl <- read_parquet(blood_eqtl_file)
  cat(sprintf("  eQTLs 文件维度：%d x %d\n", nrow(blood_eqtl), ncol(blood_eqtl)))
  cat("  列名:\n")
  print(names(blood_eqtl))
  cat("\n")
  cat("  前 3 行:\n")
  print(head(blood_eqtl, 3))
  cat("\n")
}

# 检查脑组织 eQTL
cat("2. 检查脑组织 eQTL 数据...\n")
brain_genes_file <- file.path(GTEx_DIR, "Brain_Cortex.v11.eGenes.txt")
brain_eqtl_file <- file.path(GTEx_DIR, "Brain_Cortex.v11.eQTLs.signif_pairs.parquet")

if (file.exists(brain_genes_file)) {
  cat("  读取 eGenes 文件...\n")
  brain_genes <- fread(brain_genes_file, sep = "\t", header = TRUE)
  cat(sprintf("  eGenes 文件维度：%d x %d\n", nrow(brain_genes), ncol(brain_genes)))
  cat("  列名:\n")
  print(names(brain_genes))
  cat("\n")
}

if (file.exists(brain_eqtl_file)) {
  cat("  读取 eQTLs parquet 文件...\n")
  brain_eqtl <- read_parquet(brain_eqtl_file)
  cat(sprintf("  eQTLs 文件维度：%d x %d\n", nrow(brain_eqtl), ncol(brain_eqtl)))
  cat("  列名:\n")
  print(names(brain_eqtl))
  cat("\n")
  cat("  前 3 行:\n")
  print(head(brain_eqtl, 3))
  cat("\n")
}

# 检查 CSV 格式
cat("3. 检查 CSV 格式数据...\n")
blood_csv_file <- file.path(GTEx_DIR, "blood_eqtl.csv")
brain_csv_file <- file.path(GTEx_DIR, "brain_eqtl.csv")

if (file.exists(blood_csv_file)) {
  cat("  读取全血 CSV...\n")
  blood_csv <- fread(blood_csv_file, header = TRUE)
  cat(sprintf("  维度：%d x %d\n", nrow(blood_csv), ncol(blood_csv)))
  cat("  列名:\n")
  print(names(blood_csv))
  cat("\n")
}

if (file.exists(brain_csv_file)) {
  cat("  读取脑组织 CSV...\n")
  brain_csv <- fread(brain_csv_file, header = TRUE)
  cat(sprintf("  维度：%d x %d\n", nrow(brain_csv), ncol(brain_csv)))
  cat("  列名:\n")
  print(names(brain_csv))
  cat("\n")
}

cat("\n完成！\n")
