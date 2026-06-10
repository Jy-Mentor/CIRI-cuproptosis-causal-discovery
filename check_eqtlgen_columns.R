#!/usr/bin/env Rscript
# 检查 eQTLGen 数据的列名

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
})

cat("检查 eQTLGen 数据的列名\n")
cat("=========================================\n\n")

# 加载 eQTLGen 数据
eqtlgen_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-8_kb_10000_r2_0.001.xlsx"
eqtlgen_data <- read_excel(eqtlgen_file)

cat("列名:\n")
print(names(eqtlgen_data))
cat("\n")

cat("前 5 行数据:\n")
print(head(eqtlgen_data, 5))
cat("\n")

cat("PRKCQ 基因的数据:\n")
prkcq_snps <- eqtlgen_data[eqtlgen_data$gene == "PRKCQ" | eqtlgen_data$gene == "ENSG00000065675", ]
print(prkcq_snps)
cat("\n")

cat("SNP 列的类:\n")
print(class(prkcq_snps$SNP))
print(str(prkcq_snps$SNP))
