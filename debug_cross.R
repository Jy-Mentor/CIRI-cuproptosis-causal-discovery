#!/usr/bin/env Rscript
# 调试：检查scRNA-seq数据

suppressPackageStartupMessages(library(openxlsx))

results_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery/results/L1_phenotype_anchoring"

scrna_excel <- file.path(results_dir, "L1_scRNA_GSE174574_Summary.xlsx")
scrna_raw <- read.xlsx(scrna_excel, sheet = "Cuproptosis_Genes", startRow = 3, colNames = TRUE)

cat("scRNA-seq Cuproptosis_Genes 列名:\n")
print(colnames(scrna_raw))

cat("\n前10行:\n")
print(head(scrna_raw, 10))

cat("\nlog2FC统计:\n")
print(summary(as.numeric(scrna_raw[["log2FC"]])))

cat("\nNA数量:", sum(is.na(as.numeric(scrna_raw[["log2FC"]]))), "\n")
