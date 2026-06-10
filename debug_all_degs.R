#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(openxlsx))

results_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery/results/L1_phenotype_anchoring"

# 读取All_DEGs
all_degs <- read.xlsx(file.path(results_dir, "L1_scRNA_GSE174574_Summary.xlsx"), 
                      sheet = "All_DEGs", startRow = 3, colNames = TRUE)

cat("All_DEGs 总行数:", nrow(all_degs), "\n")
cat("列名:", paste(colnames(all_degs), collapse=", "), "\n")

# 铜死亡基因列表
cuproptosis_genes <- c("FDX1", "LIAS", "DLD", "DLAT", "DLST", "PDHA1", "PDHB", "GLS", 
  "GCSH", "LIPT1", "LIPT2", "CDKN2A", "NFE2L2", "NLRP3",
  "SLC31A1", "SLC31A2", "SLC11A2", "STEAP3", "ATP7A", "ATP7B",
  "ATOX1", "CCS", "COX17", "COX11", "SCO1", "SCO2",
  "MT1A", "MT2A", "ALB", "CP", "SOD1", "SOD3",
  "COMMD1", "MTF1")

# 查找铜死亡基因
found <- all_degs[toupper(all_degs[[1]]) %in% cuproptosis_genes, ]
cat("\n在All_DEGs中找到的铜死亡基因:", nrow(found), "\n")
print(found)
