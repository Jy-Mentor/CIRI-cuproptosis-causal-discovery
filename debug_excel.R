#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(openxlsx))
scrna <- read.xlsx("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery/results/L1_phenotype_anchoring/L1_scRNA_GSE174574_Summary.xlsx", sheet = "Cuproptosis_Genes")
cat("列名:\n")
for (i in seq_len(ncol(scrna))) {
  cat("  [", i, "] ", colnames(scrna)[i], "\n", sep="")
}
cat("\n前5行:\n")
print(head(scrna, 5))
