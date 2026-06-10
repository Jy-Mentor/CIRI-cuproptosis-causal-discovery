#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(openxlsx))
scrna_raw <- read.xlsx("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery/results/L1_phenotype_anchoring/L1_scRNA_GSE174574_Summary.xlsx", sheet = "Cuproptosis_Genes", startRow = 2, colNames = FALSE)
cat("X8列(上调/下调)前10行:\n")
for (i in 1:10) {
  val <- scrna_raw[[8]][i]
  cat("  [", i, "] '", val, "' (nchar=", nchar(as.character(val)), ")\n", sep="")
}
