#!/usr/bin/env Rscript

# 安装必要的包
cat("安装必要的R包...\n")

# 安装CRAN包
cran_packages <- c("Seurat", "data.table", "dplyr", "tidyr", "pcalg", "igraph", "ggplot2")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    cat(sprintf("安装 %s...\n", pkg))
    install.packages(pkg, dependencies = TRUE, repos = "https://cran.r-project.org/")
  } else {
    cat(sprintf("%s 已安装\n", pkg))
  }
}

# 安装Bioconductor包
cat("安装Bioconductor包...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.r-project.org/")
}

if (!require("Rgraphviz", character.only = TRUE)) {
  cat("安装 Rgraphviz...\n")
  BiocManager::install("Rgraphviz")
} else {
  cat("Rgraphviz 已安装\n")
}

cat("所有包安装完成！\n")
