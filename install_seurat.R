#!/usr/bin/env Rscript

# 手动安装Seurat包
cat("安装Seurat包...\n")

# 尝试从CRAN安装
install.packages("Seurat", dependencies = TRUE, repos = "https://cran.r-project.org/")

# 检查是否安装成功
if (require("Seurat", character.only = TRUE)) {
  cat("Seurat安装成功！\n")
} else {
  cat("Seurat安装失败，尝试从GitHub安装...\n")
  # 尝试从GitHub安装
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("satijalab/seurat")
  
  if (require("Seurat", character.only = TRUE)) {
    cat("Seurat安装成功！\n")
  } else {
    cat("Seurat安装失败，请检查网络连接或手动安装。\n")
  }
}
