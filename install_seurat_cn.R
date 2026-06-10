#!/usr/bin/env Rscript

# 使用国内CRAN镜像安装Seurat包
cat("使用国内CRAN镜像安装Seurat包...\n")

# 设置国内CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 安装依赖包
cat("安装依赖包...\n")
install.packages(c("dplyr", "ggplot2", "igraph", "data.table", "tidyr"), dependencies = TRUE)

# 安装Seurat
cat("安装Seurat...\n")
tryCatch({
  install.packages("Seurat", dependencies = TRUE)
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
}, error = function(e) {
  cat("安装过程出错:", e$message, "\n")
})

# 安装其他必要的包
cat("安装其他必要的包...\n")
install.packages(c("pcalg", "igraph"))

# 安装Bioconductor包
cat("安装Bioconductor包...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("Rgraphviz", dependencies = TRUE)

cat("所有包安装完成！\n")
