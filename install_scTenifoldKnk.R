# scTenifoldKnk 安装脚本 - Windows优化版
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

# 1. 检查并配置Rtools
cat("=== 配置Rtools ===\n")
if(!require(pkgbuild, quietly = TRUE)) install.packages('pkgbuild')
library(pkgbuild)

# 设置Rtools路径
rtools_path <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/AI 代码编写/rtools45"
if(dir.exists(rtools_path)) {
  Sys.setenv(PATH = paste(rtools_path, "usr/bin", Sys.getenv("PATH"), sep = ";"))
  Sys.setenv(BINPREF = file.path(rtools_path, "mingw64/bin/"))
  cat("Rtools路径已设置:", rtools_path, "\n")
}

# 检查编译器
cat("检查编译器...\n")
tryCatch({
  find_rtools()
  cat("[OK] Rtools已正确配置\n")
}, error = function(e) {
  cat("[WARN] Rtools配置可能有问题:", conditionMessage(e), "\n")
})

# 2. 安装关键依赖
cat("\n=== 安装编译依赖 ===\n")
dependencies <- c("Rcpp", "RcppArmadillo", "Matrix", "RSpectra", 
                  "irlba", "rsvd", "doParallel", "foreach",
                  "Seurat", "SeuratObject", "pbapply")

for(pkg in dependencies) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("安装 %s...\n", pkg))
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(sprintf("[OK] %s 已安装\n", pkg))
  }
}

# 3. 安装remotes包用于从GitHub安装
cat("\n=== 安装remotes包 ===\n")
if(!require(remotes, quietly = TRUE)) {
  install.packages('remotes')
}
library(remotes)

# 4. 安装scTenifoldNet (先安装，作为基础)
cat("\n=== 安装scTenifoldNet ===\n")
if(!require(scTenifoldNet, quietly = TRUE)) {
  install_github('dosorio/scTenifoldNet')
}

# 5. 尝试安装scTenifoldKnk
cat("\n=== 安装scTenifoldKnk ===\n")
cat("尝试从GitHub安装...\n")

tryCatch({
  # 先尝试直接安装
  remotes::install_github('dosorio/scTenifoldKnk', 
                          dependencies = TRUE,
                          upgrade = "never",
                          build_vignettes = FALSE,
                          verbose = TRUE)
  cat("[OK] scTenifoldKnk 安装成功!\n")
}, error = function(e) {
  cat("[ERROR] 安装失败:", conditionMessage(e), "\n")
  cat("\n尝试备用方法...\n")
  
  # 备用方法：手动下载并安装
  tryCatch({
    cat("下载源码包...\n")
    url <- "https://github.com/dosorio/scTenifoldKnk/archive/refs/heads/main.zip"
    destfile <- tempfile(fileext = ".zip")
    download.file(url, destfile, mode = "wb")
    
    cat("解压并安装...\n")
    unzip_dir <- tempfile()
    unzip(destfile, exdir = unzip_dir)
    
    pkg_dir <- file.path(unzip_dir, "scTenifoldKnk-main")
    install.packages(pkg_dir, repos = NULL, type = "source")
    
    cat("[OK] 备用方法安装成功!\n")
  }, error = function(e2) {
    cat("[ERROR] 备用方法也失败:", conditionMessage(e2), "\n")
  })
})

# 6. 验证安装
cat("\n=== 验证安装 ===\n")
if(require(scTenifoldKnk, quietly = TRUE)) {
  cat("[SUCCESS] scTenifoldKnk 安装成功!\n")
  cat("版本:", packageVersion("scTenifoldKnk"), "\n")
  cat("可用函数:\n")
  print(ls("package:scTenifoldKnk"))
} else {
  cat("[FAIL] 安装失败，请检查错误信息\n")
}
