# scTenifoldKnk 安装脚本 - Windows版 (路径修正)
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

cat("=== scTenifoldKnk 安装脚本 ===\n\n")

# 1. 安装并加载关键依赖
cat("1. 安装依赖包...\n")
deps <- c("remotes", "Rcpp", "RcppArmadillo", "Matrix", "RSpectra", 
          "irlba", "rsvd", "doParallel", "foreach", "pbapply",
          "Seurat", "SeuratObject")

for(pkg in deps) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  安装 %s...\n", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
}

# 2. 尝试安装scTenifoldKnk
cat("\n2. 安装 scTenifoldKnk...\n")

# 先尝试安装scTenifoldNet (依赖包)
cat("  先安装 scTenifoldNet...\n")
tryCatch({
  remotes::install_github('dosorio/scTenifoldNet', upgrade = "never")
  cat("  [OK] scTenifoldNet 安装成功\n")
}, error = function(e) {
  cat("  [WARN] scTenifoldNet 安装可能有问题，继续尝试...\n")
})

# 安装scTenifoldKnk
cat("  安装 scTenifoldKnk...\n")
tryCatch({
  remotes::install_github('dosorio/scTenifoldKnk', 
                          dependencies = TRUE,
                          upgrade = "never",
                          build_vignettes = FALSE,
                          force = TRUE)
  cat("  [OK] scTenifoldKnk 安装成功!\n")
}, error = function(e) {
  cat("  [ERROR] 安装失败:", conditionMessage(e), "\n")
  cat("\n尝试备用方法...\n")
  
  # 备用：安装开发工具包
  if(!require(devtools, quietly = TRUE)) {
    install.packages('devtools')
  }
  
  tryCatch({
    devtools::install_github('dosorio/scTenifoldKnk', force = TRUE)
  }, error = function(e2) {
    cat("  [ERROR] 备用方法也失败:\n")
    cat(conditionMessage(e2), "\n")
  })
})

# 3. 验证安装
cat("\n3. 验证安装...\n")
if(require(scTenifoldKnk, quietly = TRUE)) {
  cat("[SUCCESS] scTenifoldKnk 安装成功!\n")
  cat("版本:", packageVersion("scTenifoldKnk"), "\n")
  cat("\n可用函数:\n")
  print(ls("package:scTenifoldKnk"))
} else {
  cat("[FAIL] 安装未成功\n")
  cat("可能原因:\n")
  cat("  1. Rtools未正确配置\n")
  cat("  2. 编译器权限问题\n")
  cat("  3. 网络连接问题\n")
}
