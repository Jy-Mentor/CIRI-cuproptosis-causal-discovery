# scTenifoldKnk 安装脚本 - 用户目录版
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

# 设置用户库路径
user_lib <- Sys.getenv("R_LIBS_USER")
if(user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "win-library", "4.5")
}
if(!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))
cat("使用库路径:", user_lib, "\n")

cat("=== scTenifoldKnk 安装脚本 ===\n\n")

# 1. 安装并加载关键依赖
cat("1. 安装依赖包...\n")
deps <- c("remotes", "Rcpp", "RcppArmadillo", "Matrix", "RSpectra", 
          "irlba", "rsvd", "doParallel", "foreach", "pbapply",
          "Seurat", "SeuratObject")

for(pkg in deps) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  安装 %s...\n", pkg))
    install.packages(pkg, dependencies = TRUE, lib = user_lib)
  } else {
    cat(sprintf("  [OK] %s 已安装\n", pkg))
  }
}

# 2. 尝试安装scTenifoldKnk
cat("\n2. 安装 scTenifoldKnk...\n")

# 先尝试安装scTenifoldNet (依赖包)
cat("  先安装 scTenifoldNet...\n")
tryCatch({
  remotes::install_github('dosorio/scTenifoldNet', upgrade = "never", lib = user_lib)
  cat("  [OK] scTenifoldNet 安装成功\n")
}, error = function(e) {
  cat("  [WARN] scTenifoldNet:", conditionMessage(e), "\n")
})

# 安装scTenifoldKnk
cat("  安装 scTenifoldKnk...\n")
tryCatch({
  remotes::install_github('dosorio/scTenifoldKnk', 
                          dependencies = TRUE,
                          upgrade = "never",
                          build_vignettes = FALSE,
                          force = TRUE,
                          lib = user_lib)
  cat("  [OK] scTenifoldKnk 安装成功!\n")
}, error = function(e) {
  cat("  [ERROR]", conditionMessage(e), "\n")
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
}
