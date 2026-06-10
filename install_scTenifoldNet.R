# 安装scTenifoldNet
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))

user_lib <- Sys.getenv("R_LIBS_USER")
if(user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "win-library", "4.5")
}
if(!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

cat("=== 安装scTenifoldNet ===\n\n")

# 安装依赖
cat("1. 安装依赖...\n")
deps <- c("Rcpp", "RcppArmadillo", "Matrix", "RSpectra", "irlba", "rsvd", 
          "doParallel", "foreach", "pbapply")
for(pkg in deps) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, lib = user_lib)
  }
}

# 从GitHub安装
cat("\n2. 从GitHub安装scTenifoldNet...\n")
if(!require(remotes, quietly = TRUE)) install.packages('remotes', lib = user_lib)

tryCatch({
  remotes::install_github('dosorio/scTenifoldNet', 
                          upgrade = "never", 
                          lib = user_lib,
                          dependencies = TRUE)
  cat("\n✅ scTenifoldNet 安装成功!\n")
}, error = function(e) {
  cat("\n❌ 安装失败:", conditionMessage(e), "\n")
})

# 验证
if(require(scTenifoldNet, quietly = TRUE)) {
  cat("✅ 验证成功，版本:", packageVersion("scTenifoldNet"), "\n")
}
