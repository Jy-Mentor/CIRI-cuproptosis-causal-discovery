# 安装本地scTenifoldKnk
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))

# 设置用户库路径
user_lib <- Sys.getenv("R_LIBS_USER")
if(user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "win-library", "4.5")
}
if(!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

cat("=== 安装本地scTenifoldKnk ===\n")
cat("库路径:", user_lib, "\n\n")

# 安装依赖
cat("1. 安装依赖...\n")
deps <- c("Rcpp", "RcppArmadillo", "Matrix", "RSpectra", "irlba", "rsvd")
for(pkg in deps) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, lib = user_lib)
  }
}

# 安装本地包
cat("\n2. 安装本地scTenifoldKnk...\n")
pkg_path <- "C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/scripts/scTenifoldKnk/scTenifoldKnk-master"

if(dir.exists(pkg_path)) {
  cat("包路径:", pkg_path, "\n")
  install.packages(pkg_path, repos = NULL, type = "source", lib = user_lib)
  cat("\n✅ 安装完成!\n")
} else {
  cat("❌ 包路径不存在:\n", pkg_path, "\n")
}

# 验证
cat("\n3. 验证安装...\n")
if(require(scTenifoldKnk, quietly = TRUE)) {
  cat("✅ scTenifoldKnk 安装成功!\n")
  cat("版本:", packageVersion("scTenifoldKnk"), "\n")
} else {
  cat("❌ 安装失败\n")
}
