# 安装scTenifoldKnk (依赖已满足)
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))

user_lib <- Sys.getenv("R_LIBS_USER")
if(user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "win-library", "4.5")
}
.libPaths(c(user_lib, .libPaths()))

cat("=== 安装 scTenifoldKnk ===\n\n")

# 加载scTenifoldNet
cat("1. 加载 scTenifoldNet...\n")
library(scTenifoldNet)
cat("✅ scTenifoldNet 已加载\n")

# 安装scTenifoldKnk
cat("\n2. 安装 scTenifoldKnk...\n")
pkg_path <- "C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/scripts/scTenifoldKnk/scTenifoldKnk-master"

if(dir.exists(pkg_path)) {
  cat("包路径:", pkg_path, "\n")
  tryCatch({
    install.packages(pkg_path, repos = NULL, type = "source", lib = user_lib)
    cat("✅ scTenifoldKnk 安装成功!\n")
  }, error = function(e) {
    cat("❌ 安装失败:\n")
    cat(conditionMessage(e), "\n")
  })
} else {
  cat("❌ 包路径不存在:\n", pkg_path, "\n")
}

# 验证
cat("\n3. 验证 scTenifoldKnk...\n")
if(require(scTenifoldKnk, quietly = TRUE)) {
  cat("✅✅✅ scTenifoldKnk 安装成功! ✅✅✅\n")
  cat("\n可用函数:\n")
  print(ls("package:scTenifoldKnk"))
} else {
  cat("❌ scTenifoldKnk 未能加载\n")
}

cat("\n=== 完成 ===\n")
