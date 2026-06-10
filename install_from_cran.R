# 从CRAN安装scTenifoldNet和scTenifoldKnk
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))

user_lib <- Sys.getenv("R_LIBS_USER")
if(user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "win-library", "4.5")
}
if(!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

cat("=== 从CRAN安装 ===\n\n")

# 1. 安装scTenifoldNet
cat("1. 安装 scTenifoldNet...\n")
tryCatch({
  install.packages("scTenifoldNet", lib = user_lib, dependencies = TRUE)
  cat("✅ scTenifoldNet 安装成功!\n")
}, error = function(e) {
  cat("❌ 安装失败:", conditionMessage(e), "\n")
})

# 2. 验证scTenifoldNet
cat("\n2. 验证 scTenifoldNet...\n")
if(require(scTenifoldNet, quietly = TRUE)) {
  cat("✅ scTenifoldNet 加载成功!\n")
  cat("版本:", packageVersion("scTenifoldNet"), "\n")
} else {
  cat("❌ scTenifoldNet 未能加载\n")
}

# 3. 安装本地scTenifoldKnk
cat("\n3. 安装 scTenifoldKnk...\n")
pkg_path <- "C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/scripts/scTenifoldKnk/scTenifoldKnk-master"

if(dir.exists(pkg_path) && require(scTenifoldNet, quietly = TRUE)) {
  tryCatch({
    install.packages(pkg_path, repos = NULL, type = "source", lib = user_lib)
    cat("✅ scTenifoldKnk 安装成功!\n")
  }, error = function(e) {
    cat("❌ scTenifoldKnk 安装失败:", conditionMessage(e), "\n")
  })
} else {
  cat("⚠️ 跳过scTenifoldKnk安装 (依赖不满足或路径不存在)\n")
}

# 4. 验证scTenifoldKnk
cat("\n4. 验证 scTenifoldKnk...\n")
if(require(scTenifoldKnk, quietly = TRUE)) {
  cat("✅ scTenifoldKnk 加载成功!\n")
  cat("版本:", packageVersion("scTenifoldKnk"), "\n")
  cat("\n可用函数:\n")
  print(ls("package:scTenifoldKnk"))
} else {
  cat("❌ scTenifoldKnk 未能加载\n")
}

cat("\n=== 完成 ===\n")
