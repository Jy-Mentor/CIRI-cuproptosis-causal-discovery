# 检查GSVA包的版本和方法签名
library(GSVA)

# 检查GSVA版本
cat("GSVA版本:", packageVersion("GSVA"), "\n")

# 查看gsva函数的帮助
cat("\nGSVA函数帮助:\n")
try({
  help(gsva)
}, silent = TRUE)

# 查看gsvaParam函数的帮助
cat("\nGSVAParam函数帮助:\n")
try({
  help(gsvaParam)
}, silent = TRUE)

# 查看可用的方法签名
cat("\nGSVA可用方法签名:\n")
try({
  showMethods(gsva)
}, silent = TRUE)

cat("\n检查完成\n")