# 测试GSVA包基本功能
library(GSVA)

# 创建测试数据
test_expr <- matrix(rnorm(1000), nrow=100, ncol=10, dimnames=list(paste0("gene", 1:100), paste0("sample", 1:10)))
test_genesets <- list(
  set1 = c("gene1", "gene2", "gene3"),
  set2 = c("gene4", "gene5", "gene6")
)

# 尝试不同的GSVA调用方式
cat("尝试调用方式1: 基本位置参数\n")
try({
  result1 <- gsva(test_expr, test_genesets, method = "ssgsea", verbose = FALSE)
  cat("成功! 结果维度:", dim(result1), "\n")
}, silent = FALSE)

cat("\n尝试调用方式2: 只使用两个参数\n")
try({
  result2 <- gsva(test_expr, test_genesets)
  cat("成功! 结果维度:", dim(result2), "\n")
}, silent = FALSE)

cat("\n测试完成\n")