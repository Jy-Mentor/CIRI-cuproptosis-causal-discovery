# 最简单的GSVA测试脚本
library(GSVA)
library(GSEABase)

# 创建测试数据
test_expr <- matrix(rnorm(1000), nrow=100, ncol=10, dimnames=list(paste0("gene", 1:100), paste0("sample", 1:10)))

# 创建基因集
test_genesets <- list(
  set1 = c("gene1", "gene2", "gene3"),
  set2 = c("gene4", "gene5", "gene6")
)

# 尝试运行GSVA
cat("开始运行GSVA...\n")
tryCatch({
  result <- gsva(test_expr, test_genesets)
  cat("成功! 结果维度:", dim(result), "\n")
  cat("前几个结果:\n")
  print(head(result))
}, error = function(e) {
  cat("错误:", e$message, "\n")
})

cat("测试完成\n")