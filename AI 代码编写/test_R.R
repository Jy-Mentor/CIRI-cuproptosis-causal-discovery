# 测试R环境
print("R环境测试成功！")
print(paste("R版本:", R.version$version.string))

# 测试读取文件
try {
  genes <- read.table("intersection_genes.tsv", header = TRUE, sep = "\t")
  print(paste("成功读取交集基因文件，共", nrow(genes), "个基因"))
} catch (e) {
  print("读取文件失败:")
  print(e)
}

print("测试完成！")