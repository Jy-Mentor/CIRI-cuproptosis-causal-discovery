# 基因映射：将大鼠基因映射到人类基因

# 读取DEGs数据
degs <- read.table("DEGs.tsv", header = TRUE, sep = "\t")

# 读取基因映射库
mapping <- read.table("大鼠 小鼠 人类映射库.txt", header = TRUE, sep = "\t", comment.char = "#")

# 处理基因符号，提取第一个基因符号（如果有多个）
degs$Gene.symbol <- sapply(degs$Gene.symbol, function(x) {
  if (grepl("///", x)) {
    return(strsplit(x, "///")[[1]][1])
  } else {
    return(x)
  }
})

# 映射基因
degs_mapped <- merge(degs, mapping, by.x = "Gene.symbol", by.y = "RAT_GENE_SYMBOL", all.x = TRUE)

# 筛选有人类同源基因的DEGs
degs_mapped <- degs_mapped[!is.na(degs_mapped$HUMAN_ORTHOLOG_SYMBOL), ]

# 保存映射结果
write.table(degs_mapped, "DEGs_mapped.tsv", sep = "\t", row.names = FALSE)

# 统计结果
cat("基因映射完成！\n")
cat("原始DEGs数量:", nrow(degs), "\n")
cat("成功映射到人类基因的DEGs数量:", nrow(degs_mapped), "\n")
