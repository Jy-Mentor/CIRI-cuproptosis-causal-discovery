# 构建Drug-Disease-DEGs交集网络
setwd("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写")

# 安装必要的包
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(VennDiagram)
library(ggplot2)

# 读取药物相关基因
drug_genes <- readLines("石竹烯 人.txt")
drug_genes <- drug_genes[drug_genes != ""]

# 读取疾病相关基因
disease_genes <- readLines("脑缺血 人.txt")
disease_genes <- disease_genes[disease_genes != ""]

# 读取映射后的DEGs
degs_mapped <- read.table("DEGs_mapped.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
deg_genes <- unique(degs_mapped$HUMAN_ORTHOLOG_SYMBOL)
deg_genes <- deg_genes[!is.na(deg_genes)]

# 计算交集
drug_disease_intersect <- intersect(drug_genes, disease_genes)
drug_deg_intersect <- intersect(drug_genes, deg_genes)
disease_deg_intersect <- intersect(disease_genes, deg_genes)
triple_intersect <- intersect(drug_disease_intersect, deg_genes)

# 打印交集信息
cat("药物相关基因数量:", length(drug_genes), "\n")
cat("疾病相关基因数量:", length(disease_genes), "\n")
cat("DEGs数量:", length(deg_genes), "\n")
cat("药物-疾病交集基因数量:", length(drug_disease_intersect), "\n")
cat("药物-DEGs交集基因数量:", length(drug_deg_intersect), "\n")
cat("疾病-DEGs交集基因数量:", length(disease_deg_intersect), "\n")
cat("药物-疾病-DEGs交集基因数量:", length(triple_intersect), "\n")

# 保存交集基因
write.table(data.frame(Gene = triple_intersect), "intersection_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# 绘制Venn图
venn.plot <- venn.diagram(
  x = list(
    "Drug" = drug_genes,
    "Disease" = disease_genes,
    "DEGs" = deg_genes
  ),
  filename = "Drug_Disease_DEGs_Venn.png",
  col = "black",
  fill = c("#FF9999", "#66B2FF", "#99FF99"),
  alpha = 0.5,
  label.col = "black",
  cex = 1.5,
  fontfamily = "sans",
  main = "Drug-Disease-DEGs 交集网络",
  main.cex = 2,
  category.names = c("药物相关基因", "疾病相关基因", "差异表达基因")
)

# 构建交集基因与表达数据的关联
intersection_expr <- degs_mapped[degs_mapped$HUMAN_ORTHOLOG_SYMBOL %in% triple_intersect, ]
write.table(intersection_expr, "intersection_genes_expression.tsv", sep = "\t", row.names = FALSE)

cat("交集网络构建完成！\n")
cat("交集基因已保存到 intersection_genes.tsv\n")
cat("Venn图已保存到 Drug_Disease_DEGs_Venn.png\n")
cat("交集基因表达数据已保存到 intersection_genes_expression.tsv\n")
