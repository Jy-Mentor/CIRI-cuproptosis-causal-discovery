# 交集基因差异分析与功能富集
setwd("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写")

# 安装必要的包
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("ComplexHeatmap")
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("circlize", quietly = TRUE)) {
  BiocManager::install("circlize")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(ComplexHeatmap)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(circlize)
library(ggplot2)
library(tidyverse)

# 读取交集基因
intersection_genes <- read.table("intersection_genes.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
intersection_genes <- intersection_genes$Gene

# 读取映射后的DEGs
degs_mapped <- read.table("DEGs_mapped.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 匹配交集基因的表达数据
intersection_expr <- degs_mapped[degs_mapped$HUMAN_ORTHOLOG_SYMBOL %in% intersection_genes, ]

# 构建表达矩阵
# 由于我们只有logFC和P值，我们将使用logFC作为表达值
# 为了构建热图，我们需要整理数据
heatmap_data <- intersection_expr %>%
  select(HUMAN_ORTHOLOG_SYMBOL, logFC, adj.P.Val)

# 计算基因的平均logFC（如果有多个探针）
heatmap_data <- heatmap_data %>%
  group_by(HUMAN_ORTHOLOG_SYMBOL) %>%
  summarise(mean_logFC = mean(logFC), mean_pval = mean(adj.P.Val))

# 构建矩阵
heatmap_matrix <- matrix(heatmap_data$mean_logFC, nrow = nrow(heatmap_data), dimnames = list(heatmap_data$HUMAN_ORTHOLOG_SYMBOL, "logFC"))

# 绘制热图
# 先检查是否能创建文件
tryCatch({
  pdf("intersection_genes_heatmap.pdf", width = 10, height = 15)
  draw(Heatmap(
    heatmap_matrix,
    name = "logFC",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 12),
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  ), main_heatmap = 1, row_title = "Intersection Genes", column_title = "Expression Value")
  dev.off()
  cat("热图绘制成功\n")
}, error = function(e) {
  cat("热图绘制失败，可能是文件权限问题: ", e$message, "\n")
})

# 基因ID转换
gene_list <- bitr(intersection_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO富集分析
GO_enrichment <- enrichGO(
  gene = gene_list$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# 保存GO富集结果
write.table(as.data.frame(GO_enrichment), "intersection_GO_enrichment.tsv", sep = "\t", row.names = FALSE)

# KEGG富集分析
try({
  KEGG_enrichment <- enrichKEGG(
    gene = gene_list$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  # 保存KEGG富集结果
  write.table(as.data.frame(KEGG_enrichment), "intersection_KEGG_enrichment.tsv", sep = "\t", row.names = FALSE)
}, silent = TRUE)

cat("KEGG富集分析完成（如果遇到网络超时，结果可能不完整）\n")

# 绘制GO富集柱状图
tryCatch({
  pdf("intersection_GO_barplot.pdf", width = 12, height = 10)
  goplot <- barplot(GO_enrichment, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free_y") +
    ggtitle("Intersection Genes GO Enrichment Analysis")
  print(goplot)
  dev.off()
  cat("GO富集柱状图绘制成功\n")
}, error = function(e) {
  cat("GO富集柱状图绘制失败，可能是文件权限问题: ", e$message, "\n")
})

# 绘制KEGG富集气泡图
tryCatch({
  if (exists("KEGG_enrichment")) {
    pdf("intersection_KEGG_bubble.pdf", width = 12, height = 10)
    keggplot <- dotplot(KEGG_enrichment, showCategory = 20) +
      ggtitle("Intersection Genes KEGG Enrichment Analysis")
    print(keggplot)
    dev.off()
    cat("KEGG富集气泡图绘制成功\n")
  }
}, error = function(e) {
  cat("KEGG富集气泡图绘制失败，可能是文件权限问题: ", e$message, "\n")
})

# 绘制GO富集圈图
tryCatch({
  if (nrow(GO_enrichment) > 0) {
    pdf("intersection_GO_cnetplot.pdf", width = 15, height = 12)
    cnetplot(GO_enrichment, showCategory = 10, foldChange = gene_list$ENTREZID)
    dev.off()
    cat("GO富集圈图绘制成功\n")
  }
}, error = function(e) {
  cat("GO富集圈图绘制失败，可能是文件权限问题: ", e$message, "\n")
})

# 绘制GO富集弦图
tryCatch({
  if (nrow(GO_enrichment) > 0) {
    pdf("intersection_GO_upsetplot.pdf", width = 12, height = 10)
    upsetplot(GO_enrichment, n = 10)
    dev.off()
    cat("GO富集弦图绘制成功\n")
  }
}, error = function(e) {
  cat("GO富集弦图绘制失败，可能是文件权限问题: ", e$message, "\n")
})

cat("交集基因差异分析与功能富集完成！\n")
cat("热图已保存到 intersection_genes_heatmap.pdf\n")
cat("GO富集结果已保存到 intersection_GO_enrichment.tsv\n")
cat("KEGG富集结果已保存到 intersection_KEGG_enrichment.tsv\n")
cat("GO富集柱状图已保存到 intersection_GO_barplot.pdf\n")
cat("KEGG富集气泡图已保存到 intersection_KEGG_bubble.pdf\n")
