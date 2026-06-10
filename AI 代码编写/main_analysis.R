# 主分析脚本：系统化生物信息学分析流程

# 设置工作目录
setwd("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写")

# 1. 差异表达基因分析
cat("===== 1. 差异表达基因分析 =====\n")

# 读取数据
data <- read.table("GSE61616.top.table (1).tsv", header = TRUE, sep = "\t")

# 筛选差异表达基因（绝对对数FC≥0.2，校正后P值<0.05）
degs <- data[data$adj.P.Val < 0.05 & abs(data$logFC) >= 0.2, ]

# 保存DEGs
write.table(degs, "DEGs.tsv", sep = "\t", row.names = FALSE)

# 绘制火山图
library(ggplot2)
volcano_plot <- ggplot(data, aes(x = logFC, y = -log10(adj.P.Val)))
volcano_plot <- volcano_plot +
  geom_point(aes(color = ifelse(adj.P.Val < 0.05 & abs(logFC) >= 0.2, "显著差异", "无显著差异")), 
             alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("显著差异" = "red", "无显著差异" = "gray")) +
  labs(title = "差异表达基因火山图",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "表达差异") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "black")

# 保存火山图
ggsave("volcano_plot.png", volcano_plot, width = 10, height = 8, dpi = 300)

# 统计结果
cat("总基因数:", nrow(data), "\n")
cat("显著差异表达基因数:", nrow(degs), "\n")
cat("其中上调基因数:", sum(degs$logFC > 0), "\n")
cat("其中下调基因数:", sum(degs$logFC < 0), "\n")

# 2. 基因映射
cat("\n===== 2. 基因映射 =====\n")

# 读取基因映射库
library(stringr)
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
cat("原始DEGs数量:", nrow(degs), "\n")
cat("成功映射到人类基因的DEGs数量:", nrow(degs_mapped), "\n")

# 3. 基因集富集分析
cat("\n===== 3. 基因集富集分析 =====\n")

# 加载必要的包
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# 提取人类基因符号
human_genes <- unique(degs_mapped$HUMAN_ORTHOLOG_SYMBOL)

# 转换为ENTREZ ID
gene_list <- bitr(human_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO富集分析
GO_enrichment <- enrichGO(gene = gene_list$ENTREZID, 
                         OrgDb = org.Hs.eg.db, 
                         ont = "ALL", 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.2)

# KEGG富集分析
KEGG_enrichment <- enrichKEGG(gene = gene_list$ENTREZID, 
                             organism = "hsa", 
                             pAdjustMethod = "BH", 
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 0.2)

# 保存富集分析结果
write.table(GO_enrichment@result, "GO_enrichment.tsv", sep = "\t", row.names = FALSE)
write.table(KEGG_enrichment@result, "KEGG_enrichment.tsv", sep = "\t", row.names = FALSE)

# 可视化

# GO富集分析柱状图
go_bar <- barplot(GO_enrichment, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  ggtitle("GO富集分析")
ggsave("GO_barplot.png", go_bar, width = 12, height = 10, dpi = 300)

# KEGG富集分析气泡图
kegg_bubble <- dotplot(KEGG_enrichment, showCategory = 20)
ggsave("KEGG_bubble.png", kegg_bubble, width = 10, height = 8, dpi = 300)

# 统计结果
cat("GO富集分析结果数量:", nrow(GO_enrichment@result), "\n")
cat("KEGG富集分析结果数量:", nrow(KEGG_enrichment@result), "\n")

cat("\n===== 分析完成！ =====\n")
