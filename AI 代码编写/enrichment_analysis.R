# 基因集富集分析与可视化

# 加载必要的包
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# 读取映射后的DEGs数据
degs_mapped <- read.table("DEGs_mapped.tsv", header = TRUE, sep = "\t")

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
cat("基因集富集分析完成！\n")
cat("GO富集分析结果数量:", nrow(GO_enrichment@result), "\n")
cat("KEGG富集分析结果数量:", nrow(KEGG_enrichment@result), "\n")
