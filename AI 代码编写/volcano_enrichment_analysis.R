# 差异火山图绘制及基因集富集分析可视化
setwd("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写")

# 安装必要的包
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot")
}
if (!requireNamespace("DOSE", quietly = TRUE)) {
  BiocManager::install("DOSE")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

# 加载包
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(readxl)

# 1. 读取差异表达基因数据
# 使用GSE61616.top.table (1).tsv文件
data <- read.table("GSE61616.top.table (1).tsv", header = TRUE, sep = "\t")

# 2. 绘制火山图
# 计算-log10(P值)
data$logP <- -log10(data$adj.P.Val)

# 设置显著差异基因的标记
data$significance <- ifelse(data$adj.P.Val < 0.05 & abs(data$logFC) >= 0.2, 
                           ifelse(data$logFC > 0, "Up-regulated", "Down-regulated"), 
                           "Not significant")

# 绘制火山图
pdf("volcano_plot.pdf", width = 12, height = 8)
ggplot(data, aes(x = logFC, y = logP, color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Down-regulated" = "blue", "Not significant" = "gray", "Up-regulated" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "top",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()

cat("火山图绘制完成，保存为 volcano_plot.pdf\n")

# 3. 基因集富集分析
# 提取显著差异基因
DEGs <- data[data$adj.P.Val < 0.05 & abs(data$logFC) >= 0.2, ]

# 提取基因符号
gene_list <- DEGs$Gene.symbol

# 移除空值
gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]

# 转换基因符号为ENTREZ ID
gene_entrez <- mapIds(org.Hs.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
gene_entrez <- gene_entrez[!is.na(gene_entrez)]

if (length(gene_entrez) > 0) {
  # GO富集分析
  tryCatch({
    go_enrich <- enrichGO(gene = gene_entrez, 
                         OrgDb = org.Hs.eg.db, 
                         keyType = "ENTREZID", 
                         ont = "ALL", 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.1)
    
    # 保存GO富集结果
    write.table(go_enrich@result, "go_enrichment_results.txt", sep = "\t", row.names = FALSE)
    
    # GO富集分析气泡图
    pdf("go_bubble_plot.pdf", width = 12, height = 10)
    dotplot(go_enrich, showCategory = 30) +
      labs(title = "GO Enrichment Analysis") +
      theme(plot.title = element_text(hjust = 0.5, size = 16))
    dev.off()
    
    # GO富集分析网络图
    if (nrow(go_enrich) > 0) {
      pdf("go_network_plot.pdf", width = 15, height = 12)
      tryCatch({
        cnetplot(go_enrich, showCategory = 10, foldChange = gene_entrez)
      }, error = function(e) {
        cat("GO网络绘制失败：", e$message, "\n")
      })
      dev.off()
    }
    
    cat("GO基因集富集分析完成，结果已保存\n")
    cat("GO富集分析结果：go_enrichment_results.txt\n")
    cat("GO气泡图：go_bubble_plot.pdf\n")
    cat("GO网络图：go_network_plot.pdf\n")
  }, error = function(e) {
    cat("GO富集分析失败：", e$message, "\n")
  })
  
  # KEGG富集分析（可能会因为网络问题失败）
  tryCatch({
    kegg_enrich <- enrichKEGG(gene = gene_entrez, 
                            organism = "hsa", 
                            pAdjustMethod = "BH", 
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.1)
    
    # 保存KEGG富集结果
    write.table(kegg_enrich@result, "kegg_enrichment_results.txt", sep = "\t", row.names = FALSE)
    
    # KEGG富集分析气泡图
    pdf("kegg_bubble_plot.pdf", width = 12, height = 10)
    dotplot(kegg_enrich, showCategory = 30) +
      labs(title = "KEGG Pathway Enrichment Analysis") +
      theme(plot.title = element_text(hjust = 0.5, size = 16))
    dev.off()
    
    # KEGG富集分析网络图
    if (nrow(kegg_enrich) > 0) {
      pdf("kegg_network_plot.pdf", width = 15, height = 12)
      tryCatch({
        cnetplot(kegg_enrich, showCategory = 10, foldChange = gene_entrez)
      }, error = function(e) {
        cat("KEGG网络绘制失败：", e$message, "\n")
      })
      dev.off()
    }
    
    cat("KEGG基因集富集分析完成，结果已保存\n")
    cat("KEGG富集分析结果：kegg_enrichment_results.txt\n")
    cat("KEGG气泡图：kegg_bubble_plot.pdf\n")
    cat("KEGG网络图：kegg_network_plot.pdf\n")
  }, error = function(e) {
    cat("KEGG富集分析失败（可能是网络连接问题）：", e$message, "\n")
  })
} else {
  cat("没有找到有效的基因ID，无法进行富集分析\n")
}

cat("分析完成！\n")
