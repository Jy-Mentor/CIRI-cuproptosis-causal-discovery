# GO和KEGG富集分析脚本
# 分析用户提供的基因列表

# 设置CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 安装必要的R包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ggplot2", "enrichplot", "ggpubr", "viridis"))

# 加载所需包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(ggpubr)
library(viridis)

# 定义用户提供的基因列表
gene_list <- c(
  "TP53", "IL1B", "IL6", "TNF", "STAT3", "BCL2", "NFKB1", "PTGS2", "TLR4", "SRC",
  "STAT1", "RELA", "ICAM1", "CCL2", "CCL5", "CASP8", "VCAM1", "TGFB1", "PTPRC", "IKBKB",
  "STAT5A", "CCND1", "HMOX1", "TIMP1", "NLRP3", "CDK4", "PARP1", "CCR5", "FAS", "MAPK9",
  "NFE2L2", "SREBF1", "IRF1", "IL10RA", "CXCR3", "PGR", "BID", "EGR1", "F3", "AIF1",
  "CTSS", "PTGS1", "IRAK4", "LYN", "SREBF2", "TOP2A", "GFAP", "CCNA2", "PTGES", "PTPN2",
  "ERBB4", "CTSD", "CTSB", "C3", "SQLE", "HMGCR", "LSS", "CYP51A1"
)

# 创建输出目录
dir.create("enrichment_results_user_genes", recursive = TRUE, showWarnings = FALSE)

# 转换基因符号为ENTREZ ID
gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

cat("基因映射结果:\n")
print(head(gene_entrez))
cat("成功映射的基因数:", nrow(gene_entrez), "\n")
cat("未映射的基因:", setdiff(gene_list, gene_entrez$SYMBOL), "\n")

# GO富集分析
cat("\n===== GO富集分析 =====\n")

# BP分析
GO_BP <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# CC分析
GO_CC <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# MF分析
GO_MF <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# KEGG富集分析
cat("\n===== KEGG富集分析 =====\n")
KEGG <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# 保存富集结果
write.csv(as.data.frame(GO_BP), "enrichment_results_user_genes/GO_BP_results.csv", row.names = FALSE)
write.csv(as.data.frame(GO_CC), "enrichment_results_user_genes/GO_CC_results.csv", row.names = FALSE)
write.csv(as.data.frame(GO_MF), "enrichment_results_user_genes/GO_MF_results.csv", row.names = FALSE)
write.csv(as.data.frame(KEGG), "enrichment_results_user_genes/KEGG_results.csv", row.names = FALSE)

# 可视化结果
cat("\n===== 生成可视化结果 =====\n")

# GO-BP柱状图（Top 15）
pdf("enrichment_results_user_genes/GO_BP_barplot.pdf", width = 12, height = 10)
if (nrow(as.data.frame(GO_BP)) > 0) {
  go_bp_data <- as.data.frame(GO_BP)
  go_bp_top15 <- go_bp_data[order(go_bp_data$pvalue)[1:min(15, nrow(go_bp_data))], ]
  go_bp_top15$Description <- factor(go_bp_top15$Description, levels = rev(go_bp_top15$Description))
  
  print(ggplot(go_bp_top15, aes(x = Description, y = -log10(pvalue)))
    + geom_bar(stat = "identity", fill = "#377EB8")
    + coord_flip()
    + labs(title = "Top 15 GO-BP Enrichment", x = "GO Term", y = "-log10(p-value)")
    + theme_bw()
    + theme(axis.text.y = element_text(size = 8)))
} else {
  plot.new()
  text(0.5, 0.5, "No GO-BP Enrichment Results")
}
dev.off()

# KEGG气泡图（Top 15）
pdf("enrichment_results_user_genes/KEGG_bubbleplot.pdf", width = 12, height = 10)
if (nrow(as.data.frame(KEGG)) > 0) {
  kegg_data <- as.data.frame(KEGG)
  kegg_top15 <- kegg_data[order(kegg_data$pvalue)[1:min(15, nrow(kegg_data))], ]
  
  print(ggplot(kegg_top15, aes(x = GeneRatio, y = Description, size = Count, color = pvalue))
    + geom_point()
    + scale_color_viridis(option = "viridis", direction = -1)
    + labs(title = "Top 15 KEGG Enrichment", x = "Gene Ratio", y = "Pathway")
    + theme_bw()
    + theme(axis.text.y = element_text(size = 8)))
} else {
  plot.new()
  text(0.5, 0.5, "No KEGG Enrichment Results")
}
dev.off()

# GO富集网络图
pdf("enrichment_results_user_genes/GO_network.pdf", width = 12, height = 10)
if (nrow(as.data.frame(GO_BP)) > 0) {
  tryCatch({
    cnetplot(GO_BP, showCategory = 10, foldChange = NULL)
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, "Network Plot Failed")
  })
} else {
  plot.new()
  text(0.5, 0.5, "No GO-BP Enrichment Results")
}
dev.off()

# KEGG富集网络图
pdf("enrichment_results_user_genes/KEGG_network.pdf", width = 12, height = 10)
if (nrow(as.data.frame(KEGG)) > 0) {
  tryCatch({
    cnetplot(KEGG, showCategory = 10, foldChange = NULL)
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, "Network Plot Failed")
  })
} else {
  plot.new()
  text(0.5, 0.5, "No KEGG Enrichment Results")
}
dev.off()

# 生成结果摘要
cat("\n===== 分析结果摘要 =====\n")
cat("\n1. 基因映射:")
cat("\n   - 输入基因数:", length(gene_list))
cat("\n   - 成功映射基因数:", nrow(gene_entrez))

cat("\n\n2. GO富集分析:")
cat("\n   - BP显著富集条目数:", ifelse(is.null(GO_BP), 0, nrow(as.data.frame(GO_BP))))
cat("\n   - CC显著富集条目数:", ifelse(is.null(GO_CC), 0, nrow(as.data.frame(GO_CC))))
cat("\n   - MF显著富集条目数:", ifelse(is.null(GO_MF), 0, nrow(as.data.frame(GO_MF))))

cat("\n\n3. KEGG富集分析:")
cat("\n   - 显著富集通路数:", ifelse(is.null(KEGG), 0, nrow(as.data.frame(KEGG))))

# 显示Top 5 KEGG通路
if (!is.null(KEGG) && nrow(as.data.frame(KEGG)) > 0) {
  cat("\n\n4. Top 5 KEGG通路:")
  kegg_top5 <- as.data.frame(KEGG)[order(as.data.frame(KEGG)$pvalue)[1:min(5, nrow(as.data.frame(KEGG)))], ]
  for (i in 1:nrow(kegg_top5)) {
    cat(paste0("\n   ", i, ". ", kegg_top5$Description[i], " (p-value: ", round(kegg_top5$pvalue[i], 4), ")"))
  }
}

cat("\n\n分析完成！结果已保存到 enrichment_results_user_genes 目录。\n")
