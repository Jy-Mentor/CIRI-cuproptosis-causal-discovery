# GO与KEGG富集分析

# 设置CRAN镜像
options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

# 安装和加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
if (!require("enrichplot")) BiocManager::install("enrichplot")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("tidyverse")) install.packages("tidyverse")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(tidyverse)

# 设置工作目录
setwd("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写\\GSE58294 人全血")

# 核心靶点列表（24H时间点）
core_genes <- c(
  "STAT3", "MAPK14", "SMARCA4", "ALOX5", "CASP9", "PTPRC", "MDM2", "CTSB", 
  "AIF1", "NFE2L2", "ESR1", "DDIT3", "FAS", "CDC42", "PTGS2", "NFKB1", 
  "XBP1", "CCR2", "ACTA2", "JAK2", "MYBPC3", "ZEB2", "TACR1", "C5", 
  "IL1B", "PPARA", "EPHX2", "HTR2A", "FOXP2", "MMP9", "PPARG", "IDO1", 
  "CASP8", "DPP4", "PTGS1"
)

cat("核心靶点数量：", length(core_genes), "\n")

# 转换基因符号为ENTREZID
gene_list <- bitr(core_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("成功转换的基因数：", nrow(gene_list), "\n")

# 1. GO功能富集分析
cat("\n开始GO富集分析...\n")
go_enrich <- enrichGO(
  gene = gene_list$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

# 2. KEGG通路富集分析
cat("\n开始KEGG富集分析...\n")
kegg_enrich <- enrichKEGG(
  gene = gene_list$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

# 3. 保存分析结果
output_dir <- "output/enrichment_analysis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveRDS(list(
  go_enrich = go_enrich,
  kegg_enrich = kegg_enrich,
  core_genes = core_genes,
  gene_list = gene_list
), file.path(output_dir, "enrichment_results.rds"))

# 4. 可视化结果
# 4.1 GO富集分析气泡图
png(file.path(output_dir, "go_bubble_plot.png"), width = 1200, height = 800)
dotplot(go_enrich, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 4.2 KEGG富集分析气泡图
png(file.path(output_dir, "kegg_bubble_plot.png"), width = 1000, height = 800)
dotplot(kegg_enrich) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 4.3 KEGG富集分析条形图
png(file.path(output_dir, "kegg_bar_plot.png"), width = 1000, height = 800)
barplot(kegg_enrich, showCategory = 20) +
  theme(axis.text.y = element_text(size = 10))
dev.off()

# 4.4 GO富集分析网络图
png(file.path(output_dir, "go_network_plot.png"), width = 1200, height = 1000)
if (nrow(go_enrich) > 0) {
  cnetplot(go_enrich, categorySize = "pvalue", foldChange = NULL, showCategory = 10)
}
dev.off()

# 5. 生成分析报告
report_file <- file.path(output_dir, "enrichment_analysis_report.txt")
sink(report_file)

cat("GO与KEGG富集分析报告\n")
cat("======================\n\n")

# 报告GO富集分析结果
cat("1. GO功能富集分析结果\n")
cat("====================\n")
if (nrow(go_enrich) > 0) {
  go_results <- as.data.frame(go_enrich)
  cat("总富集条目数：", nrow(go_results), "\n\n")
  
  # 按 ontology 分类
  go_bp <- go_results[go_results$ONTOLOGY == "BP", ]
  go_cc <- go_results[go_results$ONTOLOGY == "CC", ]
  go_mf <- go_results[go_results$ONTOLOGY == "MF", ]
  
  cat("生物过程（BP）富集条目数：", nrow(go_bp), "\n")
  if (nrow(go_bp) > 0) {
    top_bp <- head(go_bp, 10)
    cat("Top 10 BP条目：\n")
    for (i in 1:nrow(top_bp)) {
      cat(sprintf("%d. %s (p-value: %.4f)\n", i, top_bp$Description[i], top_bp$pvalue[i]))
    }
  }
  cat("\n")
  
  cat("细胞组分（CC）富集条目数：", nrow(go_cc), "\n")
  if (nrow(go_cc) > 0) {
    top_cc <- head(go_cc, 10)
    cat("Top 10 CC条目：\n")
    for (i in 1:nrow(top_cc)) {
      cat(sprintf("%d. %s (p-value: %.4f)\n", i, top_cc$Description[i], top_cc$pvalue[i]))
    }
  }
  cat("\n")
  
  cat("分子功能（MF）富集条目数：", nrow(go_mf), "\n")
  if (nrow(go_mf) > 0) {
    top_mf <- head(go_mf, 10)
    cat("Top 10 MF条目：\n")
    for (i in 1:nrow(top_mf)) {
      cat(sprintf("%d. %s (p-value: %.4f)\n", i, top_mf$Description[i], top_mf$pvalue[i]))
    }
  }
} else {
  cat("未获得显著富集的GO条目\n")
}
cat("\n")

# 报告KEGG富集分析结果
cat("2. KEGG通路富集分析结果\n")
cat("========================\n")
if (nrow(kegg_enrich) > 0) {
  kegg_results <- as.data.frame(kegg_enrich)
  cat("总富集通路数：", nrow(kegg_results), "\n\n")
  
  top_kegg <- head(kegg_results, 15)
  cat("Top 15 KEGG通路：\n")
  for (i in 1:nrow(top_kegg)) {
    cat(sprintf("%d. %s (p-value: %.4f)\n", i, top_kegg$Description[i], top_kegg$pvalue[i]))
  }
} else {
  cat("未获得显著富集的KEGG通路\n")
}
cat("\n")

# 报告总结
cat("3. 分析总结\n")
cat("============\n")
cat("分析的核心靶点数量：", length(core_genes), "\n")
cat("成功转换为ENTREZID的基因数：", nrow(gene_list), "\n")
cat("GO富集分析完成，结果保存在：", output_dir, "\n")
cat("KEGG富集分析完成，结果保存在：", output_dir, "\n")
cat("\n可视化结果包括：\n")
cat("- go_bubble_plot.png: GO富集分析气泡图\n")
cat("- kegg_bubble_plot.png: KEGG富集分析气泡图\n")
cat("- kegg_bar_plot.png: KEGG富集分析条形图\n")
cat("- go_network_plot.png: GO富集分析网络图\n")

sink()

print("富集分析完成！结果保存在output/enrichment_analysis目录中。")
