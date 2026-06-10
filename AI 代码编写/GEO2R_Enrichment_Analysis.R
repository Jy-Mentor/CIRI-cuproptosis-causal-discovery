# ============================================================================= 
# GEO2R 差异基因结果处理与富集分析完整脚本 
# 输入：GEO2R 下载的 .tsv 文件（通常包含 logFC、adj.P.Val 等） 
# 输出： 
#   - significant_genes.csv：显著差异基因列表 
#   - GO_BP_results.csv：GO 生物过程富集结果 
#   - GO_CC_results.csv：GO 细胞组分富集结果 
#   - GO_MF_results.csv：GO 分子功能富集结果 
#   - KEGG_results.csv：KEGG 通路富集结果 
#   - unmapped_genes.txt：未映射基因列表
#   - 可视化图表（气泡图、柱状图） 
# ============================================================================= 

# 0. 设置工作目录（可选，请改为您自己的路径） 
setwd("C:/Users/Jy-Mentor-7/Desktop/大创") 

# 1. 安装并加载必要的 R 包 
packages <- c("BiocManager", "clusterProfiler", "org.Mm.eg.db", "enrichplot", 
              "ggplot2", "ggrepel", "pheatmap") 

# 智能加载包函数
load_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("⚠  需要安装包:", pkg))
    if (pkg == "BiocManager") {
      install.packages("BiocManager")
    } else {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
  message(paste("✓ 已加载:", pkg))
}

for (pkg in packages) {
  load_package(pkg)
}

# 2. 读取 GEO2R 结果文件（请修改为您的文件路径） 
file_path <- "C:/Users/Jy-Mentor-7/Desktop/大创/GSE61616.top.table (1).tsv" 
deg <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE) 

# 查看数据结构和列名 
cat("数据维度：", dim(deg), "\n") 
cat("列名：", colnames(deg), "\n") 
head(deg) 

# 3. 确定关键列名（根据实际列名修改） 
# 通常 GEO2R 结果包含： 
#   - Gene.symbol（基因符号） 
#   - logFC（log2 倍数变化） 
#   - adj.P.Val（校正后 P 值） 
# 如果列名不同，请手动修改以下变量 
gene_col <- "Gene.symbol"    # 基因符号所在列名 
logFC_col <- "logFC"         # logFC 所在列名 
pval_col <- "adj.P.Val"      # 校正后 P 值所在列名 

# 检查列是否存在 
if (!all(c(gene_col, logFC_col, pval_col) %in% colnames(deg))) {
  stop("请根据实际列名修改 gene_col、logFC_col、pval_col 变量") 
}

# 4. 筛选显著差异基因（阈值可根据需要调整） 
logFC_cutoff <- 1            # |logFC| > 1 
pval_cutoff <- 0.05          # adj.P.Val < 0.05 

sig_genes <- subset(deg, abs(deg[[logFC_col]]) > logFC_cutoff & deg[[pval_col]] < pval_cutoff) 

cat("显著差异基因数量：", nrow(sig_genes), "\n") 

# 保存显著基因列表 
write.csv(sig_genes, "significant_genes.csv", row.names = FALSE) 

# 5. 提取基因符号列表并处理可能的多基因条目（如 "TP53 /// ABC"） 
gene_symbols <- sig_genes[[gene_col]] 

# 去除空值或 NA 
gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""] 

# 拆分多基因条目（按 " /// " 分隔），并取唯一值 
gene_symbols <- unique(unlist(strsplit(gene_symbols, " /// "))) 

cat("处理后的唯一基因符号数量：", length(gene_symbols), "\n") 

# 6. 将基因符号转换为 ENTREZ ID（用于 KEGG 富集） 
# 首先使用 bitr 进行转换
gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db) 

# 处理未映射的基因
unmapped_genes <- setdiff(gene_symbols, gene_entrez$SYMBOL)

# 尝试使用基因映射库进行补充映射
if (length(unmapped_genes) > 0) {
  cat("未映射的基因数量：", length(unmapped_genes), "\n")
  
  # 读取基因映射库
  mapping_file <- "C:/Users/Jy-Mentor-7/Desktop/大创/大鼠 小鼠 人类映射库.txt"
  if (file.exists(mapping_file)) {
    cat("正在使用基因映射库进行补充映射...\n")
    mapping_data <- read.delim(mapping_file, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
    
    # 筛选小鼠基因映射
    mouse_mapping <- mapping_data[!is.na(mapping_data$MOUSE_ORTHOLOG_SYMBOL), ]
    
    # 尝试为未映射基因找到映射
    additional_mappings <- data.frame()
    for (gene in unmapped_genes) {
      # 查找小鼠基因的人类同源基因
      match <- mouse_mapping[mouse_mapping$MOUSE_ORTHOLOG_SYMBOL == gene, ]
      if (nrow(match) > 0 && !is.na(match$HUMAN_ORTHOLOG_NCBI_GENE_ID)) {
        # 尝试使用人类同源基因的 ENTREZ ID
        human_entrez <- match$HUMAN_ORTHOLOG_NCBI_GENE_ID[1]
        additional_mappings <- rbind(additional_mappings, 
                                   data.frame(SYMBOL = gene, 
                                              ENTREZID = human_entrez))
      }
    }
    
    if (nrow(additional_mappings) > 0) {
      cat("通过映射库补充映射的基因数量：", nrow(additional_mappings), "\n")
      gene_entrez <- rbind(gene_entrez, additional_mappings)
      # 更新未映射基因列表
      unmapped_genes <- setdiff(unmapped_genes, additional_mappings$SYMBOL)
    }
  }
  
  if (length(unmapped_genes) > 0) {
    writeLines(unmapped_genes, "unmapped_genes.txt")
    cat("最终未映射基因已保存到 unmapped_genes.txt\n")
  } else {
    cat("所有基因都成功映射\n")
  }
} else {
  cat("所有基因都成功映射\n")
}

# 查看转换成功的比例 
cat("成功转换为 ENTREZ ID 的基因数：", nrow(gene_entrez), "\n") 
if (nrow(gene_entrez) == 0) {
  stop("没有基因能成功转换为 ENTREZ ID，请检查基因符号格式") 
}

# 7. GO 富集分析（分别进行 BP、CC、MF） 
go_bp <- enrichGO(gene = gene_entrez$ENTREZID, 
                   OrgDb = org.Mm.eg.db, 
                   ont = "BP", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.2, 
                   readable = TRUE) 

go_cc <- enrichGO(gene = gene_entrez$ENTREZID, 
                   OrgDb = org.Mm.eg.db, 
                   ont = "CC", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.2, 
                   readable = TRUE) 

go_mf <- enrichGO(gene = gene_entrez$ENTREZID, 
                   OrgDb = org.Mm.eg.db, 
                   ont = "MF", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.2, 
                   readable = TRUE) 

# 8. KEGG 通路富集分析 
tryCatch({
  kegg <- enrichKEGG(gene = gene_entrez$ENTREZID, 
                    organism = "mmu", 
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2)
  cat("KEGG 富集分析成功完成\n")
}, error = function(e) {
  cat("KEGG 富集分析失败（可能是网络连接问题）：", e$message, "\n")
  kegg <- NULL
}) 

# 9. 查看富集结果概要 
cat("\nGO-BP 富集条目数：", nrow(go_bp), "\n") 
cat("GO-CC 富集条目数：", nrow(go_cc), "\n") 
cat("GO-MF 富集条目数：", nrow(go_mf), "\n") 
cat("KEGG 富集条目数：", nrow(kegg), "\n") 

# 10. 保存富集结果为 CSV 文件 
write.csv(as.data.frame(go_bp), "GO_BP_results.csv", row.names = FALSE) 
write.csv(as.data.frame(go_cc), "GO_CC_results.csv", row.names = FALSE) 
write.csv(as.data.frame(go_mf), "GO_MF_results.csv", row.names = FALSE) 
write.csv(as.data.frame(kegg), "KEGG_results.csv", row.names = FALSE) 

# 11. 可视化 
# 11.1 KEGG 气泡图（Top 15） 
if (nrow(kegg) > 0) {
  p1 <- dotplot(kegg, showCategory = 15) + 
        ggtitle("KEGG Pathway Enrichment") + 
        theme(plot.title = element_text(hjust = 0.5)) 
  ggsave("KEGG_dotplot.pdf", p1, width = 10, height = 8) 
  print(p1) 
}

# 11.2 GO-BP 柱状图（Top 15） 
if (nrow(go_bp) > 0) {
  p2 <- barplot(go_bp, showCategory = 15) + 
        ggtitle("GO Biological Process Enrichment") + 
        theme(plot.title = element_text(hjust = 0.5)) 
  ggsave("GO_BP_barplot.pdf", p2, width = 10, height = 8) 
  print(p2) 
}

# 11.3 GO 网络图（需要 enrichplot） 
if (nrow(go_bp) > 0) {
  p3 <- cnetplot(go_bp, showCategory = 5, foldChange = NULL) + 
        ggtitle("GO-BP Network") 
  ggsave("GO_BP_network.pdf", p3, width = 12, height = 10) 
  print(p3) 
}

# 12. 生成分析报告（可选） 
sink("analysis_report.txt") 
cat("========================================\n") 
cat("GEO2R 差异基因富集分析报告\n") 
cat("========================================\n\n") 
cat("分析时间：", date(), "\n") 
cat("输入文件：", file_path, "\n") 
cat("筛选阈值：|logFC| >", logFC_cutoff, "且 adj.P.Val <", pval_cutoff, "\n") 
cat("显著基因数量：", nrow(sig_genes), "\n") 
cat("成功转换 ENTREZ ID 基因数：", nrow(gene_entrez), "\n")
if (length(unmapped_genes) > 0) {
  cat("未映射基因数：", length(unmapped_genes), "\n")
}
cat("\n") 
cat("GO-BP 富集条目数：", nrow(go_bp), "\n") 
cat("GO-CC 富集条目数：", nrow(go_cc), "\n") 
cat("GO-MF 富集条目数：", nrow(go_mf), "\n") 
cat("KEGG 富集条目数：", nrow(kegg), "\n\n") 
cat("Top 10 KEGG 通路：\n") 
if (nrow(kegg) > 0) {
  print(head(as.data.frame(kegg)[, c("ID", "Description", "pvalue", "p.adjust", "Count")], 10)) 
}
if (length(unmapped_genes) > 0) {
  cat("\n未映射基因（前 10 个）：\n")
  print(head(unmapped_genes, 10))
}
sink() 

cat("\n所有分析完成！结果已保存到当前工作目录：", getwd(), "\n")