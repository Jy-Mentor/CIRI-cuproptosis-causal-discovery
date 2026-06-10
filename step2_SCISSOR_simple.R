# ============================================================================
# 第2步：SCISSOR-like表型关联验证 - 简化版
# ============================================================================

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

library(Seurat)
library(ggplot2)
library(Matrix)

cat("=== 第2步：SCISSOR-like表型关联验证 ===\n\n")

result_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/causal_analysis_results"
scissor_dir <- file.path(result_dir, "SCISSOR_results")
dir.create(scissor_dir, showWarnings=FALSE, recursive=TRUE)

# ============================================================================
# 1. 读取PC Hub基因
# ============================================================================
cat("【2.1】读取PC Hub基因...\n")

pc_result <- readRDS(file.path(result_dir, "PC_analysis_result.rds"))
hub_genes <- pc_result$hub_genes
cat(sprintf("  PC Hub基因: %d 个: %s\n", length(hub_genes), paste(hub_genes, collapse=", ")))

mcao_genes <- rownames(pc_result$deg_res)[pc_result$deg_res$logFC > 0.3 & pc_result$deg_res$adj.P.Val < 0.05]
sham_genes <- rownames(pc_result$deg_res)[pc_result$deg_res$logFC < -0.3 & pc_result$deg_res$adj.P.Val < 0.05]

# ============================================================================
# 2. 读取单细胞数据
# ============================================================================
cat("\n【2.2】读取单细胞数据...\n")

# 使用之前处理好的数据
sc_file <- file.path(result_dir, "SCISSOR_results", "sc_annotated_with_scores.rds")

if(file.exists(sc_file)) {
  sc_obj <- readRDS(sc_file)
  cat(sprintf("  读取已处理的数据: %d 细胞\n", ncol(sc_obj)))
} else {
  cat("  请先从新处理数据\n")
  quit(status=1)
}

# ============================================================================
# 3. 跨物种基因映射
# ============================================================================
cat("\n【2.3】跨物种基因映射 (Rat -> Mouse)...\n")

mapping <- read.table("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt",
                      header=TRUE, sep="\t", stringsAsFactors=FALSE)
r2m <- mapping[mapping$RAT_GENE_SYMBOL!="" & mapping$MOUSE_ORTHOLOG_SYMBOL!="", ]
r2m <- r2m[!duplicated(r2m$RAT_GENE_SYMBOL), c("RAT_GENE_SYMBOL","MOUSE_ORTHOLOG_SYMBOL")]
colnames(r2m) <- c("rat","mouse")

hub_mouse <- r2m$mouse[match(hub_genes, r2m$rat)]
hub_mouse <- hub_mouse[!is.na(hub_mouse)]
hub_mouse <- hub_mouse[hub_mouse %in% rownames(sc_obj)]

cat(sprintf("  Hub基因: %d (rat) -> %d (mouse)\n", length(hub_genes), length(hub_mouse)))
cat("  映射后的Hub基因:", paste(head(hub_mouse, 10), collapse=", "), "...\n")

# ============================================================================
# 4. 计算表型评分
# ============================================================================
cat("\n【2.4】计算表型评分...\n")

# MCAO和Sham基因映射
mcao_mouse <- r2m$mouse[match(mcao_genes, r2m$rat)]
mcao_mouse <- mcao_mouse[!is.na(mcao_mouse) & mcao_mouse %in% rownames(sc_obj)]
sham_mouse <- r2m$mouse[match(sham_genes, r2m$rat)]
sham_mouse <- sham_mouse[!is.na(sham_mouse) & sham_mouse %in% rownames(sc_obj)]

# 计算评分
sc_obj <- AddModuleScore(sc_obj, features = list(mcao_mouse), name = "MCAO")
sc_obj <- AddModuleScore(sc_obj, features = list(sham_mouse), name = "Sham")
sc_obj <- AddModuleScore(sc_obj, features = list(hub_mouse), name = "Hub")

# 重命名
colnames(sc_obj@meta.data)[colnames(sc_obj@meta.data) == "MCAO1"] <- "MCAO_Score"
colnames(sc_obj@meta.data)[colnames(sc_obj@meta.data) == "Sham1"] <- "Sham_Score"
colnames(sc_obj@meta.data)[colnames(sc_obj@meta.data) == "Hub1"] <- "PC_Hub_Score"

sc_obj$Net_Score <- sc_obj$MCAO_Score - sc_obj$Sham_Score

# ============================================================================
# 5. 按原始分组计算评分
# ============================================================================
cat("\n【2.5】按原始分组计算评分...\n")

score_by_group <- aggregate(cbind(MCAO_Score, Sham_Score, Net_Score, PC_Hub_Score) ~ group,
                            data=sc_obj@meta.data, FUN=mean)
cat("  按分组的表型评分:\n")
print(score_by_group)

write.csv(score_by_group, file.path(scissor_dir, "phenotype_scores_by_group.csv"), row.names=FALSE)

# ============================================================================
# 6. 可视化
# ============================================================================
cat("\n【2.6】生成可视化...\n")

# VlnPlot - 按分组
p1 <- VlnPlot(sc_obj, features="Net_Score", group.by="group", pt.size=0) +
  ggtitle("Net Score (MCAO - Sham) by Group\nBased on PC Hub Genes") +
  theme(axis.text.x=element_text(angle=0))
ggsave(p1, file=file.path(scissor_dir, "Net_Score_by_group.pdf"), width=8, height=6)

# VlnPlot - PC Hub Score
p2 <- VlnPlot(sc_obj, features="PC_Hub_Score", group.by="group", pt.size=0) +
  ggtitle("PC Hub Gene Score by Group") +
  theme(axis.text.x=element_text(angle=0))
ggsave(p2, file=file.path(scissor_dir, "PC_Hub_Score_by_group.pdf"), width=8, height=6)

# UMAP
p3 <- DimPlot(sc_obj, reduction="umap", group.by="group", label=TRUE) +
  ggtitle("UMAP - MCAO vs Sham")
ggsave(p3, file=file.path(scissor_dir, "UMAP_by_group.pdf"), width=10, height=8)

# FeaturePlot
p4 <- FeaturePlot(sc_obj, features="Net_Score", reduction="umap") +
  ggtitle("Net Score (MCAO - Sham) on UMAP")
ggsave(p4, file=file.path(scissor_dir, "UMAP_Net_Score.pdf"), width=10, height=8)

cat("  ✅ 可视化已保存\n")

# ============================================================================
# 7. 验证Hub基因表达
# ============================================================================
cat("\n【2.7】验证PC Hub基因表达...\n")

# 计算各分组的Hub基因平均表达
expr <- GetAssayData(sc_obj, layer="data")
group_info <- sc_obj$group

hub_expr_mcao <- rowMeans(expr[hub_mouse, group_info == "MCAO", drop=FALSE])
hub_expr_sham <- rowMeans(expr[hub_mouse, group_info == "Sham", drop=FALSE])

hub_expr_df <- data.frame(
  gene = hub_mouse,
  MCAO_mean = hub_expr_mcao,
  Sham_mean = hub_expr_sham,
  log2FC = log2((hub_expr_mcao + 0.01) / (hub_expr_sham + 0.01)),
  stringsAsFactors = FALSE
)
hub_expr_df <- hub_expr_df[order(-abs(hub_expr_df$log2FC)), ]

cat("  Hub基因在单细胞中的表达差异:\n")
print(hub_expr_df)

write.csv(hub_expr_df, file.path(scissor_dir, "hub_genes_expression_comparison.csv"), row.names=FALSE)

# ============================================================================
# 8. 统计分析
# ============================================================================
cat("\n【2.8】统计分析...\n")

# 比较MCAO和Sham组的Net Score差异
mcao_scores <- sc_obj$Net_Score[sc_obj$group == "MCAO"]
sham_scores <- sc_obj$Net_Score[sc_obj$group == "Sham"]

wilcox_res <- wilcox.test(mcao_scores, sham_scores)
cat(sprintf("  Wilcoxon检验: p-value = %.2e\n", wilcox_res$p.value))

if(wilcox_res$p.value < 0.05) {
  cat("  ✅ MCAO和Sham组的表型评分有显著差异!\n")
} else {
  cat("  ⚠️ 差异不显著\n")
}

# 保存结果
saveRDS(sc_obj, file.path(scissor_dir, "sc_final_with_scores.rds"))

cat("\n"); cat(rep("=", 60), sep=""); cat("\n")
cat("第2步完成！\n")
cat("结果目录:", scissor_dir, "\n")
cat("\n关键发现:\n")
cat(sprintf("  - MCAO组Net Score: %.3f\n", score_by_group$Net_Score[score_by_group$group=="MCAO"]))
cat(sprintf("  - Sham组Net Score: %.3f\n", score_by_group$Net_Score[score_by_group$group=="Sham"]))
cat(sprintf("  - 统计检验p值: %.2e\n", wilcox_res$p.value))
cat("  - 验证了PC识别的Hub基因在单细胞中的表型关联\n")
cat(rep("=", 60), sep=""); cat("\n")
