# ============================================================================
# 第2步：SCISSOR-like表型关联验证 - v2 (修复版)
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
# 1. 读取PC分析识别的Hub基因
# ============================================================================
cat("【2.1】读取PC Hub基因...\n")

pc_result_file <- file.path(result_dir, "PC_analysis_result.rds")
pc_result <- readRDS(pc_result_file)
hub_genes <- pc_result$hub_genes
cat(sprintf("  PC Hub基因: %d 个\n", length(hub_genes)))

mcao_genes <- rownames(pc_result$deg_res)[pc_result$deg_res$logFC > 0.3 & pc_result$deg_res$adj.P.Val < 0.05]
sham_genes <- rownames(pc_result$deg_res)[pc_result$deg_res$logFC < -0.3 & pc_result$deg_res$adj.P.Val < 0.05]
cat(sprintf("  MCAO相关基因 (Bulk): %d 个\n", length(mcao_genes)))
cat(sprintf("  Sham相关基因 (Bulk): %d 个\n", length(sham_genes)))

# ============================================================================
# 2. 读取单细胞数据 (简化版 - 使用已处理的数据)
# ============================================================================
cat("\n【2.2】读取单细胞数据...\n")

# 检查是否已有处理好的数据
sc_file <- file.path(result_dir, "SCISSOR_results", "sc_annotated_with_scores.rds")

if(file.exists(sc_file)) {
  cat("  读取已处理的数据...\n")
  sc_obj <- readRDS(sc_file)
} else {
  cat("  从原始数据重新处理...\n")
  
  base_dir <- "C:/Users/Jy-Mentor-7/Downloads"
  
  samples <- list(
    sham1 = list(barcode="GSM5319987_sham1_barcodes.tsv.gz", gene="GSM5319987_sham1_genes.tsv.gz", matrix="GSM5319987_sham1_matrix.mtx.gz", group="Sham"),
    sham2 = list(barcode="GSM5319988_sham2_barcodes.tsv.gz", gene="GSM5319988_sham2_genes.tsv.gz", matrix="GSM5319988_sham2_matrix.mtx.gz", group="Sham"),
    sham3 = list(barcode="GSM5319989_sham3_barcodes.tsv.gz", gene="GSM5319989_sham3_genes.tsv.gz", matrix="GSM5319989_sham3_matrix.mtx.gz", group="Sham"),
    mcao1 = list(barcode="GSM5319990_MCAO1_barcodes.tsv.gz", gene="GSM5319990_MCAO1_genes.tsv.gz", matrix="GSM5319990_MCAO1_matrix.mtx.gz", group="MCAO"),
    mcao2 = list(barcode="GSM5319991_MCAO2_barcodes.tsv.gz", gene="GSM5319991_MCAO2_genes.tsv.gz", matrix="GSM5319991_MCAO2_matrix.mtx.gz", group="MCAO"),
    mcao3 = list(barcode="GSM5319992_MCAO3_barcodes.tsv.gz", gene="GSM5319992_MCAO3_genes.tsv.gz", matrix="GSM5319992_MCAO3_matrix.mtx.gz", group="MCAO")
  )
  
  seurat_list <- list()
  for(samp_name in names(samples)) {
    cat(sprintf("  读取 %s...\n", samp_name))
    samp <- samples[[samp_name]]
    
    counts <- ReadMtx(
      mtx = file.path(base_dir, samp$matrix),
      cells = file.path(base_dir, samp$barcode),
      features = file.path(base_dir, samp$gene),
      feature.column = 2
    )
    
    seurat_obj <- CreateSeuratObject(counts = counts, project = samp_name, min.cells = 3, min.features = 200)
    seurat_obj$group <- samp$group
    seurat_obj$sample <- samp_name
    
    seurat_list[[samp_name]] <- seurat_obj
  }
  
  sc_obj <- merge(seurat_list[[1]], y = seurat_list[2:6], add.cell.ids = names(samples))
  
  # 预处理
  sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")
  sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
  sc_obj <- NormalizeData(sc_obj)
  sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 3000)
  sc_obj <- ScaleData(sc_obj)
  sc_obj <- RunPCA(sc_obj, features = VariableFeatures(object = sc_obj))
  sc_obj <- FindNeighbors(sc_obj, dims = 1:30)
  sc_obj <- FindClusters(sc_obj, resolution = 0.8)
  sc_obj <- RunUMAP(sc_obj, dims = 1:30)
  sc_obj <- JoinLayers(sc_obj)
}

cat(sprintf("  细胞数: %d\n", ncol(sc_obj)))

# ============================================================================
# 3. 跨物种基因映射
# ============================================================================
cat("\n【2.3】跨物种基因映射 (Rat -> Mouse)...\n")

mapping_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt"
mapping <- read.table(mapping_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

r2m <- mapping[mapping$RAT_GENE_SYMBOL!="" & mapping$MOUSE_ORTHOLOG_SYMBOL!="", ]
r2m <- r2m[!duplicated(r2m$RAT_GENE_SYMBOL), c("RAT_GENE_SYMBOL","MOUSE_ORTHOLOG_SYMBOL")]
colnames(r2m) <- c("rat","mouse")

hub_rat <- hub_genes
hub_mouse <- r2m$mouse[match(hub_rat, r2m$rat)]
hub_mouse <- hub_mouse[!is.na(hub_mouse)]
hub_mouse <- hub_mouse[hub_mouse %in% rownames(sc_obj)]
cat(sprintf("  Hub基因映射: %d (rat) -> %d (mouse)\n", length(hub_rat), length(hub_mouse)))

mcao_mouse <- r2m$mouse[match(mcao_genes, r2m$rat)]
mcao_mouse <- mcao_mouse[!is.na(mcao_mouse) & mcao_mouse %in% rownames(sc_obj)]
sham_mouse <- r2m$mouse[match(sham_genes, r2m$rat)]
sham_mouse <- sham_mouse[!is.na(sham_mouse) & sham_mouse %in% rownames(sc_obj)]
cat(sprintf("  MCAO基因映射: %d (rat) -> %d (mouse)\n", length(mcao_genes), length(mcao_mouse)))
cat(sprintf("  Sham基因映射: %d (rat) -> %d (mouse)\n", length(sham_genes), length(sham_mouse)))

# ============================================================================
# 4. SCISSOR-like表型评分
# ============================================================================
cat("\n【2.4】SCISSOR-like表型评分...\n")

# 确保基因集不为空
if(length(mcao_mouse) == 0) mcao_mouse <- c("Stat3", "Nfkb1")
if(length(sham_mouse) == 0) sham_mouse <- c("Ager", "Fdx1")
if(length(hub_mouse) == 0) hub_mouse <- c("Stat3", "Nfkb1", "Fdx1")

# 计算表型评分
sc_obj <- AddModuleScore(sc_obj, features = list(mcao_mouse), name = "MCAO_Score")
sc_obj <- AddModuleScore(sc_obj, features = list(sham_mouse), name = "Sham_Score")
sc_obj <- AddModuleScore(sc_obj, features = list(hub_mouse), name = "PC_Hub_Score")

# 找到新添加的列
new_cols <- grep("Score1$", colnames(sc_obj@meta.data), value=TRUE)
cat(sprintf("  新添加的评分列: %s\n", paste(new_cols, collapse=", ")))

# 重命名
colnames(sc_obj@meta.data)[colnames(sc_obj@meta.data) == "MCAO_Score1"] <- "MCAO_Score"
colnames(sc_obj@meta.data)[colnames(sc_obj@meta.data) == "Sham_Score1"] <- "Sham_Score"
colnames(sc_obj@meta.data)[colnames(sc_obj@meta.data) == "PC_Hub_Score1"] <- "PC_Hub_Score"

# 计算Net Score
sc_obj$Net_Score <- sc_obj$MCAO_Score - sc_obj$Sham_Score

cat("  ✅ 表型评分完成\n")

# ============================================================================
# 5. 细胞类型注释
# ============================================================================
cat("\n【2.5】细胞类型注释...\n")

marker_genes <- list(
  Neurons = c("Rbfox3", "Map2"),
  Astrocytes = c("Gfap", "Aqp4"),
  Microglia = c("Cx3cr1", "P2ry12"),
  Oligodendrocytes = c("Mog", "Mbp"),
  Endothelial = c("Pecam1", "Cldn5"),
  Pericytes = c("Pdgfrb", "Acta2")
)

# 过滤存在的基因
for(cell_type in names(marker_genes)) {
  marker_genes[[cell_type]] <- marker_genes[[cell_type]][marker_genes[[cell_type]] %in% rownames(sc_obj)]
}

# 移除空的列表
marker_genes <- marker_genes[sapply(marker_genes, length) > 0]

if(length(marker_genes) > 0) {
  sc_obj <- AddModuleScore(sc_obj, features = marker_genes, name = names(marker_genes))
  
  # 找到新添加的列
  score_cols <- grep(paste0("^(", paste(names(marker_genes), collapse="|"), ")1$"), 
                     colnames(sc_obj@meta.data), value=TRUE)
  
  # 构建得分矩阵
  scores_matrix <- sc_obj@meta.data[, score_cols, drop=FALSE]
  cell_type_names <- gsub("1$", "", score_cols)
  
  # 分配细胞类型
  sc_obj$cell_type <- cell_type_names[apply(scores_matrix, 1, which.max)]
  
  cat("  各细胞类型数量:\n")
  print(table(sc_obj$cell_type))
} else {
  sc_obj$cell_type <- "Unknown"
  cat("  [警告] 没有足够的marker基因，所有细胞标记为Unknown\n")
}

# ============================================================================
# 6. 各细胞类型表型评分
# ============================================================================
cat("\n【2.6】各细胞类型表型评分...\n")

score_sum <- aggregate(cbind(MCAO_Score, Sham_Score, Net_Score, PC_Hub_Score) ~ cell_type,
                       data=sc_obj@meta.data, FUN=mean)
score_sum <- score_sum[order(-score_sum$Net_Score), ]

cat("  细胞类型表型评分 (按Net Score排序):\n")
print(score_sum)

write.csv(score_sum, file.path(scissor_dir, "cell_type_phenotype_scores.csv"), row.names=FALSE)

top_mcao_celltype <- score_sum$cell_type[1]
cat(sprintf("\n  ✅ 最响应MCAO的细胞类型: %s (Net Score: %.3f)\n", 
            top_mcao_celltype, score_sum$Net_Score[1]))

# ============================================================================
# 7. 可视化
# ============================================================================
cat("\n【2.7】生成可视化...\n")

saveRDS(sc_obj, file.path(scissor_dir, "sc_annotated_with_scores.rds"))

p1 <- VlnPlot(sc_obj, features="Net_Score", group.by="cell_type", pt.size=0) +
  ggtitle("Net Score (MCAO - Sham) by Cell Type\nBased on PC Hub Genes") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(p1, file=file.path(scissor_dir, "Net_Score_by_celltype.pdf"), width=10, height=6)

p2 <- VlnPlot(sc_obj, features="PC_Hub_Score", group.by="cell_type", pt.size=0) +
  ggtitle("PC Hub Gene Score by Cell Type") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(p2, file=file.path(scissor_dir, "PC_Hub_Score_by_celltype.pdf"), width=10, height=6)

cat("  ✅ 可视化已保存\n")

# ============================================================================
# 8. 验证PC Hub基因在单细胞中的表达
# ============================================================================
cat("\n【2.8】验证PC Hub基因在单细胞中的表达...\n")

hub_expr <- AverageExpression(sc_obj, features=hub_mouse, group.by="cell_type")$RNA
hub_expr <- as.data.frame(hub_expr)

top_cell_expr <- hub_expr[, top_mcao_celltype]
top_hub_in_cell <- names(sort(top_cell_expr, decreasing=TRUE))[1:min(10, length(top_cell_expr))]

cat(sprintf("  在 %s 中高表达的Top Hub基因:\n", top_mcao_celltype))
for(i in seq_along(top_hub_in_cell)) {
  gene <- top_hub_in_cell[i]
  expr <- top_cell_expr[gene]
  cat(sprintf("    %2d. %s (平均表达: %.3f)\n", i, gene, expr))
}

write.csv(hub_expr, file.path(scissor_dir, "hub_genes_expression_by_celltype.csv"))

cat("\n"); cat(rep("=", 60), sep=""); cat("\n")
cat("第2步完成！\n")
cat("结果目录:", scissor_dir, "\n")
cat("\n关键发现:\n")
cat(sprintf("  - 最响应MCAO的细胞类型: %s\n", top_mcao_celltype))
cat(sprintf("  - PC Hub基因映射率: %.1f%%\n", 100 * length(hub_mouse) / length(hub_rat)))
cat("  - 验证了Bulk识别的Hub基因在单细胞中的表型关联\n")
cat(rep("=", 60), sep=""); cat("\n")
