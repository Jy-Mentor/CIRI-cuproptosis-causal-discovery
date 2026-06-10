options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

options(future.globals.maxSize = 5000 * 1024^2)

if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"dplyr" %in% installed.packages()){install.packages('dplyr')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
if(!"patchwork" %in% installed.packages()){install.packages('patchwork')}
if(!"openxlsx" %in% installed.packages()){install.packages('openxlsx')}

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(openxlsx)

dir_save <- "./result"
sc_obj <- readRDS(paste0(dir_save, "/sc_merged_QC.rds"))

cat("=== Step 1: Normalization and Variable Features ===\n")
sc_obj <- NormalizeData(sc_obj, verbose = FALSE)
sc_obj <- FindVariableFeatures(sc_obj, nfeatures = 2000, verbose = FALSE)
sc_obj <- ScaleData(sc_obj, verbose = FALSE)

cat("=== Step 2: PCA and UMAP ===\n")
sc_obj <- RunPCA(sc_obj, npcs = 50, verbose = FALSE)
sc_obj <- RunUMAP(sc_obj, dims = 1:50, verbose = FALSE)

cat("=== Step 3: Clustering ===\n")
sc_obj <- FindNeighbors(sc_obj, dims = 1:50, verbose = FALSE)
sc_obj <- FindClusters(sc_obj, resolution = 0.5, verbose = FALSE)
sc_obj <- JoinLayers(sc_obj)

p1 <- DimPlot(sc_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  labs(title = "UMAP by Condition") +
  DimPlot(sc_obj, reduction = "umap", group.by = "condition", pt.size = 0.5)
ggsave(p1, filename = paste0(dir_save, "/UMAP_condition.png"), width = 14, height = 6)

p2 <- DimPlot(sc_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  labs(title = "Clusters")
ggsave(p2, filename = paste0(dir_save, "/UMAP_clusters.png"), width = 8, height = 6)

cat("=== Step 4: Manual Cell Type Annotation using Markers ===\n")
cluster_markers <- list(
  Excitatory_Neurons = c("Slc17a7", "Satb2", "Cux2", "Rorb"),
  Inhibitory_Neurons = c("Gad1", "Gad2", "Pvalb", "Sst"),
  Astrocytes = c("Gfap", "Aqp4", "Aldoc", "Gja1"),
  Oligodendrocytes = c("Mbp", "Plp1", "Cnp", "Mog"),
  Microglia = c("Cx3cr1", "P2ry12", "Tmem119", "Hexb"),
  Endothelial = c("Cldn5", "Flt1", "Pecam1", "Kdr"),
  Pericytes = c("Pdgfrb", "Notch3", "Col1a1"),
  Smooth_Muscle = c("Acta2", "Myh11", "Tagln")
)

DefaultAssay(sc_obj) <- "RNA"
sc_obj@meta.data$cell_type <- "Unknown"

for (cluster_id in unique(sc_obj@meta.data$seurat_clusters)) {
  cat(sprintf("分析 Cluster %s...\n", cluster_id))
  cells_cluster <- subset(sc_obj, seurat_clusters == cluster_id)
  avg_exp <- rowMeans(AverageExpression(cells_cluster, features = unlist(cluster_markers), verbose = FALSE)$RNA)

  best_match <- "Unknown"
  best_score <- 0
  for (cell_type in names(cluster_markers)) {
    markers_set <- cluster_markers[[cell_type]]
    markers_found <- markers_set[markers_set %in% names(avg_exp)]
    if (length(markers_found) > 0) {
      score <- mean(avg_exp[markers_found], na.rm = TRUE)
      if (score > best_score) {
        best_score <- score
        best_match <- cell_type
      }
    }
  }
  sc_obj@meta.data$cell_type[sc_obj@meta.data$seurat_clusters == cluster_id] <- best_match
  cat(sprintf("Cluster %s: %s (score=%.3f)\n", cluster_id, best_match, best_score))
}

sc_obj$cell_type <- as.factor(sc_obj$cell_type)

p3 <- DimPlot(sc_obj, reduction = "umap", group.by = "cell_type", label = TRUE, pt.size = 0.5) +
  labs(title = "Cell Type Annotation") +
  theme(legend.position = "right")
ggsave(p3, filename = paste0(dir_save, "/UMAP_celltype.png"), width = 10, height = 6)

saveRDS(sc_obj, file = paste0(dir_save, "/sc_annotated.rds"))

cat("=== Step 5: Key Genes Expression ===\n")
key_genes <- c("Ager", "Nfkb1", "Fdx1", "Tlr4", "Stat1", "Stat3",
               "Tgfb1", "Nfe2l2", "Jak1", "Ccl2", "Icam1", "Hmox1")
key_genes <- key_genes[key_genes %in% rownames(sc_obj)]

for (gene in key_genes) {
  cat(sprintf("%s: 在数据中\n", gene))
}

p5 <- VlnPlot(sc_obj, features = key_genes, group.by = "cell_type", pt.size = 0, ncol = 4) &
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p5, filename = paste0(dir_save, "/key_genes_violin.png"), width = 14, height = 10)

cat("=== Step 6: DEGs between Sham and MCAO ===\n")
deg_results <- list()

for (cell_type_id in levels(sc_obj$cell_type)) {
  cat(sprintf("分析 %s 的差异基因...\n", cell_type_id))
  cells_ct <- subset(sc_obj, cell_type == cell_type_id)
  if (length(unique(cells_ct$condition)) < 2) {
    cat(sprintf("%s 缺少两个条件，跳过\n", cell_type_id))
    next
  }
  deg <- FindMarkers(cells_ct, group.by = "condition",
                     ident.1 = "MCAO", ident.2 = "Sham",
                     test.use = "wilcox", verbose = FALSE)
  if (nrow(deg) > 0) {
    deg$gene <- rownames(deg)
    deg$cell_type <- cell_type_id
    deg_results[[cell_type_id]] <- deg
  }
}

deg_all <- do.call(rbind, deg_results)
if (!is.null(deg_all) && nrow(deg_all) > 0) {
  if (!"p_val_adj" %in% colnames(deg_all)) {
    deg_all$p_val_adj <- deg_all$p_val
  }
  deg_sig <- deg_all %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

  write.xlsx(deg_all, paste0(dir_save, "/DEG_all.xlsx"))
  write.xlsx(deg_sig, paste0(dir_save, "/DEG_significant.xlsx"))

  cat(sprintf("总DEG数: %d, 显著DEG数: %d\n", nrow(deg_all), nrow(deg_sig)))

  if (nrow(deg_sig) > 0) {
    p6 <- ggplot(deg_sig, aes(x = avg_log2FC, y = -log10(p_val_adj), color = cell_type)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      labs(title = "DEGs: MCAO vs Sham", x = "log2FC", y = "-log10(FDR)") +
      theme_bw()
    ggsave(p6, filename = paste0(dir_save, "/DEG_volcano.png"), width = 10, height = 8)
  }
}

cat("=== Step 7: Focus on BCP Axis Genes ===\n")
axis_genes <- c("Ager", "Nfkb1", "Fdx1")
if (!is.null(deg_all) && nrow(deg_all) > 0) {
  axis_deg <- deg_all[deg_all$gene %in% axis_genes, ]
  write.xlsx(axis_deg, paste0(dir_save, "/BCP_axis_DEG.xlsx"))
  cat("BCP轴基因DEG结果:\n")
  print(axis_deg)
}

cat("=== 分析完成！输出文件：===\n")
cat("1. sc_annotated.rds - 注释后的Seurat对象\n")
cat("2. UMAP图 - 条件、聚类、细胞类型\n")
cat("3. DEG_all.xlsx / DEG_significant.xlsx - 差异基因\n")
cat("4. DEG_volcano.png - 火山图\n")
cat("5. key_genes_violin.png - 关键基因表达图\n")
cat("6. BCP_axis_DEG.xlsx - BCP轴基因差异\n")
