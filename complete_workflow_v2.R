# 完整流程 v2 - 修复版
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

library(Seurat)
library(Matrix)

base_dir <- "C:/Users/Jy-Mentor-7/Downloads"
result_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/analysis_results"
dir.create(result_dir, showWarnings=FALSE, recursive=TRUE)

cat("=== 完整分析流程 v2 ===\n\n")

# ============================================================================
# 步骤1-3: 数据读取和标准化
# ============================================================================
cat("【步骤1-3】读取数据并标准化...\n")

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
cat(sprintf("  总细胞数: %d\n", ncol(sc_obj)))

# QC和标准化
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
cat(sprintf("  质控后细胞数: %d\n", ncol(sc_obj)))

sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 3000)
sc_obj <- ScaleData(sc_obj)
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(object = sc_obj))
sc_obj <- FindNeighbors(sc_obj, dims = 1:30)
sc_obj <- FindClusters(sc_obj, resolution = 0.8)
sc_obj <- RunUMAP(sc_obj, dims = 1:30)

# ============================================================================
# 步骤4: 细胞注释 (修复版)
# ============================================================================
cat("\n【步骤4】细胞类型注释...\n")

# 使用AddModuleScore计算细胞类型得分
markers <- list(
  Neurons = c("Rbfox3", "Map2", "Syn1"),
  Astrocytes = c("Gfap", "Aqp4", "Slc1a3"),
  Microglia = c("Cx3cr1", "P2ry12", "Trem2"),
  Oligodendrocytes = c("Mog", "Mbp", "Plp1"),
  Endothelial = c("Pecam1", "Cldn5", "Flt1"),
  Pericytes = c("Pdgfrb", "Cspg4", "Anpep")
)

# 过滤存在的基因
for(cell_type in names(markers)) {
  markers[[cell_type]] <- markers[[cell_type]][markers[[cell_type]] %in% rownames(sc_obj)]
}

# 计算模块得分
sc_obj <- AddModuleScore(sc_obj, features = markers, name = names(markers))

# 根据最高得分分配细胞类型
score_cols <- paste0(names(markers), "1")
scores <- sc_obj@meta.data[, score_cols]
colnames(scores) <- names(markers)

sc_obj$cell_type <- colnames(scores)[apply(scores, 1, which.max)]

cat("  各细胞类型数量:\n")
print(table(sc_obj$cell_type))

# ============================================================================
# 步骤5: scTenifoldKnk虚拟敲除
# ============================================================================
cat("\n【步骤5】scTenifoldKnk虚拟敲除...\n")

library(scTenifoldNet)
library(scTenifoldKnk)

# 提取count矩阵并转为dense matrix (scTenifoldKnk需要)
cat("  准备count矩阵...\n")
sc_obj <- JoinLayers(sc_obj)
count_matrix <- as.matrix(GetAssayData(sc_obj, layer="counts"))

# 过滤低表达基因
keep_genes <- rowSums(count_matrix > 0) >= 50
count_matrix <- count_matrix[keep_genes, ]
cat(sprintf("  Count矩阵: %d 基因 x %d 细胞\n", nrow(count_matrix), ncol(count_matrix)))

# 保存矩阵
saveRDS(count_matrix, file.path(result_dir, "sc_count_matrix.rds"))

# 目标基因
target_genes <- c("Nfkb1", "Fdx1", "Tlr4")

for(gene in target_genes) {
  if(gene %in% rownames(count_matrix)) {
    cat(sprintf("\n  虚拟敲除 %s...\n", gene))

    set.seed(666)
    knk_res <- scTenifoldKnk(
      countMatrix = count_matrix,
      gKO = gene,
      qc = TRUE,
      qc_mtThreshold = 0.1,
      qc_minLSize = 1000,
      nc_lambda = 0,
      nc_nNet = 10,
      nc_nCells = 500,
      nc_nComp = 3,
      nc_q = 0.9,
      td_K = 3,
      ma_nDim = 2,
      nCores = 2
    )

    # 保存差异基因
    degs <- knk_res$diffRegulation %>%
      filter(p.value < 0.05) %>%
      arrange(p.value)

    write.csv(degs, file.path(result_dir, paste0("knockout_", gene, "_DEGs.csv")), row.names=FALSE)
    cat(sprintf("    差异基因数: %d\n", nrow(degs)))

    # 显示top10
    cat("    Top 10差异基因:\n")
    for(i in 1:min(10, nrow(degs))) {
      cat(sprintf("      %s: logFC=%.3f, p=%.2e\n", degs$gene[i], log2(degs$FC[i]), degs$p.value[i]))
    }

    # 保存完整结果
    saveRDS(knk_res, file.path(result_dir, paste0("knockout_", gene, "_result.rds")))
  } else {
    cat(sprintf("  [跳过] %s 不在数据中\n", gene))
  }
}

# 保存Seurat对象
saveRDS(sc_obj, file.path(result_dir, "sc_annotated.rds"))

cat("\n=== 完成 ===\n")
cat("结果保存到:", result_dir, "\n")
