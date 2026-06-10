# 完整流程: 从原始数据到scTenifoldKnk虚拟敲除
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

# 加载包
if(!"Seurat" %in% installed.packages()){install.packages('Seurat', repos='https://cran.rstudio.com')}
if(!"Matrix" %in% installed.packages()){install.packages('Matrix')}
library(Seurat)
library(Matrix)

# 设置路径
base_dir <- "C:/Users/Jy-Mentor-7/Downloads"
result_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/analysis_results"
dir.create(result_dir, showWarnings=FALSE, recursive=TRUE)

cat("=== 完整分析流程 ===\n\n")

# ============================================================================
# 步骤1: 读取单细胞数据
# ============================================================================
cat("【步骤1】读取GSE174574单细胞数据...\n")

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

# 合并
cat("  合并6个样本...\n")
sc_obj <- merge(seurat_list[[1]], y = seurat_list[2:6], add.cell.ids = names(samples))
cat(sprintf("  总细胞数: %d\n", ncol(sc_obj)))

# ============================================================================
# 步骤2: QC和标准化
# ============================================================================
cat("\n【步骤2】质控和标准化...\n")
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
cat(sprintf("  质控后细胞数: %d\n", ncol(sc_obj)))

sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 3000)

# ============================================================================
# 步骤3: 降维和聚类
# ============================================================================
cat("\n【步骤3】降维和聚类...\n")
sc_obj <- ScaleData(sc_obj)
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(object = sc_obj))
sc_obj <- FindNeighbors(sc_obj, dims = 1:30)
sc_obj <- FindClusters(sc_obj, resolution = 0.8)
sc_obj <- RunUMAP(sc_obj, dims = 1:30)

# ============================================================================
# 步骤4: 细胞注释
# ============================================================================
cat("\n【步骤4】细胞类型注释...\n")
markers <- list(
  Neurons = c("Rbfox3", "Map2", "Syn1"),
  Astrocytes = c("Gfap", "Aqp4", "Slc1a3"),
  Microglia = c("Cx3cr1", "P2ry12", "Trem2"),
  Oligodendrocytes = c("Mog", "Mbp", "Plp1"),
  Endothelial = c("Pecam1", "Cldn5", "Flt1"),
  Pericytes = c("Pdgfrb", "Cspg4", "Anpep")
)

sc_obj$cell_type <- "Unknown"
for(cell_type in names(markers)) {
  genes <- markers[[cell_type]]
  genes <- genes[genes %in% rownames(sc_obj)]
  if(length(genes) > 0) {
    scores <- Matrix::rowMeans(sc_obj@assays$RNA$scale.data[genes, , drop=FALSE])
    sc_obj$cell_type[scores > 0.5] <- cell_type
  }
}

cat(sprintf("  各细胞类型数量:\n"))
print(table(sc_obj$cell_type))

# 保存
saveRDS(sc_obj, file.path(result_dir, "sc_annotated.rds"))
cat("  已保存到:", file.path(result_dir, "sc_annotated.rds"), "\n")

# ============================================================================
# 步骤5: scTenifoldKnk虚拟敲除
# ============================================================================
cat("\n【步骤5】scTenifoldKnk虚拟敲除...\n")

# 加载scTenifoldKnk
if(!require(scTenifoldKnk, quietly=TRUE)) {
  library(scTenifoldNet)
  pkg_path <- "C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/scripts/scTenifoldKnk/scTenifoldKnk-master"
  install.packages(pkg_path, repos=NULL, type="source")
}
library(scTenifoldKnk)

# 提取count矩阵
count_matrix <- GetAssayData(sc_obj, layer="counts")
cat(sprintf("  Count矩阵: %d 基因 x %d 细胞\n", nrow(count_matrix), ncol(count_matrix)))

# 目标基因
target_genes <- c("Nfkb1", "Fdx1", "Tlr4")
knockout_results <- list()

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

    knockout_results[[gene]] <- knk_res

    # 保存差异基因
    degs <- knk_res$diffRegulation %>%
      filter(p.value < 0.05) %>%
      arrange(p.value)

    write.csv(degs, file.path(result_dir, paste0("knockout_", gene, "_DEGs.csv")), row.names=FALSE)
    cat(sprintf("    差异基因数: %d\n", nrow(degs)))

    # 保存结果
    saveRDS(knk_res, file.path(result_dir, paste0("knockout_", gene, "_result.rds")))
  } else {
    cat(sprintf("  [跳过] %s 不在数据中\n", gene))
  }
}

# 保存所有结果
saveRDS(knockout_results, file.path(result_dir, "all_knockout_results.rds"))

cat("\n=== 完成 ===\n")
cat("结果保存到:", result_dir, "\n")
