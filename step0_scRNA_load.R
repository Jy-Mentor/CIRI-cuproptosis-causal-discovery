options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"BiocManager" %in% installed.packages()){install.packages('BiocManager')}
if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"SeuratObject" %in% installed.packages()){BiocManager::install('SeuratObject')}
if(!"dplyr" %in% installed.packages()){install.packages('dplyr')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
if(!"patchwork" %in% installed.packages()){install.packages('patchwork')}

library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(patchwork)

dir_data <- "./"
dir_save <- "./result"
if(!dir.exists(dir_save)){dir.create(dir_save, recursive = T)}

read_10x_mtx <- function(mtx_file, genes_file, barcodes_file) {
  genes_df <- read.table(genes_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  barcodes <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)[, 1]

  mtx_conn <- file(mtx_file, "r")
  mtx_header <- readLines(mtx_conn, n = 10)
  close(mtx_conn)

  mtx_data <- as.matrix(read.table(mtx_file, skip = 3, header = FALSE))
  ngenes <- as.integer(strsplit(mtx_header[3], "\\s+")[[1]][1])
  ncells <- as.integer(strsplit(mtx_header[3], "\\s+")[[1]][2])

  gene_ids <- genes_df[, 1]
  gene_symbols <- genes_df[, 2]

  dup_idx <- which(duplicated(gene_symbols))
  if (length(dup_idx) > 0) {
    gene_symbols[dup_idx] <- paste0(gene_symbols[dup_idx], "_dup", seq_along(dup_idx))
  }
  names(gene_symbols) <- gene_ids

  sparse_mat <- Matrix::sparseMatrix(
    i = mtx_data[, 1],
    j = mtx_data[, 2],
    x = mtx_data[, 3],
    dims = c(ngenes, ncells)
  )
  rownames(sparse_mat) <- gene_symbols
  colnames(sparse_mat) <- barcodes

  return(sparse_mat)
}

cat("=== 读取Sham样本 ===\n")
sham1 <- read_10x_mtx(file.path(dir_data, "GSM5319987_sham1_matrix.mtx"),
                      file.path(dir_data, "GSM5319987_sham1_genes.tsv"),
                      file.path(dir_data, "GSM5319987_sham1_barcodes.tsv"))
sham2 <- read_10x_mtx(file.path(dir_data, "GSM5319988_sham2_matrix.mtx"),
                      file.path(dir_data, "GSM5319988_sham2_genes.tsv"),
                      file.path(dir_data, "GSM5319988_sham2_barcodes.tsv"))
sham3 <- read_10x_mtx(file.path(dir_data, "GSM5319989_sham3_matrix.mtx"),
                      file.path(dir_data, "GSM5319989_sham3_genes.tsv"),
                      file.path(dir_data, "GSM5319989_sham3_barcodes.tsv"))

cat("=== 读取MCAO样本 ===\n")
mcao1 <- read_10x_mtx(file.path(dir_data, "GSM5319990_MCAO1_matrix.mtx"),
                      file.path(dir_data, "GSM5319990_MCAO1_genes.tsv"),
                      file.path(dir_data, "GSM5319990_MCAO1_barcodes.tsv"))
mcao2 <- read_10x_mtx(file.path(dir_data, "GSM5319991_MCAO2_matrix.mtx"),
                      file.path(dir_data, "GSM5319991_MCAO2_genes.tsv"),
                      file.path(dir_data, "GSM5319991_MCAO2_barcodes.tsv"))
mcao3 <- read_10x_mtx(file.path(dir_data, "GSM5319992_MCAO3_matrix.mtx"),
                      file.path(dir_data, "GSM5319992_MCAO3_genes.tsv"),
                      file.path(dir_data, "GSM5319992_MCAO3_barcodes.tsv"))

cat("=== 创建Seurat对象 ===\n")
sham1 <- CreateSeuratObject(sham1, project = "GSE174574", assay = "RNA")
sham2 <- CreateSeuratObject(sham2, project = "GSE174574", assay = "RNA")
sham3 <- CreateSeuratObject(sham3, project = "GSE174574", assay = "RNA")
mcao1 <- CreateSeuratObject(mcao1, project = "GSE174574", assay = "RNA")
mcao2 <- CreateSeuratObject(mcao2, project = "GSE174574", assay = "RNA")
mcao3 <- CreateSeuratObject(mcao3, project = "GSE174574", assay = "RNA")

sham1@meta.data$condition <- "Sham"
sham2@meta.data$condition <- "Sham"
sham3@meta.data$condition <- "Sham"
mcao1@meta.data$condition <- "MCAO"
mcao2@meta.data$condition <- "MCAO"
mcao3@meta.data$condition <- "MCAO"

sham1@meta.data$sample <- "sham1"
sham2@meta.data$sample <- "sham2"
sham3@meta.data$sample <- "sham3"
mcao1@meta.data$sample <- "mcao1"
mcao2@meta.data$sample <- "mcao2"
mcao3@meta.data$sample <- "mcao3"

cat(sprintf("Sham1: %d cells, Sham2: %d cells, Sham3: %d cells\n",
            ncol(sham1), ncol(sham2), ncol(sham3)))
cat(sprintf("MCAO1: %d cells, MCAO2: %d cells, MCAO3: %d cells\n",
            ncol(mcao1), ncol(mcao2), ncol(mcao3)))

cat("=== 合并样本 ===\n")
sc_obj <- merge(sham1, y = c(sham2, sham3, mcao1, mcao2, mcao3),
                add.cell.ids = c("sham1", "sham2", "sham3", "mcao1", "mcao2", "mcao3"))
sc_obj@meta.data$cell_id <- rownames(sc_obj@meta.data)
cat(sprintf("合并后总细胞数: %d\n", ncol(sc_obj)))

cat("=== 计算质控指标 ===\n")
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")
sc_obj[["percent.rb"]] <- PercentageFeatureSet(sc_obj, pattern = "^Rp[sl]")

cat("QC指标统计:\n")
print(summary(sc_obj@meta.data$nFeature_RNA))
print(summary(sc_obj@meta.data$nCount_RNA))
print(summary(sc_obj@meta.data$percent.mt))

p1 <- VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "condition", pt.size = 0.1) &
  theme(legend.position = "none")
ggsave(p1, filename = paste0(dir_save, "/QC_violin.png"), width = 10, height = 6)

cat("=== 质控过滤 ===\n")
sc_obj <- subset(sc_obj,
                 nFeature_RNA > 200 &
                 nFeature_RNA < 5000 &
                 percent.mt < 20)
cat(sprintf("QC后细胞数: %d\n", ncol(sc_obj)))

saveRDS(sc_obj, file = paste0(dir_save, "/sc_merged_QC.rds"))
cat("合并和QC后的数据已保存\n")

target_genes <- c("RAGE", "NFKB1", "FDX1", "AGER", "TLR4", "STAT1", "TGFB1", "NFE2L2",
                  "JAK1", "STAT3", "CCL2", "ICAM1", "HMOX1")
gene_check <- rownames(sc_obj)[toupper(rownames(sc_obj)) %in% toupper(target_genes)]
cat("目标基因检查:\n")
print(gene_check)
