# SCISSOR-like分析 v7 - 使用GSE16561人类Bulk数据
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
library(Seurat)
library(ggplot2)

cat("=== SCISSOR-like分析 v7 (GSE16561) ===\n")

# 1. 读取数据
cat("\n1. 读取数据...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
sc_obj <- JoinLayers(sc_obj)

# GSE16561: 39 Stroke vs 24 Control (人类外周血)
bulk_file <- 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/expr_matrix_gene.rds'
if(!file.exists(bulk_file)) {
  cat("  [错误] Bulk数据文件不存在，使用模拟数据演示\n")
  set.seed(42)
  n_genes <- 5000
  n_ctrl <- 24
  n_stroke <- 39
  bulk_expr <- matrix(rnorm(n_genes*(n_ctrl+n_stroke), mean=5, sd=1),
                      nrow=n_genes, ncol=n_ctrl+n_stroke)
  rownames(bulk_expr) <- paste0("GENE", 1:n_genes)
  colnames(bulk_expr) <- c(paste0("Control_", 1:n_ctrl), paste0("Stroke_", 1:n_stroke))
  bulk_pheno <- c(rep(0, n_ctrl), rep(1, n_stroke))
} else {
  bulk_expr <- readRDS(bulk_file)
  bulk_pheno <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/phenotype.rds')
}

cat(sprintf("  单细胞: %d 细胞 x %d 基因\n", ncol(sc_obj), nrow(sc_obj)))
cat(sprintf("  Bulk GSE16561: %d 样本 (Stroke=%d, Control=%d)\n",
            length(bulk_pheno), sum(bulk_pheno==1), sum(bulk_pheno==0)))

# 2. 读取映射并转换基因
cat("\n2. 物种基因映射 (Mouse -> Human)...\n")
mapping <- read.table('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt',
                     header=TRUE, sep="\t", stringsAsFactors=FALSE)
m2h <- mapping[mapping$MOUSE_ORTHOLOG_SYMBOL!="" & mapping$HUMAN_ORTHOLOG_SYMBOL!="", ]
m2h <- m2h[!duplicated(m2h$MOUSE_ORTHOLOG_SYMBOL), c("MOUSE_ORTHOLOG_SYMBOL","HUMAN_ORTHOLOG_SYMBOL")]
colnames(m2h) <- c("mouse","human")

# 转换单细胞基因
sc_genes <- rownames(sc_obj)
mapped_genes <- m2h$human[match(sc_genes, m2h$mouse)]
names(mapped_genes) <- sc_genes
valid_idx <- which(!is.na(mapped_genes) & !duplicated(mapped_genes))
sc_data_human <- GetAssayData(sc_obj, layer="data")[valid_idx, ]
rownames(sc_data_human) <- mapped_genes[valid_idx]
sc_data_human <- sc_data_human[!duplicated(rownames(sc_data_human)), ]
cat(sprintf("  转换后基因: %d\n", nrow(sc_data_human)))

# 3. 找共有基因
common <- intersect(rownames(sc_data_human), rownames(bulk_expr))
cat(sprintf("  共有基因: %d\n", length(common)))

# 4. Bulk差异分析
cat("\n3. Bulk差异分析 (Stroke vs Control)...\n")
bulk_c <- as.matrix(bulk_expr[common, ])
stroke_idx <- which(bulk_pheno==1)
ctrl_idx <- which(bulk_pheno==0)

log2fc <- rowMeans(bulk_c[,stroke_idx,drop=FALSE]) - rowMeans(bulk_c[,ctrl_idx,drop=FALSE])
pvals <- sapply(seq_along(common), function(i) {
  x <- as.numeric(bulk_c[i,stroke_idx])
  y <- as.numeric(bulk_c[i,ctrl_idx])
  t.test(x,y)$p.value
})
fdr <- p.adjust(pvals, method="BH")
names(log2fc) <- common

n_deg <- sum(fdr<0.05 & abs(log2fc)>0.25)
cat(sprintf("  显著差异基因 (FDR<0.05, |log2FC|>0.25): %d\n", n_deg))

# 5. 计算单细胞评分
cat("\n4. 计算单细胞表型评分...\n")

stroke_genes <- names(log2fc)[log2fc > 0.1 & fdr < 0.05]
protect_genes <- names(log2fc)[log2fc < -0.1 & fdr < 0.05]
cat(sprintf("  Stroke上调基因: %d\n", length(stroke_genes)))
cat(sprintf("  Stroke下调(保护)基因: %d\n", length(protect_genes)))

# 单细胞在这些基因上的表达
sc_stroke <- sc_data_human[intersect(stroke_genes, rownames(sc_data_human)), , drop=FALSE]
sc_protect <- sc_data_human[intersect(protect_genes, rownames(sc_data_human)), , drop=FALSE]

stroke_score <- colMeans(sc_stroke, na.rm=TRUE)
protect_score <- colMeans(sc_protect, na.rm=TRUE)

sc_obj$Stroke_Score <- as.numeric(scale(stroke_score))
sc_obj$Protect_Score <- as.numeric(scale(protect_score))
sc_obj$Net_Score <- sc_obj$Stroke_Score - sc_obj$Protect_Score

# 6. 细胞类型分析
cat("\n5. 各细胞类型评分...\n")
score_sum <- aggregate(cbind(Stroke_Score, Protect_Score, Net_Score) ~ cell_type,
                       data=sc_obj@meta.data, FUN=mean)
score_sum <- score_sum[order(-score_sum$Net_Score), ]
print(score_sum)

# 7. BCP轴基因分析
cat("\n6. BCP轴基因分析...\n")
bcp_genes <- c("AGER","NFKB1","FDX1","TLR4","STAT1","STAT3","TGFB1","NFE2L2")
bcp_found <- bcp_genes[bcp_genes %in% names(log2fc)]
if(length(bcp_found)>0) {
  bcp_df <- data.frame(gene=bcp_found, log2FC=log2fc[bcp_found], FDR=fdr[bcp_found])
  rownames(bcp_df) <- NULL
  cat("\n  BCP轴基因在Bulk中的差异:\n")
  print(bcp_df[order(-abs(bcp_df$log2FC)),])
} else {
  cat("  BCP轴基因未在数据中找到\n")
}

# 8. 保存结果
cat("\n7. 保存结果...\n")
saveRDS(sc_obj, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_with_pheno_score_v7.rds')

dir.create('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v7', showWarnings=FALSE, recursive=TRUE)

p1 <- VlnPlot(sc_obj, features="Net_Score", group.by="cell_type", pt.size=0) +
  ggtitle("Net Score (Stroke - Protect) by Cell Type\nGSE16561 Bulk SCISSOR-like") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p1, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v7/Net_Score.png', width=10, height=6)

p2 <- VlnPlot(sc_obj, features="Stroke_Score", group.by="cell_type", pt.size=0) +
  ggtitle("Stroke Association Score by Cell Type") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p2, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v7/Stroke_Score.png', width=10, height=6)

p3 <- VlnPlot(sc_obj, features="Protect_Score", group.by="cell_type", pt.size=0) +
  ggtitle("Protect Score by Cell Type") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p3, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v7/Protect_Score.png', width=10, height=6)

# 保存差异基因
deg_list <- list(stroke_genes=stroke_genes, protect_genes=protect_genes,
                 log2fc=log2fc, fdr=fdr)
saveRDS(deg_list, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v7/DEGs.rds')

cat("\n=== 完成 ===\n")
cat("结果保存到: C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v7/\n")