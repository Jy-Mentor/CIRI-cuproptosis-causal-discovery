# 简化版SCISSOR分析
# 使用相关性方法识别与卒中表型相关的细胞类型
# 不依赖SCISSOR包，直接使用单细胞和Bulk数据的共有基因

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"openxlsx" %in% installed.packages()){install.packages('openxlsx')}

library(Seurat)
library(openxlsx)

cat("=== 简化版SCISSOR分析 ===\n")

# 1. 读取单细胞数据
cat("\n1. 读取单细胞数据...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
sc_obj <- JoinLayers(sc_obj)
cat(sprintf("  单细胞: %d 细胞 x %d 基因\n", ncol(sc_obj), nrow(sc_obj)))
cat("  细胞类型分布:\n")
print(table(sc_obj$cell_type))

# 2. 读取Bulk表达矩阵
cat("\n2. 读取Bulk表达矩阵...\n")
bulk_expr <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/expr_matrix_gene.rds')
bulk_pheno <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/phenotype.rds')
cat(sprintf("  Bulk: %d 基因 x %d 样本\n", nrow(bulk_expr), ncol(bulk_expr)))
cat(sprintf("  表型: Control=%d, Stroke=%d\n", sum(bulk_pheno==0), sum(bulk_pheno==1)))

# 3. 匹配基因
cat("\n3. 匹配单细胞和Bulk数据的基因...\n")
sc_genes <- rownames(sc_obj)
bulk_genes <- rownames(bulk_expr)
common_genes <- intersect(sc_genes, bulk_genes)
cat(sprintf("  共有基因数: %d\n", length(common_genes)))

# 4. 计算Bulk数据的表型相关基因
cat("\n4. 计算Bulk数据中与卒中表型相关的基因...\n")
stroke_samples <- names(bulk_pheno)[bulk_pheno == 1]
control_samples <- names(bulk_pheno)[bulk_pheno == 0]

bulk_expr_df <- as.data.frame(bulk_expr[common_genes, ])
bulk_expr_df$Stroke_mean <- rowMeans(bulk_expr_df[, stroke_samples])
bulk_expr_df$Control_mean <- rowMeans(bulk_expr_df[, control_samples])
bulk_expr_df$log2FC <- bulk_expr_df$Stroke_mean - bulk_expr_df$Control_mean
bulk_expr_df$gene <- rownames(bulk_expr_df)

# 计算p-value (简单的t-test)
compute_ttest_p <- function(gene_expr, group1_samples, group2_samples) {
  x <- as.numeric(gene_expr[group1_samples])
  y <- as.numeric(gene_expr[group2_samples])
  if(sd(x) == 0 & sd(y) == 0) return(1)
  t.test(x, y)$p.value
}

cat("  计算差异表达...\n")
p_values <- sapply(1:nrow(bulk_expr_df), function(i) {
  compute_ttest_p(bulk_expr_df[i, ], stroke_samples, control_samples)
})
bulk_expr_df$p_value <- p_values
bulk_expr_df$FDR <- p.adjust(bulk_expr_df$p_value, method = "BH")

# 筛选显著基因
sig_genes <- bulk_expr_df[bulk_expr_df$FDR < 0.05 & abs(bulk_expr_df$log2FC) > 0.25, ]
cat(sprintf("  显著差异基因数 (FDR<0.05, |log2FC|>0.25): %d\n", nrow(sig_genes)))

# 5. 对单细胞数据计算表型相关评分
cat("\n5. 计算单细胞表型相关评分...\n")

# 方法：计算每个细胞的表达与Bulk表型相关基因的相关性
# 为每个细胞计算一个"Stroke association score"

# 获取单细胞的归一化数据
sc_data <- as.matrix(GetAssayData(sc_obj, layer = "data")[common_genes, ])

# 计算Bulk表型相关向量
phenotype_vector <- as.numeric(bulk_pheno)
names(phenotype_vector) <- names(bulk_pheno)

# 计算每个基因与表型的相关性
gene_pheno_cor <- sapply(1:nrow(sc_data), function(i) {
  gene_expr <- as.numeric(sc_data[i, names(phenotype_vector)])
  cor(gene_expr, phenotype_vector, method = "spearman", use = "complete.obs")
})
names(gene_pheno_cor) <- rownames(sc_data)

# 为每个细胞计算加权评分
cat("  计算细胞评分...\n")
stroke_associated_genes <- names(gene_pheno_cor)[gene_pheno_cor > 0]
protec_associated_genes <- names(gene_pheno_cor)[gene_pheno_cor < 0]

# 计算每个细胞的stroke score
stroke_score <- colSums(sc_data[stroke_associated_genes, ] * gene_pheno_cor[stroke_associated_genes], na.rm = TRUE)
protect_score <- colSums(sc_data[protec_associated_genes, ] * abs(gene_pheno_cor[protec_associated_genes]), na.rm = TRUE)

sc_obj$Stroke_Score <- scale(stroke_score)
sc_obj$Protect_Score <- scale(protect_score)
sc_obj$Net_Score <- sc_obj$Stroke_Score - sc_obj$Protect_Score

# 6. 按细胞类型分析
cat("\n6. 各细胞类型的表型关联评分...\n")
score_by_celltype <- aggregate(cbind(Stroke_Score, Protect_Score, Net_Score) ~ cell_type,
                              data = sc_obj@meta.data, FUN = mean)
score_by_celltype <- score_by_celltype[order(-score_by_celltype$Net_Score), ]
print(score_by_celltype)

# 7. 识别与BCP轴相关的细胞
cat("\n7. BCP轴基因在各细胞类型中的表达...\n")
bcp_genes <- c("Ager", "Nfkb1", "Fdx1", "Tlr4", "Stat1", "Stat3", "Tgfb1", "Nfe2l2")
bcp_genes_found <- bcp_genes[bcp_genes %in% rownames(sc_obj)]

for(g in bcp_genes_found) {
  cat(sprintf("\n  %s:\n", g))
  expr_by_ct <- aggregate(as.formula(paste0(g, " ~ cell_type")),
                          data = sc_obj@meta.data, FUN = mean)
  expr_by_ct <- expr_by_ct[order(-expr_by_ct[, 2]), ]
  for(i in 1:min(3, nrow(expr_by_ct))) {
    cat(sprintf("    %s: %.3f\n", expr_by_ct$cell_type[i], expr_by_ct[i, 2]))
  }
}

# 8. 保存结果
cat("\n8. 保存结果...\n")
saveRDS(sc_obj, file = 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_with_pheno_score.rds')

score_summary <- list(
  celltype_scores = score_by_celltype,
  bulk_DEG = sig_genes[, c("gene", "log2FC", "p_value", "FDR")],
  gene_pheno_cor = gene_pheno_cor,
  bcp_genes = bcp_genes_found
)
saveRDS(score_summary, file = 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/phenotype_association_summary.rds')

# 9. 可视化
cat("\n9. 生成可视化...\n")
dir_plot <- 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots'
if(!dir.exists(dir_plot)){dir.create(dir_plot)}

p1 <- VlnPlot(sc_obj, features = "Stroke_Score", group.by = "cell_type", pt.size = 0) +
  labs(title = "Stroke Association Score by Cell Type") +
  theme(legend.position = "none")
ggsave(p1, file = paste0(dir_plot, "/Stroke_Score_by_celltype.png"), width = 10, height = 6)

p2 <- VlnPlot(sc_obj, features = "Net_Score", group.by = "cell_type", pt.size = 0) +
  labs(title = "Net Phenotype Score (Stroke - Protective) by Cell Type") +
  theme(legend.position = "none")
ggsave(p2, file = paste0(dir_plot, "/Net_Score_by_celltype.png"), width = 10, height = 6)

# 热图：BCP基因在各细胞类型的表达
bcp_expr <- GetAssayData(sc_obj, layer = "data")[bcp_genes_found, ]
bcp_expr_ct <- aggregate(t(as.matrix(bcp_expr)), by = list(cell_type = sc_obj$cell_type), mean)
bcp_expr_matrix <- as.matrix(bcp_expr_ct[, -1])
rownames(bcp_expr_matrix) <- bcp_expr_ct$cell_type
bcp_expr_matrix <- t(scale(bcp_expr_matrix))

p3 <- pheatmap::pheatmap(bcp_expr_matrix,
                         main = "BCP Axis Genes Expression by Cell Type",
                         color = colorRampPalette(c("blue", "white", "red"))(100))
ggsave(plot = p3, file = paste0(dir_plot, "/BCP_genes_heatmap.png"), width = 8, height = 6)

cat("\n=== 分析完成 ===\n")
cat("结果保存在:\n")
cat("  - sc_with_pheno_score.rds: 带评分的Seurat对象\n")
cat("  - phenotype_association_summary.rds: 分析摘要\n")
cat("  - scissor_plots/: 可视化图\n")
