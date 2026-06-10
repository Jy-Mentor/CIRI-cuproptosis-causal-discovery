# 简化版SCISSOR分析 v3
# 解决物种问题：小鼠单细胞 vs 人类Bulk数据

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"openxlsx" %in% installed.packages()){install.packages('openxlsx')}

library(Seurat)
library(openxlsx)

cat("=== 简化版SCISSOR分析 (跨物种) ===\n")

# 1. 读取数据
cat("\n1. 读取数据...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
sc_obj <- JoinLayers(sc_obj)
cat(sprintf("  单细胞: %d 细胞\n", ncol(sc_obj)))

bulk_expr <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/expr_matrix_gene.rds')
bulk_pheno <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/phenotype.rds')
cat(sprintf("  Bulk: %d 基因 x %d 样本\n", nrow(bulk_expr), ncol(bulk_expr)))

# 2. 读取映射表
cat("\n2. 读取物种映射...\n")
mapping <- read.table('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt',
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mouse_to_human <- mapping[mapping$MOUSE_ORTHOLOG_SYMBOL != "" & mapping$HUMAN_ORTHOLOG_SYMBOL != "", ]
mouse_to_human <- mouse_to_human[!duplicated(mouse_to_human$MOUSE_ORTHOLOG_SYMBOL), ]
mouse_to_human <- mouse_to_human[, c("MOUSE_ORTHOLOG_SYMBOL", "HUMAN_ORTHOLOG_SYMBOL")]
colnames(mouse_to_human) <- c("Mouse", "Human")
cat(sprintf("  映射表: %d 对\n", nrow(mouse_to_human)))

# 3. 获取单细胞基因
cat("\n3. 基因转换...\n")
sc_genes <- rownames(sc_obj)
sc_data <- GetAssayData(sc_obj, layer = "data")

# 转换小鼠基因到人类
mouse_genes <- rownames(sc_data)
human_genes <- mouse_to_human$Human[match(mouse_genes, mouse_to_human$Mouse)]
names(human_genes) <- mouse_genes
human_genes <- human_genes[!is.na(human_genes) & !duplicated(human_genes)]

# 创建转换后的表达矩阵
sc_data_human <- sc_data[names(human_genes), ]
rownames(sc_data_human) <- human_genes[names(human_genes)]
sc_data_human <- sc_data_human[!duplicated(rownames(sc_data_human)), ]
cat(sprintf("  转换后单细胞基因: %d\n", nrow(sc_data_human)))

# 4. 共有基因
common_genes <- intersect(rownames(sc_data_human), rownames(bulk_expr))
cat(sprintf("  共有基因: %d\n", length(common_genes)))

# 5. Bulk差异基因
cat("\n4. Bulk差异分析...\n")
stroke_idx <- which(bulk_pheno == 1)
control_idx <- which(bulk_pheno == 0)
bulk_common <- bulk_expr[common_genes, ]
stroke_means <- rowMeans(bulk_common[, stroke_idx, drop=FALSE])
control_means <- rowMeans(bulk_common[, control_idx, drop=FALSE])
log2fc <- stroke_means - control_means

# t-test
pvals <- sapply(1:nrow(bulk_common), function(i) {
  x <- as.numeric(bulk_common[i, stroke_idx])
  y <- as.numeric(bulk_common[i, control_idx])
  t.test(x, y)$p.value
})
fdr <- p.adjust(pvals, method = "BH")

deg_df <- data.frame(gene=common_genes, log2FC=log2fc, p_value=pvals, FDR=fdr)
deg_sig <- deg_df[deg_df$FDR < 0.05 & abs(deg_df$log2FC) > 0.25, ]
cat(sprintf("  显著差异基因: %d\n", nrow(deg_sig)))

# 6. 计算细胞评分
cat("\n5. 计算细胞表型评分...\n")
phenotype_vec <- as.numeric(bulk_pheno)
names(phenotype_vec) <- colnames(bulk_expr)

# 基因-表型相关性
sc_common <- sc_data_human[common_genes, names(phenotype_vec)]
gene_cor <- sapply(1:nrow(sc_common), function(i) {
  cor(as.numeric(sc_common[i,]), phenotype_vec, method="spearman", use="complete.obs")
})
names(gene_cor) <- rownames(sc_common)

# 上下调基因
pos_genes <- names(gene_cor)[gene_cor > 0.05]
neg_genes <- names(gene_cor)[gene_cor < -0.05]
cat(sprintf("  正相关基因: %d, 负相关基因: %d\n", length(pos_genes), length(neg_genes)))

# 细胞评分
pos_score <- colSums(sc_data_human[pos_genes[pos_genes %in% rownames(sc_data_human)], ] *
                      gene_cor[pos_genes[pos_genes %in% names(gene_cor)]], na.rm = TRUE)
neg_score <- colSums(sc_data_human[neg_genes[neg_genes %in% rownames(sc_data_human)], ] *
                      abs(gene_cor[neg_genes[neg_genes %in% names(gene_cor)]]), na.rm = TRUE)

sc_obj$Stroke_Score <- as.numeric(scale(pos_score))
sc_obj$Protect_Score <- as.numeric(scale(neg_score))
sc_obj$Net_Score <- sc_obj$Stroke_Score - sc_obj$Protect_Score

# 7. 按细胞类型汇总
cat("\n6. 各细胞类型评分...\n")
score_summary <- aggregate(cbind(Stroke_Score, Protect_Score, Net_Score) ~ cell_type,
                          data = sc_obj@meta.data, FUN = mean)
score_summary <- score_summary[order(-score_summary$Net_Score), ]
print(score_summary)

# 8. BCP轴基因
cat("\n7. BCP轴基因...\n")
bcp_genes_human <- c("AGER", "NFKB1", "FDX1", "TLR4", "STAT1", "STAT3", "TGFB1", "NFE2L2")
bcp_found <- bcp_genes_human[bcp_genes_human %in% rownames(deg_df)]
cat(sprintf("  BCP轴基因在Bulk中: %s\n", paste(bcp_found, collapse=", ")))
bcp_deg <- deg_df[deg_df$gene %in% bcp_found, c("gene", "log2FC", "FDR")]
bcp_deg <- bcp_deg[order(-abs(bcp_deg$log2FC)), ]
print(bcp_deg)

# 9. 保存
cat("\n8. 保存结果...\n")
saveRDS(sc_obj, file = 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_with_pheno_score.rds')

dir_plot <- 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots'
if(!dir.exists(dir_plot)){dir.create(dir_plot)}

p1 <- VlnPlot(sc_obj, features = "Net_Score", group.by = "cell_type", pt.size = 0) +
  labs(title = "Net Phenotype Score by Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p1, file = paste0(dir_plot, "/Net_Score.png"), width = 10, height = 6)

cat("\n=== 完成 ===\n")
