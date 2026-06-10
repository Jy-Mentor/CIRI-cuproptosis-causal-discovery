# 简化版SCISSOR分析 v2
# 解决物种问题：小鼠单细胞 vs 人类Bulk数据
# 使用基因ID转换：mouse gene → human ortholog

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"openxlsx" %in% installed.packages()){install.packages('openxlsx')}

library(Seurat)
library(openxlsx)

cat("=== 简化版SCISSOR分析 (跨物种) ===\n")

# 1. 读取单细胞数据 (小鼠)
cat("\n1. 读取单细胞数据 (小鼠)...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
sc_obj <- JoinLayers(sc_obj)
cat(sprintf("  单细胞(小鼠): %d 细胞 x %d 基因\n", ncol(sc_obj), nrow(sc_obj)))

# 2. 读取Bulk表达矩阵 (人类)
cat("\n2. 读取Bulk表达矩阵 (人类)...\n")
bulk_expr <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/expr_matrix_gene.rds')
bulk_pheno <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/phenotype.rds')
cat(sprintf("  Bulk(人类): %d 基因 x %d 样本\n", nrow(bulk_expr), ncol(bulk_expr)))

# 3. 读取物种映射表
cat("\n3. 读取物种映射表...\n")
mapping_file <- 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt'
if(file.exists(mapping_file)) {
  mapping <- read.table(mapping_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat(sprintf("  映射表: %d 行\n", nrow(mapping)))
  head(mapping)
} else {
  cat("  映射文件不存在!\n")
  stop("请确保映射文件存在")
}

# 4. 创建小鼠→人类基因映射
cat("\n4. 创建小鼠→人类基因映射...\n")
# 映射表中小鼠基因在 MOUSE_ORTHOLOG_SYMBOL 列，人类基因在 HUMAN_ORTHOLOG_SYMBOL 列
mouse_to_human <- mapping[mapping$MOUSE_ORTHOLOG_SYMBOL != "" & mapping$HUMAN_ORTHOLOG_SYMBOL != "", ]
mouse_to_human <- mouse_to_human[!duplicated(mouse_to_human$MOUSE_ORTHOLOG_SYMBOL), ]
mouse_to_human <- mouse_to_human[, c("MOUSE_ORTHOLOG_SYMBOL", "HUMAN_ORTHOLOG_SYMBOL")]
colnames(mouse_to_human) <- c("Mouse", "Human")
cat(sprintf("  小鼠→人类映射基因数: %d\n", nrow(mouse_to_human)))

# 5. 将单细胞数据的小鼠基因转换为人类同源基因
cat("\n5. 转换单细胞基因名 (Mouse → Human)...\n")
sc_genes <- rownames(sc_obj)
names(sc_genes) <- sc_genes

# 查找小鼠基因在映射表中的对应
mouse_genes_in_mapping <- sc_genes[sc_genes %in% mouse_to_human$Mouse]
human_genes_for_sc <- mouse_to_human$Human[match(mouse_genes_in_mapping, mouse_to_human$Mouse)]
names(human_genes_for_sc) <- mouse_genes_in_mapping

cat(sprintf("  单细胞基因中可转换的: %d/%d\n", length(mouse_genes_in_mapping), length(sc_genes)))

# 6. 构建转换后的单细胞表达矩阵
cat("\n6. 构建转换后的表达矩阵...\n")
sc_data_raw <- GetAssayData(sc_obj, layer = "data")

# 只保留可以转换为人类基因的小鼠基因
convertible_genes <- names(mouse_genes_in_mapping)
sc_data_converted <- sc_data_raw[convertible_genes, ]
rownames(sc_data_converted) <- human_genes_for_sc[convertible_genes]

# 去除重复的人类基因名（同一基因多个小鼠基因对应）
sc_data_converted <- sc_data_converted[!duplicated(rownames(sc_data_converted)), ]

# 7. 获取Bulk数据的共有基因
common_genes <- intersect(rownames(sc_data_converted), rownames(bulk_expr))
cat(sprintf("  共有基因数: %d\n", length(common_genes)))

# 8. 计算Bulk数据中与卒中表型相关的基因
cat("\n7. 计算Bulk数据中与卒中表型相关的基因...\n")
stroke_samples <- names(bulk_pheno)[bulk_pheno == 1]
control_samples <- names(bulk_pheno)[bulk_pheno == 0]

bulk_common <- bulk_expr[common_genes, ]
bulk_common$Stroke_mean <- rowMeans(bulk_common[, stroke_samples])
bulk_common$Control_mean <- rowMeans(bulk_common[, control_samples])
bulk_common$log2FC <- bulk_common$Stroke_mean - bulk_common$Control_mean

# 简单t-test
compute_ttest_p <- function(expr_values, idx1, idx2) {
  x <- as.numeric(expr_values[idx1])
  y <- as.numeric(expr_values[idx2])
  if(sd(x) == 0 & sd(y) == 0) return(1)
  t.test(x, y)$p.value
}

cat("  计算差异...\n")
p_values <- sapply(1:nrow(bulk_common), function(i) {
  compute_ttest_p(bulk_common[i, 1:ncol(bulk_expr)], stroke_samples, control_samples)
})
bulk_common$p_value <- p_values
bulk_common$FDR <- p.adjust(bulk_common$p_value, method = "BH")
bulk_common$gene <- rownames(bulk_common)

# 9. 计算单细胞表型相关评分
cat("\n8. 计算单细胞表型相关评分...\n")

# 表型向量
phenotype_vector <- as.numeric(bulk_pheno)
names(phenotype_vector) <- names(bulk_pheno)

# 计算每个基因与表型的相关性
sc_data_common <- sc_data_converted[common_genes, names(phenotype_vector)]

gene_pheno_cor <- sapply(1:nrow(sc_data_common), function(i) {
  gene_expr <- as.numeric(sc_data_common[i, ])
  cor(gene_expr, phenotype_vector, method = "spearman", use = "complete.obs")
})
names(gene_pheno_cor) <- rownames(sc_data_common)

# 为每个细胞计算评分
positive_genes <- names(gene_pheno_cor)[gene_pheno_cor > 0.1]
negative_genes <- names(gene_pheno_cor)[gene_pheno_cor < -0.1]

cat(sprintf("  正相关基因数 (rho>0.1): %d\n", length(positive_genes)))
cat(sprintf("  负相关基因数 (rho<-0.1): %d\n", length(negative_genes)))

# 计算Stroke score
stroke_score <- colSums(sc_data_converted[positive_genes, ] * gene_pheno_cor[positive_genes], na.rm = TRUE)
protect_score <- colSums(sc_data_converted[negative_genes, ] * abs(gene_pheno_cor[negative_genes]), na.rm = TRUE)

sc_obj$Stroke_Score <- as.numeric(scale(stroke_score))
sc_obj$Protect_Score <- as.numeric(scale(protect_score))
sc_obj$Net_Score <- sc_obj$Stroke_Score - sc_obj$Protect_Score

# 10. 按细胞类型分析
cat("\n9. 各细胞类型的表型关联评分...\n")
score_summary <- aggregate(cbind(Stroke_Score, Protect_Score, Net_Score) ~ cell_type,
                          data = sc_obj@meta.data, FUN = mean)
score_summary <- score_summary[order(-score_summary$Net_Score), ]
print(score_summary)

# 11. BCP轴基因分析
cat("\n10. BCP轴基因在各细胞类型中的表达...\n")
bcp_genes_human <- c("AGER", "NFKB1", "FDX1", "TLR4", "STAT1", "STAT3", "TGFB1", "NFE2L2")
bcp_in_data <- bcp_genes_human[bcp_genes_human %in% rownames(sc_data_converted)]
cat(sprintf("  BCP轴基因在数据中: %s\n", paste(bcp_in_data, collapse = ", ")))

# 检查Bulk中这些基因的差异
cat("\n  BCP基因在Bulk(Stroke vs Control)中的差异:\n")
bcp_bulk <- bulk_common[bulk_common$gene %in% bcp_in_data, c("gene", "log2FC", "p_value", "FDR")]
bcp_bulk <- bcp_bulk[order(abs(bcp_bulk$log2FC), decreasing = TRUE), ]
print(bcp_bulk)

# 12. 保存结果
cat("\n11. 保存结果...\n")
saveRDS(sc_obj, file = 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_with_pheno_score.rds')

score_results <- list(
  celltype_scores = score_summary,
  bulk_DEG = bulk_common[, c("gene", "log2FC", "p_value", "FDR")],
  bcp_bulk = bcp_bulk
)
saveRDS(score_results, file = 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/phenotype_score_results.rds')

# 12. 可视化
cat("\n12. 生成可视化...\n")
dir_plot <- 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots'
if(!dir.exists(dir_plot)){dir.create(dir_plot)}

p1 <- VlnPlot(sc_obj, features = "Stroke_Score", group.by = "cell_type", pt.size = 0) +
  labs(title = "Stroke Association Score by Cell Type") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p1, file = paste0(dir_plot, "/Stroke_Score_by_celltype.png"), width = 10, height = 6)

p2 <- VlnPlot(sc_obj, features = "Net_Score", group.by = "cell_type", pt.size = 0) +
  labs(title = "Net Score (Stroke - Protective) by Cell Type") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p2, file = paste0(dir_plot, "/Net_Score_by_celltype.png"), width = 10, height = 6)

cat("\n=== 分析完成 ===\n")
cat("结果保存在:\n")
cat("  - sc_with_pheno_score.rds\n")
cat("  - phenotype_score_results.rds\n")
cat("  - scissor_plots/\n")
