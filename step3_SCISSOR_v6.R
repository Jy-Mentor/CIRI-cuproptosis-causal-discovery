# 简化版SCISSOR分析 v6
# 核心逻辑：用Bulk数据的表型相关基因来评估单细胞
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
library(Seurat)
library(ggplot2)

cat("=== SCISSOR-like分析 v6 ===\n")

# 1. 读取数据
cat("\n1. 读取数据...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
sc_obj <- JoinLayers(sc_obj)
bulk_expr <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/expr_matrix_gene.rds')
bulk_pheno <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/phenotype.rds')

cat(sprintf("  单细胞: %d 细胞 x %d 基因\n", ncol(sc_obj), nrow(sc_obj)))

# 2. 读取映射并转换基因
cat("\n2. 物种基因映射...\n")
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
cat("\n3. Bulk差异分析...\n")
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
cat(sprintf("  显著差异基因: %d\n", sum(fdr<0.05 & abs(log2fc)>0.25)))

# 5. 计算单细胞评分
cat("\n4. 计算单细胞表型评分...\n")

# 取Bulk中与表型相关的基因
stroke_genes <- names(log2fc)[log2fc > 0.1 & fdr < 0.05]
protect_genes <- names(log2fc)[log2fc < -0.1 & fdr < 0.05]
cat(sprintf("  Stroke相关基因: %d\n", length(stroke_genes)))
cat(sprintf("  Protect相关基因: %d\n", length(protect_genes)))

# 单细胞在这些基因上的表达
sc_stroke <- sc_data_human[intersect(stroke_genes, rownames(sc_data_human)), , drop=FALSE]
sc_protect <- sc_data_human[intersect(protect_genes, rownames(sc_data_human)), , drop=FALSE]

# 细胞评分 = 基因表达量的加权和
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

# 7. BCP轴基因
cat("\n6. BCP轴基因Bulk差异...\n")
bcp <- c("AGER","NFKB1","FDX1","TLR4","STAT1","STAT3","TGFB1","NFE2L2")
bcp_found <- bcp[bcp %in% names(log2fc)]
if(length(bcp_found)>0) {
  bcp_df <- data.frame(gene=bcp_found, log2FC=log2fc[bcp_found], FDR=fdr[bcp_found])
  print(bcp_df[order(-abs(bcp_df$log2FC)),])
}

# 8. 保存
cat("\n7. 保存结果...\n")
saveRDS(sc_obj, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_with_pheno_score.rds')

dir.create('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots', showWarnings=FALSE)
p1 <- VlnPlot(sc_obj, features="Net_Score", group.by="cell_type", pt.size=0) +
  ggtitle("Net Phenotype Score (Stroke - Protective) by Cell Type") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p1, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots/Net_Score.png', width=10, height=6)

p2 <- VlnPlot(sc_obj, features="Stroke_Score", group.by="cell_type", pt.size=0) +
  ggtitle("Stroke Association Score by Cell Type") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p2, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots/Stroke_Score.png', width=10, height=6)

cat("\n=== 完成 ===\n")
