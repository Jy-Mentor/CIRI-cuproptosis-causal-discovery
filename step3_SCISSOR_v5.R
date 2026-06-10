# 简化版SCISSOR分析 v5
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
library(Seurat)

cat("=== SCISSOR-like分析 ===\n")

# 1. 读取数据
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
sc_obj <- JoinLayers(sc_obj)
bulk_expr <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/expr_matrix_gene.rds')
bulk_pheno <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/phenotype.rds')

cat(sprintf("单细胞: %d细胞, Bulk: %d基因x%d样本\n", ncol(sc_obj), nrow(bulk_expr), ncol(bulk_expr)))
cat(sprintf("Bulk样本名: %s\n", paste(head(colnames(bulk_expr),3), collapse=",")))
cat(sprintf("Bulk表型名: %s\n", paste(head(names(bulk_pheno),3), collapse=",")))

# 2. 读取映射
mapping <- read.table('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt',
                     header=TRUE, sep="\t", stringsAsFactors=FALSE)
m2h <- mapping[mapping$MOUSE_ORTHOLOG_SYMBOL!="" & mapping$HUMAN_ORTHOLOG_SYMBOL!="", ]
m2h <- m2h[!duplicated(m2h$MOUSE_ORTHOLOG_SYMBOL), c("MOUSE_ORTHOLOG_SYMBOL","HUMAN_ORTHOLOG_SYMBOL")]
colnames(m2h) <- c("mouse","human")
cat(sprintf("映射: %d对\n", nrow(m2h)))

# 3. 转换基因为人类
sc_genes <- rownames(sc_obj)
sc_data <- as.matrix(GetAssayData(sc_obj, layer="data"))
mapped_genes <- m2h$human[match(sc_genes, m2h$mouse)]
names(mapped_genes) <- sc_genes
valid_idx <- which(!is.na(mapped_genes) & !duplicated(mapped_genes))
sc_data_human <- sc_data[valid_idx, ]
rownames(sc_data_human) <- mapped_genes[valid_idx]
sc_data_human <- sc_data_human[!duplicated(rownames(sc_data_human)), ]
cat(sprintf("转换后基因: %d\n", nrow(sc_data_human)))

# 4. 共有基因
common <- intersect(rownames(sc_data_human), rownames(bulk_expr))
cat(sprintf("共有基因: %d\n", length(common)))

# 5. Bulk差异
stroke_idx <- which(bulk_pheno==1)
ctrl_idx <- which(bulk_pheno==0)
bulk_c <- as.matrix(bulk_expr[common,,drop=FALSE])
log2fc <- rowMeans(bulk_c[,stroke_idx,drop=FALSE]) - rowMeans(bulk_c[,ctrl_idx,drop=FALSE])
pvals <- sapply(seq_along(common), function(i) {
  x <- as.numeric(bulk_c[i,stroke_idx])
  y <- as.numeric(bulk_c[i,ctrl_idx])
  t.test(x,y)$p.value
})
fdr <- p.adjust(pvals, method="BH")
names(log2fc) <- common
cat(sprintf("Bulk显著基因(FDR<0.05,|log2FC|>0.25): %d\n", sum(fdr<0.05 & abs(log2fc)>0.25)))

# 6. 计算细胞评分
phenovec <- as.numeric(bulk_pheno)
names(phenovec) <- colnames(bulk_expr)

# 只保留Bulk样本名匹配的列
sc_cols <- intersect(colnames(sc_data_human), names(phenovec))
cat(sprintf("单细胞样本匹配: %d/%d\n", length(sc_cols), ncol(sc_data_human)))

sc_c <- sc_data_human[, sc_cols, drop=FALSE]
phenovec_matched <- phenovec[sc_cols]

# 基因-表型相关性
gene_cor <- sapply(seq_len(nrow(sc_c)), function(i) {
  cor(as.numeric(sc_c[i,]), phenovec_matched, method="spearman", use="complete.obs")
})
names(gene_cor) <- rownames(sc_c)

# 正负基因
pos_g <- names(gene_cor)[gene_cor > 0.05]
neg_g <- names(gene_cor)[gene_cor < -0.05]
cat(sprintf("正相关基因: %d, 负相关基因: %d\n", length(pos_g), length(neg_g)))

# 细胞评分
stroke_score <- colSums(sc_data_human[intersect(pos_g,rownames(sc_data_human)),,drop=FALSE] *
                        gene_cor[intersect(pos_g,names(gene_cor))], na.rm=TRUE)
protect_score <- colSums(sc_data_human[intersect(neg_g,rownames(sc_data_human)),,drop=FALSE] *
                        abs(gene_cor[intersect(neg_g,names(gene_cor))]), na.rm=TRUE)
sc_obj$Stroke_Score <- as.numeric(scale(stroke_score))
sc_obj$Protect_Score <- as.numeric(scale(protect_score))
sc_obj$Net_Score <- sc_obj$Stroke_Score - sc_obj$Protect_Score

# 7. 细胞类型评分
cat("\n各细胞类型Net Score:\n")
score_sum <- aggregate(Net_Score ~ cell_type, data=sc_obj@meta.data, FUN=mean)
print(score_sum[order(-score_sum$Net_Score),])

# 8. BCP轴
cat("\nBCP轴基因Bulk差异:\n")
bcp <- c("AGER","NFKB1","FDX1","TLR4","STAT1","STAT3","TGFB1","NFE2L2")
bcp_found <- bcp[bcp %in% names(log2fc)]
if(length(bcp_found)>0) {
  bcp_df <- data.frame(gene=bcp_found, log2FC=log2fc[bcp_found], FDR=fdr[bcp_found])
  print(bcp_df[order(-abs(bcp_df$log2FC)),])
}

# 9. 保存
saveRDS(sc_obj, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_with_pheno_score.rds')
dir.create('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots', showWarnings=FALSE)
p1 <- VlnPlot(sc_obj, features="Net_Score", group.by="cell_type", pt.size=0) +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p1, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots/Net_Score.png', width=10, height=6)
cat("\n=== 完成 ===\n")
