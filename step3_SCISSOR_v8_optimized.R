# SCISSOR-like分析 v8 - 优化版 (limma + AddModuleScore)
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
if(!"limma" %in% installed.packages()){BiocManager::install('limma')}
library(Seurat)
library(ggplot2)
library(limma)

cat("=== SCISSOR-like分析 v8 (优化版) ===\n")

# 1. 读取数据
cat("\n1. 读取数据...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
sc_obj <- JoinLayers(sc_obj)

bulk_file <- 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/expr_matrix_gene.rds'
bulk_expr <- readRDS(bulk_file)
bulk_pheno <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/phenotype.rds')

cat(sprintf("  单细胞: %d 细胞 x %d 基因\n", ncol(sc_obj), nrow(sc_obj)))
cat(sprintf("  Bulk数据: %d 样本 (Stroke=%d, Control=%d)\n",
            length(bulk_pheno), sum(bulk_pheno==1), sum(bulk_pheno==0)))

# 2. 读取映射并转换基因
cat("\n2. 物种基因映射 (Mouse -> Human)...\n")
mapping <- read.table('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt',
                     header=TRUE, sep="\t", stringsAsFactors=FALSE)
m2h <- mapping[mapping$MOUSE_ORTHOLOG_SYMBOL!="" & mapping$HUMAN_ORTHOLOG_SYMBOL!="", ]
m2h <- m2h[!duplicated(m2h$MOUSE_ORTHOLOG_SYMBOL), c("MOUSE_ORTHOLOG_SYMBOL","HUMAN_ORTHOLOG_SYMBOL")]
colnames(m2h) <- c("mouse","human")

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

# 4. Bulk差异分析 (limma)
cat("\n3. Bulk差异分析 (limma)...\n")
bulk_c <- as.matrix(bulk_expr[common, ])
stroke_idx <- which(bulk_pheno==1)
ctrl_idx <- which(bulk_pheno==0)

design <- model.matrix(~0 + factor(bulk_pheno))
colnames(design) <- c("Control", "Stroke")
contrast.matrix <- makeContrasts(Stroke - Control, levels=design)

fit <- lmFit(bulk_c, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg_res <- topTable(fit2, adjust="BH", number=Inf)

log2fc <- deg_res$logFC
pvals <- deg_res$P.Value
fdr <- deg_res$adj.P.Val
names(log2fc) <- rownames(deg_res)

n_deg <- sum(fdr<0.05 & abs(log2fc)>0.25)
cat(sprintf("  显著差异基因 (FDR<0.05, |log2FC|>0.25): %d\n", n_deg))

# 5. 计算单细胞评分 (AddModuleScore)
cat("\n4. 计算单细胞表型评分 (AddModuleScore)...\n")

stroke_genes <- names(log2fc)[log2fc > 0.1 & fdr < 0.05]
protect_genes <- names(log2fc)[log2fc < -0.1 & fdr < 0.05]
cat(sprintf("  Stroke上调基因: %d\n", length(stroke_genes)))
cat(sprintf("  Stroke下调(保护)基因: %d\n", length(protect_genes)))

sc_obj <- AddModuleScore(
  sc_obj,
  features = list(stroke_genes, protect_genes),
  name = c("Stroke_Score", "Protect_Score")
)
colnames(sc_obj@meta.data)[colnames(sc_obj@meta.data) == "Stroke_Score1"] <- "Stroke_Score"
colnames(sc_obj@meta.data)[colnames(sc_obj@meta.data) == "Protect_Score2"] <- "Protect_Score"
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
}

# 8. 保存结果
cat("\n7. 保存结果...\n")
saveRDS(sc_obj, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_with_pheno_score_v8.rds')

dir.create('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v8', showWarnings=FALSE, recursive=TRUE)

p1 <- VlnPlot(sc_obj, features="Net_Score", group.by="cell_type", pt.size=0) +
  ggtitle("Net Score (Stroke - Protect) by Cell Type\n(limma + AddModuleScore)") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p1, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v8/Net_Score.pdf', width=10, height=6)

p2 <- VlnPlot(sc_obj, features="Stroke_Score", group.by="cell_type", pt.size=0) +
  ggtitle("Stroke Association Score by Cell Type") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p2, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v8/Stroke_Score.pdf', width=10, height=6)

p3 <- VlnPlot(sc_obj, features="Protect_Score", group.by="cell_type", pt.size=0) +
  ggtitle("Protect Score by Cell Type") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p3, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v8/Protect_Score.pdf', width=10, height=6)

deg_list <- list(stroke_genes=stroke_genes, protect_genes=protect_genes,
                 log2fc=log2fc, fdr=fdr, deg_res=deg_res)
saveRDS(deg_list, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v8/DEGs.rds')

cat("\n=== 完成 ===\n")
cat("结果保存到: C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v8/\n")