# SCISSOR-like分析 v9 - 正确版 (同物种: 大鼠MCAO模型)
# Bulk: GSE163614 (大鼠MCAO Bulk)
# Single-cell: GSE174574 (大鼠MCAO scRNA-seq)
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"Seurat" %in% installed.packages()){install.packages('Seurat', repos='https://cran.rstudio.com')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
if(!"limma" %in% installed.packages()){BiocManager::install('limma')}
if(!"readxl" %in% installed.packages()){install.packages('readxl')}
library(Seurat)
library(ggplot2)
library(limma)
library(readxl)

cat("=== SCISSOR-like分析 v9 (正确版: 同物种大鼠MCAO) ===\n")
cat("Bulk: GSE163614 | Single-cell: GSE174574\n\n")

# 1. 读取Bulk数据 (GSE163614)
cat("1. 读取Bulk数据 (GSE163614 大鼠MCAO)...\n")
bulk_df <- read_excel('C:/Users/Jy-Mentor-7/Downloads/GSE163614_mRNA_Expression_Profiling.xlsx', skip=5)
colnames(bulk_df)[1] <- "gene_short_name"
colnames(bulk_df)[7:12] <- c("MCAO1","MCAO2","MCAO3","Sham1","Sham2","Sham3")
bulk_expr <- as.matrix(bulk_df[,7:12])
rownames(bulk_expr) <- bulk_df$gene_short_name
bulk_expr <- apply(bulk_expr, 2, as.numeric)
rownames(bulk_expr) <- bulk_df$gene_short_name
cat(sprintf("  Bulk数据: %d 基因 x %d 样本\n", nrow(bulk_expr), ncol(bulk_expr)))

bulk_pheno <- c(rep(1,3), rep(0,3))
names(bulk_pheno) <- colnames(bulk_expr)

# 2. 读取单细胞数据
cat("\n2. 读取单细胞数据 (GSE174574)...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
sc_obj <- JoinLayers(sc_obj)
cat(sprintf("  单细胞: %d 细胞 x %d 基因\n", ncol(sc_obj), nrow(sc_obj)))

# 3. 找共有基因 (无需物种转换，同为大鼠)
cat("\n3. 物种基因匹配...\n")
sc_genes <- rownames(sc_obj)
bulk_genes <- rownames(bulk_expr)
common_genes <- intersect(sc_genes, bulk_genes)
cat(sprintf("  单细胞基因: %d | Bulk基因: %d | 共有基因: %d\n",
            length(sc_genes), length(bulk_genes), length(common_genes)))

# 4. Bulk差异分析 (limma)
cat("\n4. Bulk差异分析 (limma - MCAO vs Sham)...\n")
bulk_c <- as.matrix(bulk_expr[common_genes, ])

design <- model.matrix(~0 + factor(bulk_pheno))
colnames(design) <- c("MCAO", "Sham")
contrast.matrix <- makeContrasts(MCAO - Sham, levels=design)

fit <- lmFit(bulk_c, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg_res <- topTable(fit2, adjust="BH", number=Inf)

log2fc <- deg_res$logFC
pvals <- deg_res$P.Value
fdr <- deg_res$adj.P.Val
names(log2fc) <- rownames(deg_res)

n_deg <- sum(fdr<0.05 & abs(log2fc)>0.5)
cat(sprintf("  显著差异基因 (FDR<0.05, |log2FC|>0.5): %d\n", n_deg))

# 5. 计算单细胞评分 (AddModuleScore)
cat("\n5. 计算单细胞表型评分 (AddModuleScore)...\n")

mcao_genes <- names(log2fc)[log2fc > 0.5 & fdr < 0.05]
sham_genes <- names(log2fc)[log2fc < -0.5 & fdr < 0.05]
cat(sprintf("  MCAO上调基因: %d\n", length(mcao_genes)))
cat(sprintf("  MCAO下调(保护)基因: %d\n", length(sham_genes)))

sc_obj <- AddModuleScore(
  sc_obj,
  features = list(mcao_genes, sham_genes),
  name = c("MCAO_Score", "Sham_Score")
)
colnames(sc_obj@meta.data)[colnames(sc_obj@meta.data) == "MCAO_Score1"] <- "MCAO_Score"
colnames(sc_obj@meta.data)[colnames(sc_obj@meta.data) == "Sham_Score2"] <- "Sham_Score"
sc_obj$Net_Score <- sc_obj$MCAO_Score - sc_obj$Sham_Score

# 6. 细胞类型分析
cat("\n6. 各细胞类型Net Score...\n")
score_sum <- aggregate(cbind(MCAO_Score, Sham_Score, Net_Score) ~ cell_type,
                       data=sc_obj@meta.data, FUN=mean)
score_sum <- score_sum[order(-score_sum$Net_Score), ]
print(score_sum)

# 7. BCP轴+铜死亡基因分析
cat("\n7. BCP轴与铜死亡核心基因...\n")
cuproptosis_genes <- c("SLC31A1","SLC31A2","ATP7A","ATP7B","ATOX1","LIAS","LIPT1","DLAT","PDHA1","PDHB","GCSH","DLD")
bcp_genes <- c("AGER","NFKB1","FDX1","TLR4","STAT1","STAT3","TGFB1","NFE2L2","SOD1","CAT")
all_target_genes <- c(cuproptosis_genes, bcp_genes)
target_found <- all_target_genes[all_target_genes %in% names(log2fc)]

if(length(target_found)>0) {
  target_df <- data.frame(
    gene=target_found,
    log2FC=log2fc[target_found],
    FDR=fdr[target_found],
    Category=ifelse(target_found %in% cuproptosis_genes, "Cuproptosis", "BCP_Axis")
  )
  rownames(target_df) <- NULL
  cat("\n  铜死亡与BCP轴基因差异:\n")
  print(target_df[order(-abs(target_df$log2FC)),])
}

# 8. 保存结果
cat("\n8. 保存结果...\n")
saveRDS(sc_obj, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_with_pheno_score_v9.rds')

dir.create('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v9', showWarnings=FALSE, recursive=TRUE)

p1 <- VlnPlot(sc_obj, features="Net_Score", group.by="cell_type", pt.size=0.1) +
  ggtitle("Net Score (MCAO - Sham) by Cell Type\nGSE163614 Bulk + GSE174574 sc") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p1, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v9/Net_Score.pdf', width=10, height=6)

p2 <- VlnPlot(sc_obj, features="MCAO_Score", group.by="cell_type", pt.size=0.1) +
  ggtitle("MCAO Association Score by Cell Type") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p2, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v9/MCAO_Score.pdf', width=10, height=6)

p3 <- VlnPlot(sc_obj, features="Sham_Score", group.by="cell_type", pt.size=0.1) +
  ggtitle("Sham (Protect) Score by Cell Type") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p3, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v9/Sham_Score.pdf', width=10, height=6)

deg_list <- list(mcao_genes=mcao_genes, sham_genes=sham_genes,
                 log2fc=log2fc, fdr=fdr, deg_res=deg_res)
saveRDS(deg_list, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v9/DEGs.rds')

cat("\n=== 完成 ===\n")
cat("结果保存到: C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v9/\n")