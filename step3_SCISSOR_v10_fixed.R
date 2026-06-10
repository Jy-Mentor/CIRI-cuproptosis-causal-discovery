# SCISSOR-like分析 v10 - 修正版
# Bulk: GSE163614 (大鼠MCAO Bulk)
# Single-cell: GSE174574 (小鼠MCAO scRNA-seq)
# 修正: 1)基因名大小写 2)limma对比方向 3)跨物种映射检查
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

cat("=== SCISSOR-like分析 v10 (修正版) ===\n")
cat("Bulk: GSE163614 (Rat) | Single-cell: GSE174574 (Mouse)\n\n")

# 1. 读取Bulk数据 (GSE163614 大鼠)
cat("1. 读取Bulk数据 (GSE163614 大鼠MCAO)...\n")
bulk_df <- read_excel('C:/Users/Jy-Mentor-7/Downloads/GSE163614_mRNA_Expression_Profiling.xlsx', skip=5)
colnames(bulk_df)[1] <- "gene_short_name"
colnames(bulk_df)[7:12] <- c("MCAO1","MCAO2","MCAO3","Sham1","Sham2","Sham3")
bulk_expr <- as.matrix(bulk_df[,7:12])
rownames(bulk_expr) <- bulk_df$gene_short_name
bulk_expr <- apply(bulk_expr, 2, as.numeric)
rownames(bulk_expr) <- bulk_df$gene_short_name
cat(sprintf("  Bulk数据: %d 基因 x %d 样本\n", nrow(bulk_expr), ncol(bulk_expr)))
cat(sprintf("  样本名: %s\n", paste(colnames(bulk_expr), collapse=", ")))

# 2. 读取单细胞数据 (GSE174574 小鼠)
cat("\n2. 读取单细胞数据 (GSE174574 小鼠MCAO)...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
sc_obj <- JoinLayers(sc_obj)
cat(sprintf("  单细胞: %d 细胞 x %d 基因\n", ncol(sc_obj), nrow(sc_obj)))

# 3. 跨物种基因映射 (小鼠 -> 大鼠)
cat("\n3. 跨物种基因映射 (Mouse -> Rat)...\n")
mapping <- read.table('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt',
                     header=TRUE, sep="\t", stringsAsFactors=FALSE)
m2r <- mapping[mapping$MOUSE_ORTHOLOG_SYMBOL!="" & mapping$RAT_GENE_SYMBOL!="", ]
m2r <- m2r[!duplicated(m2r$MOUSE_ORTHOLOG_SYMBOL), c("MOUSE_ORTHOLOG_SYMBOL","RAT_GENE_SYMBOL")]
colnames(m2r) <- c("mouse","rat")
cat(sprintf("  映射库: %d 个小鼠->大鼠映射\n", nrow(m2r)))

sc_genes <- rownames(sc_obj)
mapped_genes <- m2r$rat[match(sc_genes, m2r$mouse)]
names(mapped_genes) <- sc_genes
valid_idx <- which(!is.na(mapped_genes) & !duplicated(mapped_genes))
sc_data_rat <- GetAssayData(sc_obj, layer="data")[valid_idx, ]
rownames(sc_data_rat) <- mapped_genes[valid_idx]
sc_data_rat <- sc_data_rat[!duplicated(rownames(sc_data_rat)), ]
cat(sprintf("  转换后单细胞基因: %d\n", nrow(sc_data_rat)))

# 4. 找共有基因
common_genes <- intersect(rownames(sc_data_rat), rownames(bulk_expr))
cat(sprintf("  共有基因(大鼠): %d\n", length(common_genes)))

# 5. Bulk差异分析 (limma - 修正对比方向)
cat("\n4. Bulk差异分析 (limma - MCAO vs Sham)...\n")
bulk_c <- as.matrix(bulk_expr[common_genes, ])

bulk_pheno <- c(rep("MCAO",3), rep("Sham",3))
design <- model.matrix(~0 + factor(bulk_pheno, levels=c("Sham","MCAO")))
colnames(design) <- c("Sham", "MCAO")
cat("  Design matrix:\n")
print(design)

contrast.matrix <- makeContrasts(MCAO - Sham, levels=design)
cat("  Contrast: MCAO - Sham\n")

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

# 6. 计算单细胞评分 (AddModuleScore)
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

# 7. 细胞类型分析
cat("\n6. 各细胞类型Net Score...\n")
score_sum <- aggregate(cbind(MCAO_Score, Sham_Score, Net_Score) ~ cell_type,
                       data=sc_obj@meta.data, FUN=mean)
score_sum <- score_sum[order(-score_sum$Net_Score), ]
print(score_sum)

# 8. BCP轴+铜死亡基因分析 (使用小鼠基因名)
cat("\n7. BCP轴与铜死亡核心基因 (小鼠基因名)...\n")
cuproptosis_genes <- c("Slc31a1","Slc31a2","Atp7a","Atp7b","Atox1","Lias","Lipt1","Dlat","Pdha1","Pdhb","Gcsh","Dld")
bcp_genes <- c("Ager","Nfkb1","Fdx1","Tlr4","Stat1","Stat3","Tgfbr1","Nfe2l2","Sod1","Cat")
all_target_genes <- c(cuproptosis_genes, bcp_genes)

mouse_in_sc <- rownames(sc_obj)[rownames(sc_obj) %in% all_target_genes]
cat(sprintf("  目标基因在单细胞中: %d/%d\n", length(mouse_in_sc), length(all_target_genes)))

mouse_found_df <- data.frame(
  gene=mouse_in_sc,
  in_sc=TRUE,
  Category=ifelse(mouse_in_sc %in% cuproptosis_genes, "Cuproptosis", "BCP_Axis")
)
print(mouse_found_df)

# 9. 保存结果
cat("\n8. 保存结果...\n")
saveRDS(sc_obj, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_with_pheno_score_v10.rds')

dir.create('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v10', showWarnings=FALSE, recursive=TRUE)

p1 <- VlnPlot(sc_obj, features="Net_Score", group.by="cell_type", pt.size=0) +
  ggtitle("Net Score (MCAO - Sham) by Cell Type\nv10: Mouse sc -> Rat Bulk") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(p1, file='C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v10/Net_Score.pdf', width=10, height=6)

deg_list <- list(mcao_genes=mcao_genes, sham_genes=sham_genes,
                 log2fc=log2fc, fdr=fdr, deg_res=deg_res)
saveRDS(deg_list, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v10/DEGs.rds')

cat("\n=== 完成 ===\n")
cat("结果保存到: C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v10/\n")