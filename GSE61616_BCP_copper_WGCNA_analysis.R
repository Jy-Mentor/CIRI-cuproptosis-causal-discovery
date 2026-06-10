# ============================================
# GSE61616 BCP靶点-铜死亡-WGCNA整合分析
# 7d MCAO vs Sham | 药物-疾病交集 | WGCNA模块 | PPI Hub基因
# ============================================

# ==================== 1. 包安装与加载 ====================
cat("正在检查和安装必要的R包...\n")

packages <- c("limma", "ggplot2", "dplyr", "tidyr", "GEOquery", 
              "preprocessCore", "ComplexHeatmap", "circlize", 
              "WGCNA", "igraph", "ggraph", "tidygraph", "readr")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Installing package: ", pkg, "\n"))
    if (pkg %in% c("GEOquery", "limma", "preprocessCore", "ComplexHeatmap", "circlize", "WGCNA")) {
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org/")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
    library(pkg, character.only = TRUE)
  }
}

for (pkg in packages) {
  tryCatch({
    install_if_missing(pkg)
  }, error = function(e) {
    cat(paste0("警告: 无法安装/加载包 ", pkg, ": ", e$message, "\n"))
  })
}

cat("包加载完成!\n\n")

# ==================== 2. 设置工作目录和文件路径 ====================
work_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
setwd(work_dir)

series_matrix_file <- "C:/Users/Jy-Mentor-7/Downloads/GSE61616_series_matrix (1).txt.gz"
platform_file <- "C:/Users/Jy-Mentor-7/Downloads/GPL1355-10794 (1).txt"
bcp_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/石竹烯 人.txt"
copper_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/铜死亡 15 特异性.txt"

stopifnot(file.exists(series_matrix_file))
stopifnot(file.exists(platform_file))
stopifnot(file.exists(bcp_file))
stopifnot(file.exists(copper_file))

output_dir <- file.path(work_dir, "GSE61616_BCP_copper_WGCNA_results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

set.seed(42)

# ==================== 3. 读取系列矩阵数据 ====================
cat("读取 GSE61616 系列矩阵数据...\n")

tryCatch({
  gse <- getGEO(filename = series_matrix_file, GSEMatrix = TRUE)
  if (is.list(gse)) {
    gse <- gse[[1]]
  }
  exprs <- exprs(gse)
  cat(paste0("表达矩阵维度: ", nrow(exprs), " probes x ", ncol(exprs), " samples\n"))
}, error = function(e) {
  stop(paste0("读取GEO数据失败: ", e$message))
})

stopifnot(nrow(exprs) > 0, ncol(exprs) > 0)

# ==================== 4. 定义样本分组（7d MCAO vs Sham） ====================
cat("定义样本分组信息...\n")

# GSE61616: 7d MCAO (n=3) vs Sham (n=3)
sample_info <- data.frame(
  geo_accession = colnames(exprs),
  group = NA_character_,
  stringsAsFactors = FALSE
)

# 从phenoData提取分组信息
pd <- pData(phenoData(gse))
cat("样本特征表列名:\n")
print(names(pd))

# 尝试自动识别分组
title_col <- grep("title|characteristics|source", names(pd), ignore.case = TRUE, value = TRUE)
if (length(title_col) > 0) {
  cat("\n样本标题/特征:\n")
  print(pd[, title_col[1], drop = FALSE])
}

# GSE61616实际样本: Sham_1-5, Model_1-5, XST_1-5
# 根据实际数据: Sham = Sham_1-5, MCAO = Model_1-5
mcaO_samples <- grep("Model", pd[[title_col[1]]], ignore.case = TRUE)
sham_samples <- grep("Sham", pd[[title_col[1]]], ignore.case = TRUE)

if (length(mcaO_samples) == 0 || length(sham_samples) == 0) {
  cat("无法自动识别分组，使用硬编码样本ID...\n")
  # 根据实际样本名称
  mcaO_ids <- c("GSM1509427", "GSM1509428", "GSM1509429", "GSM1509430", "GSM1509431")
  sham_ids <- c("GSM1509422", "GSM1509423", "GSM1509424", "GSM1509425", "GSM1509426")
  
  sample_info$group <- ifelse(sample_info$geo_accession %in% mcaO_ids, "MCAO_7d",
                              ifelse(sample_info$geo_accession %in% sham_ids, "Sham_7d", NA))
} else {
  sample_info$group[mcaO_samples] <- "MCAO_7d"
  sample_info$group[sham_samples] <- "Sham_7d"
}

sample_info <- sample_info[!is.na(sample_info$group), ]
cat("样本分组情况:\n")
print(table(sample_info$group))

stopifnot(nrow(sample_info) > 0)

# ==================== 5. 提取表达矩阵子集 ====================
exprs_subset <- exprs[, sample_info$geo_accession, drop = FALSE]
cat(paste0("子集表达矩阵维度: ", nrow(exprs_subset), " probes x ", ncol(exprs_subset), " samples\n\n"))

# ==================== 6. 数据标准化 ====================
cat("检查表达值范围...\n")
data_range <- range(exprs_subset, na.rm = TRUE)
cat(paste0("表达值范围: ", round(data_range[1], 2), " - ", round(data_range[2], 2), "\n"))

if (data_range[2] > 100) {
  cat("执行log2转换...\n")
  exprs_subset <- log2(exprs_subset + 1)
} else {
  cat("数据已经是log2转换后的值\n")
}

cat("执行quantile标准化...\n")
library(preprocessCore)
exprs_norm <- normalize.quantiles(exprs_subset)
rownames(exprs_norm) <- rownames(exprs_subset)
colnames(exprs_norm) <- colnames(exprs_subset)

# ==================== 7. 读取平台注释信息 ====================
cat("读取 GPL1355 平台注释...\n")

plat_lines <- readLines(platform_file, warn = FALSE)
header_line <- which(grepl("^ID", plat_lines))[1]
stopifnot(length(header_line) > 0, !is.na(header_line))

platform_annot <- read.delim(platform_file, 
                             skip = header_line - 1, 
                             header = TRUE,
                             check.names = FALSE,
                             stringsAsFactors = FALSE,
                             quote = "",
                             fill = TRUE)

cat(paste0("平台注释探针数: ", nrow(platform_annot), "\n"))

annot_table <- data.frame(
  ID = platform_annot$ID,
  Gene_Symbol = if("Gene Symbol" %in% names(platform_annot)) platform_annot[["Gene Symbol"]] else NA,
  Gene_Title = if("Gene Title" %in% names(platform_annot)) platform_annot[["Gene Title"]] else NA,
  Entrez_ID = if("ENTREZ_GENE_ID" %in% names(platform_annot)) platform_annot[["ENTREZ_GENE_ID"]] else NA,
  stringsAsFactors = FALSE
)

# ==================== 8. limma差异表达分析（7d MCAO vs Sham） ====================
cat("\n开始进行 limma 差异表达分析 (7d MCAO vs Sham)...\n")

sample_info$group <- factor(sample_info$group, levels = c("Sham_7d", "MCAO_7d"))
design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- gsub("group", "", colnames(design))

cat("设计矩阵:\n")
print(design)

fit <- lmFit(exprs_norm, design)
contrast.matrix <- makeContrasts(MCAO_7d_vs_Sham_7d = MCAO_7d - Sham_7d, levels = design)

cat("\n对比矩阵:\n")
print(contrast.matrix)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 提取所有基因的差异表达结果
all_results <- topTable(fit2, coef = "MCAO_7d_vs_Sham_7d", number = Inf, adjust.method = "fdr", sort.by = "p")
all_results$ProbeID <- rownames(all_results)
all_results$Gene <- annot_table$Gene_Symbol[match(all_results$ProbeID, annot_table$ID)]

cat(paste0("差异表达分析完成，共 ", nrow(all_results), " 个探针\n"))

# 筛选显著差异基因 |logFC| > 1, FDR < 0.05
sig_threshold <- 0.05
fc_threshold <- 1
sig_degs <- all_results[!is.na(all_results$logFC) & abs(all_results$logFC) > fc_threshold & 
                        !is.na(all_results$adj.P.Val) & all_results$adj.P.Val < sig_threshold, ]
cat(paste0("显著差异基因 (|logFC| > 1, FDR < 0.05): ", nrow(sig_degs), "\n"))

# 保存差异表达结果
write.table(all_results, file = file.path(output_dir, "GSE61616_all_DEG_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig_degs, file = file.path(output_dir, "GSE61616_sig_DEGs_logFC1_FDR0.05.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# ==================== 9. 读取BCP靶点和铜死亡基因 ====================
cat("\n读取BCP靶点基因和铜死亡核心基因...\n")

bcp_targets <- readLines(bcp_file, warn = FALSE)
bcp_targets <- trimws(bcp_targets)
bcp_targets <- bcp_targets[bcp_targets != ""]
bcp_targets <- unique(bcp_targets)
cat(paste0("BCP靶点基因数: ", length(bcp_targets), "\n"))

copper_genes <- readLines(copper_file, warn = FALSE)
copper_genes <- trimws(copper_genes)
copper_genes <- copper_genes[copper_genes != ""]
copper_genes <- unique(copper_genes)
cat(paste0("铜死亡核心基因数: ", length(copper_genes), "\n"))

# ==================== 10. 药物-疾病交集基因 ====================
cat("\n计算药物-疾病交集基因...\n")

# BCP靶点与差异基因取交集
sig_gene_symbols <- unique(sig_degs$Gene)
sig_gene_symbols <- sig_gene_symbols[!is.na(sig_gene_symbols)]

drug_disease_genes <- intersect(bcp_targets, sig_gene_symbols)
cat(paste0("药物-疾病交集基因 (BCP靶点 ∩ 差异基因): ", length(drug_disease_genes), "\n"))
print(drug_disease_genes)

# 保存交集基因
write.table(drug_disease_genes, file = file.path(output_dir, "drug_disease_intersection_genes.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==================== 11. 扩展基因集：交集基因 ∪ 铜死亡核心基因 ====================
cat("\n构建扩展基因集...\n")

extended_genes <- union(drug_disease_genes, copper_genes)
cat(paste0("扩展基因集 (药物-疾病交集 ∪ 铜死亡核心): ", length(extended_genes), "\n"))

write.table(extended_genes, file = file.path(output_dir, "extended_gene_set.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==================== 12. WGCNA分析 ====================
cat("\n开始 WGCNA 分析...\n")

# 准备表达矩阵（探针→基因转换，取每个基因表达量最高的探针）
exprs_df <- as.data.frame(exprs_norm)
exprs_df$ProbeID <- rownames(exprs_df)
exprs_df$Gene <- annot_table$Gene_Symbol[match(exprs_df$ProbeID, annot_table$ID)]
exprs_df <- exprs_df[!is.na(exprs_df$Gene) & exprs_df$Gene != "", ]

# 对每个基因保留方差最大的探针
exprs_df$var <- apply(exprs_df[, 1:(ncol(exprs_df)-2)], 1, var)
exprs_df <- exprs_df[order(exprs_df$Gene, -exprs_df$var), ]
exprs_df <- exprs_df[!duplicated(exprs_df$Gene), ]

gene_exprs <- exprs_df[, c("Gene", sample_info$geo_accession)]
rownames(gene_exprs) <- gene_exprs$Gene
gene_exprs$Gene <- NULL
gene_exprs <- as.matrix(gene_exprs)

cat(paste0("基因级表达矩阵维度: ", nrow(gene_exprs), " genes x ", ncol(gene_exprs), " samples\n"))

# WGCNA数据预处理
WGCNA::allowWGCNAThreads()

datExpr0 <- t(gene_exprs)

gsg <- WGCNA::goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# 选择软阈值
powers <- c(1:20)
sft <- WGCNA::pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# 绘制软阈值选择图
pdf(file.path(output_dir, "WGCNA_soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

# 使用推荐的软阈值
softPower <- sft$powerEstimate
if (is.na(softPower) || softPower < 1) {
  softPower <- 6
  cat(paste0("无法自动确定软阈值，使用默认值: ", softPower, "\n"))
} else {
  cat(paste0("选择的软阈值: ", softPower, "\n"))
}

# 构建邻接矩阵和TOM
adjacency <- WGCNA::adjacency(datExpr0, power = softPower)
TOM <- WGCNA::TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# 层次聚类
geneTree <- hclust(as.dist(dissTOM), method = "average")

# 动态树切割识别模块
minModuleSize <- 30
library(dynamicTreeCut)
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                    deepSplit = 2, pamRespectsDendro = FALSE,
                                    minClusterSize = minModuleSize)
dynamicColors <- WGCNA::labels2colors(dynamicMods)

# 模块合并
MEList <- WGCNA::moduleEigengenes(datExpr0, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

mergeCutHeight <- 0.25
merge <- WGCNA::mergeCloseModules(datExpr0, dynamicColors, cutHeight = mergeCutHeight, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

cat(paste0("WGCNA模块数: ", length(unique(mergedColors)), "\n"))
cat("模块颜色:\n")
print(table(mergedColors))

# 保存模块分配
module_df <- data.frame(Gene = colnames(datExpr0), Module = mergedColors)
write.table(module_df, file = file.path(output_dir, "WGCNA_gene_module_assignment.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 绘制模块聚类图
pdf(file.path(output_dir, "WGCNA_module_clustering.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
WGCNA::plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                           c("Dynamic Tree Cut", "Merged dynamic"),
                           dendroLabels = FALSE, hang = 0.03,
                           addGuide = TRUE, guideHang = 0.05)
dev.off()

# ==================== 13. 模块与性状关联 ====================
cat("\n计算模块与MCAO性状的关联...\n")

# 定义性状矩阵（MCAO=1, Sham=0）
traitData <- data.frame(MCAO = ifelse(sample_info$group == "MCAO_7d", 1, 0))
rownames(traitData) <- sample_info$geo_accession

# 计算模块-性状关联
moduleTraitCor <- cor(mergedMEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr0))

# 保存模块-性状关联
cor_df <- data.frame(Module = rownames(moduleTraitCor), 
                     Correlation = moduleTraitCor[, 1], 
                     Pvalue = moduleTraitPvalue[, 1])
write.table(cor_df, file = file.path(output_dir, "WGCNA_module_trait_correlation.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 绘制模块-性状关联热图
pdf(file.path(output_dir, "WGCNA_module_trait_heatmap.pdf"), width = 6, height = 8)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                      xLabels = names(traitData),
                      yLabels = rownames(moduleTraitCor),
                      ySymbols = rownames(moduleTraitCor),
                      colorLabels = FALSE,
                      colors = WGCNA::blueWhiteRed(50),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      cex.text = 0.5,
                      zlim = c(-1, 1),
                      main = "Module-trait relationships (GSE61616 MCAO 7d)")
dev.off()

# 识别与MCAO最相关的模块
best_module_idx <- which.max(abs(moduleTraitCor))
best_module <- rownames(moduleTraitCor)[best_module_idx]
best_module_color <- gsub("ME", "", best_module)
best_cor <- moduleTraitCor[best_module_idx, 1]
cat(paste0("与MCAO最相关的模块: ", best_module, " (颜色: ", best_module_color, ", r=", round(best_cor, 3), ")\n"))

# ==================== 14. 模块富集验证：BCP靶点 + 铜死亡基因覆盖度 ====================
cat("\n进行模块富集验证...\n")

# 获取最佳模块中的基因（使用模块颜色名称匹配）
best_module_genes <- module_df$Gene[module_df$Module == best_module_color]
cat(paste0(best_module, " (", best_module_color, ") 模块基因数: ", length(best_module_genes), "\n"))

# 计算BCP靶点在模块中的富集
bcp_in_module <- intersect(bcp_targets, best_module_genes)
bcp_coverage <- length(bcp_in_module) / length(bcp_targets)
cat(paste0("BCP靶点在", best_module, "模块中的覆盖度: ", length(bcp_in_module), "/", length(bcp_targets), 
           " (", round(bcp_coverage * 100, 1), "%)\n"))

# 计算铜死亡基因在模块中的富集
copper_in_module <- intersect(copper_genes, best_module_genes)
copper_coverage <- length(copper_in_module) / length(copper_genes)
cat(paste0("铜死亡基因在", best_module, "模块中的覆盖度: ", length(copper_in_module), "/", length(copper_genes), 
           " (", round(copper_coverage * 100, 1), "%)\n"))

# 计算扩展基因集在模块中的覆盖
extended_in_module <- intersect(extended_genes, best_module_genes)
extended_coverage <- length(extended_in_module) / length(extended_genes)
cat(paste0("扩展基因集在", best_module, "模块中的覆盖度: ", length(extended_in_module), "/", length(extended_genes), 
           " (", round(extended_coverage * 100, 1), "%)\n"))

# 超几何检验
module_enrichment_test <- function(target_genes, module_genes, all_genes) {
  x <- length(intersect(target_genes, module_genes))
  m <- length(target_genes)
  n <- length(all_genes) - m
  k <- length(module_genes)
  pval <- phyper(x - 1, m, n, k, lower.tail = FALSE)
  return(list(hits = x, expected = k * m / length(all_genes), pvalue = pval))
}

all_gene_symbols <- rownames(gene_exprs)

bcp_enrich <- module_enrichment_test(bcp_targets, best_module_genes, all_gene_symbols)
copper_enrich <- module_enrichment_test(copper_genes, best_module_genes, all_gene_symbols)
extended_enrich <- module_enrichment_test(extended_genes, best_module_genes, all_gene_symbols)

enrich_df <- data.frame(
  GeneSet = c("BCP_targets", "Copper_death", "Extended"),
  Module = best_module,
  Hits = c(bcp_enrich$hits, copper_enrich$hits, extended_enrich$hits),
  Total = c(length(bcp_targets), length(copper_genes), length(extended_genes)),
  Coverage = c(bcp_coverage, copper_coverage, extended_coverage),
  Expected = c(bcp_enrich$expected, copper_enrich$expected, extended_enrich$expected),
  Pvalue = c(bcp_enrich$pvalue, copper_enrich$pvalue, extended_enrich$pvalue)
)
write.table(enrich_df, file = file.path(output_dir, "module_enrichment_validation.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n模块富集验证结果:\n")
print(enrich_df)

# 保存模块中的关键基因
write.table(bcp_in_module, file = file.path(output_dir, paste0("BCP_targets_in_", best_module, "_module.txt")), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(copper_in_module, file = file.path(output_dir, paste0("Copper_genes_in_", best_module, "_module.txt")), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==================== 15. 提取最佳模块基因表达矩阵用于PPI ====================
cat("\n提取最佳模块基因用于PPI网络分析...\n")

module_expr <- gene_exprs[best_module_genes, ]
module_expr_df <- as.data.frame(module_expr)
module_expr_df$Gene <- rownames(module_expr_df)

# 计算模块内基因与MCAO性状的关联（GS）和模块隶属度（MM）
ME_best <- mergedMEs[, best_module]
GS <- cor(t(module_expr), traitData$MCAO, use = "p")
MM <- cor(t(module_expr), ME_best, use = "p")

module_hub_df <- data.frame(
  Gene = best_module_genes,
  GS = as.numeric(GS),
  GS_pvalue = corPvalueStudent(as.numeric(GS), ncol(module_expr)),
  MM = as.numeric(MM),
  MM_pvalue = corPvalueStudent(as.numeric(MM), ncol(module_expr))
)
module_hub_df <- module_hub_df[order(-abs(module_hub_df$MM)), ]
write.table(module_hub_df, file = file.path(output_dir, paste0(best_module, "_module_hub_genes.txt")), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# ==================== 16. 准备String PPI网络输入 ====================
cat("\n准备String PPI网络输入文件...\n")

# 筛选高连接度基因作为String输入（模块内top基因 + 扩展基因集交集）
module_top_genes <- module_hub_df$Gene[1:min(200, nrow(module_hub_df))]
ppi_input_genes <- union(module_top_genes, extended_in_module)
ppi_input_genes <- intersect(ppi_input_genes, best_module_genes)

cat(paste0("PPI网络输入基因数: ", length(ppi_input_genes), "\n"))

# 保存String输入
write.table(ppi_input_genes, file = file.path(output_dir, "String_PPI_input_genes.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# 同时保存扩展基因集的String输入
write.table(extended_genes, file = file.path(output_dir, "String_PPI_extended_genes_input.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==================== 17. 可视化：热图、火山图、模块图 ====================
cat("\n生成可视化图表...\n")

# 扩展基因集表达热图
extended_expr <- gene_exprs[intersect(extended_genes, rownames(gene_exprs)), ]
if (nrow(extended_expr) > 0) {
  sample_order <- order(sample_info$group)
  extended_expr_ordered <- extended_expr[, sample_order, drop = FALSE]
  extended_z <- t(scale(t(extended_expr_ordered)))
  
  group_colors <- c("MCAO_7d" = "#E74C3C", "Sham_7d" = "#2ECC71")
  ha <- HeatmapAnnotation(
    Group = sample_info$group[sample_order],
    col = list(Group = group_colors)
  )
  
  ht <- Heatmap(extended_z,
                name = "Z-score",
                top_annotation = ha,
                cluster_columns = FALSE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 8),
                column_title = "GSE61616 Extended Gene Set Expression (MCAO 7d vs Sham)")
  
  pdf(file.path(output_dir, "extended_gene_set_heatmap.pdf"), width = 10, height = 12)
  print(ht)
  dev.off()
  
  png(file.path(output_dir, "extended_gene_set_heatmap.png"), width = 1000, height = 1200, res = 150)
  print(ht)
  dev.off()
}

# 火山图（所有基因，标记扩展基因集）
all_results$logP <- -log10(all_results$P.Value)
all_results$color <- "Other"
all_results$color[all_results$Gene %in% extended_genes & all_results$adj.P.Val < 0.05 & abs(all_results$logFC) > 1] <- "Extended_Sig"
all_results$color[all_results$Gene %in% copper_genes & all_results$adj.P.Val < 0.05 & abs(all_results$logFC) > 1] <- "Copper_Sig"
all_results$color[all_results$Gene %in% drug_disease_genes & all_results$adj.P.Val < 0.05 & abs(all_results$logFC) > 1] <- "DrugDisease_Sig"

colors_volcano <- c("Other" = "grey80", "Extended_Sig" = "#9B59B6", 
                    "Copper_Sig" = "#E74C3C", "DrugDisease_Sig" = "#3498DB")

p_volcano <- ggplot(all_results, aes(x = logFC, y = logP, color = color)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = colors_volcano, name = "Gene Set") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title = "GSE61616: MCAO 7d vs Sham (Extended Gene Sets)",
       x = expression(paste("log"[2], " Fold Change")),
       y = expression(paste("-log"[10], " P-value"))) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "volcano_extended_genes.png"), plot = p_volcano, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "volcano_extended_genes.pdf"), plot = p_volcano, width = 8, height = 6)

# 模块基因GS-MM散点图
p_gsmm <- ggplot(module_hub_df, aes(x = GS, y = MM)) +
  geom_point(aes(color = abs(GS) > 0.5 & abs(MM) > 0.8), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "grey60"), name = "Hub-like") +
  geom_hline(yintercept = c(-0.8, 0.8), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  labs(title = paste0(best_module, " Module: Gene Significance vs Module Membership"),
       x = "Gene Significance (correlation with MCAO)",
       y = "Module Membership") +
  theme_bw(base_size = 12)

ggsave(file.path(output_dir, "GS_MM_scatter.png"), plot = p_gsmm, width = 8, height = 6, dpi = 300)

# ==================== 18. 总结报告 ====================
cat("\n========================================\n")
cat("         分析完成总结报告\n")
cat("========================================\n")
cat(paste0("数据集: GSE61616 (7d MCAO vs Sham)\n"))
cat(paste0("平台: GPL1355\n"))
cat(paste0("样本: MCAO 7d (n=3), Sham 7d (n=3)\n\n"))

cat(paste0("差异表达基因 (|logFC|>1, FDR<0.05): ", nrow(sig_degs), "\n"))
cat(paste0("药物-疾病交集基因 (BCP ∩ DEGs): ", length(drug_disease_genes), "\n"))
cat(paste0("扩展基因集 (药物-疾病 ∪ 铜死亡): ", length(extended_genes), "\n\n"))

cat(paste0("WGCNA最佳模块: ", best_module, "\n"))
cat(paste0("  模块基因数: ", length(best_module_genes), "\n"))
cat(paste0("  与MCAO相关性: r=", round(best_cor, 3), "\n"))
cat(paste0("  BCP靶点覆盖: ", length(bcp_in_module), "/", length(bcp_targets), 
           " (", round(bcp_coverage*100, 1), "%)\n"))
cat(paste0("  铜死亡基因覆盖: ", length(copper_in_module), "/", length(copper_genes), 
           " (", round(copper_coverage*100, 1), "%)\n\n"))

cat(paste0("String PPI输入基因: ", length(ppi_input_genes), "\n"))
cat("  请将这些基因提交至 String DB (https://string-db.org) 进行PPI分析\n")
cat("  然后在Cytoscape中使用CytoHubba插件计算MCC得分筛选Hub基因\n\n")

cat(paste0("输出目录: ", output_dir, "\n"))
cat("生成文件列表:\n")
for (f in list.files(output_dir)) {
  cat(paste0("  - ", f, "\n"))
}
cat("========================================\n")
