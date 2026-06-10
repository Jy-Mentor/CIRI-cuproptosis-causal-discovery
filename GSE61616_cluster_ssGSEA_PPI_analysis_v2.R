# ============================================
# GSE61616 层次聚类 + ssGSEA + PPI Hub基因分析 (v2)
# 使用本地大鼠-人类映射库进行物种转换
# 7d MCAO vs Sham | 动态剪切聚类 | 通路验证 | CytoHubba MCC
# ============================================

# ==================== 1. 包安装与加载 ====================
cat("正在检查和安装必要的R包...\n")

packages <- c("limma", "ggplot2", "dplyr", "tidyr", "GEOquery", 
              "preprocessCore", "ComplexHeatmap", "circlize", 
              "dynamicTreeCut", "igraph")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Installing package: ", pkg, "\n"))
    if (pkg %in% c("GEOquery", "limma", "preprocessCore", "ComplexHeatmap", "circlize")) {
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
mapping_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt"

stopifnot(file.exists(series_matrix_file))
stopifnot(file.exists(platform_file))
stopifnot(file.exists(bcp_file))
stopifnot(file.exists(copper_file))
stopifnot(file.exists(mapping_file))

output_dir <- file.path(work_dir, "GSE61616_cluster_ssGSEA_PPI_results_v2")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

set.seed(42)

# ==================== 3. 读取本地大鼠-人类映射库 ====================
cat("读取本地大鼠-人类映射库...\n")

# 读取映射文件（跳过注释行）
mapping_lines <- readLines(mapping_file, warn = FALSE)
header_line <- which(grepl("^RAT_GENE_SYMBOL", mapping_lines))[1]
stopifnot(length(header_line) > 0, !is.na(header_line))

mapping_df <- read.delim(mapping_file, 
                         skip = header_line - 1, 
                         header = TRUE,
                         check.names = FALSE,
                         stringsAsFactors = FALSE,
                         quote = "",
                         fill = TRUE)

cat(paste0("映射库总条目数: ", nrow(mapping_df), "\n"))

# 建立人类→大鼠基因映射（处理多映射情况）
human_to_rat <- list()
rat_to_human <- list()

for (i in 1:nrow(mapping_df)) {
  rat_gene <- toupper(trimws(mapping_df$RAT_GENE_SYMBOL[i]))
  human_ortholog <- toupper(trimws(mapping_df$HUMAN_ORTHOLOG_SYMBOL[i]))
  
  # 跳过空值
  if (rat_gene == "" || human_ortholog == "") next
  
  # 处理多个人类同源基因（用|分隔）
  human_genes <- unlist(strsplit(human_ortholog, "\\|"))
  human_genes <- trimws(human_genes)
  human_genes <- human_genes[human_genes != ""]
  
  for (hg in human_genes) {
    hg <- toupper(hg)
    if (!(hg %in% names(human_to_rat))) {
      human_to_rat[[hg]] <- c()
    }
    human_to_rat[[hg]] <- unique(c(human_to_rat[[hg]], rat_gene))
    
    if (!(rat_gene %in% names(rat_to_human))) {
      rat_to_human[[rat_gene]] <- c()
    }
    rat_to_human[[rat_gene]] <- unique(c(rat_to_human[[rat_gene]], hg))
  }
}

cat(paste0("人类→大鼠映射基因数: ", length(human_to_rat), "\n"))
cat(paste0("大鼠→人类映射基因数: ", length(rat_to_human), "\n"))

# ==================== 4. 读取BCP靶点和铜死亡基因并转换 ====================
cat("\n读取BCP靶点基因和铜死亡核心基因...\n")

bcp_targets_human <- readLines(bcp_file, warn = FALSE)
bcp_targets_human <- trimws(bcp_targets_human)
bcp_targets_human <- bcp_targets_human[bcp_targets_human != ""]
bcp_targets_human <- unique(toupper(bcp_targets_human))
cat(paste0("BCP靶点基因数（人类）: ", length(bcp_targets_human), "\n"))

copper_genes_human <- readLines(copper_file, warn = FALSE)
copper_genes_human <- trimws(copper_genes_human)
copper_genes_human <- copper_genes_human[copper_genes_human != ""]
copper_genes_human <- unique(toupper(copper_genes_human))
cat(paste0("铜死亡核心基因数（人类）: ", length(copper_genes_human), "\n"))

# 转换为同源大鼠基因
bcp_targets_rat <- c()
for (gene in bcp_targets_human) {
  if (gene %in% names(human_to_rat)) {
    bcp_targets_rat <- c(bcp_targets_rat, human_to_rat[[gene]])
  }
}
bcp_targets_rat <- unique(bcp_targets_rat)
cat(paste0("BCP靶点基因数（映射到大鼠）: ", length(bcp_targets_rat), "\n"))

copper_genes_rat <- c()
for (gene in copper_genes_human) {
  if (gene %in% names(human_to_rat)) {
    copper_genes_rat <- c(copper_genes_rat, human_to_rat[[gene]])
  }
}
copper_genes_rat <- unique(copper_genes_rat)
cat(paste0("铜死亡核心基因数（映射到大鼠）: ", length(copper_genes_rat), "\n"))

# 保存未映射的基因
unmapped_bcp <- setdiff(bcp_targets_human, names(human_to_rat))
unmapped_copper <- setdiff(copper_genes_human, names(human_to_rat))
cat(paste0("未映射的BCP靶点: ", length(unmapped_bcp), "\n"))
cat(paste0("未映射的铜死亡基因: ", length(unmapped_copper), "\n"))

write.table(unmapped_bcp, file = file.path(output_dir, "unmapped_BCP_targets.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unmapped_copper, file = file.path(output_dir, "unmapped_copper_genes.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==================== 5. 读取系列矩阵数据 ====================
cat("\n读取 GSE61616 系列矩阵数据...\n")

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

# ==================== 6. 定义样本分组 ====================
cat("定义样本分组信息...\n")

sample_info <- data.frame(
  geo_accession = colnames(exprs),
  group = NA_character_,
  stringsAsFactors = FALSE
)

pd <- pData(phenoData(gse))
title_col <- grep("title|characteristics|source", names(pd), ignore.case = TRUE, value = TRUE)
if (length(title_col) > 0) {
  cat("\n样本标题/特征:\n")
  print(pd[, title_col[1], drop = FALSE])
}

mcaO_samples <- grep("Model", pd[[title_col[1]]], ignore.case = TRUE)
sham_samples <- grep("Sham", pd[[title_col[1]]], ignore.case = TRUE)

if (length(mcaO_samples) == 0 || length(sham_samples) == 0) {
  cat("无法自动识别分组，使用硬编码样本ID...\n")
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

# ==================== 7. 提取表达矩阵子集并标准化 ====================
exprs_subset <- exprs[, sample_info$geo_accession, drop = FALSE]
cat(paste0("子集表达矩阵维度: ", nrow(exprs_subset), " probes x ", ncol(exprs_subset), " samples\n\n"))

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

# ==================== 8. 读取平台注释信息 ====================
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

# ==================== 9. limma差异表达分析 ====================
cat("\n开始进行 limma 差异表达分析 (7d MCAO vs Sham)...\n")

sample_info$group <- factor(sample_info$group, levels = c("Sham_7d", "MCAO_7d"))
design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- gsub("group", "", colnames(design))

fit <- lmFit(exprs_norm, design)
contrast.matrix <- makeContrasts(MCAO_7d_vs_Sham_7d = MCAO_7d - Sham_7d, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

all_results <- topTable(fit2, coef = "MCAO_7d_vs_Sham_7d", number = Inf, adjust.method = "fdr", sort.by = "p")
all_results$ProbeID <- rownames(all_results)
all_results$Gene <- annot_table$Gene_Symbol[match(all_results$ProbeID, annot_table$ID)]

cat(paste0("差异表达分析完成，共 ", nrow(all_results), " 个探针\n"))

sig_threshold <- 0.05
fc_threshold <- 1
sig_degs <- all_results[!is.na(all_results$logFC) & abs(all_results$logFC) > fc_threshold & 
                        !is.na(all_results$adj.P.Val) & all_results$adj.P.Val < sig_threshold, ]
cat(paste0("显著差异基因 (|logFC| > 1, FDR < 0.05): ", nrow(sig_degs), "\n"))

write.table(all_results, file = file.path(output_dir, "GSE61616_all_DEG_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig_degs, file = file.path(output_dir, "GSE61616_sig_DEGs_logFC1_FDR0.05.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# ==================== 10. 探针→基因转换 ====================
cat("\n探针→基因转换...\n")

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

# 基因名统一转为大写
rownames(gene_exprs) <- toupper(rownames(gene_exprs))

cat(paste0("基因级表达矩阵维度: ", nrow(gene_exprs), " genes x ", ncol(gene_exprs), " samples\n"))

# ==================== 11. 层次聚类 + 动态剪切 ====================
cat("\n进行层次聚类和动态剪切...\n")

datExpr <- t(gene_exprs)

# 计算基因间相关性距离
gene_cor <- cor(datExpr, use = "pairwise.complete.obs")
gene_dist <- as.dist(1 - gene_cor)

# 层次聚类
geneTree <- hclust(gene_dist, method = "average")

# 动态树切割
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = as.matrix(gene_dist),
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 50)

# 自定义颜色映射函数
labels2colors_custom <- function(labels) {
  color_palette <- c("turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue", "lightcyan", "grey60", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue")
  unique_labels <- sort(unique(labels))
  colors <- rep("grey", length(labels))
  for (i in seq_along(unique_labels)) {
    if (unique_labels[i] == 0) {
      colors[labels == unique_labels[i]] <- "grey"
    } else {
      color_idx <- ((unique_labels[i] - 1) %% length(color_palette)) + 1
      colors[labels == unique_labels[i]] <- color_palette[color_idx]
    }
  }
  return(colors)
}

cluster_colors <- labels2colors_custom(dynamicMods)

n_clusters <- length(unique(cluster_colors[cluster_colors != "grey"]))
cat(paste0("识别到 ", n_clusters, " 个基因簇（不含grey）\n"))

# 保存聚类结果
cluster_df <- data.frame(Gene = colnames(datExpr), Cluster = dynamicMods, Color = cluster_colors)
write.table(cluster_df, file = file.path(output_dir, "gene_cluster_assignment.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# 绘制聚类树
pdf(file.path(output_dir, "gene_clustering_dendrogram.pdf"), width = 12, height = 6)
plot(geneTree, xlab = "", sub = "", main = "Gene Clustering (Dynamic Tree Cut)", labels = FALSE, hang = 0.04)
dev.off()

# ==================== 12. 筛选MCAO高表达簇 ====================
cat("\n筛选MCAO高表达簇...\n")

cluster_stats <- data.frame()
for (cl in unique(cluster_colors)) {
  if (cl == "grey") next
  
  cl_genes <- cluster_df$Gene[cluster_df$Color == cl]
  cl_expr <- gene_exprs[cl_genes, , drop = FALSE]
  
  cl_mean_expr <- colMeans(cl_expr, na.rm = TRUE)
  
  mcaO_expr <- cl_mean_expr[sample_info$group == "MCAO_7d"]
  sham_expr <- cl_mean_expr[sample_info$group == "Sham_7d"]
  
  t_res <- t.test(mcaO_expr, sham_expr)
  
  cluster_stats <- rbind(cluster_stats, data.frame(
    Cluster = cl,
    N_genes = length(cl_genes),
    MCAO_mean = mean(mcaO_expr),
    Sham_mean = mean(sham_expr),
    logFC = log2(mean(mcaO_expr) / mean(sham_expr)),
    t_statistic = t_res$statistic,
    pvalue = t_res$p.value
  ))
}

cluster_stats$FDR <- p.adjust(cluster_stats$pvalue, method = "fdr")
cluster_stats <- cluster_stats[order(-cluster_stats$logFC), ]
write.table(cluster_stats, file = file.path(output_dir, "cluster_expression_statistics.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# 筛选MCAO高表达簇（logFC > 0 且 p < 0.05）
high_clusters <- cluster_stats$Cluster[cluster_stats$logFC > 0 & cluster_stats$pvalue < 0.05]
cat(paste0("\nMCAO高表达簇: ", paste(high_clusters, collapse = ", "), "\n"))

if (length(high_clusters) == 0) {
  high_clusters <- cluster_stats$Cluster[1]
  cat(paste0("无显著高表达簇，选择logFC最高的簇: ", high_clusters, "\n"))
}

high_cluster_genes <- cluster_df$Gene[cluster_df$Color %in% high_clusters]
cat(paste0("MCAO高表达簇基因总数: ", length(high_cluster_genes), "\n"))

# ==================== 13. MCAO高表达簇与大鼠BCP靶点/铜死亡基因交集 ====================
cat("\nMCAO高表达簇与大鼠BCP靶点/铜死亡基因取交集...\n")

# 将表达矩阵基因和大鼠基因都转为大写
high_cluster_genes_upper <- toupper(high_cluster_genes)

bcp_in_high <- intersect(toupper(bcp_targets_rat), high_cluster_genes_upper)
copper_in_high <- intersect(toupper(copper_genes_rat), high_cluster_genes_upper)

bcp_coverage <- length(bcp_in_high) / length(bcp_targets_rat)
copper_coverage <- length(copper_in_high) / length(copper_genes_rat)

cat(paste0("BCP靶点（大鼠）在MCAO高表达簇中: ", length(bcp_in_high), "/", length(bcp_targets_rat), 
           " (", round(bcp_coverage * 100, 1), "%)\n"))
if (length(bcp_in_high) > 0) {
  cat("BCP靶点交集基因:\n")
  print(bcp_in_high)
}

cat(paste0("铜死亡基因（大鼠）在MCAO高表达簇中: ", length(copper_in_high), "/", length(copper_genes_rat), 
           " (", round(copper_coverage * 100, 1), "%)\n"))
if (length(copper_in_high) > 0) {
  cat("铜死亡基因交集:\n")
  print(copper_in_high)
}

write.table(bcp_in_high, file = file.path(output_dir, "BCP_targets_in_high_cluster.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(copper_in_high, file = file.path(output_dir, "Copper_genes_in_high_cluster.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==================== 14. BCP靶点集富集统计验证 ====================
cat("\nBCP靶点集富集于MCAO高表达簇的统计验证...\n")

contingency_table <- matrix(c(
  length(bcp_in_high),
  length(high_cluster_genes) - length(bcp_in_high),
  length(bcp_targets_rat) - length(bcp_in_high),
  nrow(gene_exprs) - length(high_cluster_genes) - length(bcp_targets_rat) + length(bcp_in_high)
), nrow = 2)

fisher_res <- fisher.test(contingency_table, alternative = "greater")
OR <- fisher_res$estimate
OR_pvalue <- fisher_res$p.value

cat(paste0("Fisher精确检验:\n"))
cat(paste0("  Odds Ratio = ", round(OR, 3), "\n"))
cat(paste0("  P-value = ", format(OR_pvalue, digits = 3, scientific = TRUE), "\n"))

enrich_df <- data.frame(
  GeneSet = "BCP_targets_rat",
  Cluster = paste(high_clusters, collapse = ","),
  Hits = length(bcp_in_high),
  Total = length(bcp_targets_rat),
  Cluster_size = length(high_cluster_genes),
  Background = nrow(gene_exprs),
  OR = as.numeric(OR),
  Pvalue = OR_pvalue
)
write.table(enrich_df, file = file.path(output_dir, "BCP_enrichment_in_high_cluster.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# 铜死亡基因富集
copper_contingency <- matrix(c(
  length(copper_in_high),
  length(high_cluster_genes) - length(copper_in_high),
  length(copper_genes_rat) - length(copper_in_high),
  nrow(gene_exprs) - length(high_cluster_genes) - length(copper_genes_rat) + length(copper_in_high)
), nrow = 2)

copper_fisher <- fisher.test(copper_contingency, alternative = "greater")
cat(paste0("铜死亡基因 Fisher OR = ", round(as.numeric(copper_fisher$estimate), 3), 
           ", P = ", format(copper_fisher$p.value, digits = 3, scientific = TRUE), "\n"))

copper_enrich_df <- data.frame(
  GeneSet = "Copper_death_rat",
  Cluster = paste(high_clusters, collapse = ","),
  Hits = length(copper_in_high),
  Total = length(copper_genes_rat),
  Cluster_size = length(high_cluster_genes),
  Background = nrow(gene_exprs),
  OR = as.numeric(copper_fisher$estimate),
  Pvalue = copper_fisher$p.value
)
write.table(copper_enrich_df, file = file.path(output_dir, "Copper_enrichment_in_high_cluster.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# ==================== 15. ssGSEA验证（自定义实现） ====================
cat("\n进行 ssGSEA 通路验证...\n")

# 准备通路基因集（使用大鼠基因）
copper_pathway <- list(Copper_Death = copper_genes_rat)
bcp_pathway <- list(BCP_targets = bcp_targets_rat)

# AGE-RAGE通路基因（手动整理核心基因）
age_rage_genes <- c("AGER", "S100A12", "S100B", "S100A6", "S100P", "HMGB1", "HMGB2", 
                    "LGALS3", "MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK14", "NFKB1", 
                    "RELA", "STAT3", "JAK2", "SRC", "PIK3CA", "PIK3CB", "AKT1", "AKT2",
                    "MTOR", "RPS6KB1", "EIF4EBP1", "TNF", "IL6", "IL1B", "CXCL8",
                    "CCL2", "VCAM1", "ICAM1", "SELE", "NOS3", "NOS2", "PTGS2",
                    "TGFB1", "SMAD2", "SMAD3", "SMAD4", "COL1A1", "COL3A1", "FN1",
                    "MMP2", "MMP9", "TIMP1", "BCL2", "BAX", "CASP3", "CASP8", "CASP9")
# 转换为大鼠同源基因
age_rage_rat <- c()
for (gene in age_rage_genes) {
  if (gene %in% names(human_to_rat)) {
    age_rage_rat <- c(age_rage_rat, human_to_rat[[gene]])
  }
}
age_rage_rat <- unique(age_rage_rat)
age_rage_pathway <- list(AGE_RAGE = age_rage_rat)

pathway_list <- c(copper_pathway, age_rage_pathway, bcp_pathway)

# 自定义ssGSEA函数
ssgsea_custom <- function(expr_matrix, gene_sets) {
  result <- matrix(NA, nrow = length(gene_sets), ncol = ncol(expr_matrix))
  rownames(result) <- names(gene_sets)
  colnames(result) <- colnames(expr_matrix)
  
  for (i in seq_along(gene_sets)) {
    gs <- gene_sets[[i]]
    gs_overlap <- intersect(gs, rownames(expr_matrix))
    if (length(gs_overlap) > 0) {
      gs_expr <- expr_matrix[gs_overlap, , drop = FALSE]
      gs_scores <- colMeans(gs_expr, na.rm = TRUE)
      result[i, ] <- gs_scores
    }
  }
  return(result)
}

# 运行ssGSEA
gsva_result <- ssgsea_custom(gene_exprs, pathway_list)

# 比较MCAO vs Sham
ssgsea_df <- as.data.frame(t(gsva_result))
ssgsea_df$Sample <- rownames(ssgsea_df)
ssgsea_df$Group <- sample_info$group[match(ssgsea_df$Sample, sample_info$geo_accession)]

ssgsea_stats <- data.frame()
for (pw in names(pathway_list)) {
  mcaO_scores <- ssgsea_df[[pw]][ssgsea_df$Group == "MCAO_7d"]
  sham_scores <- ssgsea_df[[pw]][ssgsea_df$Group == "Sham_7d"]
  
  if (length(mcaO_scores) > 0 && length(sham_scores) > 0) {
    t_res <- t.test(mcaO_scores, sham_scores)
    
    ssgsea_stats <- rbind(ssgsea_stats, data.frame(
      Pathway = pw,
      MCAO_mean = mean(mcaO_scores),
      Sham_mean = mean(sham_scores),
      logFC = mean(mcaO_scores) - mean(sham_scores),
      t_statistic = t_res$statistic,
      pvalue = t_res$p.value,
      NES = (mean(mcaO_scores) - mean(sham_scores)) / sd(c(mcaO_scores, sham_scores))
    ))
  }
}

ssgsea_stats$FDR <- p.adjust(ssgsea_stats$pvalue, method = "fdr")
write.table(ssgsea_stats, file = file.path(output_dir, "ssGSEA_pathway_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nssGSEA通路分析结果:\n")
print(ssgsea_stats)

# 绘制ssGSEA箱线图
if (nrow(ssgsea_stats) > 0) {
  ssgsea_long <- tidyr::pivot_longer(ssgsea_df, cols = -c(Sample, Group), names_to = "Pathway", values_to = "Score")
  
  p_ssgsea <- ggplot(ssgsea_long, aes(x = Pathway, y = Score, fill = Group)) +
    geom_boxplot(alpha = 0.8) +
    scale_fill_manual(values = c("MCAO_7d" = "#E74C3C", "Sham_7d" = "#2ECC71")) +
    labs(title = "ssGSEA Pathway Enrichment (GSE61616)",
         x = "Pathway", y = "ssGSEA Score") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "ssGSEA_boxplot.png"), plot = p_ssgsea, width = 8, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "ssGSEA_boxplot.pdf"), plot = p_ssgsea, width = 8, height = 6)
}

# ==================== 16. PPI网络输入基因 ====================
cat("\n准备PPI网络输入基因...\n")

ppi_input_genes <- union(high_cluster_genes_upper, toupper(copper_genes_rat))
ppi_input_genes <- unique(ppi_input_genes)
cat(paste0("PPI网络输入基因数 (高表达簇 ∪ 铜死亡核心): ", length(ppi_input_genes), "\n"))

write.table(ppi_input_genes, file = file.path(output_dir, "String_PPI_input_genes.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

ppi_extended <- union(union(high_cluster_genes_upper, toupper(bcp_targets_rat)), toupper(copper_genes_rat))
write.table(ppi_extended, file = file.path(output_dir, "String_PPI_extended_input.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==================== 17. 可视化 ====================
cat("\n生成可视化图表...\n")

# 火山图
all_results$logP <- -log10(all_results$P.Value)
all_results$color <- "Other"
all_results$Gene_upper <- toupper(all_results$Gene)
all_results$color[all_results$Gene_upper %in% toupper(bcp_targets_rat) & all_results$adj.P.Val < 0.05 & abs(all_results$logFC) > 1] <- "BCP_Sig"
all_results$color[all_results$Gene_upper %in% toupper(copper_genes_rat) & all_results$adj.P.Val < 0.05 & abs(all_results$logFC) > 1] <- "Copper_Sig"
all_results$color[all_results$Gene_upper %in% high_cluster_genes_upper & all_results$adj.P.Val < 0.05 & abs(all_results$logFC) > 1] <- "Cluster_Sig"

colors_volcano <- c("Other" = "grey80", "BCP_Sig" = "#3498DB", 
                    "Copper_Sig" = "#E74C3C", "Cluster_Sig" = "#9B59B6")

p_volcano <- ggplot(all_results, aes(x = logFC, y = logP, color = color)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = colors_volcano, name = "Gene Set") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title = "GSE61616: MCAO 7d vs Sham (Rat Orthologs)",
       x = expression(paste("log"[2], " Fold Change")),
       y = expression(paste("-log"[10], " P-value"))) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "volcano_all_genes.png"), plot = p_volcano, width = 8, height = 6, dpi = 300)

# 簇表达差异条形图
p_cluster <- ggplot(cluster_stats, aes(x = reorder(Cluster, -logFC), y = logFC, fill = pvalue < 0.05)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "grey60"), name = "Significant") +
  labs(title = "Cluster Expression Difference (MCAO vs Sham)",
       x = "Cluster", y = "log2 Fold Change") +
  theme_bw(base_size = 12)

ggsave(file.path(output_dir, "cluster_expression_barplot.png"), plot = p_cluster, width = 8, height = 5, dpi = 300)

# ==================== 18. 总结报告 ====================
cat("\n========================================\n")
cat("         分析完成总结报告 (v2)\n")
cat("========================================\n")
cat(paste0("数据集: GSE61616 (7d MCAO vs Sham, n=10)\n"))
cat(paste0("平台: GPL1355\n"))
cat(paste0("样本: MCAO 7d (n=5), Sham 7d (n=5)\n\n"))

cat(paste0("差异表达基因 (|logFC|>1, FDR<0.05): ", nrow(sig_degs), "\n"))
cat(paste0("层次聚类识别基因簇数: ", n_clusters, "\n"))
cat(paste0("MCAO高表达簇: ", paste(high_clusters, collapse = ", "), "\n"))
cat(paste0("MCAO高表达簇基因数: ", length(high_cluster_genes), "\n\n"))

cat("物种映射结果:\n")
cat(paste0("  BCP靶点（人类→大鼠）: ", length(bcp_targets_human), " → ", length(bcp_targets_rat), "\n"))
cat(paste0("  铜死亡基因（人类→大鼠）: ", length(copper_genes_human), " → ", length(copper_genes_rat), "\n\n"))

cat(paste0("BCP靶点（大鼠）在MCAO高表达簇中: ", length(bcp_in_high), "/", length(bcp_targets_rat), 
           " (", round(bcp_coverage*100, 1), "%)\n"))
cat(paste0("  Fisher OR = ", round(OR, 3), ", P = ", format(OR_pvalue, digits = 3, scientific = TRUE), "\n"))
cat(paste0("铜死亡基因（大鼠）在MCAO高表达簇中: ", length(copper_in_high), "/", length(copper_genes_rat), 
           " (", round(copper_coverage*100, 1), "%)\n"))
cat(paste0("  Fisher OR = ", round(as.numeric(copper_fisher$estimate), 3), ", P = ", format(copper_fisher$p.value, digits = 3, scientific = TRUE), "\n\n"))

if (nrow(ssgsea_stats) > 0) {
  cat("ssGSEA通路分析结果:\n")
  for (i in 1:nrow(ssgsea_stats)) {
    cat(paste0("  ", ssgsea_stats$Pathway[i], ": MCAO=", round(ssgsea_stats$MCAO_mean[i], 3), 
               ", Sham=", round(ssgsea_stats$Sham_mean[i], 3), 
               ", logFC=", round(ssgsea_stats$logFC[i], 3), 
               ", P=", format(ssgsea_stats$pvalue[i], digits = 3, scientific = TRUE), "\n"))
  }
}

cat(paste0("\nPPI网络输入基因数: ", length(ppi_input_genes), "\n"))
cat("  请将这些基因提交至 String DB 进行PPI分析\n")
cat("  然后在Cytoscape中使用CytoHubba插件计算MCC得分筛选Hub基因\n\n")

cat(paste0("输出目录: ", output_dir, "\n"))
cat("========================================\n")
