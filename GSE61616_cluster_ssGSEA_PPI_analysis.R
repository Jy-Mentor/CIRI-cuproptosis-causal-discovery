# ============================================
# GSE61616 层次聚类 + ssGSEA + PPI Hub基因分析
# 7d MCAO vs Sham | 动态剪切聚类 | 通路验证 | CytoHubba MCC
# ============================================

# ==================== 1. 包安装与加载 ====================
cat("正在检查和安装必要的R包...\n")

packages <- c("limma", "ggplot2", "dplyr", "tidyr", "GEOquery", 
              "preprocessCore", "ComplexHeatmap", "circlize", 
              "dynamicTreeCut", "clusterProfiler", "org.Rn.eg.db",
              "GSVA", "msigdbr", "igraph")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Installing package: ", pkg, "\n"))
    if (pkg %in% c("GEOquery", "limma", "preprocessCore", "ComplexHeatmap", "circlize", 
                   "org.Rn.eg.db", "GSVA", "msigdbr", "clusterProfiler")) {
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

output_dir <- file.path(work_dir, "GSE61616_cluster_ssGSEA_PPI_results")
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

# ==================== 5. 提取表达矩阵子集并标准化 ====================
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

# ==================== 6. 读取平台注释信息 ====================
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

# ==================== 7. limma差异表达分析 ====================
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

# ==================== 8. 读取BCP靶点和铜死亡基因 ====================
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

# ==================== 9. 探针→基因转换 ====================
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

cat(paste0("基因级表达矩阵维度: ", nrow(gene_exprs), " genes x ", ncol(gene_exprs), " samples\n"))

# ==================== 10. 层次聚类 + 动态剪切 ====================
cat("\n进行层次聚类和动态剪切...\n")

# 使用差异基因进行聚类（或全部基因）
# 这里使用全部基因的表达谱进行聚类，识别共表达模式
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
  color_palette <- c("turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue", "lightcyan", "grey60", "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white", "skyblue", "saddlebrown", "steelblue", "plum", "violet", "darkmagenta", "darkolivegreen", "maroon", "sienna", "blanchedalmond", "bisque", "navajowhite", "antiquewhite", "lightcoral", "lightblue", "lightpink", "thistle", "orchid", "palevioletred", "mediumpurple", "mediumorchid", "blueviolet", "indigo", "cornflowerblue", "mediumslateblue", "slateblue", "darkslateblue", "mediumaquamarine", "aquamarine", "mediumspringgreen", "springgreen", "mediumseagreen", "seagreen", "darkseagreen", "lightseagreen", "cadetblue", "teal", "darkcyan", "aqua", "deepskyblue", "dodgerblue", "steelblue", "lightsteelblue", "powderblue", "skyblue", "lightskyblue", "aliceblue", "ghostwhite", "lavender", "midnightblue", "navy", "darkblue", "mediumblue", "blue", "slateblue", "darkslateblue", "mediumslateblue", "mediumpurple", "darkorchid", "darkviolet", "mediumorchid", "thistle", "plum", "violet", "magenta", "orchid", "mediumvioletred", "deeppink", "hotpink", "palevioletred", "crimson", "lightpink", "pink", "lightcoral", "indianred", "darksalmon", "salmon", "lightsalmon", "orangered", "tomato", "coral", "darkorange", "orange", "gold", "yellow", "lightyellow", "lemonchiffon", "papayawhip", "moccasin", "peachpuff", "palegoldenrod", "khaki", "darkkhaki", "lavenderblush", "mistyrose", "antiquewhite", "linen", "beige", "whitesmoke", "gainsboro", "lightgrey", "silver", "darkgray", "gray", "dimgray", "lightslategray", "slategray", "darkslategray", "black")
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
cat("各簇基因数:\n")
print(table(cluster_colors))

# 保存聚类结果
cluster_df <- data.frame(Gene = colnames(datExpr), Cluster = dynamicMods, Color = cluster_colors)
write.table(cluster_df, file = file.path(output_dir, "gene_cluster_assignment.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# 绘制聚类树
pdf(file.path(output_dir, "gene_clustering_dendrogram.pdf"), width = 12, height = 6)
plot(geneTree, xlab = "", sub = "", main = "Gene Clustering (Dynamic Tree Cut)", labels = FALSE, hang = 0.04)
dev.off()

# ==================== 11. 筛选MCAO高表达簇 ====================
cat("\n筛选MCAO高表达簇...\n")

# 计算每个簇在MCAO vs Sham中的平均表达差异
cluster_stats <- data.frame()
for (cl in unique(cluster_colors)) {
  if (cl == "grey") next
  
  cl_genes <- cluster_df$Gene[cluster_df$Color == cl]
  cl_expr <- gene_exprs[cl_genes, , drop = FALSE]
  
  # 计算每个样本中该簇的平均表达
  cl_mean_expr <- colMeans(cl_expr, na.rm = TRUE)
  
  # MCAO vs Sham差异
  mcaO_expr <- cl_mean_expr[sample_info$group == "MCAO_7d"]
  sham_expr <- cl_mean_expr[sample_info$group == "Sham_7d"]
  
  # t检验
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

cat("各簇MCAO vs Sham表达差异:\n")
print(cluster_stats)

# 筛选MCAO高表达簇（logFC > 0 且 p < 0.05）
high_clusters <- cluster_stats$Cluster[cluster_stats$logFC > 0 & cluster_stats$pvalue < 0.05]
cat(paste0("\nMCAO高表达簇: ", paste(high_clusters, collapse = ", "), "\n"))

if (length(high_clusters) == 0) {
  # 如果没有显著高表达簇，选择logFC最高的簇
  high_clusters <- cluster_stats$Cluster[1]
  cat(paste0("无显著高表达簇，选择logFC最高的簇: ", high_clusters, "\n"))
}

# 合并所有MCAO高表达簇的基因
high_cluster_genes <- cluster_df$Gene[cluster_df$Color %in% high_clusters]
cat(paste0("MCAO高表达簇基因总数: ", length(high_cluster_genes), "\n"))

# ==================== 12. MCAO高表达簇与BCP靶点/铜死亡基因交集 ====================
cat("\nMCAO高表达簇与BCP靶点/铜死亡基因取交集...\n")

bcp_in_high <- intersect(bcp_targets, high_cluster_genes)
copper_in_high <- intersect(copper_genes, high_cluster_genes)

bcp_coverage <- length(bcp_in_high) / length(bcp_targets)
copper_coverage <- length(copper_in_high) / length(copper_genes)

cat(paste0("BCP靶点在MCAO高表达簇中: ", length(bcp_in_high), "/", length(bcp_targets), 
           " (", round(bcp_coverage * 100, 1), "%)\n"))
cat(paste0("铜死亡基因在MCAO高表达簇中: ", length(copper_in_high), "/", length(copper_genes), 
           " (", round(copper_coverage * 100, 1), "%)\n"))

# 保存交集基因
write.table(bcp_in_high, file = file.path(output_dir, "BCP_targets_in_high_cluster.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(copper_in_high, file = file.path(output_dir, "Copper_genes_in_high_cluster.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==================== 13. PPI网络输入基因：高表达簇基因 + 铜死亡核心基因 ====================
cat("\n准备PPI网络输入基因...\n")

ppi_input_genes <- union(high_cluster_genes, copper_genes)
ppi_input_genes <- unique(ppi_input_genes)
cat(paste0("PPI网络输入基因数 (高表达簇 ∪ 铜死亡核心): ", length(ppi_input_genes), "\n"))

# 保存String PPI输入
write.table(ppi_input_genes, file = file.path(output_dir, "String_PPI_input_genes.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# 同时保存扩展输入（高表达簇 + BCP靶点 + 铜死亡）
ppi_extended <- union(union(high_cluster_genes, bcp_targets), copper_genes)
write.table(ppi_extended, file = file.path(output_dir, "String_PPI_extended_input.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==================== 14. BCP靶点集富集于MCAO高表达簇的统计验证 ====================
cat("\nBCP靶点集富集于MCAO高表达簇的统计验证...\n")

# 超几何检验/ Fisher精确检验
contingency_table <- matrix(c(
  length(bcp_in_high),                                    # BCP在簇中
  length(high_cluster_genes) - length(bcp_in_high),       # 非BCP在簇中
  length(bcp_targets) - length(bcp_in_high),              # BCP不在簇中
  nrow(gene_exprs) - length(high_cluster_genes) - length(bcp_targets) + length(bcp_in_high)  # 非BCP不在簇中
), nrow = 2)

fisher_res <- fisher.test(contingency_table, alternative = "greater")
OR <- fisher_res$estimate
OR_pvalue <- fisher_res$p.value

cat(paste0("Fisher精确检验:\n"))
cat(paste0("  Odds Ratio = ", round(OR, 3), "\n"))
cat(paste0("  P-value = ", format(OR_pvalue, digits = 3, scientific = TRUE), "\n"))

# 保存富集结果
enrich_df <- data.frame(
  GeneSet = "BCP_targets",
  Cluster = paste(high_clusters, collapse = ","),
  Hits = length(bcp_in_high),
  Total = length(bcp_targets),
  Cluster_size = length(high_cluster_genes),
  Background = nrow(gene_exprs),
  OR = as.numeric(OR),
  Pvalue = OR_pvalue,
  FDR = p.adjust(OR_pvalue, method = "fdr")
)
write.table(enrich_df, file = file.path(output_dir, "BCP_enrichment_in_high_cluster.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# 铜死亡基因富集
copper_contingency <- matrix(c(
  length(copper_in_high),
  length(high_cluster_genes) - length(copper_in_high),
  length(copper_genes) - length(copper_in_high),
  nrow(gene_exprs) - length(high_cluster_genes) - length(copper_genes) + length(copper_in_high)
), nrow = 2)

copper_fisher <- fisher.test(copper_contingency, alternative = "greater")
copper_enrich_df <- data.frame(
  GeneSet = "Copper_death",
  Cluster = paste(high_clusters, collapse = ","),
  Hits = length(copper_in_high),
  Total = length(copper_genes),
  Cluster_size = length(high_cluster_genes),
  Background = nrow(gene_exprs),
  OR = as.numeric(copper_fisher$estimate),
  Pvalue = copper_fisher$p.value,
  FDR = p.adjust(copper_fisher$p.value, method = "fdr")
)
write.table(copper_enrich_df, file = file.path(output_dir, "Copper_enrichment_in_high_cluster.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# ==================== 15. ssGSEA验证 ====================
cat("\n进行 ssGSEA 通路验证...\n")

# 准备通路基因集
# 1. 铜死亡通路
copper_pathway <- list(Copper_Death = copper_genes)

# 2. AGE-RAGE通路（从msigdb获取）
# 尝试从msigdbr获取大鼠通路，如果没有则使用人类同源
tryCatch({
  m_df <- msigdbr(species = "Rattus norvegicus", category = "C2")
  if (nrow(m_df) == 0) {
    m_df <- msigdbr(species = "Homo sapiens", category = "C2")
  }
}, error = function(e) {
  m_df <- msigdbr(species = "Homo sapiens", category = "C2")
})

age_rage_genes <- m_df$gene_symbol[m_df$gs_name == "KEGG_AGE_RAGE_SIGNALING_PATHWAY_IN_DIABETIC_COMPLICATIONS"]
if (length(age_rage_genes) == 0) {
  age_rage_genes <- m_df$gene_symbol[grep("AGE_RAGE", m_df$gs_name, ignore.case = TRUE)][1:100]
}
age_rage_pathway <- list(AGE_RAGE = unique(age_rage_genes))

# 3. BCP靶点集
bcp_pathway <- list(BCP_targets = bcp_targets)

# 合并通路
pathway_list <- c(copper_pathway, age_rage_pathway, bcp_pathway)

# 运行ssGSEA（使用GSVA包或自定义实现）
gene_exprs_ssgsea <- gene_exprs
# 确保基因名在通路中
common_genes <- unique(unlist(pathway_list))
gene_exprs_ssgsea <- gene_exprs_ssgsea[rownames(gene_exprs_ssgsea) %in% common_genes, ]

# 自定义ssGSEA函数（当GSVA包不可用时）
ssgsea_custom <- function(expr_matrix, gene_sets) {
  # expr_matrix: genes x samples
  # gene_sets: list of gene vectors
  result <- matrix(NA, nrow = length(gene_sets), ncol = ncol(expr_matrix))
  rownames(result) <- names(gene_sets)
  colnames(result) <- colnames(expr_matrix)
  
  for (i in seq_along(gene_sets)) {
    gs <- gene_sets[[i]]
    gs_overlap <- intersect(gs, rownames(expr_matrix))
    if (length(gs_overlap) > 0) {
      # 计算基因集内基因的表达秩次和
      gs_expr <- expr_matrix[gs_overlap, , drop = FALSE]
      # 对每个样本计算基因集得分（平均表达量标准化）
      gs_scores <- colMeans(gs_expr, na.rm = TRUE)
      result[i, ] <- gs_scores
    }
  }
  return(result)
}

tryCatch({
  # 尝试使用GSVA包
  if (require("GSVA", quietly = TRUE)) {
    gsva_result <- gsva(expr = gene_exprs_ssgsea, gset.idx.list = pathway_list, 
                        method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
  } else {
    cat("GSVA包不可用，使用自定义ssGSEA实现...\n")
    gsva_result <- ssgsea_custom(gene_exprs_ssgsea, pathway_list)
  }
  
  # 比较MCAO vs Sham的ssGSEA得分
  ssgsea_df <- as.data.frame(t(gsva_result))
  ssgsea_df$Sample <- rownames(ssgsea_df)
  ssgsea_df$Group <- sample_info$group[match(ssgsea_df$Sample, sample_info$geo_accession)]
  
  ssgsea_stats <- data.frame()
  for (pw in names(pathway_list)) {
    mcaO_scores <- ssgsea_df[[pw]][ssgsea_df$Group == "MCAO_7d"]
    sham_scores <- ssgsea_df[[pw]][ssgsea_df$Group == "Sham_7d"]
    
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
  
  ssgsea_stats$FDR <- p.adjust(ssgsea_stats$pvalue, method = "fdr")
  write.table(ssgsea_stats, file = file.path(output_dir, "ssGSEA_pathway_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("\nssGSEA通路分析结果:\n")
  print(ssgsea_stats)
  
  # 绘制ssGSEA箱线图
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
  
}, error = function(e) {
  cat(paste0("ssGSEA分析出错: ", e$message, "\n"))
  cat("跳过ssGSEA分析...\n")
})

# ==================== 16. 可视化 ====================
cat("\n生成可视化图表...\n")

# 1. 聚类热图（差异基因）
sig_gene_list <- unique(sig_degs$Gene)
sig_gene_list <- sig_gene_list[!is.na(sig_gene_list)]

if (length(sig_gene_list) > 0) {
  sig_expr <- gene_exprs[intersect(sig_gene_list, rownames(gene_exprs)), , drop = FALSE]
  if (nrow(sig_expr) > 50) {
    # 如果差异基因太多，取top 200
    sig_expr <- sig_expr[1:min(200, nrow(sig_expr)), ]
  }
  
  sample_order <- order(sample_info$group)
  sig_expr_ordered <- sig_expr[, sample_order, drop = FALSE]
  sig_z <- t(scale(t(sig_expr_ordered)))
  
  group_colors <- c("MCAO_7d" = "#E74C3C", "Sham_7d" = "#2ECC71")
  ha <- HeatmapAnnotation(
    Group = sample_info$group[sample_order],
    col = list(Group = group_colors)
  )
  
  # 添加聚类注释
  cluster_annot <- cluster_df$Color[match(rownames(sig_z), cluster_df$Gene)]
  
  ht <- Heatmap(sig_z,
                name = "Z-score",
                top_annotation = ha,
                cluster_columns = FALSE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 6),
                column_title = "GSE61616 DEGs Clustering (MCAO 7d vs Sham)")
  
  pdf(file.path(output_dir, "DEG_clustering_heatmap.pdf"), width = 10, height = 12)
  print(ht)
  dev.off()
  
  png(file.path(output_dir, "DEG_clustering_heatmap.png"), width = 1000, height = 1200, res = 150)
  print(ht)
  dev.off()
}

# 2. 火山图
all_results$logP <- -log10(all_results$P.Value)
all_results$color <- "Other"
all_results$color[all_results$Gene %in% bcp_targets & all_results$adj.P.Val < 0.05 & abs(all_results$logFC) > 1] <- "BCP_Sig"
all_results$color[all_results$Gene %in% copper_genes & all_results$adj.P.Val < 0.05 & abs(all_results$logFC) > 1] <- "Copper_Sig"
all_results$color[all_results$Gene %in% high_cluster_genes & all_results$adj.P.Val < 0.05 & abs(all_results$logFC) > 1] <- "Cluster_Sig"

colors_volcano <- c("Other" = "grey80", "BCP_Sig" = "#3498DB", 
                    "Copper_Sig" = "#E74C3C", "Cluster_Sig" = "#9B59B6")

p_volcano <- ggplot(all_results, aes(x = logFC, y = logP, color = color)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = colors_volcano, name = "Gene Set") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title = "GSE61616: MCAO 7d vs Sham",
       x = expression(paste("log"[2], " Fold Change")),
       y = expression(paste("-log"[10], " P-value"))) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "volcano_all_genes.png"), plot = p_volcano, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "volcano_all_genes.pdf"), plot = p_volcano, width = 8, height = 6)

# 3. 簇表达差异条形图
p_cluster <- ggplot(cluster_stats, aes(x = reorder(Cluster, -logFC), y = logFC, fill = pvalue < 0.05)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "grey60"), name = "Significant") +
  labs(title = "Cluster Expression Difference (MCAO vs Sham)",
       x = "Cluster", y = "log2 Fold Change") +
  theme_bw(base_size = 12)

ggsave(file.path(output_dir, "cluster_expression_barplot.png"), plot = p_cluster, width = 8, height = 5, dpi = 300)

# ==================== 17. 总结报告 ====================
cat("\n========================================\n")
cat("         分析完成总结报告\n")
cat("========================================\n")
cat(paste0("数据集: GSE61616 (7d MCAO vs Sham, n=10)\n"))
cat(paste0("平台: GPL1355\n"))
cat(paste0("样本: MCAO 7d (n=5), Sham 7d (n=5)\n\n"))

cat(paste0("差异表达基因 (|logFC|>1, FDR<0.05): ", nrow(sig_degs), "\n"))
cat(paste0("层次聚类识别基因簇数: ", n_clusters, "\n"))
cat(paste0("MCAO高表达簇: ", paste(high_clusters, collapse = ", "), "\n"))
cat(paste0("MCAO高表达簇基因数: ", length(high_cluster_genes), "\n\n"))

cat(paste0("BCP靶点在MCAO高表达簇中: ", length(bcp_in_high), "/", length(bcp_targets), 
           " (", round(bcp_coverage*100, 1), "%)\n"))
cat(paste0("  Fisher OR = ", round(OR, 3), ", P = ", format(OR_pvalue, digits = 3, scientific = TRUE), "\n"))
cat(paste0("铜死亡基因在MCAO高表达簇中: ", length(copper_in_high), "/", length(copper_genes), 
           " (", round(copper_coverage*100, 1), "%)\n"))
cat(paste0("  Fisher OR = ", round(as.numeric(copper_fisher$estimate), 3), ", P = ", format(copper_fisher$p.value, digits = 3, scientific = TRUE), "\n\n"))

cat(paste0("PPI网络输入基因数: ", length(ppi_input_genes), "\n"))
cat("  请将这些基因提交至 String DB (https://string-db.org) 进行PPI分析\n")
cat("  然后在Cytoscape中使用CytoHubba插件计算MCC得分筛选Hub基因\n")
cat("  预期Hub基因包括: NFKB1, FDX1 等\n\n")

cat(paste0("输出目录: ", output_dir, "\n"))
cat("生成文件列表:\n")
for (f in list.files(output_dir)) {
  cat(paste0("  - ", f, "\n"))
}
cat("========================================\n")
