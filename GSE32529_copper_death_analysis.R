# ============================================
# GSE32529 铜死亡差异基因分析 (CIRI模型)
# 筛选条件: |logFC| > 1, FDR < 0.05
# ============================================

# ==================== 1. 包安装与加载 ====================
cat("正在检查和安装必要的R包...\n")

# 基础包列表
packages <- c("limma", "ggplot2", "dplyr", "tidyr", "GEOquery", 
              "preprocessCore", "ComplexHeatmap", "circlize")

# 智能安装函数
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

# 加载包
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

# 输入文件路径
series_matrix_file <- "C:/Users/Jy-Mentor-7/Downloads/GSE32529_series_matrix (1).txt.gz"
platform_file <- "C:/Users/Jy-Mentor-7/Downloads/GPL1261-56135.txt"

# 验证输入文件存在
if (!file.exists(series_matrix_file)) {
  stop("系列矩阵文件不存在!")
}
if (!file.exists(platform_file)) {
  stop("平台注释文件不存在!")
}

cat(paste0("系列矩阵文件: ", series_matrix_file, "\n"))
cat(paste0("平台注释文件: ", platform_file, "\n\n"))

# 输出目录
output_dir <- file.path(work_dir, "GSE32529_copper_death_results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==================== 3. 设置随机种子 ====================
set.seed(42)

# ==================== 4. 读取系列矩阵数据 ====================
cat("读取 GSE32529 系列矩阵数据...\n")

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

# 验证表达矩阵
if (nrow(exprs) == 0 || ncol(exprs) == 0) {
  stop("表达矩阵为空!")
}

# ==================== 5. 定义样本分组信息 ====================
cat("定义样本分组信息...\n")

sample_info <- data.frame(
  geo_accession = c(
    paste0("GSM805", 728:731),  # CIRI 3h
    paste0("GSM805", 732:735),  # CIRI 24h
    paste0("GSM805", 704:707),  # Sham 3h
    paste0("GSM805", 716:719)   # Sham 24h
  ),
  group = c(
    rep("CIRI_3h", 4), rep("CIRI_24h", 4),
    rep("Sham_3h", 4), rep("Sham_24h", 4)
  ),
  treatment = c(rep("CIRI", 8), rep("Sham", 8)),
  time = c(rep("3h", 4), rep("24h", 4), rep("3h", 4), rep("24h", 4)),
  stringsAsFactors = FALSE
)

# 验证所有样本存在于表达矩阵中
missing_samples <- setdiff(sample_info$geo_accession, colnames(exprs))
if (length(missing_samples) > 0) {
  cat("警告: 以下样本在表达矩阵中不存在，尝试使用实际列名...\n")
  print(missing_samples)
  # 使用实际的样本名称
  actual_samples <- colnames(exprs)
  cat("实际样本名称:\n")
  print(head(actual_samples))
  
  # 尝试匹配样本
  if (length(actual_samples) >= 16) {
    sample_info$geo_accession <- actual_samples[1:16]
    cat("使用前16个样本进行分析\n")
  }
}

cat("样本分组情况:\n")
print(table(sample_info$group))

# ==================== 6. 提取表达矩阵子集 ====================
exprs_subset <- exprs[, sample_info$geo_accession, drop = FALSE]
cat(paste0("子集表达矩阵维度: ", nrow(exprs_subset), " probes x ", ncol(exprs_subset), " samples\n\n"))

# ==================== 7. 检查数据范围并标准化 ====================
cat("检查表达值范围...\n")
data_range <- range(exprs_subset, na.rm = TRUE)
cat(paste0("表达值范围: ", round(data_range[1], 2), " - ", round(data_range[2], 2), "\n"))

# 如果范围在0-20之间，说明已经是log2转换后的数据
if (data_range[2] > 100) {
  cat("执行log2转换...\n")
  exprs_subset <- log2(exprs_subset + 1)
} else {
  cat("数据已经是log2转换后的值\n")
}

# quantile标准化
cat("执行quantile标准化...\n")
library(preprocessCore)
exprs_norm <- normalize.quantiles(exprs_subset)
rownames(exprs_norm) <- rownames(exprs_subset)
colnames(exprs_norm) <- colnames(exprs_subset)

# ==================== 8. 读取平台注释信息 ====================
cat("读取 GPL1261 平台注释...\n")

plat_lines <- readLines(platform_file, warn = FALSE)
header_line <- which(grepl("^ID", plat_lines))[1]
if (length(header_line) == 0 || is.na(header_line)) {
  stop("无法找到平台文件头部行")
}

platform_annot <- read.delim(platform_file, 
                             skip = header_line - 1, 
                             header = TRUE,
                             check.names = FALSE,
                             stringsAsFactors = FALSE,
                             quote = "",
                             fill = TRUE)

cat(paste0("平台注释探针数: ", nrow(platform_annot), "\n"))

# 提取关键列
annot_table <- data.frame(
  ID = platform_annot$ID,
  Gene_Symbol = if("Gene Symbol" %in% names(platform_annot)) platform_annot[["Gene Symbol"]] else NA,
  Gene_Title = if("Gene Title" %in% names(platform_annot)) platform_annot[["Gene Title"]] else NA,
  Entrez_ID = if("ENTREZ_GENE_ID" %in% names(platform_annot)) platform_annot[["ENTREZ_GENE_ID"]] else NA,
  stringsAsFactors = FALSE
)

# ==================== 9. limma差异表达分析 ====================
cat("\n开始进行 limma 差异表达分析...\n")

# 创建设计矩阵
sample_info$group <- factor(sample_info$group, 
                            levels = c("Sham_3h", "CIRI_3h", "Sham_24h", "CIRI_24h"))
design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- gsub("group", "", colnames(design))

cat("设计矩阵:\n")
print(head(design))

# 线性模型拟合
fit <- lmFit(exprs_norm, design)

# 定义对比矩阵
contrast.matrix <- makeContrasts(
  CIRI_3h_vs_Sham_3h = CIRI_3h - Sham_3h,
  CIRI_24h_vs_Sham_24h = CIRI_24h - Sham_24h,
  levels = design
)

cat("\n对比矩阵:\n")
print(contrast.matrix)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# ==================== 10. 定义铜死亡基因列表 ====================
copper_genes <- c("FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX",
                  "SLC31A1", "ATP7A", "ATP7B", "ATOX1",
                  "NFE2L2", "HIF1A", "MTOR", "NFKB1", "GPX4")

cat(paste0("\n铜死亡基因列表 (n=", length(copper_genes), "):\n"))
print(copper_genes)

# ==================== 11. 查找铜死亡基因探针 ====================
cat("\n查找铜死亡基因对应的探针...\n")

# 使用正则表达式匹配基因符号（支持多种分隔符）
copper_probes_idx <- c()
for (gene in copper_genes) {
  # 精确匹配或作为第一个基因符号
  pattern <- paste0("(^|/// )", gene, "($| ///)")
  idx <- grep(pattern, annot_table$Gene_Symbol, ignore.case = TRUE)
  if (length(idx) > 0) {
    copper_probes_idx <- c(copper_probes_idx, idx)
    cat(paste0(gene, ": 找到 ", length(idx), " 个探针\n"))
  }
}

copper_probes <- unique(annot_table$ID[copper_probes_idx])
copper_genes_found <- unique(annot_table$Gene_Symbol[copper_probes_idx])

# 过滤出存在于表达矩阵中的探针
copper_probes <- intersect(copper_probes, rownames(exprs_norm))

cat(paste0("\n找到 ", length(copper_probes), " 个铜死亡基因探针\n"))
cat(paste0("匹配的基因数: ", length(copper_genes_found), "\n"))

# 输出未映射的基因
unmapped_genes <- setdiff(copper_genes, copper_genes_found)
if (length(unmapped_genes) > 0) {
  cat(paste0("未找到探针的基因: ", paste(unmapped_genes, collapse = ", "), "\n"))
  write.table(unmapped_genes, 
              file = file.path(output_dir, "unmapped_copper_genes.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# ==================== 12. 提取铜死亡基因统计结果 ====================
cat("\n提取铜死亡基因统计结果...\n")

# 3h对比结果
results_3h <- topTable(fit2, coef = "CIRI_3h_vs_Sham_3h",
                       number = Inf, adjust.method = "fdr", sort.by = "p")
results_3h$ProbeID <- rownames(results_3h)
results_3h <- results_3h[results_3h$ProbeID %in% copper_probes, ]

# 24h对比结果  
results_24h <- topTable(fit2, coef = "CIRI_24h_vs_Sham_24h",
                        number = Inf, adjust.method = "fdr", sort.by = "p")
results_24h$ProbeID <- rownames(results_24h)
results_24h <- results_24h[results_24h$ProbeID %in% copper_probes, ]

cat(paste0("3h对比: ", nrow(results_3h), " 个铜死亡基因探针\n"))
cat(paste0("24h对比: ", nrow(results_24h), " 个铜死亡基因探针\n\n"))

# ==================== 13. 合并结果并添加基因注释 ====================
cat("合并两个时间点结果...\n")

# 合并结果
common_cols <- c("ProbeID", "logFC", "P.Value", "adj.P.Val")
if (nrow(results_3h) > 0) {
  results_3h_sel <- results_3h[, common_cols]
  names(results_3h_sel)[2:4] <- c("logFC_3h", "P.Value_3h", "adj.P.Val_3h")
} else {
  results_3h_sel <- data.frame(ProbeID = character(0), logFC_3h = numeric(0), 
                               P.Value_3h = numeric(0), adj.P.Val_3h = numeric(0))
}

if (nrow(results_24h) > 0) {
  results_24h_sel <- results_24h[, common_cols]
  names(results_24h_sel)[2:4] <- c("logFC_24h", "P.Value_24h", "adj.P.Val_24h")
} else {
  results_24h_sel <- data.frame(ProbeID = character(0), logFC_24h = numeric(0),
                                P.Value_24h = numeric(0), adj.P.Val_24h = numeric(0))
}

# 合并两个结果
copper_results <- merge(results_3h_sel, results_24h_sel, by = "ProbeID", all = TRUE)

# 添加基因注释
copper_results$Gene <- annot_table$Gene_Symbol[match(copper_results$ProbeID, annot_table$ID)]

# 排序
copper_results <- copper_results[order(copper_results$Gene), ]

# ==================== 14. 筛选显著差异的铜死亡基因 ====================
cat("\n筛选显著差异的铜死亡基因...\n")

# 使用用户指定的筛选条件: |logFC| > 1, FDR < 0.05
sig_threshold <- 0.05
fc_threshold <- 1  # logFC阈值

sig_3h <- !is.na(copper_results$logFC_3h) & abs(copper_results$logFC_3h) > fc_threshold & 
          !is.na(copper_results$adj.P.Val_3h) & copper_results$adj.P.Val_3h < sig_threshold
sig_24h <- !is.na(copper_results$logFC_24h) & abs(copper_results$logFC_24h) > fc_threshold &
           !is.na(copper_results$adj.P.Val_24h) & copper_results$adj.P.Val_24h < sig_threshold

copper_sig <- copper_results[sig_3h | sig_24h, ]

cat(paste0("显著差异的铜死亡基因 (|logFC| > 1, FDR < 0.05): ", nrow(copper_sig), "\n"))

# 统计上下调
cat("\n===== 铜死亡基因差异表达汇总 =====\n")
cat(paste0("3h对比 (CIRI vs Sham):\n"))
cat(paste0("  上调: ", sum(results_3h$logFC > fc_threshold & results_3h$adj.P.Val < sig_threshold, na.rm = TRUE), "\n"))
cat(paste0("  下调: ", sum(results_3h$logFC < -fc_threshold & results_3h$adj.P.Val < sig_threshold, na.rm = TRUE), "\n"))

cat(paste0("\n24h对比 (CIRI vs Sham):\n"))
cat(paste0("  上调: ", sum(results_24h$logFC > fc_threshold & results_24h$adj.P.Val < sig_threshold, na.rm = TRUE), "\n"))
cat(paste0("  下调: ", sum(results_24h$logFC < -fc_threshold & results_24h$adj.P.Val < sig_threshold, na.rm = TRUE), "\n"))
cat("================================\n\n")

# ==================== 15. 保存结果文件 ====================
cat("保存结果文件...\n")

# 保存所有铜死亡基因结果
write.table(copper_results,
            file = file.path(output_dir, "GSE32529_copper_death_all_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 保存显著差异的铜死亡基因
if (nrow(copper_sig) > 0) {
  write.table(copper_sig,
              file = file.path(output_dir, "GSE32529_copper_death_DEGs_logFC1_FDR0.05.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# 保存3h对比结果
write.table(results_3h,
            file = file.path(output_dir, "GSE32529_copper_death_3h_results.txt"),
            sep = "\t", quote = FALSE)

# 保存24h对比结果
write.table(results_24h,
            file = file.path(output_dir, "GSE32529_copper_death_24h_results.txt"),
            sep = "\t", quote = FALSE)

# 保存铜死亡基因表达矩阵
if (length(copper_probes) > 0) {
  copper_expr <- exprs_norm[copper_probes, , drop = FALSE]
  write.table(copper_expr,
              file = file.path(output_dir, "GSE32529_copper_death_expression_matrix.txt"),
              sep = "\t", quote = FALSE)
}

cat(paste0("结果文件已保存到: ", output_dir, "\n\n"))

# ==================== 16. 生成可视化 ====================
cat("生成可视化图表...\n")

# 准备热图数据
if (length(copper_probes) > 0) {
  copper_expr <- exprs_norm[copper_probes, , drop = FALSE]
  
  # 添加基因符号作为行名
  gene_symbols <- annot_table$Gene_Symbol[match(rownames(copper_expr), annot_table$ID)]
  rownames(copper_expr) <- ifelse(is.na(gene_symbols), rownames(copper_expr), gene_symbols)
  
  # 去除NA行
  copper_expr <- copper_expr[!is.na(rownames(copper_expr)), , drop = FALSE]
  
  if (nrow(copper_expr) > 0) {
    # 按分组排序样本
    sample_order <- order(sample_info$group)
    copper_expr_ordered <- copper_expr[, sample_order, drop = FALSE]
    
    # Z-score标准化
    copper_expr_z <- t(scale(t(copper_expr_ordered)))
    
    # 创建分组颜色注释
    group_colors <- c("CIRI_3h" = "#E74C3C", "CIRI_24h" = "#C0392B",
                      "Sham_3h" = "#2ECC71", "Sham_24h" = "#27AE60")
    
    library(ComplexHeatmap)
    library(circlize)
    
    ha <- HeatmapAnnotation(
      Group = sample_info$group[sample_order],
      col = list(Group = group_colors)
    )
    
    # 绘制热图
    ht <- Heatmap(copper_expr_z,
                  name = "Z-score",
                  top_annotation = ha,
                  cluster_columns = FALSE,
                  show_column_names = FALSE,
                  row_names_gp = gpar(fontsize = 10),
                  column_title = "GSE32529 Copper Death Genes Expression (CIRI vs Sham)")
    
    # 保存热图
    pdf(file.path(output_dir, "GSE32529_copper_death_heatmap.pdf"), width = 10, height = 8)
    print(ht)
    dev.off()
    
    png(file.path(output_dir, "GSE32529_copper_death_heatmap.png"), width = 1000, height = 800, res = 150)
    print(ht)
    dev.off()
    
    cat("热图已保存!\n")
  }
}

# ==================== 17. 生成火山图 ====================
cat("生成火山图...\n")

# 3h火山图
if (nrow(results_3h) > 0) {
  results_3h$logP <- -log10(results_3h$P.Value)
  results_3h$color <- ifelse(results_3h$adj.P.Val < 0.05 & results_3h$logFC > 1, "Upregulated",
                             ifelse(results_3h$adj.P.Val < 0.05 & results_3h$logFC < -1, "Downregulated", 
                                    "Not significant"))
  
  colors_volcano <- c("Upregulated" = "#FF4040", "Downregulated" = "#4169E1", "Not significant" = "grey80")
  
  p1 <- ggplot(results_3h, aes(x = logFC, y = logP, color = color)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = colors_volcano, name = "Significance") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(title = "GSE32529: CIRI 3h vs Sham 3h (Copper Death Genes)",
         x = expression(paste("log"[2], " Fold Change")),
         y = expression(paste("-log"[10], " P-value"))) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
  
  ggsave(file.path(output_dir, "GSE32529_copper_death_volcano_3h.png"), plot = p1, width = 8, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "GSE32529_copper_death_volcano_3h.pdf"), plot = p1, width = 8, height = 6)
}

# 24h火山图
if (nrow(results_24h) > 0) {
  results_24h$logP <- -log10(results_24h$P.Value)
  results_24h$color <- ifelse(results_24h$adj.P.Val < 0.05 & results_24h$logFC > 1, "Upregulated",
                              ifelse(results_24h$adj.P.Val < 0.05 & results_24h$logFC < -1, "Downregulated", 
                                     "Not significant"))
  
  p2 <- ggplot(results_24h, aes(x = logFC, y = logP, color = color)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = colors_volcano, name = "Significance") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(title = "GSE32529: CIRI 24h vs Sham 24h (Copper Death Genes)",
         x = expression(paste("log"[2], " Fold Change")),
         y = expression(paste("-log"[10], " P-value"))) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
  
  ggsave(file.path(output_dir, "GSE32529_copper_death_volcano_24h.png"), plot = p2, width = 8, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "GSE32529_copper_death_volcano_24h.pdf"), plot = p2, width = 8, height = 6)
}

cat("火山图已保存!\n")

# ==================== 18. 显示结果 ====================
cat("\n===== 铜死亡基因差异表达结果 =====\n")
print(copper_results[, c("Gene", "logFC_3h", "adj.P.Val_3h", "logFC_24h", "adj.P.Val_24h")], row.names = FALSE)

if (nrow(copper_sig) > 0) {
  cat("\n===== 显著差异的铜死亡基因 =====\n")
  print(copper_sig[, c("Gene", "logFC_3h", "adj.P.Val_3h", "logFC_24h", "adj.P.Val_24h")], row.names = FALSE)
} else {
  cat("\n===== 显著差异的铜死亡基因 =====\n")
  cat("无显著差异的铜死亡基因 (|logFC| > 1, FDR < 0.05)\n")
}

# ==================== 19. 总结报告 ====================
cat("\n===== 分析完成 =====\n")
cat(paste0("输出目录: ", output_dir, "\n"))
cat("\n生成文件:\n")
cat("1. GSE32529_copper_death_all_results.txt - 所有铜死亡基因结果\n")
cat("2. GSE32529_copper_death_DEGs_logFC1_FDR0.05.txt - 显著差异铜死亡基因\n")
cat("3. GSE32529_copper_death_3h_results.txt - 3h对比结果\n")
cat("4. GSE32529_copper_death_24h_results.txt - 24h对比结果\n")
cat("5. GSE32529_copper_death_expression_matrix.txt - 表达矩阵\n")
cat("6. GSE32529_copper_death_heatmap.png/pdf - 热图\n")
cat("7. GSE32529_copper_death_volcano_3h.png/pdf - 3h火山图\n")
cat("8. GSE32529_copper_death_volcano_24h.png/pdf - 24h火山图\n")
cat("9. unmapped_copper_genes.txt - 未找到探针的基因\n")
cat("\n分析参数:\n")
cat("- 数据集: GSE32529 (CIRI模型)\n")
cat("- 平台: GPL1261\n")
cat("- 样本: CIRI 3h (n=4), CIRI 24h (n=4), Sham 3h (n=4), Sham 24h (n=4)\n")
cat("- 筛选条件: |logFC| > 1, FDR < 0.05\n")
cat("- 铜死亡基因: ", paste(copper_genes, collapse = ", "), "\n")
cat("==================\n")
