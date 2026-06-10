# 大鼠脑缺血再灌注模型免疫浸润分析脚本
# 基于 GPL1355-10794 平台完整注释与 Human→Rat 跨物种映射

# 加载必要的包
required_packages <- c("data.table", "GSVA", "GSEABase", "limma", "pheatmap", "ggpubr", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 输入文件路径
human_genes_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/78_genes.txt"
gpl_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GPL1355-10794.txt"
series_matrix_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/nb/GSE61616_series_matrix.txt"
mapping_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt"
output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/"

# 步骤 1：GPL1355-10794.txt 智能解析（40MB 大文件修复）
cat("=== 步骤 1：GPL1355-10794.txt 智能解析 ===\n")
# GPL1355 40MB PARSER: 动态检测分隔符与表头，确保读取31,099探针（修复前次57行错误）

# 读取前50行检测结构
lines <- readLines(gpl_file, n = 50)

# 检测分隔符
test_line <- lines[grep("^ID", lines)[1]]
if (grepl("\t", test_line)) {
  detected_sep <- "\t"
} else if (grepl(",", test_line)) {
  detected_sep <- ","
} else {
  detected_sep <- "\t"  # 默认使用制表符
}
cat("检测到分隔符：", detected_sep, "\n")

# 定位表头位置
header_line <- which(grepl("^ID\\s*\\t|^ID\\s*,", lines))[1]
if (is.na(header_line)) {
  # 尝试其他可能的表头
  header_line <- which(grepl("Probe Set ID|Gene Symbol", lines))[1]
}
cat("表头位置：第", header_line, "行\n")

# 读取完整注释表
cat("读取完整 GPL1355 注释文件...\n")
gpl_df <- fread(gpl_file, sep = detected_sep, skip = header_line - 1, stringsAsFactors = FALSE)

# 验证文件完整性
if (nrow(gpl_df) < 30000) {
  stop("错误：GPL1355 文件损坏，行数不足 30000")
}
cat("GPL1355 注释文件维度：", nrow(gpl_df), "x", ncol(gpl_df), "\n")

# 识别 ID 列和 Gene Symbol 列
id_col <- grep("ID|Probe Set ID", colnames(gpl_df), ignore.case = TRUE)[1]
gene_col <- grep("Gene Symbol|GeneSymbol", colnames(gpl_df), ignore.case = TRUE)[1]

if (is.na(id_col) || is.na(gene_col)) {
  stop("未找到 ID 列或 Gene Symbol 列")
}

cat("ID 列：", colnames(gpl_df)[id_col], "\n")
cat("Gene Symbol 列：", colnames(gpl_df)[gene_col], "\n")

# 提取探针-基因映射
probe_to_gene <- gpl_df[, c(id_col, gene_col), with = FALSE]
colnames(probe_to_gene) <- c("ID", "GeneSymbol")

# 处理多基因情况
probe_to_gene$GeneSymbol <- sapply(probe_to_gene$GeneSymbol, function(x) {
  if (is.na(x) || x == "") return(NA)
  symbols <- unlist(strsplit(x, "///"))
  symbols <- trimws(symbols)
  symbols <- symbols[symbols != ""]
  if (length(symbols) > 0) symbols[1] else NA
})

# 移除 NA 值
probe_to_gene <- probe_to_gene[!is.na(probe_to_gene$GeneSymbol), ]

# 创建映射向量
probe_to_gene_vec <- setNames(probe_to_gene$GeneSymbol, probe_to_gene$ID)
cat("探针-基因映射表维度：", nrow(probe_to_gene), "x", ncol(probe_to_gene), "\n")

# 检查前几个探针 ID
cat("GPL 文件前 5 个探针 ID：", paste(head(probe_to_gene$ID, 5), collapse = ", "), "\n")

# 步骤 2：读取 78 个交集基因并跨物种映射（Human→Rat）
cat("\n=== 步骤 2：Human→Rat 跨物种映射 ===\n")
# HUMAN→RAT MAPPING: 使用本地映射库RAT列，基于78_genes.txt路径

# 读取 78 个交集基因
human_genes <- readLines(human_genes_file)
human_genes <- trimws(human_genes)
human_genes <- human_genes[human_genes != ""]

if (length(human_genes) != 78) {
  warning("警告：78_genes.txt 文件中的基因数量不是 78，实际为：", length(human_genes))
}
cat("读取到人类基因数量：", length(human_genes), "\n")

# 读取本地映射库
mapping_df <- fread(mapping_file, sep = "\t", stringsAsFactors = FALSE, skip = "#")

# 识别列名
human_col <- grep("HUMAN_ORTHOLOG_SYMBOL", colnames(mapping_df), ignore.case = TRUE)[1]
rat_col <- grep("RAT_GENE_SYMBOL|RAT_ORTHOLOG_SYMBOL", colnames(mapping_df), ignore.case = TRUE)[1]

if (is.na(human_col) || is.na(rat_col)) {
  stop("未找到人类或大鼠基因列")
}

cat("人类基因列：", colnames(mapping_df)[human_col], "\n")
cat("大鼠基因列：", colnames(mapping_df)[rat_col], "\n")

# 执行映射
ortholog_mapping <- data.frame(
  Human_Symbol = human_genes,
  Rat_Symbol = NA,
  Match_Status = "No Match"
)

for (i in 1:length(human_genes)) {
  gene <- human_genes[i]
  match_idx <- which(toupper(mapping_df[[human_col]]) == toupper(gene))
  if (length(match_idx) > 0) {
    rat_gene <- mapping_df[[rat_col]][match_idx[1]]
    if (!is.na(rat_gene) && rat_gene != "") {
      ortholog_mapping$Rat_Symbol[i] <- rat_gene
      ortholog_mapping$Match_Status[i] <- "Match"
    }
  }
}

# 保存映射结果
output_file <- paste0(output_dir, "ortholog_mapping_rat.csv")
write.csv(ortholog_mapping, output_file, row.names = FALSE)
cat("映射成功数量：", sum(ortholog_mapping$Match_Status == "Match"), "/", length(human_genes), "\n")

if (sum(ortholog_mapping$Match_Status == "Match") < 70) {
  warning("警告：映射成功数量少于 70，可能影响后续分析")
}

# 步骤 3：表达矩阵处理（探针→基因→过滤）
cat("\n=== 步骤 3：表达矩阵处理 ===\n")

# 读取 GSE61616_series_matrix.txt
# 先读取所有行，然后过滤出表达数据
lines <- readLines(series_matrix_file)
# 找到数据开始的行（以引号和数字开头的行）
data_start <- which(grepl('^"[0-9]+_at"', lines))[1]
# 读取数据部分
data_lines <- lines[data_start:length(lines)]
# 移除可能的页脚
 data_lines <- data_lines[!grepl('^!series_matrix_table_end', data_lines)]
# 读取数据
gse_data <- fread(text = data_lines, sep = "\t", stringsAsFactors = FALSE)

# 提取表达矩阵
if (nrow(gse_data) > 0) {
  # 检查第一列的内容
  first_col <- gse_data[[1]]
  cat("第一列前 5 个值：", paste(head(first_col, 5), collapse = ", "), "\n")
  
  # 直接使用 data.frame 来处理行名
  gse_df <- as.data.frame(gse_data)
  rownames(gse_df) <- gse_df[[1]]
  expr_data <- gse_df[, -1]
  # 移除行名中的引号
  rownames(expr_data) <- gsub('"', '', rownames(expr_data))
  
  # 检查行名
  cat("表达矩阵前 5 个行名：", paste(head(rownames(expr_data), 5), collapse = ", "), "\n")
} else {
  stop("错误：无法读取表达数据")
}

# 转换为矩阵
expr_matrix <- as.matrix(expr_data)
cat("原始表达矩阵维度：", nrow(expr_matrix), "x", ncol(expr_matrix), "\n")

# 转换探针 ID 为基因符号
cat("转换探针 ID 为基因符号...\n")
# 获取探针 ID
clean_probe_ids <- rownames(expr_matrix)

# 检查前几个探针 ID
cat("表达矩阵前 5 个探针 ID：", paste(head(clean_probe_ids, 5), collapse = ", "), "\n")

# 检查探针 ID 匹配情况
matching_probes <- clean_probe_ids %in% names(probe_to_gene_vec)
cat("匹配的探针数量：", sum(matching_probes), "\n")

# 转换探针 ID 为基因符号
new_rownames <- probe_to_gene_vec[clean_probe_ids]

# 移除未注释探针
temp_matrix <- expr_matrix
rownames(temp_matrix) <- new_rownames
temp_matrix <- temp_matrix[!is.na(rownames(temp_matrix)), ]
temp_matrix <- temp_matrix[rownames(temp_matrix) != "", ]

cat("转换后表达矩阵维度：", nrow(temp_matrix), "x", ncol(temp_matrix), "\n")

# 决定使用哪个矩阵
if (nrow(temp_matrix) > 0) {
  expr_matrix <- temp_matrix
} else {
  cat("警告：探针 ID 转换失败，使用原始表达矩阵\n")
}

# 基因去重
cat("基因去重...\n")
if (nrow(expr_matrix) > 0) {
  expr_df <- as.data.frame(expr_matrix)
  expr_df$GeneSymbol <- rownames(expr_matrix)
  # 使用 aggregate 进行基因去重
  expr_gene_level <- aggregate(. ~ GeneSymbol, data = expr_df, FUN = mean)
  rownames(expr_gene_level) <- expr_gene_level$GeneSymbol
  expr_gene_level <- expr_gene_level[, -1]
  expr_gene_level <- as.matrix(expr_gene_level)
} else {
  stop("错误：转换后表达矩阵为空")
}

cat("基因级表达矩阵维度：", nrow(expr_gene_level), "x", ncol(expr_gene_level), "\n")

# 步骤 4：样本筛选（剔除 XST 组）
cat("\n=== 步骤 4：样本筛选 ===\n")
# XST EXCLUDED: 剔除药物组，仅Sham vs MCAO/R（6样本）

# 提取分组信息
gse_header <- readLines(series_matrix_file)

# 从 !Sample_title 行提取样本标题
sample_title_line <- gse_header[grep("^!Sample_title", gse_header)][1]
# 解析样本标题
if (!is.na(sample_title_line)) {
  # 移除前缀和引号
  sample_title_line <- sub("^!Sample_title\\s+", "", sample_title_line)
  # 分割样本标题
  sample_titles <- unlist(strsplit(sample_title_line, "\t"))
  # 移除引号
  sample_titles <- gsub('"', '', sample_titles)
  # 移除空字符串
  sample_titles <- sample_titles[sample_titles != ""]
  
  # 从样本标题中提取分组信息
  groups <- sapply(sample_titles, function(x) {
    x <- tolower(x)
    if (grepl("sham", x)) return("Sham")
    if (grepl("model", x)) return("MCAO_R")
    if (grepl("xst", x)) return("XST")
    return("unknown")
  })
} else {
  # 从表达矩阵列名中提取分组信息
  sample_names <- colnames(expr_gene_level)
  groups <- sapply(sample_names, function(x) {
    x <- tolower(x)
    if (grepl("sham", x)) return("Sham")
    if (grepl("model|mcao|stroke", x)) return("MCAO_R")
    if (grepl("xst", x)) return("XST")
    return("unknown")
  })
}

# 确保groups长度与表达矩阵列数一致
if (length(groups) != ncol(expr_gene_level)) {
  # 直接从表达矩阵列名提取
  sample_names <- colnames(expr_gene_level)
  groups <- sapply(sample_names, function(x) {
    x <- tolower(x)
    if (grepl("sham", x)) return("Sham")
    if (grepl("model|mcao|stroke", x)) return("MCAO_R")
    if (grepl("xst", x)) return("XST")
    return("unknown")
  })
}

# 筛选样本（剔除 XST 组，仅保留前3个Sham和前3个Model样本）
keep_samples <- groups != "XST" & groups != "unknown"
expr_data_filtered <- expr_gene_level[, keep_samples]
groups_filtered <- groups[keep_samples]

# 只保留前3个Sham和前3个Model样本
sham_indices <- which(groups_filtered == "Sham")
model_indices <- which(groups_filtered == "MCAO_R")

if (length(sham_indices) >= 3 && length(model_indices) >= 3) {
  selected_indices <- c(sham_indices[1:3], model_indices[1:3])
  expr_data_filtered <- expr_data_filtered[, selected_indices]
  groups_filtered <- groups_filtered[selected_indices]
} else {
  # 如果样本数量不足，使用所有可用样本
  cat("警告：样本数量不足，使用所有可用样本\n")
}

groups_filtered <- factor(groups_filtered)

cat("筛选后样本数量：", ncol(expr_data_filtered), "\n")
cat("筛选后矩阵维度：", nrow(expr_data_filtered), "x", ncol(expr_data_filtered), "\n")
cat("分组情况：", table(groups_filtered), "\n")

# 步骤 5：交集基因过滤与标准化
cat("\n=== 步骤 5：交集基因过滤与标准化 ===\n")

# 提取映射成功的大鼠基因
mapped_rat_genes <- ortholog_mapping$Rat_Symbol[ortholog_mapping$Match_Status == "Match"]
mapped_rat_genes <- mapped_rat_genes[!is.na(mapped_rat_genes)]

# 大小写统一
expr_rownames <- rownames(expr_data_filtered)
expr_rownames_std <- stringr::str_to_title(expr_rownames)
mapped_rat_genes_std <- stringr::str_to_title(mapped_rat_genes)

# 过滤表达矩阵
expr_subset <- expr_data_filtered[expr_rownames_std %in% mapped_rat_genes_std, ]

# 检查存在率
 existence_rate <- nrow(expr_subset) / length(mapped_rat_genes) * 100
cat("基因存在率：", round(existence_rate, 2), "%\n")
cat("过滤后矩阵维度：", nrow(expr_subset), "x", ncol(expr_subset), "\n")

if (nrow(expr_subset) < 60) {
  warning("警告：匹配基因数量少于 60，尝试其他大小写转换")
  # 尝试全小写匹配
  expr_rownames_lower <- tolower(expr_rownames)
  mapped_rat_genes_lower <- tolower(mapped_rat_genes)
  expr_subset <- expr_data_filtered[expr_rownames_lower %in% mapped_rat_genes_lower, ]
  cat("全小写匹配后维度：", nrow(expr_subset), "x", ncol(expr_subset), "\n")
}

# Z-score 标准化
expr_scaled <- t(scale(t(expr_subset)))

# 保存过滤后的表达矩阵
output_file <- paste0(output_dir, "gene_level_expression_78genes_6samples.csv")
write.csv(expr_subset, output_file, row.names = TRUE)

# 步骤 6：真实 ssGSEA 分析（非模拟）
cat("\n=== 步骤 6：真实 ssGSEA 分析 ===\n")
# REAL ssGSEA: 基于78个映射大鼠基因的真实表达计算，非模拟

# 28 种免疫细胞基因集（大鼠基因符号）
immune_cell_genesets <- list(
  B.cells = c("Cd19", "Cd79a", "Ms4a1"),
  T.cells = c("Cd3e", "Cd3d", "Cd3g"),
  CD4.T.cells = c("Cd4", "Foxp3"),
  CD8.T.cells = c("Cd8a", "Cd8b1"),
  NK.cells = c("Ncr1", "Klrb1c", "Klrb1b"),
  Monocytes = c("Cd14", "Ly6c1", "Ly6c2"),
  Macrophages = c("Cd68", "Adgre1", "Mrc1"),
  Macrophages.M1 = c("Il6", "Tnf", "Nos2"),
  Macrophages.M2 = c("Arg1", "Mrc1", "Cd206"),
  Neutrophils = c("Ly6g", "S100a8", "S100a9"),
  Dendritic.cells = c("Cd11c", "Itgax", "H2-Aa"),
  Mast.cells = c("Cpa3", "Tpsb2", "Kit"),
  Eosinophils = c("Siglecf", "Prg2", "Ear1"),
  Basophils = c("Mcpt8", "Cd203c", "Fcer1a"),
  Endothelial.cells = c("Cdh5", "Vwf", "Pecam1"),
  Fibroblasts = c("Col1a1", "Col3a1", "Fap"),
  Astrocytes = c("Gfap", "Aqp4", "Aldh1l1"),
  Microglia = c("Iba1", "Aif1", "Tmem119"),
  Oligodendrocytes = c("Mbp", "Plp1", "Olig2"),
  Neurons = c("Snap25", "Syn1", "Map2"),
  Tregs = c("Foxp3", "Il2ra", "Ctla4"),
  Th1.cells = c("Tbx21", "Ifng", "Stat4"),
  Th2.cells = c("Gata3", "Il4", "Il13"),
  Th17.cells = c("Rorc", "Il17a", "Il17f"),
  Tfh.cells = c("Bcl6", "Cxcr5", "Il21"),
  gamma.delta.T.cells = c("Trdc", "Trdv2", "Trdv3"),
  NKT.cells = c("Cd1d1", "Klrk1", "Cd3e"),
  Plasma.cells = c("Ighg1", "Jchain", "Sdc1")
)

# 过滤基因集，只保留在表达矩阵中存在的基因
filtered_genesets <- lapply(immune_cell_genesets, function(genes) {
  genes[genes %in% rownames(expr_scaled)]
})

# 移除空基因集
filtered_genesets <- filtered_genesets[sapply(filtered_genesets, length) > 0]

if (length(filtered_genesets) == 0) {
  stop("错误：没有找到匹配的免疫细胞基因集")
}

# 运行 ssGSEA
expr_matrix <- as.matrix(expr_scaled)
# 创建ExpressionSet对象
library(Biobase)
expr_eset <- ExpressionSet(assayData = expr_matrix)
# 创建gsvaParam对象
param <- gsvaParam(exprData = expr_eset, geneSets = filtered_genesets)
# 运行GSVA
result <- gsva(param)
ssgsea_matrix <- as.matrix(result)

# 保存 ssGSEA 结果
output_file <- paste0(output_dir, "ssgsea_rat_real.csv")
write.csv(ssgsea_matrix, output_file, row.names = TRUE)
cat("ssGSEA 分析完成，保存结果到", output_file, "\n")

# 步骤 7：差异分析与可视化
cat("\n=== 步骤 7：差异分析与可视化 ===\n")

# 设计矩阵
design <- model.matrix(~0 + groups_filtered)
colnames(design) <- levels(groups_filtered)

# 添加调试信息
cat("ssgsea_matrix 维度:", dim(ssgsea_matrix), "\n")
cat("t(ssgsea_matrix) 维度:", dim(t(ssgsea_matrix)), "\n")
cat("design 矩阵维度:", dim(design), "\n")
cat("groups_filtered 长度:", length(groups_filtered), "\n")
cat("groups_filtered:", groups_filtered, "\n")

# 对比系数
contrast.matrix <- makeContrasts(MCAO_R_vs_Sham = MCAO_R - Sham, levels = design)

# 线性模型
fit <- lmFit(ssgsea_matrix, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 提取结果
diff_result <- topTable(fit2, adjust.method = "fdr", sort.by = "P", n = Inf)

# 保存差异分析结果
output_file <- paste0(output_dir, "Differential_Immune_Cells_Rat_Real.csv")
write.csv(diff_result, output_file, row.names = TRUE)
cat("差异分析完成，保存结果到", output_file, "\n")

# 筛选显著差异细胞
significant_cells <- rownames(diff_result[diff_result$adj.P.Val < 0.05 & abs(diff_result$logFC) > 0.1, ])
cat("显著差异免疫细胞数量：", length(significant_cells), "\n")

# 可视化：热图
output_file <- paste0(output_dir, "Immune_Infiltration_Heatmap_Rat_Real.pdf")
pdf(output_file, width = 10, height = 8, pointsize = 12)
try({
  annotation_col <- data.frame(Group = groups_filtered)
  rownames(annotation_col) <- colnames(ssgsea_matrix)
  annotation_colors <- list(Group = c(Sham = "#2E86AB", MCAO_R = "#F24236"))
  pheatmap(ssgsea_matrix, 
           annotation_col = annotation_col, 
           annotation_colors = annotation_colors,
           cluster_cols = FALSE, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           main = "Immune Infiltration Heatmap (Rat)")
}, silent = TRUE)
dev.off()

# 可视化：小提琴图（Top 6 显著差异细胞）
top_cells <- head(rownames(diff_result), 6)
if (length(top_cells) > 0) {
  output_file <- paste0(output_dir, "Immune_Cell_Difference_Violin_Rat_Real.pdf")
  pdf(output_file, width = 12, height = 6, pointsize = 12)
  try({
    # 准备数据
    plot_data <- data.frame()
    for (cell in top_cells) {
      if (cell %in% rownames(ssgsea_matrix)) {
        cell_data <- data.frame(
          Cell = cell,
          Score = as.numeric(ssgsea_matrix[cell, ]),
          Group = as.character(groups_filtered)
        )
        plot_data <- rbind(plot_data, cell_data)
      }
    }
    
    if (nrow(plot_data) > 0) {
      p <- ggviolin(plot_data, x = "Group", y = "Score", fill = "Group", 
                    palette = c("#2E86AB", "#F24236"),
                    add = "boxplot", add.params = list(fill = "white"))
      p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif")
      p <- p + facet_wrap(~Cell, scales = "free_y", nrow = 2)
      p <- p + ggtitle("Top 6 Significant Immune Cell Differences")
      print(p)
    }
  }, silent = TRUE)
  dev.off()
}

cat("\n=== 分析完成 ===\n")
cat("输出文件：\n")
cat("1. ortholog_mapping_rat.csv\n")
cat("2. gene_level_expression_78genes_6samples.csv\n")
cat("3. ssgsea_rat_real.csv\n")
cat("4. Differential_Immune_Cells_Rat_Real.csv\n")
cat("5. Immune_Infiltration_Heatmap_Rat_Real.pdf\n")
cat("6. Immune_Cell_Difference_Violin_Rat_Real.pdf\n")
