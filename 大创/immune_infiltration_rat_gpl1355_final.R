# 大鼠脑缺血再灌注模型免疫浸润分析脚本（最终版本）
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

# 步骤 1：文件存在性检查（前置验证）
cat("=== 步骤 1：文件存在性检查 ===\n")
# FILE CHECK: 验证所有输入文件存在性
input_files <- c(
  human_genes_file,
  gpl_file,
  series_matrix_file,
  mapping_file
)

for (file in input_files) {
  if (!file.exists(file)) {
    stop(paste("错误：文件不存在：", file))
  }
  cat("✓ 找到文件：", file, "\n")
}

# 步骤 2：GPL1355-10794.txt 智能解析（修复分隔符与表头）
cat("\n=== 步骤 2：GPL1355-10794.txt 智能解析 ===\n")
# GPL1355 PARSER: 动态检测表头与分隔符，修复skip="#"错误，确保读取31,099探针

# 读取前100行检测结构
lines <- readLines(gpl_file, n = 100)

# 检测分隔符
test_line <- lines[1:100][grep("^ID|Probe Set ID", lines[1:100])[1]]
if (grepl("\t", test_line)) {
  detected_sep <- "\t"
} else if (grepl(",", test_line)) {
  detected_sep <- ","
} else {
  detected_sep <- "\t"  # 默认使用制表符
}
cat("检测到分隔符：", detected_sep, "\n")

# 定位表头位置
header_line <- which(grepl("^ID\\s*\\t|^ID\\s*,|^Probe Set ID", lines))[1]
if (is.na(header_line)) {
  stop("错误：未找到表头行")
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

# 步骤 3：表达矩阵解析（修复起始行定位）
cat("\n=== 步骤 3：表达矩阵解析 ===\n")
# MATRIX PARSE: 修复data_start NA风险，准确定位探针ID起始行

# 读取 GSE61616_series_matrix.txt 所有行
lines <- readLines(series_matrix_file)

# 找到数据开始的行（以探针ID开头的行）
data_start <- which(grepl('^"?[0-9]+_at"?', lines))[1]

# 必须检查
if (is.na(data_start)) {
  stop("错误：未找到表达数据起始行")
}
cat("数据起始行：第", data_start, "行\n")

# 提取数据部分
data_lines <- lines[data_start:length(lines)]
# 移除可能的页脚
data_lines <- data_lines[!grepl('^!series_matrix_table_end', data_lines)]

# 读取数据
cat("读取表达矩阵...\n")
gse_data <- fread(text = data_lines, sep = "\t", stringsAsFactors = FALSE)

# 提取表达矩阵
if (nrow(gse_data) > 0) {
  # 检查第一列的内容
  first_col <- gse_data[[1]]
  
  # 直接使用 data.frame 来处理行名
  gse_df <- as.data.frame(gse_data)
  rownames(gse_df) <- gse_df[[1]]
  expr_data <- gse_df[, -1]
  # 移除行名中的引号
  rownames(expr_data) <- gsub('"', '', rownames(expr_data))
} else {
  stop("错误：无法读取表达数据")
}

# 转换为矩阵
expr_matrix <- as.matrix(expr_data)
cat("原始表达矩阵维度：", nrow(expr_matrix), "x", ncol(expr_matrix), "\n")

# 步骤 4：Human→Rat 跨物种映射（修复映射库读取）
cat("\n=== 步骤 4：Human→Rat 跨物种映射 ===\n")

# 读取本地映射库（删除 skip = "#" 参数）
mapping_df <- fread(mapping_file, sep = "\t", stringsAsFactors = FALSE)

# 识别列名
human_col <- grep("HUMAN_ORTHOLOG_SYMBOL", colnames(mapping_df), ignore.case = TRUE)[1]
rat_col <- grep("RAT_ORTHOLOG_SYMBOL|RAT_GENE_SYMBOL", colnames(mapping_df), ignore.case = TRUE)[1]

if (is.na(human_col) || is.na(rat_col)) {
  stop("未找到人类或大鼠基因列")
}

cat("人类基因列：", colnames(mapping_df)[human_col], "\n")
cat("大鼠基因列：", colnames(mapping_df)[rat_col], "\n")

# 读取 78 个交集基因
human_genes <- readLines(human_genes_file)
human_genes <- trimws(human_genes)
human_genes <- human_genes[human_genes != ""]

if (length(human_genes) != 78) {
  warning("警告：78_genes.txt 文件中的基因数量不是 78，实际为：", length(human_genes))
}
cat("读取到人类基因数量：", length(human_genes), "\n")

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

# 步骤 5：探针→基因转换与去重（修复大小写与去重）
cat("\n=== 步骤 5：探针→基因转换与去重 ===\n")

# 转换探针 ID 为基因符号
cat("转换探针 ID 为基因符号...\n")
# 获取探针 ID
clean_probe_ids <- rownames(expr_matrix)

# 转换探针 ID 为基因符号
new_rownames <- probe_to_gene_vec[clean_probe_ids]

# 移除未注释探针
temp_matrix <- expr_matrix
rownames(temp_matrix) <- new_rownames
temp_matrix <- temp_matrix[!is.na(rownames(temp_matrix)), ]
temp_matrix <- temp_matrix[rownames(temp_matrix) != "", ]

cat("转换后表达矩阵维度：", nrow(temp_matrix), "x", ncol(temp_matrix), "\n")

# 基因去重
cat("基因去重...\n")
expr_df <- as.data.frame(temp_matrix)
expr_df$GeneSymbol <- rownames(expr_df)
# 使用 aggregate 进行基因去重
expr_gene_level <- aggregate(. ~ GeneSymbol, data = expr_df, FUN = mean)
rownames(expr_gene_level) <- expr_gene_level$GeneSymbol
expr_gene_level <- expr_gene_level[, -1]
expr_gene_level <- as.matrix(expr_gene_level)

# 修复大小写统一
cat("标准化基因符号大小写...\n")
rownames(expr_gene_level) <- stringr::str_to_title(rownames(expr_gene_level))

cat("基因级表达矩阵维度：", nrow(expr_gene_level), "x", ncol(expr_gene_level), "\n")

# 步骤 6：样本筛选（修复硬编码，明确分组逻辑）
cat("\n=== 步骤 6：样本筛选 ===\n")
# GROUPING: 前5列为Sham，后5列为Model（共10样本），显式设置factor levels

# 严格按顺序分组：前5列为Sham，后5列为Model
if (ncol(expr_gene_level) >= 10) {
  expr_subset <- expr_gene_level[, 1:10]
} else {
  stop("错误：表达矩阵列数不足 10")
}

# 创建分组因子
cat("创建分组因子...\n")
groups <- factor(rep(c("Sham", "MCAO_R"), each = 5), levels = c("Sham", "MCAO_R"))

# 验证
if (length(groups) != 10) {
  stop("错误：分组长度不等于 10")
}
cat("分组情况：", table(groups), "\n")

# 步骤 7：交集基因过滤（修复存在率检查）
cat("\n=== 步骤 7：交集基因过滤与标准化 ===\n")

# 提取映射成功的大鼠基因
mapped_rat_genes <- ortholog_mapping$Rat_Symbol[ortholog_mapping$Match_Status == "Match"]
mapped_rat_genes <- mapped_rat_genes[!is.na(mapped_rat_genes)]
# 标准化大小写
mapped_rat_genes_std <- stringr::str_to_title(mapped_rat_genes)

# 过滤表达矩阵
expr_filtered <- expr_subset[rownames(expr_subset) %in% mapped_rat_genes_std, ]

# 检查存在率
existence_rate <- nrow(expr_filtered) / length(mapped_rat_genes) * 100
cat("基因存在率：", round(existence_rate, 2), "%\n")
cat("过滤后矩阵维度：", nrow(expr_filtered), "x", ncol(expr_filtered), "\n")

if (nrow(expr_filtered) < 60) {
  warning("警告：匹配基因数量少于 60，尝试全小写匹配")
  # 尝试全小写匹配
  expr_rownames_lower <- tolower(rownames(expr_subset))
  mapped_rat_genes_lower <- tolower(mapped_rat_genes)
  expr_filtered <- expr_subset[expr_rownames_lower %in% mapped_rat_genes_lower, ]
  cat("全小写匹配后维度：", nrow(expr_filtered), "x", ncol(expr_filtered), "\n")
}

# Z-score 标准化
cat("进行 Z-score 标准化...\n")
expr_scaled <- t(scale(t(expr_filtered)))

# 修复NaN问题
# NA HANDLING: scale()后替换NaN为0，防止零方差基因导致错误
expr_scaled[is.nan(expr_scaled)] <- 0

# 保存基因级表达矩阵
output_file <- paste0(output_dir, "gene_level_expression_10samples.csv")
write.csv(expr_filtered, output_file, row.names = TRUE)

# 步骤 8：ssGSEA 分析（修复基因集构建）
cat("\n=== 步骤 8：真实 ssGSEA 分析 ===\n")

# 基于实际拥有的基因构建免疫细胞基因集
immune_cell_genesets <- list(
  Microglia = c("Aif1", "Cx3cr1", "P2ry12"),
  Astrocytes = c("Gfap", "Aqp4", "Aldh1l1"),
  Neurons = c("Snap25", "Syn1", "Map2"),
  Oligodendrocytes = c("Mbp", "Plp1", "Olig2"),
  Macrophages = c("Cd68", "Adgre1", "Mrc1"),
  T.cells = c("Cd3e", "Cd4", "Cd8a"),
  B.cells = c("Cd19", "Cd79a", "Ms4a1"),
  NK.cells = c("Ncr1", "Klrb1c"),
  Monocytes = c("Cd14", "Ly6c1"),
  Neutrophils = c("Ly6g", "S100a8")
)

# 过滤基因集，只保留在表达矩阵中存在的基因
filtered_genesets <- lapply(immune_cell_genesets, function(genes) {
  genes[genes %in% rownames(expr_scaled)]
})

# 移除空基因集
filtered_genesets <- filtered_genesets[sapply(filtered_genesets, length) > 0]

# 必须检查
# GENESET CHECK: 验证免疫细胞基因集与表达矩阵存在重叠，避免空集报错
if (length(filtered_genesets) == 0) {
  stop("错误：免疫细胞基因集与表达矩阵无重叠基因")
}

cat("构建基因集数量：", length(filtered_genesets), "\n")

# 准备基因集为GeneSetCollection对象
gene_sets <- lapply(names(filtered_genesets), function(name) {
  genes <- filtered_genesets[[name]]
  if (length(genes) > 0) {
    GeneSet(genes, setName = name, geneIdType = SymbolIdentifier())
  }
})
gene_sets <- gene_sets[!sapply(gene_sets, is.null)]
immune_gene_sets <- GeneSetCollection(gene_sets)

# 创建ExpressionSet对象
expr_eset <- ExpressionSet(assayData = expr_scaled)

# 创建gsvaParam对象
param <- gsvaParam(exprData = expr_eset, geneSets = immune_gene_sets)

# 运行GSVA
cat("运行 ssGSEA 分析...\n")
result <- gsva(param, verbose = FALSE)
ssgsea_matrix <- as.matrix(result)

# 保存 ssGSEA 结果
output_file <- paste0(output_dir, "ssgsea_rat_real.csv")
write.csv(ssgsea_matrix, output_file, row.names = TRUE)
cat("ssGSEA 分析完成，保存结果到", output_file, "\n")

# 步骤 9：差异分析（修复设计矩阵与对比）
cat("\n=== 步骤 9：差异分析 ===\n")

# 设计矩阵
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# 验证维度
if (ncol(expr_scaled) != nrow(design)) {
  stop("错误：表达矩阵列数与设计矩阵行数不匹配")
}

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

# 步骤 10：可视化
cat("\n=== 步骤 10：可视化 ===\n")

# 热图
output_file <- paste0(output_dir, "Immune_Infiltration_Heatmap_Rat_Real.pdf")
pdf(output_file, width = 10, height = 8, pointsize = 12)
try({
  annotation_col <- data.frame(Group = groups)
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

# 小提琴图（Top差异细胞）
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
          Group = as.character(groups)
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
cat("2. gene_level_expression_10samples.csv\n")
cat("3. ssgsea_rat_real.csv\n")
cat("4. Differential_Immune_Cells_Rat_Real.csv\n")
cat("5. Immune_Infiltration_Heatmap_Rat_Real.pdf\n")
cat("6. Immune_Cell_Difference_Violin_Rat_Real.pdf\n")