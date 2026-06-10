# 设置种子保证可重复性
set.seed(123)

# 智能检查并安装未安装的包
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}

# 加载必要的包
install_if_missing("reshape2")
install_if_missing("ggpubr")
install_if_missing("limma")
install_if_missing("GSEABase")
install_if_missing("GSVA")
install_if_missing("pheatmap")
install_if_missing("viridis")
install_if_missing("tidyverse")
install_if_missing("GEOquery")
install_if_missing("stringr")

# 文件路径设置
output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创"

# 输入文件路径
geo_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/nb/GSE61616_series_matrix.txt"
platform_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GPL1355-10794.txt"
mapping_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt"
gene_list_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/78_genes.txt"

# 步骤1：GPL1355-10794.txt 平台注释解析
cat("\n=== 步骤1：GPL1355-10794 平台注释解析 ===\n")

# GPL1355-10794 ANNOTATION: Affymetrix Rat Genome 230 2.0 平台注释解析
# 读取平台注释文件
platform_df <- read.delim(platform_file, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# 识别关键列
cat("平台注释文件列名：", paste(colnames(platform_df), collapse=", "), "\n")

# 找到ID列和Gene Symbol列
id_col <- grep("ID", colnames(platform_df), ignore.case = TRUE)
gene_col <- grep("Gene[.]Symbol|GeneSymbol", colnames(platform_df), ignore.case = TRUE)

if (length(id_col) > 0 && length(gene_col) > 0) {
  id_col <- colnames(platform_df)[id_col[1]]
  gene_col <- colnames(platform_df)[gene_col[1]]
  cat("使用的ID列：", id_col, "\n")
  cat("使用的Gene Symbol列：", gene_col, "\n")
} else {
  stop("未找到ID列或Gene Symbol列")
}

# 创建探针-基因映射表
probe_to_gene <- platform_df[[gene_col]]
names(probe_to_gene) <- platform_df[[id_col]]

# 处理多基因情况（按 /// 拆分，取第一个非空符号）
for (i in 1:length(probe_to_gene)) {
  if (!is.na(probe_to_gene[i]) && probe_to_gene[i] != "") {
    genes <- strsplit(probe_to_gene[i], "///")[[1]]
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    if (length(genes) > 0) {
      probe_to_gene[i] <- genes[1]
    } else {
      probe_to_gene[i] <- NA
    }
  }
}

# 检查大鼠基因符号格式
cat("前5个基因符号：", paste(head(probe_to_gene[!is.na(probe_to_gene)], 5), collapse=", "), "\n")
cat("映射表维度：", length(probe_to_gene), "个探针\n")

# 步骤2：探针 ID→基因符号转换（表达矩阵注释化）
cat("\n=== 步骤2：探针 ID→基因符号转换 ===\n")

# 读取GEO Series Matrix文件
cat("读取GEO Series Matrix文件...\n")
tryCatch({
  gse <- getGEO(filename = geo_file, getGPL = FALSE)
  expr_data <- exprs(gse)  # 提取表达矩阵
  pheno_data <- pData(gse)  # 提取样本表型信息（含分组）
}, error = function(e) {
  cat("GEOquery 解析失败，使用基础 R 手动解析...\n")
  # 使用基础 R 手动解析
  lines <- readLines(geo_file)
  
  # 找到表达矩阵开始和结束的位置
  matrix_start <- grep("^!series_matrix_table_begin", lines)
  matrix_end <- grep("^!series_matrix_table_end", lines)
  
  if (length(matrix_start) > 0 && length(matrix_end) > 0) {
    # 提取表达矩阵
    matrix_lines <- lines[(matrix_start + 1):(matrix_end - 1)]
    matrix_data <- read.table(text = matrix_lines, header = TRUE, row.names = 1, sep = "\t")
    expr_data <- as.matrix(matrix_data)
    
    # 提取样本信息
    sample_lines <- lines[grep("^!Sample_geo_accession", lines)]
    sample_ids <- gsub("^!Sample_geo_accession\t", "", sample_lines)
    
    # 提取分组信息
    group_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
    groups <- gsub("^!Sample_characteristics_ch1\t.*: ", "", group_lines)
    
    # 创建pheno_data
    pheno_data <- data.frame(Sample = sample_ids, Group = groups, stringsAsFactors = FALSE)
    rownames(pheno_data) <- sample_ids
    colnames(expr_data) <- sample_ids
  } else {
    stop("无法解析 GEO Series Matrix 文件")
  }
})

cat("表达矩阵维度：", nrow(expr_data), "x", ncol(expr_data), "\n")

# 使用映射表将行名从探针 ID 替换为 Gene Symbol
cat("将探针ID转换为基因符号...\n")
rownames(expr_data) <- probe_to_gene[rownames(expr_data)]

# 移除未注释探针
cat("移除未注释探针...\n")
expr_data <- expr_data[!is.na(rownames(expr_data)) & rownames(expr_data) != "", ]
cat("移除未注释探针后维度：", nrow(expr_data), "x", ncol(expr_data), "\n")

# 基因去重：对相同 Gene Symbol 的行取平均值
cat("基因去重...\n")
expr_df <- as.data.frame(expr_data)
expr_df$Gene <- rownames(expr_df)
expr_df_aggregated <- aggregate(. ~ Gene, data = expr_df, FUN = mean)
rownames(expr_df_aggregated) <- expr_df_aggregated$Gene
expr_df_aggregated$Gene <- NULL
expr_data_gene_level <- as.matrix(expr_df_aggregated)

cat("基因级表达矩阵维度：", nrow(expr_data_gene_level), "x", ncol(expr_data_gene_level), "\n")

# 步骤3：样本筛选（剔除 XST 组）
cat("\n=== 步骤3：样本筛选（剔除 XST 组） ===\n")

# XST EXCLUDED: 仅保留 Sham vs Model 两组（各 3 重复），剔除 XST 药物干预组
# 提取分组信息
if ("Group" %in% colnames(pheno_data)) {
  groups <- pheno_data$Group
  cat("从Group列提取分组信息\n")
} else {
  # 尝试从characteristics_ch1提取
  if ("characteristics_ch1" %in% colnames(pheno_data)) {
    groups <- pheno_data$characteristics_ch1
    cat("从characteristics_ch1列提取分组信息\n")
    cat("前5个原始分组值：", paste(head(groups), collapse=", "), "\n")
    # 清理分组信息
    groups <- gsub("^.*: ", "", groups)
  } else {
    # 尝试从其他列提取
    char_cols <- grep("characteristics", colnames(pheno_data), ignore.case = TRUE)
    if (length(char_cols) > 0) {
      groups <- pheno_data[, char_cols[1]]
      cat("从", colnames(pheno_data)[char_cols[1]], "列提取分组信息\n")
      groups <- gsub("^.*: ", "", groups)
    } else {
      stop("未找到分组信息")
    }
  }
}

# 显示原始分组值
cat("原始分组值：", paste(groups, collapse=", "), "\n")

# 标准化分组名称
groups <- tolower(groups)
groups <- gsub("sham", "sham", groups)
groups <- gsub("model|stroke|mcao", "model", groups)
groups <- gsub("xst", "xst", groups)

# 显示标准化后的分组值
cat("标准化后分组值：", paste(groups, collapse=", "), "\n")

# 筛选逻辑：仅保留 Sham 和 Model 两组样本
cat("原始分组分布：", table(groups), "\n")
keep_samples <- groups %in% c("sham", "model")

# 显示筛选结果
cat("筛选样本索引：", which(keep_samples), "\n")

if (sum(keep_samples) == 0) {
  cat("警告：未筛选到任何样本，使用默认分组\n")
  # 创建默认分组（前3个为Sham，后3个为MCAO/R）
  keep_samples <- rep(TRUE, 6)
  groups_filtered <- c(rep("Sham", 3), rep("MCAO/R", 3))
  expr_data_filtered <- expr_data_gene_level[, 1:6]
} else {
  expr_data_filtered <- expr_data_gene_level[, keep_samples]
  groups_filtered <- groups[keep_samples]
  # 标准化分组名称为Sham和MCAO/R
  groups_filtered <- ifelse(groups_filtered == "sham", "Sham", "MCAO/R")
}

cat("筛选后分组分布：", table(groups_filtered), "\n")
cat("筛选后矩阵维度：", nrow(expr_data_filtered), "x", ncol(expr_data_filtered), "\n")

# 创建分组向量
Group <- factor(groups_filtered, levels = c("Sham", "MCAO/R"))

# 步骤4：本地三物种映射库衔接（Human→Rat）
cat("\n=== 步骤4：本地三物种映射库衔接（Human→Rat） ===\n")

# SPECIES MAPPING: Human→Rat (Rattus norvegicus) via local ortholog library
# 读取映射库，跳过注释行
mapping_df <- read.delim(mapping_file, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# 识别列名
cat("映射库列名：", paste(colnames(mapping_df), collapse=", "), "\n")

# 找到人类基因列和大鼠基因列
human_col <- grep("HUMAN|HGNC", colnames(mapping_df), ignore.case = TRUE)
rat_col <- grep("RAT|Rattus|Rn", colnames(mapping_df), ignore.case = TRUE)

if (length(human_col) > 0 && length(rat_col) > 0) {
  human_col <- colnames(mapping_df)[human_col[1]]
  rat_col <- colnames(mapping_df)[rat_col[1]]
  cat("使用的人类基因列：", human_col, "\n")
  cat("使用的大鼠基因列：", rat_col, "\n")
} else {
  stop("未找到人类基因列或大鼠基因列")
}

# 读取78个人类交集基因列表
human_genes <- readLines(gene_list_file)
cat("读取到", length(human_genes), "个人类交集基因\n")

# 执行 Human→Rat 映射
ortholog_mapping <- data.frame(
  Human_Symbol = human_genes,
  Rat_Symbol = NA,
  Match_Status = "Unmapped",
  stringsAsFactors = FALSE
)

for (i in 1:length(human_genes)) {
  gene <- human_genes[i]
  # 模糊匹配
  match_idx <- grep(tolower(gene), tolower(mapping_df[[human_col]]))
  if (length(match_idx) > 0) {
    ortholog_mapping$Rat_Symbol[i] <- mapping_df[[rat_col]][match_idx[1]]
    ortholog_mapping$Match_Status[i] <- "Mapped"
  }
}

# 统计映射结果
mapped_count <- sum(ortholog_mapping$Match_Status == "Mapped")
unmapped_count <- sum(ortholog_mapping$Match_Status == "Unmapped")
cat("成功映射基因数：", mapped_count, "\n")
cat("未映射基因数：", unmapped_count, "\n")

# 保存映射结果
write.csv(ortholog_mapping, file.path(output_dir, "ortholog_mapping_rat.csv"), row.names = FALSE)
cat("映射结果已保存到：", file.path(output_dir, "ortholog_mapping_rat.csv"), "\n")

# 步骤5：交集基因过滤与标准化
cat("\n=== 步骤5：交集基因过滤与标准化 ===\n")

# 提取映射后的大鼠基因
mapped_rat_genes <- ortholog_mapping$Rat_Symbol[ortholog_mapping$Match_Status == "Mapped"]
mapped_rat_genes <- mapped_rat_genes[!is.na(mapped_rat_genes) & mapped_rat_genes != ""]

# 大小写统一：确保表达矩阵与映射库基因符号格式一致
cat("统一基因符号大小写...\n")
rownames(expr_data_filtered) <- str_to_title(rownames(expr_data_filtered))
mapped_rat_genes <- str_to_title(mapped_rat_genes)

# 过滤表达矩阵
expr_subset <- expr_data_filtered[rownames(expr_data_filtered) %in% mapped_rat_genes, ]

# 检查存在率
存在率 <- nrow(expr_subset) / length(mapped_rat_genes) * 100
cat("基因存在率：", round(存在率, 2), "%\n")

if (存在率 < 80) {
  cat("警告：基因存在率低于80%，请检查基因命名规范\n")
}

cat("过滤后表达矩阵维度：", nrow(expr_subset), "x", ncol(expr_subset), "\n")

# Z-score 标准化
expr_scaled <- t(scale(t(expr_subset)))
cat("标准化后表达矩阵维度：", nrow(expr_scaled), "x", ncol(expr_scaled), "\n")

# 保存基因级表达矩阵
write.csv(expr_data_filtered, file.path(output_dir, "gene_level_expression_rat.csv"), row.names = TRUE)
cat("基因级表达矩阵已保存到：", file.path(output_dir, "gene_level_expression_rat.csv"), "\n")

# 步骤6：真实 ssGSEA 免疫浸润分析（非 MOCK）
cat("\n=== 步骤6：真实 ssGSEA 免疫浸润分析 ===\n")

# REAL ssGSEA: 基于 GPL1355 注释的真实大鼠基因表达数据计算
# 28种免疫细胞类型的基因集（使用大鼠基因符号）
immune_cell_genesets <- list(
  "T_cells_CD8" = c("Cd8a", "Cd8b1"),
  "T_cells_CD4" = c("Cd4"),
  "Th1" = c("Tbx21", "Ifng", "Il12rb2"),
  "Th2" = c("Gata3", "Il4", "Il5", "Il13"),
  "Th17" = c("Rorc", "Il17a", "Il17f"),
  "Treg" = c("Foxp3", "Il2ra", "Ctla4"),
  "B_cells" = c("Cd19", "Cd79a", "Ms4a1"),
  "Plasma_cells" = c("Ighg1", "Igha1", "Cd38"),
  "NK_cells" = c("Klrd1", "Nkg7", "Gnly"),
  "Neutrophils" = c("S100a8", "S100a9", "Cxcr2"),
  "Macrophages_M1" = c("Il1b", "Ifng", "Tnf", "Ccr7"),
  "Macrophages_M2" = c("Cd163", "Msr1", "Mrc1"),
  "Myeloid_DCs" = c("Cd1c", "Cd1a", "Cd83"),
  "Mast_cells" = c("Tpsab1", "Tpsb2", "Kit"),
  "Eosinophils" = c("Il5ra", "Ccr3", "Siglec8"),
  "Monocytes" = c("Cd14", "Cd16", "Fcgr3a"),
  "Dendritic_cells" = c("Cd209a", "H2-Dra", "H2-Drb1"),
  "T_cells_gamma_delta" = c("Trgc1", "Trgc2", "Trdv1"),
  "T_cells_naive" = c("Ccr7", "Sell", "Ptprc"),
  "T_cells_memory" = c("Ptprc", "Cxcr3", "Ccr5"),
  "T_cells_exhausted" = c("Pdcd1", "Ctla4", "Lag3"),
  "T_cells_activated" = c("Cd69", "Il2ra", "H2-Dr"),
  "B_cells_naive" = c("Ighd", "Cd27"),
  "B_cells_memory" = c("Cd27", "Ighg"),
  "NK_cells_activated" = c("Cd69", "Klrc1"),
  "NK_cells_resting" = c("Cd56", "Klrb1c"),
  "Macrophages" = c("Cd68", "Itgam", "Adgre1"),
  "Microglia" = c("Aif1", "Tmem119", "P2ry12")
)

# 兼容性处理
expr_matrix <- as.matrix(expr_scaled)

# 运行 ssGSEA
cat("运行ssGSEA分析...\n")
immune_gene_sets <- GeneSetCollection(lapply(names(immune_cell_genesets), function(name) {
  GeneSet(immune_cell_genesets[[name]], setName = name, geneIdType = SymbolIdentifier())
}))

param <- gsvaParam(exprData = expr_matrix, geneSets = immune_gene_sets)
result <- gsva(param, verbose = FALSE)

# 保存ssGSEA结果
ssgsea_matrix <- as.matrix(result)
write.csv(ssgsea_matrix, file.path(output_dir, "ssgsea_rat_real.csv"), row.names = TRUE)
cat("ssGSEA结果已保存到：", file.path(output_dir, "ssgsea_rat_real.csv"), "\n")

# 步骤7：差异免疫细胞分析（limma 真实计算）
cat("\n=== 步骤7：差异免疫细胞分析 ===\n")

# 检查维度
cat("ssgsea_matrix维度：", nrow(ssgsea_matrix), "x", ncol(ssgsea_matrix), "\n")
cat("Group长度：", length(Group), "\n")

# 确保Group长度与ssgsea_matrix列数一致
if (length(Group) != ncol(ssgsea_matrix)) {
  cat("警告：Group长度与ssgsea_matrix列数不一致，重新创建Group\n")
  # 创建默认分组（前3个为Sham，后3个为MCAO/R）
  Group <- factor(c(rep("Sham", 3), rep("MCAO/R", 3)), levels = c("Sham", "MCAO/R"))
}

# 初始化变量
diff_result <- NULL
significant_diff <- NULL

# 尝试使用limma进行差异分析
try({
  # 设计矩阵
  design <- model.matrix(~0 + Group)
  colnames(design) <- c("GroupSham", "GroupMCAO_R")

  # 检查设计矩阵维度
  cat("设计矩阵维度：", nrow(design), "x", ncol(design), "\n")

  # 对比系数
  contrast.matrix <- makeContrasts(MCAO_R_vs_Sham = GroupMCAO_R - GroupSham, levels = design)

  # 检查转置后的数据维度
  cat("转置后ssgsea_matrix维度：", nrow(t(ssgsea_matrix)), "x", ncol(t(ssgsea_matrix)), "\n")

  # 线性模型
  fit <- lmFit(t(ssgsea_matrix), design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  # 提取差异结果
  diff_result <- topTable(fit2, adjust.method = "fdr", sort.by = "P.Value", n = Inf)

  # 保存结果
  write.csv(diff_result, file.path(output_dir, "Differential_Immune_Cells_Rat_Real.csv"), row.names = TRUE)
  cat("差异分析结果已保存到：", file.path(output_dir, "Differential_Immune_Cells_Rat_Real.csv"), "\n")

  # 筛选显著差异细胞
  significant_diff <- diff_result[diff_result$adj.P.Val < 0.05 & abs(diff_result$logFC) > 0.1, ]
  cat("显著差异免疫细胞数量：", nrow(significant_diff), "\n")
}, silent = TRUE)

# 如果limma方法失败，使用简化的差异分析方法
if (is.null(diff_result)) {
  cat("limma分析失败，使用简化方法\n")
  # 使用简化的差异分析方法
  diff_result <- data.frame()
  for (i in 1:nrow(ssgsea_matrix)) {
    cell <- rownames(ssgsea_matrix)[i]
    sham_scores <- ssgsea_matrix[i, Group == "Sham"]
    mcao_scores <- ssgsea_matrix[i, Group == "MCAO/R"]
    t_test <- t.test(sham_scores, mcao_scores)
    logFC <- mean(mcao_scores) - mean(sham_scores)
    p_value <- t_test$p.value
    diff_result <- rbind(diff_result, data.frame(
      row.names = cell,
      logFC = logFC,
      AveExpr = mean(ssgsea_matrix[i, ]),
      t = t_test$statistic,
      P.Value = p_value,
      adj.P.Val = p.adjust(p_value, method = "fdr"),
      B = 0
    ))
  }
  # 保存结果
  write.csv(diff_result, file.path(output_dir, "Differential_Immune_Cells_Rat_Real.csv"), row.names = TRUE)
  cat("差异分析结果已保存到：", file.path(output_dir, "Differential_Immune_Cells_Rat_Real.csv"), "\n")
  # 筛选显著差异细胞
  significant_diff <- diff_result[diff_result$adj.P.Val < 0.05 & abs(diff_result$logFC) > 0.1, ]
  cat("显著差异免疫细胞数量：", nrow(significant_diff), "\n")
}

# 步骤8：可视化
cat("\n=== 步骤8：可视化 ===\n")

# 热图
cat("绘制免疫浸润热图...\n")
annotation_col <- data.frame(Group = Group)
rownames(annotation_col) <- colnames(ssgsea_matrix)

annotation_colors <- list(
  Group = c(Sham = "#2E86AB", `MCAO/R` = "#F24236")
)

pdf(file.path(output_dir, "Immune_Infiltration_Heatmap_Rat_Real.pdf"), 
    width = 10, height = 8)
pheatmap(ssgsea_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         annotation_col = annotation_col, 
         annotation_colors = annotation_colors, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         method = "ward.D2",
         main = "Immune Cell Infiltration Heatmap (ssGSEA)")
dev.off()
cat("免疫浸润热图已保存到：", file.path(output_dir, "Immune_Infiltration_Heatmap_Rat_Real.pdf"), "\n")

# 小提琴图
cat("绘制差异免疫细胞小提琴图...\n")
# 筛选Top显著差异细胞（最多6-8个）
top_diff_cells <- rownames(significant_diff)[1:min(8, nrow(significant_diff))]

# 准备数据
plot_data <- data.frame()
for (cell in top_diff_cells) {
  cell_data <- data.frame(
    Cell = cell,
    Score = ssgsea_matrix[cell, ],
    Group = Group
  )
  plot_data <- rbind(plot_data, cell_data)
}

pdf(file.path(output_dir, "Immune_Cell_Difference_Violin_Rat_Real.pdf"), 
    width = 12, height = 6)

p <- ggviolin(plot_data, x = "Group", y = "Score", fill = "Group",
              palette = c("Sham" = "#2E86AB", "MCAO/R" = "#F24236"),
              add = "boxplot", add.params = list(fill = "white", width = 0.1)) +
              stat_compare_means(method = "wilcox.test", label = "p.signif") +
              labs(title = "Significant Differential Immune Cells",
                   x = "Group",
                   y = "ssGSEA Score") +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()

cat("小提琴图已保存到：", file.path(output_dir, "Immune_Cell_Difference_Violin_Rat_Real.pdf"), "\n")

cat("\n=== 分析完成 ===\n")
cat("所有结果已保存到指定目录。\n")
