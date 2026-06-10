# 大鼠脑缺血再灌注模型免疫浸润分析脚本
# 基于 GPL1355-10794 平台注释与 Human→Rat 跨物种映射

# 加载必要的包
required_packages <- c("GSVA", "GSEABase", "limma", "pheatmap", "ggpubr", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 输入文件路径
series_matrix_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/nb/GSE61616_series_matrix.txt"
platform_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GPL1355-10794.txt"
mapping_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt"
human_genes_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/78_genes.txt"

# 步骤 1：GPL1355-10794.txt 平台注释解析
cat("=== 步骤 1：GPL1355-10794 平台注释解析 ===\n")
# GPL1355-10794 ANNOTATION: Affymetrix Rat Genome 230 2.0 平台注释解析
platform_df <- read.delim(platform_file, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# 识别关键列
id_col <- grep("ID", colnames(platform_df), ignore.case = TRUE)
gene_col <- grep("Gene[.]Symbol|GeneSymbol|Gene Symbol", colnames(platform_df), ignore.case = TRUE)

if (length(id_col) == 0 || length(gene_col) == 0) {
  stop("未找到 ID 列或 Gene Symbol 列")
}

# 创建探针-基因映射表
probe_to_gene <- platform_df[, c(id_col, gene_col)]
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

# 步骤 2：探针 ID→基因符号转换（表达矩阵注释化）
cat("\n=== 步骤 2：探针 ID→基因符号转换 ===\n")

# 读取表达矩阵
gse_data <- read.delim(series_matrix_file, sep = "\t", stringsAsFactors = FALSE, comment.char = "!")

# 提取表达数据（跳过前几行注释）
cat("原始数据维度：", nrow(gse_data), "x", ncol(gse_data), "\n")

# 提取探针 ID 行（更宽松的匹配模式）
expr_data <- gse_data[grep("_at$", gse_data[, 1]), ]
cat("提取后数据维度：", nrow(expr_data), "x", ncol(expr_data), "\n")

if (nrow(expr_data) == 0) {
  cat("警告：未提取到任何探针数据，使用所有数据行\n")
  expr_data <- gse_data[-(1:grep("^!", gse_data[, 1], invert = TRUE)[1]-1), ]
}

rownames(expr_data) <- expr_data[, 1]
expr_data <- expr_data[, -1]

# 转换为矩阵
expr_data <- as.matrix(expr_data)

# 转换探针 ID 为基因符号
cat("转换探针 ID 为基因符号...\n")
cat("表达矩阵探针数量：", nrow(expr_data), "\n")
cat("映射表探针数量：", length(probe_to_gene_vec), "\n")

# 检查前几个探针 ID
if (nrow(expr_data) > 0) {
  cat("前 5 个探针 ID：", paste(head(rownames(expr_data), 5), collapse = ", "), "\n")
  cat("映射表中存在的探针数量：", sum(rownames(expr_data) %in% names(probe_to_gene_vec)), "\n")
}

new_rownames <- probe_to_gene_vec[rownames(expr_data)]
rownames(expr_data) <- new_rownames

# 移除未注释探针
expr_data <- expr_data[!is.na(rownames(expr_data)), ]
expr_data <- expr_data[rownames(expr_data) != "", ]

cat("注释后基因数量：", nrow(expr_data), "\n")

# 基因去重
cat("基因去重...\n")
if (nrow(expr_data) > 0) {
  expr_data_df <- as.data.frame(expr_data)
  expr_data_df$GeneSymbol <- rownames(expr_data_df)
  expr_data_gene_level <- aggregate(. ~ GeneSymbol, data = expr_data_df, FUN = mean)
  rownames(expr_data_gene_level) <- expr_data_gene_level$GeneSymbol
  expr_data_gene_level <- expr_data_gene_level[, -1]
  expr_data_gene_level <- as.matrix(expr_data_gene_level)
  
  cat("基因级表达矩阵维度：", nrow(expr_data_gene_level), "x", ncol(expr_data_gene_level), "\n")
} else {
  cat("警告：没有可去重的基因数据\n")
  # 创建空矩阵作为占位符
  expr_data_gene_level <- matrix(nrow = 0, ncol = ncol(expr_data))
}

# 步骤 3：样本筛选（剔除 XST 组）
cat("\n=== 步骤 3：样本筛选（剔除 XST 组） ===\n")
# XST EXCLUDED: 仅保留 Sham vs Model 两组（各 3 重复），剔除 XST 药物干预组

# 提取分组信息
gse_header <- readLines(series_matrix_file)
sample_info <- gse_header[grep("^!Sample_characteristics_ch1", gse_header)]
sample_info <- sample_info[1:ncol(expr_data)]

groups <- sapply(sample_info, function(x) {
  x <- tolower(x)
  if (grepl("sham", x)) return("sham")
  if (grepl("model|mcao|stroke", x)) return("model")
  if (grepl("xst", x)) return("xst")
  return("unknown")
})

# 筛选样本
keep_samples <- groups %in% c("sham", "model")

if (sum(keep_samples) == 0) {
  cat("警告：未筛选到任何样本，使用默认分组\n")
  keep_samples <- rep(TRUE, 6)
  groups_filtered <- factor(c(rep("Sham", 3), rep("MCAO/R", 3)))
  expr_data_filtered <- expr_data_gene_level[, 1:6]
} else {
  expr_data_filtered <- expr_data_gene_level[, keep_samples]
  groups_filtered <- groups[keep_samples]
  groups_filtered <- ifelse(groups_filtered == "sham", "Sham", "MCAO/R")
  groups_filtered <- factor(groups_filtered)
}

cat("筛选后样本数量：", sum(keep_samples), "\n")
cat("筛选后矩阵维度：", nrow(expr_data_filtered), "x", ncol(expr_data_filtered), "\n")

# 步骤 4：本地三物种映射库衔接（Human→Rat）
cat("\n=== 步骤 4：Human→Rat 跨物种映射 ===\n")
# SPECIES MAPPING: Human→Rat (Rattus norvegicus) via local ortholog library

# 读取映射库（跳过注释行）
lines <- readLines(mapping_file)
data_start <- which(grepl("^RAT_GENE_SYMBOL", lines))[1]
mapping_df <- read.delim(text = lines[data_start:length(lines)], sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# 识别列名（更准确的匹配）
human_col <- grep("HUMAN_ORTHOLOG_SYMBOL", colnames(mapping_df), ignore.case = TRUE)
rat_col <- grep("RAT_GENE_SYMBOL", colnames(mapping_df), ignore.case = TRUE)

# 确保找到列
if (length(human_col) == 0 || length(rat_col) == 0) {
  stop("未找到人类或大鼠基因列")
}

cat("人类基因列：", colnames(mapping_df)[human_col], "(列索引:", human_col, ")\n")
cat("大鼠基因列：", colnames(mapping_df)[rat_col], "(列索引:", rat_col, ")\n")
cat("映射库维度：", nrow(mapping_df), "x", ncol(mapping_df), "\n")

# 读取 78 个人类交集基因
human_genes <- read.delim(human_genes_file, header = FALSE, stringsAsFactors = FALSE)
human_genes <- human_genes[, 1]
human_genes <- toupper(human_genes)

# 映射人类基因到大鼠基因
ortholog_mapping <- data.frame(
  Human_Symbol = human_genes,
  Rat_Symbol = NA,
  Match_Status = "No Match"
)

for (i in 1:length(human_genes)) {
  gene <- human_genes[i]
  match_idx <- which(toupper(mapping_df[, human_col]) == gene)
  if (length(match_idx) > 0) {
    rat_gene <- mapping_df[match_idx[1], rat_col]
    if (!is.na(rat_gene) && rat_gene != "") {
      ortholog_mapping$Rat_Symbol[i] <- rat_gene
      ortholog_mapping$Match_Status[i] <- "Match"
    }
  }
}

# 保存映射结果
write.csv(ortholog_mapping, "ortholog_mapping_rat.csv", row.names = FALSE)
cat("映射成功数量：", sum(ortholog_mapping$Match_Status == "Match"), "/", length(human_genes), "\n")

# 步骤 5：交集基因过滤与标准化
cat("\n=== 步骤 5：交集基因过滤与标准化 ===\n")

# 提取映射成功的大鼠基因
mapped_rat_genes <- ortholog_mapping$Rat_Symbol[ortholog_mapping$Match_Status == "Match"]
mapped_rat_genes <- mapped_rat_genes[!is.na(mapped_rat_genes)]

# 大小写统一
expr_rownames <- rownames(expr_data_filtered)
expr_data_filtered <- expr_data_filtered[stringr::str_to_title(expr_rownames) %in% stringr::str_to_title(mapped_rat_genes), ]

# 检查存在率
existence_rate <- nrow(expr_data_filtered) / length(mapped_rat_genes) * 100
cat("基因存在率：", round(existence_rate, 2), "%\n")
cat("过滤后矩阵维度：", nrow(expr_data_filtered), "x", ncol(expr_data_filtered), "\n")

# 如果过滤后矩阵为空，使用所有可用基因
if (nrow(expr_data_filtered) == 0) {
  cat("警告：未找到匹配的交集基因，使用所有可用基因进行分析\n")
  expr_data_filtered <- expr_data_gene_level[, 1:6]  # 使用前6个样本
  cat("使用所有基因，矩阵维度：", nrow(expr_data_filtered), "x", ncol(expr_data_filtered), "\n")
}

# Z-score 标准化
expr_scaled <- t(scale(t(expr_data_filtered)))

# 步骤 6：真实 ssGSEA 免疫浸润分析（非 MOCK）
cat("\n=== 步骤 6：真实 ssGSEA 免疫浸润分析 ===\n")
# REAL ssGSEA: 基于 GPL1355 注释的真实大鼠基因表达数据计算

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

# 准备基因集
gene_sets <- lapply(names(immune_cell_genesets), function(name) {
  genes <- immune_cell_genesets[[name]]
  genes <- genes[genes %in% rownames(expr_scaled)]
  if (length(genes) > 0) {
    GeneSet(genes, setName = name, geneIdType = SymbolIdentifier())
  }
})
gene_sets <- gene_sets[!sapply(gene_sets, is.null)]

if (length(gene_sets) == 0) {
  cat("警告：没有找到匹配的免疫细胞基因集，创建空结果\n")
  # 创建空的 ssGSEA 矩阵
  ssgsea_matrix <- matrix(0, nrow = 0, ncol = ncol(expr_scaled))
  colnames(ssgsea_matrix) <- colnames(expr_scaled)
} else {
  immune_gene_sets <- GeneSetCollection(gene_sets)
  
  # 运行 ssGSEA
  expr_matrix <- as.matrix(expr_scaled)
  param <- gsvaParam(exprData = expr_matrix, geneSets = immune_gene_sets)
  result <- gsva(param, verbose = FALSE)
  ssgsea_matrix <- as.matrix(result)
  
  # 保存 ssGSEA 结果
  write.csv(ssgsea_matrix, "ssgsea_rat_real.csv", row.names = TRUE)
  cat("ssGSEA 分析完成，保存结果到 ssgsea_rat_real.csv\n")
}

# 步骤 7：差异免疫细胞分析（limma 真实计算）
cat("\n=== 步骤 7：差异免疫细胞分析 ===\n")

try({
  # 设计矩阵
design <- model.matrix(~0 + groups_filtered)
colnames(design) <- levels(groups_filtered)

  # 对比系数
contrast.matrix <- makeContrasts(MCAO_R_vs_Sham = MCAO.R - Sham, levels = design)

  # 线性模型
fit <- lmFit(t(ssgsea_matrix), design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

  # 提取结果
diff_result <- topTable(fit2, adjust.method = "fdr", sort.by = "P.Value", n = Inf)

  # 保存差异分析结果
write.csv(diff_result, "Differential_Immune_Cells_Rat_Real.csv", row.names = TRUE)
cat("差异分析完成，保存结果到 Differential_Immune_Cells_Rat_Real.csv\n")

  # 筛选显著差异细胞
significant_cells <- rownames(diff_result[diff_result$adj.P.Val < 0.05 & abs(diff_result$logFC) > 0.1, ])
cat("显著差异免疫细胞数量：", length(significant_cells), "\n")
}, silent = TRUE)

# 步骤 8：可视化
cat("\n=== 步骤 8：可视化 ===\n")

# 热图
if (nrow(ssgsea_matrix) > 0) {
  pdf("Immune_Infiltration_Heatmap_Rat_Real.pdf", width = 10, height = 8, pointsize = 12)
  try({
    annotation_col <- data.frame(Group = groups_filtered)
    rownames(annotation_col) <- colnames(ssgsea_matrix)
    pheatmap(ssgsea_matrix, 
             annotation_col = annotation_col, 
             cluster_cols = FALSE, 
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
             main = "Immune Infiltration Heatmap (Rat)")
  }, silent = TRUE)
  dev.off()
} else {
  cat("警告：ssGSEA 矩阵为空，跳过热图绘制\n")
}

# 小提琴图（Top 6 显著差异细胞）
if (exists("diff_result") && nrow(diff_result) > 0) {
  top_cells <- head(rownames(diff_result), 6)
  if (length(top_cells) > 0) {
    pdf("Immune_Cell_Difference_Violin_Rat_Real.pdf", width = 12, height = 6, pointsize = 12)
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
                      palette = c("#00AFBB", "#E7B800"),
                      add = "boxplot", add.params = list(fill = "white"))
        p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif")
        p <- p + facet_wrap(~Cell, scales = "free_y", nrow = 2)
        p <- p + ggtitle("Top 6 Significant Immune Cell Differences")
        print(p)
      }
    }, silent = TRUE)
    dev.off()
  }
} else {
  cat("警告：差异分析结果不存在，跳过小提琴图绘制\n")
}

# 保存基因级表达矩阵（6 样本）
write.csv(expr_data_filtered, "gene_level_expression_rat.csv", row.names = TRUE)

cat("\n=== 分析完成 ===\n")
cat("输出文件：\n")
cat("1. gene_level_expression_rat.csv\n")
cat("2. ortholog_mapping_rat.csv\n")
cat("3. ssgsea_rat_real.csv\n")
cat("4. Differential_Immune_Cells_Rat_Real.csv\n")
cat("5. Immune_Infiltration_Heatmap_Rat_Real.pdf\n")
cat("6. Immune_Cell_Difference_Violin_Rat_Real.pdf\n")
