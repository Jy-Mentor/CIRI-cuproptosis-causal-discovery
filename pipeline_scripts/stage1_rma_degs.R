# ============================================================
# 阶段1: GSE61616 Affymetrix .CEL文件 RMA归一化 + limma差异表达
# 数据: 大鼠MCAO模型, Sham vs Model vs XST(药物)
# 芯片: Affymetrix Rat Genome 230 2.0 Array (GPL1355)
# 输出: 表达矩阵(expr_matrix.csv) + DEGs(limma_degs.csv)
# ============================================================

suppressPackageStartupMessages({
  library(oligo)
  library(limma)
  library(rat2302.db)
  library(AnnotationDbi)
})

set.seed(123)

# ---- 路径配置 ----
base_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results", "stage1_rma_degs")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

raw_tar <- "C:/Users/Jy-Mentor-7/Downloads/GSE61616_RAW (1).tar"
cel_dir <- file.path(data_dir, "GSE61616_CEL")
dir.create(cel_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 解压CEL文件 ----
cat("[1/6] 解压CEL文件...\n")
untar(raw_tar, exdir = cel_dir)
cel_files <- list.files(cel_dir, pattern = "\\.CEL$|\\.CEL\\.gz$", 
                        full.names = TRUE, ignore.case = TRUE)
cat(sprintf("  找到 %d 个CEL文件\n", length(cel_files)))

# ---- 读取CEL文件 ----
cat("[2/6] 读取CEL文件 (oligo)...\n")
raw_data <- read.celfiles(cel_files)
cat(sprintf("  探针数: %d, 样本数: %d\n", nrow(raw_data), ncol(raw_data)))

# ---- 过滤: 仅保留Sham和Model样本 ----
cat("[2b/6] 过滤样本 (仅保留Sham和Model)...\n")
sample_names_raw <- sampleNames(raw_data)
keep_idx <- c()
for (i in seq_along(sample_names_raw)) {
  sn <- tolower(sample_names_raw[i])
  if (grepl("sham", sn) || grepl("model", sn)) {
    keep_idx <- c(keep_idx, i)
  }
}
cat(sprintf("  原始样本: %d, 保留Sham+Model: %d\n", 
            length(sample_names_raw), length(keep_idx)))
cat("  保留样本:\n")
for (i in keep_idx) {
  cat(sprintf("    %s\n", sample_names_raw[i]))
}
raw_data <- raw_data[, keep_idx]

# ---- RMA归一化 ----
cat("[3/6] RMA归一化 (仅Sham+Model)...\n")
eset_rma <- rma(raw_data)
expr_matrix <- exprs(eset_rma)
cat(sprintf("  归一化后: %d 探针 × %d 样本\n", nrow(expr_matrix), ncol(expr_matrix)))

# ---- 保存表达矩阵 ----
expr_file <- file.path(results_dir, "expr_matrix.csv")
write.csv(expr_matrix, file = expr_file, row.names = TRUE)
cat(sprintf("  表达矩阵已保存: %s\n", expr_file))

# ---- 构建分组信息 ----
cat("[4/6] 构建分组信息...\n")
sample_names <- colnames(expr_matrix)
group <- rep(NA, length(sample_names))

for (i in seq_along(sample_names)) {
  sn <- tolower(sample_names[i])
  if (grepl("sham", sn)) {
    group[i] <- "Sham"
  } else if (grepl("model", sn)) {
    group[i] <- "Model"
  } else if (grepl("xst", sn)) {
    group[i] <- "XST"
  }
}

cat("  分组:\n")
for (i in seq_along(sample_names)) {
  cat(sprintf("    %s -> %s\n", sample_names[i], group[i]))
}

# ---- limma差异表达 (Sham vs Model) ----
cat("[5/6] limma差异表达分析 (Sham vs Model)...\n")

# 只取Sham和Model组
keep_idx <- which(group %in% c("Sham", "Model"))
expr_sub <- expr_matrix[, keep_idx]
group_sub <- factor(group[keep_idx], levels = c("Sham", "Model"))

design <- model.matrix(~ group_sub)
colnames(design) <- c("Intercept", "Model_vs_Sham")

fit <- lmFit(expr_sub, design)
fit <- eBayes(fit)

degs <- topTable(fit, coef = "Model_vs_Sham", number = Inf, adjust.method = "BH")
degs$ProbeID <- rownames(degs)

cat(sprintf("  总探针: %d\n", nrow(degs)))
cat(sprintf("  显著DEGs (adj.P<0.05 & |logFC|>1): %d\n", 
            sum(degs$adj.P.Val < 0.05 & abs(degs$logFC) > 1, na.rm = TRUE)))

# ---- 探针注释 (使用rat2302.db) ----
cat("[6/6] 探针注释 (rat2302.db)...\n")

probe_ids <- rownames(degs)
anno <- tryCatch({
  select(rat2302.db, keys = probe_ids, columns = c("SYMBOL", "GENENAME"), 
         keytype = "PROBEID")
}, error = function(e) {
  cat(sprintf("  rat2302.db查询失败: %s, 尝试biomaRt...\n", e$message))
  NULL
})

if (!is.null(anno) && nrow(anno) > 0) {
  # 去重: 每个探针取第一个基因符号
  anno_unique <- anno[!duplicated(anno$PROBEID), ]
  degs$GeneSymbol <- anno_unique$SYMBOL[match(degs$ProbeID, anno_unique$PROBEID)]
  degs$GeneName <- anno_unique$GENENAME[match(degs$ProbeID, anno_unique$PROBEID)]
} else {
  degs$GeneSymbol <- ""
  degs$GeneName <- ""
}

degs$GeneSymbol[is.na(degs$GeneSymbol)] <- ""

# 保存DEGs
deg_file <- file.path(results_dir, "limma_degs.csv")
write.csv(degs, file = deg_file, row.names = FALSE)
cat(sprintf("  DEGs已保存: %s\n", deg_file))

# ---- 保存分组信息 ----
group_df <- data.frame(
  Sample = sample_names,
  Group = group,
  stringsAsFactors = FALSE
)
group_file <- file.path(results_dir, "sample_groups.csv")
write.csv(group_df, file = group_file, row.names = FALSE)

# ---- 输出摘要 ----
cat("\n========================================\n")
cat("阶段1完成! 输出文件:\n")
cat(sprintf("  表达矩阵: %s\n", expr_file))
cat(sprintf("  DEGs: %s\n", deg_file))
cat(sprintf("  分组: %s\n", group_file))
cat(sprintf("  显著DEGs: %d\n", 
            sum(degs$adj.P.Val < 0.05 & abs(degs$logFC) > 1, na.rm = TRUE)))
cat("========================================\n")