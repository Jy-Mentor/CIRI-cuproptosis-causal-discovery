# ============================================================
# 阶段4 WGCNA R脚本
# 输入: 种子池基因表达子矩阵 (CSV)
# 输出: WGCNA模块分配 + hub基因
# ============================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
})

# ---- 解析命令行参数 ----
args <- commandArgs(trailingOnly = TRUE)
expr_file <- NULL
out_dir <- "."
n_samples <- 15

for (i in seq_along(args)) {
  if (args[i] == "--expr_file" && i < length(args)) {
    expr_file <- args[i + 1]
  } else if (args[i] == "--out_dir" && i < length(args)) {
    out_dir <- args[i + 1]
  } else if (args[i] == "--n_samples" && i < length(args)) {
    n_samples <- as.integer(args[i + 1])
  }
}

cat(sprintf("WGCNA参数:\n"))
cat(sprintf("  表达文件: %s\n", expr_file))
cat(sprintf("  输出目录: %s\n", out_dir))
cat(sprintf("  样本数: %d\n", n_samples))

if (is.null(expr_file) || !file.exists(expr_file)) {
  stop("表达文件不存在")
}

# ---- 加载数据 ----
cat("[1/5] 加载表达矩阵...\n")
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
cat(sprintf("  维度: %d 基因 x %d 样本\n", nrow(expr), ncol(expr)))

# 转置: WGCNA需要 样本×基因
datExpr <- t(expr)

# ---- 检查缺失值 ----
gsg <- goodSamplesGenes(datExpr, verbose = 1)
if (!gsg$allOK) {
  datExpr <- datExpr[, gsg$goodGenes]
  cat(sprintf("  移除坏基因后: %d 基因\n", ncol(datExpr)))
}

if (ncol(datExpr) < 20) {
  cat(sprintf("  基因数不足 (%d), 跳过WGCNA\n", ncol(datExpr)))
  quit(save = "no", status = 0)
}

# ---- 选择软阈值 ----
cat("[2/5] 选择软阈值...\n")

powers <- c(1:20)
sft <- tryCatch({
  pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
}, error = function(e) {
  cat(sprintf("  pickSoftThreshold失败: %s\n", e$message))
  NULL
})

if (is.null(sft)) {
  cat("  使用默认power=6\n")
  softPower <- 6
} else {
  sft_df <- sft$fitIndices
  # 选择R² > 0.8的最小power
  suitable <- sft_df$Power[sft_df$SFT.R.sq > 0.8]
  if (length(suitable) > 0) {
    softPower <- min(suitable)
  } else {
    softPower <- 6
  }
  cat(sprintf("  选择软阈值: %d\n", softPower))
}

# ---- 构建网络 ----
cat("[3/5] 构建共表达网络...\n")

adjacency <- adjacency(datExpr, power = softPower, type = "unsigned")
TOM <- TOMsimilarity(adjacency, TOMType = "unsigned")
dissTOM <- 1 - TOM

# 层次聚类
geneTree <- hclust(as.dist(dissTOM), method = "average")

# 动态剪切树
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 10
)

moduleColors <- labels2colors(dynamicMods)
cat(sprintf("  模块数: %d\n", length(unique(moduleColors))))

# ---- 计算模块特征基因 ----
cat("[4/5] 计算模块特征基因...\n")

MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes

# ---- 识别hub基因 ----
cat("[5/5] 识别hub基因...\n")

gene_names <- colnames(datExpr)
hub_genes_all <- c()

# 计算kME
mod_kME <- signedKME(datExpr, MEs, outputColumnName = "")

# 查看实际列名
cat(sprintf("  kME列名: %s\n", paste(colnames(mod_kME), collapse=", ")))

# 使用模块颜色直接匹配列名
for (mod in unique(moduleColors)) {
  if (mod == "grey") next
  
  mod_genes <- gene_names[moduleColors == mod]
  if (length(mod_genes) < 5) next
  
  # 尝试多种可能的列名格式
  possible_cols <- c(paste0("kME.", mod), paste0("ME", mod), paste0("kME", mod), mod)
  kme_col <- NULL
  
  for (col_name in possible_cols) {
    if (col_name %in% colnames(mod_kME)) {
      kme_col <- col_name
      break
    }
  }
  
  if (is.null(kme_col)) {
    cat(sprintf("  模块 %s: 无匹配kME列 (尝试: %s)，使用模块内连通度替代\n",
                mod, paste(possible_cols, collapse=", ")))
    
    # 备选方案: 使用模块内基因间相关性作为hub评分
    mod_expr <- datExpr[, mod_genes, drop=FALSE]
    mod_cor <- cor(mod_expr, use="pairwise.complete.obs")
    mod_connectivity <- rowSums(abs(mod_cor))
    names(mod_connectivity) <- mod_genes
    
    mod_kme_vals <- mod_connectivity
  } else {
    cat(sprintf("  模块 %s: 使用kME列 '%s'\n", mod, kme_col))
    mod_kme_vals <- mod_kME[mod_genes, kme_col]
  }
  
  if (length(mod_kme_vals) == 0 || all(is.na(mod_kme_vals))) {
    cat(sprintf("  模块 %s: kME值无效，跳过\n", mod))
    next
  }
  
  mod_kme_vals <- mod_kME[mod_genes, kme_col]
  
  if (length(mod_kme_vals) == 0 || all(is.na(mod_kme_vals))) {
    cat(sprintf("  模块 %s: kME值无效，跳过\n", mod))
    next
  }
  
  # Top 10 hub
  top_n <- min(10, length(mod_genes))
  top_idx <- order(mod_kme_vals, decreasing = TRUE, na.last = TRUE)[1:top_n]
  top_genes <- mod_genes[top_idx]
  
  hub_genes_all <- c(hub_genes_all, top_genes)
  
  cat(sprintf("  模块 %s: %d 基因, hub: %s\n",
              mod, length(mod_genes),
              paste(top_genes[1:min(5, top_n)], collapse = ", ")))
}

hub_genes_all <- unique(hub_genes_all)
hub_genes_all <- as.character(hub_genes_all)
cat(sprintf("  总hub基因: %d\n", length(hub_genes_all)))

# ---- 保存结果 ----
# 模块分配
module_df <- data.frame(
  ProbeID = gene_names,
  Module = moduleColors,
  stringsAsFactors = FALSE
)
write.csv(module_df, file.path(out_dir, "wgcna_modules.csv"), row.names = FALSE)

# Hub基因
writeLines(hub_genes_all, file.path(out_dir, "wgcna_hub_genes.txt"))

# 模块特征基因
write.csv(MEs, file.path(out_dir, "wgcna_MEs.csv"))

cat("WGCNA完成!\n")