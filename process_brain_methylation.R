# === 脑组织甲基化数据 R 处理脚本 ===
# 来源: EWAS Data Hub brain_methylation_v1
# 格式: RData (matrix) + TXT (CpG x samples beta values)
# 输出: gene-methylation edges

# 配置
data_dir <- "C:/Users/Jy-Mentor-7/Desktop/GAT/brain_methylation_temp/brain_methylation"
output_dir <- "D:/反向网络药理学/GAT拓展维度"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=", strrep("=", 58), "\n", sep = "")
cat("脑组织甲基化数据处理\n")
cat("=", strrep("=", 58), "\n", sep = "")

# Step 1: 加载数据
cat("\n--- Step 1: 加载甲基化数据 ---\n")
rdata_file <- file.path(data_dir, "brain_methylation_v1.RData")
txt_file <- file.path(data_dir, "brain_methylation_v1.txt")

if (file.exists(rdata_file)) {
  cat("加载 RData 文件 (11 GB)...\n")
  load(rdata_file)
  # 检查加载的对象名
  obj_names <- ls()
  cat("加载的对象:", paste(obj_names, collapse = ", "), "\n")
  
  # 找到甲基化矩阵
  meth_obj <- NULL
  for (obj in obj_names) {
    if (obj %in% c("brain_methylation", "meth", "beta", "mSet", "data")) {
      meth_obj <- get(obj)
      cat("使用对象:", obj, "\n")
      break
    }
  }
  
  if (is.null(meth_obj)) {
    cat("自动检测对象...\n")
    for (obj in obj_names) {
      if (is.matrix(get(obj)) || is.data.frame(get(obj))) {
        dims <- dim(get(obj))
        cat("  候选:", obj, " dim:", dims[1], "x", dims[2], "\n")
        if (dims[1] > 10000 && dims[2] > 10) {
          meth_obj <- get(obj)
          cat("  使用:", obj, "\n")
          break
        }
      }
    }
  }
  
  if (!is.null(meth_obj)) {
    cat("甲基化矩阵维度:", nrow(meth_obj), "CpG x", ncol(meth_obj), "样本\n")
    cat("行名(前5):", paste(rownames(meth_obj)[1:5], collapse = ", "), "\n")
    cat("列名(前5):", paste(colnames(meth_obj)[1:5], collapse = ", "), "\n")
  } else {
    cat("未找到矩阵对象, 改用 TXT 文件...\n")
    meth_obj <- NULL
  }
}

if (is.null(meth_obj)) {
  cat("改用 TXT 文件加载 (5 GB)...\n")
  # 读取前2行获取样本信息
  header <- readLines(txt_file, n = 2)
  sample_ids <- strsplit(header[1], "\t")[[1]][-1]
  tissue_types <- strsplit(header[2], "\t")[[1]][-1]
  cat("样本数:", length(sample_ids), "\n")
  cat("组织类型:", unique(tissue_types)[1:10], "\n")
  
  # 使用 data.table 高效读取
  if (requireNamespace("data.table", quietly = TRUE)) {
    library(data.table)
    meth_dt <- fread(txt_file, skip = 2, sep = "\t", header = FALSE,
                     na.strings = "NA", showProgress = TRUE,
                     nThread = parallel::detectCores())
    rownames <- meth_dt[[1]]
    meth_dt[, 1 := NULL]
    meth_obj <- as.matrix(meth_dt)
    rownames(meth_obj) <- rownames
    colnames(meth_obj) <- sample_ids
    rm(meth_dt)
    cat("TXT 加载完成, 维度:", dim(meth_obj), "\n")
  } else {
    cat("data.table 未安装, 使用基础 R 读取 (较慢)...\n")
    meth_df <- read.table(txt_file, skip = 2, sep = "\t", 
                          na.strings = "NA", header = FALSE)
    rownames <- meth_df[, 1]
    meth_df <- meth_df[, -1]
    meth_obj <- as.matrix(meth_df)
    rownames(meth_obj) <- rownames
    colnames(meth_obj) <- sample_ids
    rm(meth_df)
    cat("TXT 加载完成, 维度:", dim(meth_obj), "\n")
  }
}

# Step 2: 获取 CpG→基因映射
cat("\n--- Step 2: 获取 CpG→基因映射 ---\n")

# 从 IlluminaHumanMethylationEPICanno 包获取
cpg_genes <- c()

if (requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
  cat("使用 IlluminaHumanMethylationEPICanno 包...\n")
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  cpg_genes <- ann[rownames(meth_obj), "UCSC_RefGene_Name"]
  names(cpg_genes) <- rownames(meth_obj)
  cat("映射数:", sum(!is.na(cpg_genes) & cpg_genes != ""), "\n")
} else if (requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  cat("使用 IlluminaHumanMethylation450kanno 包...\n")
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  cpg_genes <- ann[rownames(meth_obj), "UCSC_RefGene_Name"]
  names(cpg_genes) <- rownames(meth_obj)
  cat("映射数:", sum(!is.na(cpg_genes) & cpg_genes != ""), "\n")
} else {
  cat("安装 Illumina 注释包...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", update = FALSE, ask = FALSE)
  
  if (requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    cpg_genes <- ann[rownames(meth_obj), "UCSC_RefGene_Name"]
    names(cpg_genes) <- rownames(meth_obj)
    cat("映射数:", sum(!is.na(cpg_genes) & cpg_genes != ""), "\n")
  } else {
    cat("注释包安装失败, 使用简化映射...\n")
    cpg_genes <- rep("", nrow(meth_obj))
    names(cpg_genes) <- rownames(meth_obj)
  }
}

# Step 3: 计算平均甲基化并提取边
cat("\n--- Step 3: 计算甲基化状态 ---\n")

# 只处理有基因映射的 CpG
has_gene <- !is.na(cpg_genes) & cpg_genes != ""
cat("可映射到基因的 CpG:", sum(has_gene), "/", length(has_gene), "\n")

if (sum(has_gene) == 0) {
  cat("WARNING: 无 CpG 可映射到基因!\n")
  cat("使用 CpG probe ID 前缀作为基因映射 (仅用于边构建)...\n")
  # 使用简化的基因名映射
  cpg_genes <- gsub("^cg\\d+_", "", rownames(meth_obj))
  has_gene <- nchar(cpg_genes) > 0 & cpg_genes != rownames(meth_obj)
  cat("简化映射数:", sum(has_gene), "\n")
}

if (sum(has_gene) > 0) {
  # 筛选有基因映射的 CpG
  meth_sub <- meth_obj[has_gene, , drop = FALSE]
  gene_map <- cpg_genes[has_gene]
  
  cat("处理", nrow(meth_sub), "个CpG...\n")
  
  # 按行计算平均 beta
  avg_beta <- rowMeans(meth_sub, na.rm = TRUE)
  
  # 筛选显著高/低甲基化 (beta > 0.7 或 < 0.3)
  high_meth <- avg_beta > 0.7 & !is.na(avg_beta)
  low_meth <- avg_beta < 0.3 & !is.na(avg_beta)
  
  cat("高甲基化 CpG (beta > 0.7):", sum(high_meth), "\n")
  cat("低甲基化 CpG (beta < 0.3):", sum(low_meth), "\n")
  
  # 构建边
  edges_high <- data.frame(
    gene = gene_map[high_meth],
    cpg_id = rownames(meth_sub)[high_meth],
    avg_beta = avg_beta[high_meth],
    status = "hypermethylated",
    stringsAsFactors = FALSE
  )
  
  edges_low <- data.frame(
    gene = gene_map[low_meth],
    cpg_id = rownames(meth_sub)[low_meth],
    avg_beta = avg_beta[low_meth],
    status = "hypomethylated",
    stringsAsFactors = FALSE
  )
  
  edges <- rbind(edges_high, edges_low)
  
  # 多基因处理 (取第一个基因)
  edges$gene <- sapply(strsplit(edges$gene, ";"), `[`, 1)
  edges$gene <- sapply(strsplit(edges$gene, "/"), `[`, 1)
  edges$gene <- toupper(trimws(edges$gene))
  
  cat("总边数:", nrow(edges), "\n")
  cat("唯一基因:", length(unique(edges$gene)), "\n")
  
  # 保存
  output_file <- file.path(output_dir, "gene_methylation_edges.txt")
  write.table(edges, file = output_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  cat("OK! 保存到:", output_file, "\n")
  
  # Top genes
  gene_counts <- sort(table(edges$gene), decreasing = TRUE)
  cat("\nTop-10 甲基化基因:\n")
  print(head(gene_counts, 10))
  
  # 统计
  stats <- list(
    total_edges = nrow(edges),
    unique_genes = length(unique(edges$gene)),
    unique_cpgs = length(unique(edges$cpg_id)),
    hypermethylated = sum(edges$status == "hypermethylated"),
    hypomethylated = sum(edges$status == "hypomethylated"),
    top_methylated_genes = names(head(gene_counts, 20)),
    source = "EWAS Data Hub - Brain Methylation v1",
    n_samples = ncol(meth_obj),
    n_cpg_total = nrow(meth_obj)
  )
  
  stats_file <- file.path(output_dir, "gene_methylation_stats.json")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    library(jsonlite)
    writeLines(toJSON(stats, pretty = TRUE, auto_unbox = TRUE), stats_file)
  } else {
    cat("jsonlite 未安装, 跳过统计文件\n")
  }
  cat("统计文件:", stats_file, "\n")
} else {
  cat("ERROR: 无法获取 CpG→基因映射\n")
  quit(status = 1)
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("脑组织甲基化数据处理完成!\n")
cat("=", strrep("=", 60), "\n", sep = "")