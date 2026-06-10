#!/usr/bin/env Rscript

# 设置全局选项
options(stringsAsFactors = FALSE)
set.seed(42)

# 参数区
COR_THRESHOLD <- 0.5
ALPHA_BASE <- 0.05
SAVE_ADJACENCY <- FALSE  # 是否保存邻接矩阵
N_CORES <- 4  # 并行计算核心数

# 设置输出目录
output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# 缓存文件路径
cache_file <- file.path(output_dir, "analysis_cache.rds")

# 检查缓存是否存在
if (file.exists(cache_file)) {
  cat("发现缓存文件，尝试加载...\n")
  tryCatch({
    # 确保Seurat包已加载
    if (!require("Seurat", character.only = TRUE)) {
      install.packages("Seurat", dependencies = TRUE)
      library(Seurat)
    }
    cached_data <- readRDS(cache_file)
    # 检查缓存是否包含必要的对象
    if (all(c("neurons", "target_genes", "output_dir") %in% names(cached_data))) {
      cat("缓存加载成功，跳过步骤1-4\n")
      # 从缓存中加载对象
      neurons <- cached_data$neurons
      target_genes <- cached_data$target_genes
      output_dir <- cached_data$output_dir
      # 直接跳转到步骤5
      goto_step5 <- TRUE
    } else {
      cat("缓存文件不完整，重新运行步骤1-4\n")
      goto_step5 <- FALSE
    }
  }, error = function(e) {
    cat(sprintf("加载缓存失败: %s\n", e$message))
    goto_step5 <- FALSE
  })
} else {
  cat("未发现缓存文件，运行步骤1-4\n")
  goto_step5 <- FALSE
}

# 93个目标基因（人类Symbol）
target_genes <- c("ACTA2","ADORA1","AIF1","ALDH1A1","ALDH9A1","AOC3","ATF4","BRD4","C3","CASP8","CCL2","CCND1","CCR5","CDC42","CNDP2","CNR2","COL1A1","CP","CPT1A","CPT2","CTSB","CTSD","CTSS","CXCR3","DDIT3","EGR1","F3","FABP3","FABP5","FAS","FASN","GAD1","GFAP","GPT","HMGCR","HMOX1","HSPA5","HTR2A","ICAM1","IGF1R","IL6","IRF1","JAK1","MAOB","MAPK9","MDM2","MGLL","NFE2L2","NFKB1","NOTCH1","NR1H3","PARP1","PLA2G4A","PRKCQ","PTGES","PTGS1","PTGS2","PTPN6","PTPRC","RELA","S100A6","S1PR1","SAT1","SOD2","SREBF1","STAT1","STAT3","STAT5A","TIMP1","TGFB1","TSPO","XDH","FDX1","DLAT","DLD","LIPT1","PDHX","PDHB","SLC31A1","ATP7B","ATP7A","ATOX1","COMMD1","MT2A","TLR4")

# 步骤1：物种转换（人类→小鼠）
if (!goto_step5) {
  cat(">>>> 步骤1: 物种转换（人类→小鼠）\n")

# 读取映射文件
mapping_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt"
if (!file.exists(mapping_file)) {
  stop("物种映射文件路径错误: ", mapping_file)
}

# 跳过前36行注释，读取数据
cat("读取映射文件...\n")
mapping_data <- read.table(mapping_file, sep="\t", header=TRUE, skip=36, stringsAsFactors=FALSE)

# 检查必要的列是否存在
if (!all(c("HUMAN_ORTHOLOG_SYMBOL", "MOUSE_ORTHOLOG_SYMBOL") %in% colnames(mapping_data))) {
  stop("映射文件缺少必要的列")
}

# 处理多对多映射（拆分|符号，取第一个小鼠基因）
process_mouse_gene <- function(x) {
  if (is.na(x) || x == "") {
    return(NA)
  }
  # 拆分|符号，取第一个
  genes <- strsplit(x, "|", fixed=TRUE)[[1]]
  return(genes[1])
}

# 创建映射表
gene_mapping <- data.frame(
  human_symbol = target_genes,
  mouse_symbol = NA,
  stringsAsFactors = FALSE
)

# 进行映射
for (i in 1:length(target_genes)) {
  human_gene <- target_genes[i]
  # 找到对应的小鼠基因
  matches <- mapping_data[mapping_data$HUMAN_ORTHOLOG_SYMBOL == human_gene, "MOUSE_ORTHOLOG_SYMBOL"]
  if (length(matches) > 0) {
    # 处理多对多映射
    mouse_gene <- process_mouse_gene(matches[1])
    gene_mapping$mouse_symbol[i] <- mouse_gene
  }
}

# 统计映射结果
success_count <- sum(!is.na(gene_mapping$mouse_symbol))
failure_count <- sum(is.na(gene_mapping$mouse_symbol))
unmapped_genes <- gene_mapping$human_symbol[is.na(gene_mapping$mouse_symbol)]

# 输出统计结果
cat(sprintf("映射成功: %d, 映射失败: %d\n", success_count, failure_count))
if (failure_count > 0) {
  cat("未映射基因: ", paste(unmapped_genes, collapse=", "), "\n")
}

# 保存映射表
mapping_output <- file.path(output_dir, "gene_mapping_93.csv")
write.csv(gene_mapping, mapping_output, row.names=FALSE)
cat(sprintf("映射表已保存到: %s\n", mapping_output))

cat(">>>> 物种转换完成 | 关键统计: 成功映射数 =", success_count, "失败映射数 =", failure_count, "\n\n")

# 检查GSE目录是否存在
cat("检查GSE目录...\n")
gse_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/GSE"
if (!dir.exists(gse_dir)) {
  cat("GSE目录不存在，跳过后续步骤\n")
  quit()
}

# 检查样本目录是否存在
samples <- list(
  sham1 = list(path = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/GSE/GSM5319987_sham1", group = "control", sample_name = "sham1"),
  sham2 = list(path = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/GSE/GSM5319988_sham2", group = "control", sample_name = "sham2"),
  sham3 = list(path = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/GSE/GSM5319989_sham3", group = "control", sample_name = "sham3"),
  mcao1 = list(path = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/GSE/GSM5319990_mcao1", group = "stroke", sample_name = "mcao1"),
  mcao2 = list(path = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/GSE/GSM5319991_mcao2", group = "stroke", sample_name = "mcao2"),
  mcao3 = list(path = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/GSE/GSM5319992_mcao3", group = "stroke", sample_name = "mcao3")
)

for (sample_name in names(samples)) {
  sample_path <- samples[[sample_name]]$path
  if (dir.exists(sample_path)) {
    cat(sprintf("样本 %s 目录存在\n", sample_name))
    # 检查必要文件
    barcodes <- file.exists(file.path(sample_path, "barcodes.tsv"))
    features <- file.exists(file.path(sample_path, "features.tsv"))
    genes <- file.exists(file.path(sample_path, "genes.tsv"))
    matrix <- file.exists(file.path(sample_path, "matrix.mtx"))
    cat(sprintf("  barcodes.tsv: %s\n", ifelse(barcodes, "存在", "缺失")))
    cat(sprintf("  features.tsv: %s\n", ifelse(features, "存在", "缺失")))
    cat(sprintf("  genes.tsv: %s\n", ifelse(genes, "存在", "缺失")))
    cat(sprintf("  matrix.mtx: %s\n", ifelse(matrix, "存在", "缺失")))
  } else {
    cat(sprintf("样本 %s 目录不存在\n", sample_name))
  }
}

# 尝试加载Seurat和ggplot2
cat("\n尝试加载必要的包...\n")
if (require("Seurat", character.only = TRUE) && require("ggplot2", character.only = TRUE)) {
  cat("Seurat和ggplot2包加载成功\n")
  library(data.table)
} else {
  cat("必要的包未安装，尝试安装...\n")
  install.packages(c("Seurat", "ggplot2"), dependencies = TRUE)
  if (require("Seurat", character.only = TRUE) && require("ggplot2", character.only = TRUE)) {
    cat("包安装并加载成功\n")
    library(data.table)
  } else {
    cat("包安装失败，跳过后续步骤\n")
    quit()
  }
}

# 步骤2：读取单细胞数据
cat("\n>>>> 步骤2: 读取单细胞数据\n")

seurat_objects <- list()

# 创建临时英文路径
 temp_dir <- "C:/Users/Jy-Mentor-7/Desktop/test/GSE"
 if (!dir.exists(temp_dir)) {
   dir.create(temp_dir, recursive = TRUE)
 }

for (sample_name in names(samples)) {
  sample_info <- samples[[sample_name]]
  sample_path <- sample_info$path
  
  cat(sprintf("处理样本: %s (路径: %s)\n", sample_name, sample_path))
  
  # 检查文件存在性
  barcodes_path <- file.path(sample_path, "barcodes.tsv")
  features_path <- file.path(sample_path, "features.tsv")
  genes_path <- file.path(sample_path, "genes.tsv")
  matrix_path <- file.path(sample_path, "matrix.mtx")
  
  # 检查barcodes.tsv
  if (!file.exists(barcodes_path)) {
    stop(sprintf("样本 %s 缺失文件: %s", sample_name, barcodes_path))
  }
  
  # 检查matrix.mtx
  if (!file.exists(matrix_path)) {
    stop(sprintf("样本 %s 缺失文件: %s", sample_name, matrix_path))
  }
  
  # 检查features.tsv或genes.tsv
  if (!file.exists(features_path) && !file.exists(genes_path)) {
    stop(sprintf("样本 %s 缺失features.tsv和genes.tsv", sample_name))
  }
  
  # 创建Seurat对象
  cat(sprintf("创建Seurat对象: %s\n", sample_name))
  
  # 尝试手动读取和处理文件
  cat("尝试手动读取和处理文件...\n")
  
  # 读取barcodes.tsv（去除空行）
  barcodes_path <- file.path(sample_path, "barcodes.tsv")
  genes_path <- file.path(sample_path, "genes.tsv")
  matrix_path <- file.path(sample_path, "matrix.mtx")
  
  cat("读取barcodes.tsv...\n")
  barcodes <- readLines(barcodes_path, encoding = "UTF-8")
  barcodes <- barcodes[barcodes != ""]  # 去除空行
  cat(sprintf("读取到 %d 个barcodes\n", length(barcodes)))
  
  cat("读取genes.tsv...\n")
  features <- read.delim(genes_path, header = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("读取到 %d 个基因\n", nrow(features)))
  
  cat("读取matrix.mtx...\n")
  library(Matrix)
  counts <- readMM(matrix_path)
  cat(sprintf("矩阵维度: %d x %d\n", nrow(counts), ncol(counts)))
  
  # 确保barcodes数量与矩阵列数匹配
  if (length(barcodes) != ncol(counts)) {
    cat(sprintf("警告: barcodes数量 (%d) 与矩阵列数 (%d) 不匹配\n", length(barcodes), ncol(counts)))
    # 取最小长度
    min_length <- min(length(barcodes), ncol(counts))
    barcodes <- barcodes[1:min_length]
    counts <- counts[, 1:min_length]
    cat(sprintf("已调整为 %d 个细胞\n", min_length))
  }
  
  # 设置行名和列名
  # 关键修复：使用第2列Gene Symbol作为基因名（而非第1列Ensembl ID）
  if (ncol(features) >= 2) {
    gene_symbols <- features[, 2]  # 第二列是Symbol
  } else {
    gene_symbols <- features[, 1]
  }
  # 去除空值并设置唯一名称
  gene_symbols <- make.names(gene_symbols, unique = TRUE)
  # 全部转换为小写，统一基因名格式
  gene_symbols <- tolower(gene_symbols)
  rownames(counts) <- gene_symbols
  colnames(counts) <- barcodes
  
  # 创建Seurat对象
  cat("创建Seurat对象...\n")
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    min.cells = 3,
    min.features = 200,
    project = sample_name
  )
  cat("Seurat对象创建成功\n")
  
  # 添加metadata
  seurat_obj$group <- sample_info$group
  seurat_obj$sample <- sample_info$sample_name
  
  seurat_objects[[sample_name]] <- seurat_obj
}

# 清理临时目录
unlink(temp_dir, recursive = TRUE)

# 合并6个样本
cat("合并样本...\n")
sc_total <- merge(seurat_objects[[1]], seurat_objects[2:length(seurat_objects)])

# 输出合并结果
cat(sprintf("合并完成，总细胞数: %d\n", ncol(sc_total)))
cat(">>>> 读取单细胞数据完成 | 关键统计: 总细胞数 =", ncol(sc_total), "\n\n")

# 步骤3：Seurat标准化流程
cat(">>>> 步骤3: Seurat标准化流程\n")

# 质控：计算 percent.mt 和 percent.rb
cat("计算质控指标...\n")
sc_total <- PercentageFeatureSet(sc_total, pattern = "^mt-", col.name = "percent.mt")
sc_total <- PercentageFeatureSet(sc_total, pattern = "^Rp[sl]", col.name = "percent.rb")

# 过滤
cat("过滤细胞...\n")
sc_total <- subset(sc_total, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 15)
cat(sprintf("过滤后细胞数: %d\n", ncol(sc_total)))

# 标准化
cat("标准化数据...\n")
sc_total <- NormalizeData(sc_total)

# 找可变基因
cat("寻找可变基因...\n")
sc_total <- FindVariableFeatures(sc_total, selection.method = "vst", nfeatures = 3000)

# 缩放数据
cat("缩放数据...\n")
sc_total <- ScaleData(sc_total)

# 降维
cat("PCA降维...\n")
sc_total <- RunPCA(sc_total, dims = 1:30)

# 找邻居
cat("寻找邻居...\n")
sc_total <- FindNeighbors(sc_total, dims = 1:30)

# 聚类
cat("聚类...\n")
sc_total <- FindClusters(sc_total, resolution = 0.8)

# UMAP
cat("UMAP降维...\n")
sc_total <- RunUMAP(sc_total, dims = 1:30)

# 输出UMAP可视化
umap_output <- file.path(output_dir, "umap_overview.pdf")
pdf(umap_output, width = 10, height = 8)
DimPlot(sc_total, group.by = "group", cols = c("control" = "#2E86AB", "stroke" = "#F24236")) + 
  ggtitle("UMAP Overview (Control vs Stroke)")
dev.off()

cat(sprintf("UMAP可视化已保存到: %s\n", umap_output))
cat(">>>> Seurat标准化流程完成 | 关键统计: 过滤后细胞数 =", ncol(sc_total), "\n\n")

# 步骤4：细胞类型注释与神经元提取
cat(">>>> 步骤4: 细胞类型注释与神经元提取\n")

# 检查基因名称格式
cat("检查基因名称格式...\n")
head_genes <- head(rownames(sc_total), 10)
cat("前10个基因名:\n")
print(head_genes)

# 定义Marker基因（多种形式）
markers <- list(
  Excitatory = c("Slc17a7", "Camk2a", "slc17a7", "camk2a", "ENSMUSG00000026337", "ENSMUSG00000026341"),
  Inhibitory = c("Gad1", "Gad2", "gad1", "gad2", "ENSMUSG00000020143", "ENSMUSG00000020144"),
  Microglia = c("Cx3cr1", "P2ry12", "cx3cr1", "p2ry12", "ENSMUSG00000020670", "ENSMUSG00000032048"),
  Astrocytes = c("Aqp4", "Gfap", "aqp4", "gfap", "ENSMUSG00000001606", "ENSMUSG00000001310"),
  Oligo = c("Mbp", "Plp1", "mbp", "plp1", "ENSMUSG00000026107", "ENSMUSG00000026108")
)

# 检查并筛选存在的标记基因
for (cell_type in names(markers)) {
  present_markers <- markers[[cell_type]][markers[[cell_type]] %in% rownames(sc_total)]
  if (length(present_markers) > 0) {
    cat(sprintf("%s: 找到 %d 个标记基因\n", cell_type, length(present_markers)))
    print(present_markers)
    # 使用存在的标记基因
    markers[[cell_type]] <- present_markers
  } else {
    cat(sprintf("%s: 未找到标记基因\n", cell_type))
    # 移除没有标记基因的细胞类型
    markers[[cell_type]] <- NULL
  }
}

# 过滤掉没有标记基因的细胞类型
markers <- markers[!sapply(markers, is.null)]

# 添加ModuleScore
for (cell_type in names(markers)) {
  if (length(markers[[cell_type]]) >= 1) {
    cat(sprintf("计算%s的ModuleScore...\n", cell_type))
    sc_total <- AddModuleScore(
      sc_total,
      features = list(markers[[cell_type]]),
      name = paste0(cell_type, "_score")
    )
  } else {
    cat(sprintf("跳过%s，标记基因不足\n", cell_type))
  }
}

# 自动分配细胞类型
score_cols <- paste0(names(markers), "_score1")
cell_type_scores <- sc_total@meta.data[, score_cols]
colnames(cell_type_scores) <- names(markers)

# 每个细胞取最高评分对应的细胞类型
sc_total$cell_type <- apply(cell_type_scores, 1, function(x) names(markers)[which.max(x)])

# 提取神经元子集：Excitatory + Inhibitory
neurons <- subset(sc_total, subset = cell_type %in% c("Excitatory", "Inhibitory"))

# 报告统计结果
neuron_count <- ncol(neurons)
total_count <- ncol(sc_total)
neuron_ratio <- neuron_count / total_count * 100

cat(sprintf("神经元数量: %d, 占比: %.2f%%\n", neuron_count, neuron_ratio))

# control vs stroke 细胞数对比
control_neurons <- subset(neurons, subset = group == "control")
stroke_neurons <- subset(neurons, subset = group == "stroke")

cat(sprintf("Control组神经元数: %d\n", ncol(control_neurons)))
cat(sprintf("Stroke组神经元数: %d\n", ncol(stroke_neurons)))

# 检查神经元数量是否足够
if (neuron_count < 50) {
  cat("警告: 神经元数量过少（<50），无法进行后续分析\n")
  flush.console()
  quit()
}

cat(">>>> 细胞类型注释与神经元提取完成 | 关键统计: 神经元数 =", neuron_count, "占比 =", round(neuron_ratio, 2), "%\n\n")

  # 释放内存
  cat("释放内存...\n")
  flush.console()
  rm(sc_total)
  gc()
  
  # 保存缓存
  cat("保存缓存到文件...\n")
  flush.console()
  cached_data <- list(
    neurons = neurons,
    target_genes = target_genes,
    output_dir = output_dir
  )
  saveRDS(cached_data, cache_file)
  cat(sprintf("缓存已保存到: %s\n", cache_file))
  flush.console()
}

# 步骤5：93基因提取（Windows兼容版）
cat(">>>> 步骤5: 93基因提取（Windows兼容版）\n") 

# 首先诊断：显示Seurat对象中的实际基因名样本
cat("Seurat对象中的基因名样本（前20个）:\n")
print(head(rownames(neurons), 20))

# 铜死亡核心基因检查
cat("铜死亡核心基因检查:\n")
core_candidates <- c("fdx1", "lias", "slc31a1", "dlat", "Fdx1", "Lias", "Slc31a1", "Dlat", 
                     "FDX1", "LIAS", "SLC31A1", "DLAT")
for(g in core_candidates) {
  if(g %in% rownames(neurons)) {
    cat(sprintf("  ✓ 找到: %s\n", g))
  }
}

# 加载基因映射表，获取小鼠基因符号
cat("加载基因映射表...\n")
gene_mapping_file <- file.path(output_dir, "gene_mapping_93.csv")

# 检查文件存在性
if (!file.exists(gene_mapping_file)) {
  cat("错误: 基因映射文件不存在，请先完成步骤1\n")
  flush.console()
  quit()
}

gene_mapping <- read.csv(gene_mapping_file, stringsAsFactors = FALSE)

# 最简匹配方案（忽略大小写）
target_lower <- tolower(gene_mapping$mouse_symbol[!is.na(gene_mapping$mouse_symbol)])
available_lower <- tolower(rownames(neurons))
matched_lower <- intersect(target_lower, available_lower)
# 转回原始格式
matched_genes <- rownames(neurons)[available_lower %in% matched_lower]

# 强制确保铜死亡核心基因被包含（即使表达低）
core_priority <- c("fdx1", "lias", "slc31a1", "dlat")  # 小写格式
core_found <- intersect(tolower(core_priority), matched_lower)
if(length(core_found) < 4) {
  cat(sprintf("警告: 仅找到%d/4个铜死亡核心基因，尝试查找其他格式...\n", length(core_found)))
  # 尝试查找任何包含这些基因名片段的基因
  for(core in c("fdx", "lias", "slc31a", "dlat")) {
    partial_match <- grep(core, rownames(neurons), value=TRUE, ignore.case=TRUE)
    if(length(partial_match) > 0) {
      cat(sprintf("  发现部分匹配 '%s': %s\n", core, paste(partial_match, collapse=", ")))
      matched_genes <- unique(c(matched_genes, partial_match))
    }
  }
}

cat(sprintf("目标基因匹配: %d/%d\n", length(matched_genes), length(target_lower)))

# 检查未匹配的基因
unmatched_genes <- setdiff(target_lower, matched_lower)
if (length(unmatched_genes) > 0) {
  cat("未匹配的基因: ", paste(unmatched_genes, collapse=", "), "\n")
}

# 2. 使用Seurat v5 LayerData接口直接提取（保持稀疏格式） 
cat("提取表达矩阵...\n")
expr_sparse <- LayerData(neurons, assay="RNA", layer="counts", features=matched_genes) 
# 保持稀疏矩阵格式，避免内存爆炸
exp_matrix <- t(expr_sparse)  # 转置为细胞×基因，保持稀疏格式

# 3. 智能过滤：保留高表达基因 + 强制保留核心基因（即使低表达）
gene_counts <- colSums(exp_matrix > 0) 
high_genes <- names(gene_counts[gene_counts >= 5]) 

# 动态检测核心基因实际格式
core_candidates <- c("fdx1", "lias", "slc31a1", "dlat", "mt2a", "atox1", "nfkb1") 
core_found <- c() 

cat("动态检测核心基因格式...\n")
for(core in core_candidates) { 
  # 尝试精确匹配（忽略大小写） 
  match_idx <- grep(paste0("^", core, "$"), colnames(exp_matrix), ignore.case=TRUE) 
  if(length(match_idx) > 0) { 
    actual_name <- colnames(exp_matrix)[match_idx[1]] 
    core_found <- c(core_found, actual_name) 
    cat(sprintf("  ✓ 核心基因匹配: %s → %s\n", core, actual_name)) 
  } else { 
    # 尝试模糊匹配（包含子串） 
    partial_idx <- grep(core, colnames(exp_matrix), ignore.case=TRUE) 
    if(length(partial_idx) > 0) { 
      actual_name <- colnames(exp_matrix)[partial_idx[1]] 
      core_found <- c(core_found, actual_name) 
      cat(sprintf("  ✓ 核心基因模糊匹配: %s → %s\n", core, actual_name)) 
    } 
  } 
}

# 强制合并（使用实际检测到的格式） 
keep_genes <- unique(c(high_genes, core_found)) 
exp_matrix <- exp_matrix[, keep_genes, drop=FALSE] 

cat(sprintf("高表达基因: %d个 | 强制保留核心基因: %d个 | 最终: %d基因\n",  
             length(high_genes), length(core_found), ncol(exp_matrix)))

# 检查铜死亡核心基因是否存在
core_genes <- c("fdx1", "lias", "slc31a1", "dlat")
# 转换为实际格式进行检查
missing_core <- c()
for(core in core_genes) {
  match_idx <- grep(paste0("^", core, "$"), colnames(exp_matrix), ignore.case=TRUE)
  if(length(match_idx) == 0) {
    missing_core <- c(missing_core, core)
  }
}
if(length(missing_core) > 0) {
  cat("警告: 铜死亡核心基因缺失: ", paste(missing_core, collapse=", "), "\n")
  flush.console()
}

# 保存key_genes变量，避免后续代码出错
key_genes <- c()

cat(">>>> 93基因提取与预处理完成 | 关键统计: 最终基因数 =", ncol(exp_matrix), "\n\n")

# 步骤6：分块PC网络（一举两得版）
cat(">>>> 步骤6: 分块PC网络分析（全部93基因保留）\n")

# 安装并加载必要的包
cat("安装并加载pcalg及其依赖包...\n")
flush.console()
tryCatch({
  # 设置CRAN镜像
  options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  
  if (!require("pcalg", character.only = TRUE)) {
    # 安装BiocManager（如果未安装）
    if (!require("BiocManager", character.only = TRUE)) {
      install.packages("BiocManager")
      library(BiocManager)
    }
    # 安装RBGL（pcalg的依赖）
    BiocManager::install("RBGL")
    # 安装pcalg
    install.packages("pcalg")
    library(pcalg)
  }
  cat("pcalg包加载成功\n")
  flush.console()
}, error = function(e) {
  cat(sprintf("包安装失败: %s\n", e$message))
  flush.console()
  stop("无法安装必要的包，请手动安装pcalg及其依赖")
})

# 创建与exp_matrix行数匹配的逻辑向量
exp_cells <- rownames(exp_matrix)
cell_indices <- match(exp_cells, colnames(neurons))
cell_groups <- neurons$group[cell_indices]

# 自适应调整模块基因名格式以匹配exp_matrix
adjust_module_format <- function(module_genes, available_genes) {
  adjusted <- c()
  for(g in module_genes) {
    idx <- grep(paste0("^", g, "$"), available_genes, ignore.case=TRUE)
    if(length(idx) > 0) adjusted <- c(adjusted, available_genes[idx[1]])
  }
  return(unique(adjusted))
}

# 1. 定义功能模块（93个基因的科学分类）
gene_modules <- list(
  cuproptosis_core = adjust_module_format(
    c("fdx1", "lias", "dlat", "dld", "lipt1", "pdhx", "pdhb", "slc31a1"),
    colnames(exp_matrix)),
  copper_homeostasis = adjust_module_format(
    c("atp7b", "atp7a", "atox1", "commd1", "mt2a", "cp"),
    colnames(exp_matrix)),
  inflammation = adjust_module_format(
    c("nfkb1", "rela", "tlr4", "il6", "stat1", "stat3", "ccl2", "ptgs2"),
    colnames(exp_matrix)),
  others = setdiff(colnames(exp_matrix),
                   unlist(lapply(list(
                     c("fdx1", "lias", "dlat", "dld", "lipt1", "pdhx", "pdhb", "slc31a1"),
                     c("atp7b", "atp7a", "atox1", "commd1", "mt2a", "cp"),
                     c("nfkb1", "rela", "tlr4", "il6", "stat1", "stat3", "ccl2", "ptgs2")
                   ), function(x) adjust_module_format(x, colnames(exp_matrix))))
)

cat("模块基因格式调整完成:\n")
for(mod_name in names(gene_modules)) {
  cat(sprintf("  %s: %d个基因 (%s...)\n", mod_name, length(gene_modules[[mod_name]]), 
              paste(head(gene_modules[[mod_name]], 3), collapse=", ")))
}
flush.console()

# 辅助函数：检测跨模块边
detect_cross_module_edges <- function(exp_matrix, gene_modules, cor_threshold) {
  cat("检测跨模块边...\n")
  flush.console()
  
  cross_edges <- data.frame()
  
  # 收集所有模块基因
  all_module_genes <- c()
  module_indices <- list()
  
  if (!is.null(gene_modules$cuproptosis_core)) {
    all_module_genes <- c(all_module_genes, gene_modules$cuproptosis_core)
    module_indices$cuproptosis_core <- length(all_module_genes) - length(gene_modules$cuproptosis_core) + 1:length(gene_modules$cuproptosis_core)
  }
  
  if (!is.null(gene_modules$copper_homeostasis)) {
    all_module_genes <- c(all_module_genes, gene_modules$copper_homeostasis)
    module_indices$copper_homeostasis <- length(all_module_genes) - length(gene_modules$copper_homeostasis) + 1:length(gene_modules$copper_homeostasis)
  }
  
  if (!is.null(gene_modules$inflammation)) {
    all_module_genes <- c(all_module_genes, gene_modules$inflammation)
    module_indices$inflammation <- length(all_module_genes) - length(gene_modules$inflammation) + 1:length(gene_modules$inflammation)
  }
  
  # 去重并筛选存在的基因
  all_module_genes <- unique(all_module_genes)
  present_genes <- intersect(all_module_genes, colnames(exp_matrix))
  
  if (length(present_genes) < 2) {
    return(cross_edges)
  }
  
  # 向量化计算相关系数矩阵
  cat("计算相关系数矩阵...\n")
  flush.console()
  
  # 临时转为dense矩阵用于计算
  if (nrow(exp_matrix) > 10000) {
    cat("细胞数过多，采样10000个细胞计算相关系数...\n")
    flush.console()
    sample_idx <- sample(nrow(exp_matrix), 10000)
    exp_matrix_sample <- exp_matrix[sample_idx, present_genes, drop=FALSE]
    exp_matrix_dense <- as.matrix(exp_matrix_sample)
  } else {
    exp_matrix_dense <- as.matrix(exp_matrix[, present_genes, drop=FALSE])
  }
  
  # 计算相关系数矩阵（使用spearman方法）
  cor_matrix <- cor(exp_matrix_dense, method="spearman")
  
  # 检测cuproptosis_core与其他模块间的边
  if (!is.null(gene_modules$cuproptosis_core) && !is.null(module_indices$cuproptosis_core)) {
    cupro_genes <- intersect(gene_modules$cuproptosis_core, present_genes)
    if (length(cupro_genes) > 0) {
      # 与copper_homeostasis模块
      if (!is.null(gene_modules$copper_homeostasis)) {
        copper_genes <- intersect(gene_modules$copper_homeostasis, present_genes)
        if (length(copper_genes) > 0) {
          # 提取相关系数
          cor_subset <- cor_matrix[cupro_genes, copper_genes, drop=FALSE]
          # 筛选超过阈值的边
          high_cor <- which(abs(cor_subset) > cor_threshold, arr.ind=TRUE)
          if (nrow(high_cor) > 0) {
            for (i in 1:nrow(high_cor)) {
              g1 <- rownames(cor_subset)[high_cor[i, 1]]
              g2 <- colnames(cor_subset)[high_cor[i, 2]]
              cor_val <- cor_subset[high_cor[i, 1], high_cor[i, 2]]
              cross_edges <- rbind(cross_edges, data.frame(from=g1, to=g2, weight=cor_val))
              cat(sprintf("添加跨模块边: %s-%s (cor=%.3f)\n", g1, g2, cor_val))
              flush.console()
            }
          }
        }
      }
      
      # 与inflammation模块
      if (!is.null(gene_modules$inflammation)) {
        inflam_genes <- intersect(gene_modules$inflammation, present_genes)
        if (length(inflam_genes) > 0) {
          # 提取相关系数
          cor_subset <- cor_matrix[cupro_genes, inflam_genes, drop=FALSE]
          # 筛选超过阈值的边
          high_cor <- which(abs(cor_subset) > cor_threshold, arr.ind=TRUE)
          if (nrow(high_cor) > 0) {
            for (i in 1:nrow(high_cor)) {
              g1 <- rownames(cor_subset)[high_cor[i, 1]]
              g2 <- colnames(cor_subset)[high_cor[i, 2]]
              cor_val <- cor_subset[high_cor[i, 1], high_cor[i, 2]]
              cross_edges <- rbind(cross_edges, data.frame(from=g1, to=g2, weight=cor_val))
              cat(sprintf("添加跨模块边: %s-%s (cor=%.3f)\n", g1, g2, cor_val))
              flush.console()
            }
          }
        }
      }
    }
  }
  
  # 检测copper_homeostasis与inflammation模块间的边
  if (!is.null(gene_modules$copper_homeostasis) && !is.null(gene_modules$inflammation)) {
    copper_genes <- intersect(gene_modules$copper_homeostasis, present_genes)
    inflam_genes <- intersect(gene_modules$inflammation, present_genes)
    if (length(copper_genes) > 0 && length(inflam_genes) > 0) {
      # 提取相关系数
      cor_subset <- cor_matrix[copper_genes, inflam_genes, drop=FALSE]
      # 筛选超过阈值的边
      high_cor <- which(abs(cor_subset) > cor_threshold, arr.ind=TRUE)
      if (nrow(high_cor) > 0) {
        for (i in 1:nrow(high_cor)) {
          g1 <- rownames(cor_subset)[high_cor[i, 1]]
          g2 <- colnames(cor_subset)[high_cor[i, 2]]
          cor_val <- cor_subset[high_cor[i, 1], high_cor[i, 2]]
          cross_edges <- rbind(cross_edges, data.frame(from=g1, to=g2, weight=cor_val))
          cat(sprintf("添加跨模块边: %s-%s (cor=%.3f)\n", g1, g2, cor_val))
          flush.console()
        }
      }
    }
  }
  
  return(cross_edges)
}

# 辅助函数：创建邻接矩阵
create_adjacency_matrix <- function(edges, all_genes) {
  if (nrow(edges) == 0) {
    return(matrix(0, nrow=length(all_genes), ncol=length(all_genes), dimnames=list(all_genes, all_genes)))
  }
  
  # 使用Matrix包创建稀疏矩阵
  library(Matrix)
  i <- match(edges$from, all_genes)
  j <- match(edges$to, all_genes)
  x <- rep(1, nrow(edges))
  
  adj_matrix <- sparseMatrix(
    i = i,
    j = j,
    x = x,
    dims = c(length(all_genes), length(all_genes)),
    dimnames = list(all_genes, all_genes),
    symmetric = TRUE
  )
  
  # 转为普通矩阵
  return(as.matrix(adj_matrix))
}

# 2. 分块构建网络（每块基因数<25，计算极快）
build_module_network <- function(module_genes, module_name, data) {
  tryCatch({
    # 确保基因存在，尝试不同的大小写格式
    present_genes <- c()
    for (gene in module_genes) {
      # 尝试原始格式
      if (gene %in% colnames(data)) {
        present_genes <- c(present_genes, gene)
      } else {
        # 尝试小写
        gene_lower <- tolower(gene)
        if (gene_lower %in% colnames(data)) {
          present_genes <- c(present_genes, gene_lower)
        } else {
          # 尝试大写
          gene_upper <- toupper(gene)
          if (gene_upper %in% colnames(data)) {
            present_genes <- c(present_genes, gene_upper)
          } else {
            # 尝试首字母大写
            gene_title <- tools::toTitleCase(gene)
            if (gene_title %in% colnames(data)) {
              present_genes <- c(present_genes, gene_title)
            }
          }
        }
      }
    }
    
    # 去重
    present_genes <- unique(present_genes)
    
    if (length(present_genes) < 3) {
      cat(sprintf("%s模块: 基因数不足（%d个），跳过\n", module_name, length(present_genes)))
      flush.console()
      return(NULL)
    }
    
    module_data <- data[, present_genes, drop=FALSE]
    
    # 根据样本量动态调整alpha
    n_cells <- nrow(module_data)
    alpha <- ifelse(n_cells < 100, 0.1, ifelse(n_cells > 500, 0.01, ALPHA_BASE))
    cat(sprintf("%s模块: %d细胞，使用alpha=%.3f\n", module_name, n_cells, alpha))
    flush.console()
    
    # 临时转为dense矩阵用于计算，若细胞数过多则采样
    if (n_cells > 10000) {
      cat("细胞数过多，采样10000个细胞计算相关系数...\n")
      flush.console()
      sample_idx <- sample(n_cells, 10000)
      module_data_sample <- module_data[sample_idx, , drop=FALSE]
      module_data_dense <- as.matrix(module_data_sample)
    } else {
      module_data_dense <- as.matrix(module_data)
    }
    
    # PC算法（小数据量极快）
    library(pcalg)
    suff <- list(C=cor(module_data_dense), n=nrow(module_data_dense))
    pc.fit <- pc(suff, indepTest=gaussCItest, alpha=alpha, 
                 labels=colnames(module_data_dense), verbose=FALSE)
    
    # 提取边
    adj <- as(pc.fit@graph, "matrix")
    edges <- which(adj[upper.tri(adj)] == 1, arr.ind=TRUE)
    
    cat(sprintf("%s模块: %d基因, %d条边\n", module_name, length(present_genes), length(edges)))
    flush.console()
    return(list(adj=adj, genes=present_genes))
  }, error = function(e) {
    cat(sprintf("%s模块计算失败: %s\n", module_name, e$message))
    flush.console()
    return(NULL)
  })
}

# 3. 为Control和Stroke分别构建
all_edges_list <- list()
network_statistics <- data.frame()

# 加载必要的包
if (!require("data.table", character.only = TRUE)) {
  install.packages("data.table")
  library(data.table)
}

if (!require("parallel", character.only = TRUE)) {
  install.packages("parallel")
  library(parallel)
}

# Stroke组诊断
cat("\n==== Stroke组诊断 ====\n")
stroke_idx <- cell_groups == "stroke"
cat(sprintf("Stroke组细胞数: %d\n", sum(stroke_idx)))
if(sum(stroke_idx) > 0) {
  stroke_data <- exp_matrix[stroke_idx, , drop=FALSE]
  cat(sprintf("Stroke组表达矩阵维度: %d细胞 × %d基因\n", nrow(stroke_data), ncol(stroke_data)))
  cat("核心基因在Stroke组的表达情况:\n")
  for(g in core_found) {
    if(g %in% colnames(stroke_data)) {
      expr_cells <- sum(stroke_data[, g] > 0)
      cat(sprintf("  %s: %d/%d细胞表达 (%.1f%%)\n", 
                  g, expr_cells, nrow(stroke_data), 100*expr_cells/nrow(stroke_data)))
    }
  }
  if(sum(stroke_idx) < 30) {
    cat("警告: Stroke组细胞数<30，PC算法无法运行（样本量不足）\n")
  }
}
cat("====================\n\n")

for (group in c("control", "stroke")) {
  tryCatch({
    idx <- cell_groups == group
    group_data <- exp_matrix[idx, , drop=FALSE]
    
    cat(sprintf("\n构建%s组网络:\n", group))
    flush.console()  # 强制刷新控制台输出
    
    # 并行计算4个模块
    cat("并行计算模块网络...\n")
    flush.console()
    
    # 检测操作系统并选择并行方案
    is_windows <- Sys.info()["sysname"] == "Windows"
    
    if (!is_windows && N_CORES > 1) {
      # Mac/Linux使用mclapply
      library(parallel)
      cat(sprintf("使用mclapply进行%d核并行计算\n", N_CORES))
      networks <- mclapply(names(gene_modules), function(mod_name) {
        cat(sprintf("处理模块: %s\n", mod_name))
        build_module_network(gene_modules[[mod_name]], mod_name, group_data)
      }, mc.cores=N_CORES)
    } else if (is_windows && N_CORES > 1) {
      # Windows使用parLapply（支持多核）
      library(parallel)
      cat(sprintf("Windows系统：使用parLapply进行%d核并行计算\n", N_CORES))
      # 创建集群
      cl <- makeCluster(N_CORES)
      # 导出必要变量到集群
      clusterExport(cl, c("build_module_network", "gene_modules", "group_data", 
                          "ALPHA_BASE", "COR_THRESHOLD"), 
                    envir=environment())
      # 加载必要包到集群节点
      clusterEvalQ(cl, {
        library(pcalg)
        library(Matrix)
      })
      # 并行计算
      networks <- parLapply(cl, names(gene_modules), function(mod_name) {
        cat(sprintf("处理模块: %s\n", mod_name))
        build_module_network(gene_modules[[mod_name]], mod_name, group_data)
      })
      # 关闭集群（重要！）
      stopCluster(cl)
    } else {
      # 单核模式（Windows/Mac/Linux通用）
      cat("使用单核lapply计算\n")
      networks <- lapply(names(gene_modules), function(mod_name) {
        cat(sprintf("处理模块: %s\n", mod_name))
        build_module_network(gene_modules[[mod_name]], mod_name, group_data)
      })
    }
    
    names(networks) <- names(gene_modules)
    
    # 4. 合并所有模块的边（跨模块边用简单相关性补充）
    cat("合并模块边...\n")
    flush.console()  # 强制刷新控制台输出
    
    all_edges <- data.frame(from=character(), to=character(), weight=numeric())
    
    # 收集模块内边
    for (mod in networks) {
      if (is.null(mod)) next
      adj <- mod$adj
      genes <- mod$genes
      for (i in 1:nrow(adj)) {
        for (j in i:ncol(adj)) {
          if (adj[i,j] == 1) {
            all_edges <- rbind(all_edges, 
                             data.frame(from=genes[i], to=genes[j], weight=1))
          }
        }
      }
    }
    
    # 去重（使用data.table提高效率）
    cat("去重边列表...\n")
    flush.console()
    
    if (nrow(all_edges) > 0) {
      all_edges <- as.data.table(all_edges)
      all_edges[, edge_id := paste(pmin(from, to), pmax(from, to), sep="_")]
      all_edges <- unique(all_edges, by="edge_id")
      all_edges <- as.data.frame(all_edges[, c("from", "to", "weight")])
    }
    
    # 优化1：补充模块间强相关边（关键）
    cross_edges <- detect_cross_module_edges(group_data, gene_modules, COR_THRESHOLD)
    
    # 合并跨模块边
    if(nrow(cross_edges) > 0) {
      # 去重（使用data.table提高效率）
      cross_edges <- as.data.table(cross_edges)
      cross_edges[, edge_id := paste(pmin(from, to), pmax(from, to), sep="_")]
      
      if (nrow(all_edges) > 0) {
        all_edges_dt <- as.data.table(all_edges)
        all_edges_dt[, edge_id := paste(pmin(from, to), pmax(from, to), sep="_")]
        cross_edges <- cross_edges[!edge_id %in% all_edges_dt$edge_id, ]
      }
      
      cross_edges <- as.data.frame(cross_edges[, c("from", "to", "weight")])
      all_edges <- rbind(all_edges, cross_edges)
    }
  
  # 保存边列表（默认输出格式）
  cat("保存网络边列表...\n")
  flush.console()  # 强制刷新控制台输出
  write.csv(all_edges, file.path(output_dir, sprintf("network_%s_edges.csv", group)), row.names=FALSE)
  
  # 保存边列表（用于差异分析）
  all_edges_list[[group]] <- all_edges
  
  # 可选：保存邻接矩阵
  if (SAVE_ADJACENCY) {
    cat("保存邻接矩阵...\n")
    flush.console()
    all_genes <- unique(c(all_edges$from, all_edges$to))
    adj_matrix <- create_adjacency_matrix(all_edges, all_genes)
    write.csv(adj_matrix, file.path(output_dir, sprintf("network_%s_adj.csv", group)))
  }
  
  # 网络可视化
  if (nrow(all_edges) > 0) {
    cat("生成网络可视化...\n")
    flush.console()
    
    if (!require("igraph", character.only = TRUE)) {
      install.packages("igraph")
      library(igraph)
    }
    
    # 创建igraph对象
    g <- graph_from_data_frame(all_edges, directed=FALSE)
    
    # 按模块着色
    V(g)$color <- "gray"
    
    # 铜死亡核心基因红色
    if (!is.null(gene_modules$cuproptosis_core)) {
      cupro_genes <- intersect(gene_modules$cuproptosis_core, V(g)$name)
      if (length(cupro_genes) > 0) {
        V(g)[cupro_genes]$color <- "red"
      }
    }
    
    # 炎症基因蓝色
    if (!is.null(gene_modules$inflammation)) {
      inflam_genes <- intersect(gene_modules$inflammation, V(g)$name)
      if (length(inflam_genes) > 0) {
        V(g)[inflam_genes]$color <- "blue"
      }
    }
    
    # 铜稳态基因绿色
    if (!is.null(gene_modules$copper_homeostasis)) {
      copper_genes <- intersect(gene_modules$copper_homeostasis, V(g)$name)
      if (length(copper_genes) > 0) {
        V(g)[copper_genes]$color <- "green"
      }
    }
    
    # 绘制网络
    pdf(file.path(output_dir, sprintf("network_%s_plot.pdf", group)), width=12, height=10)
    plot(g, 
         layout=layout_with_fr(g), 
         vertex.size=5, 
         edge.width=0.5, 
         main=sprintf("%s Group Network", group),
         vertex.label.cex=0.7)
    # 添加图例
    legend("bottomright", 
           legend=c("Cuproptosis Core", "Inflammation", "Copper Homeostasis", "Others"),
           col=c("red", "blue", "green", "gray"),
           pch=19, 
           cex=0.8)
    dev.off()
    
    cat(sprintf("网络可视化已保存到: network_%s_plot.pdf\n", group))
    flush.console()
    
    # 计算网络拓扑统计
    cat("计算网络拓扑统计...\n")
    flush.console()
    
    n_nodes <- vcount(g)
    n_edges <- ecount(g)
    density <- graph.density(g)
    avg_degree <- mean(degree(g))
    transitivity <- transitivity(g, type="global")
    
    # 保存统计信息
    network_statistics <- rbind(network_statistics, data.frame(
      group=group,
      nodes=n_nodes,
      edges=n_edges,
      density=density,
      avg_degree=avg_degree,
      transitivity=ifelse(is.na(transitivity), 0, transitivity)
    ))
    
    # 输出统计信息
    cat(sprintf("%s组网络统计: %d节点, %d边, 密度=%.4f, 平均度=%.2f, 聚类系数=%.4f\n", 
                group, n_nodes, n_edges, density, avg_degree, 
                ifelse(is.na(transitivity), 0, transitivity)))
    flush.console()
  } else {
    # 空网络的统计信息
    network_statistics <- rbind(network_statistics, data.frame(
      group=group,
      nodes=0,
      edges=0,
      density=0,
      avg_degree=0,
      transitivity=0
    ))
    
    cat(sprintf("%s组网络构建完成，无边\n", group))
    flush.console()
  }
  
  }, error = function(e) {
    cat(sprintf("%s组网络构建失败: %s\n", group, e$message))
    flush.console()
    # 保存空的边列表
    all_edges <- data.frame(from=character(), to=character(), weight=numeric())
    write.csv(all_edges, file.path(output_dir, sprintf("network_%s_edges.csv", group)), row.names=FALSE)
    all_edges_list[[group]] <- all_edges
    # 保存空的统计信息
    network_statistics <- rbind(network_statistics, data.frame(
      group=group,
      nodes=0,
      edges=0,
      density=0,
      avg_degree=0,
      transitivity=0
    ))
  })
}

cat(">>>> PC因果网络分析完成\n\n")

# 步骤7：差异网络分析
cat(">>>> 步骤7: 差异网络分析\n")

# 检查两组网络是否都构建成功
if (!is.null(all_edges_list$control) && !is.null(all_edges_list$stroke)) {
  # 提取两组的边
  control_edges <- all_edges_list$control
  stroke_edges <- all_edges_list$stroke
  
  # 创建边的唯一标识符
  control_edges$edge_id <- apply(control_edges[, c("from", "to")], 1, function(x) paste(sort(x), collapse="_"))
  stroke_edges$edge_id <- apply(stroke_edges[, c("from", "to")], 1, function(x) paste(sort(x), collapse="_"))
  
  # 优化3：使用dplyr::anti_join（比%in%逻辑快10倍）
  cat("使用anti_join识别差异边...\n")
  flush.console()
  
  # 加载dplyr包
  if (!require("dplyr", character.only = TRUE)) {
    install.packages("dplyr")
    library(dplyr)
  }
  
  # 识别Stroke特异性新边（在stroke中存在，在control中不存在）
  stroke_specific_edges <- stroke_edges %>% 
    anti_join(control_edges, by = "edge_id")
  
  if (nrow(stroke_specific_edges) > 0) {
    # 移除edge_id列，重命名列名
    stroke_specific_edges <- stroke_specific_edges[, c("from", "to", "weight")]
    colnames(stroke_specific_edges) <- c("From", "To", "weight")
    stroke_specific_edges$Type <- "Stroke_Specific"
    
    # 输出novel_stroke_edges.csv
    edges_output <- file.path(output_dir, "novel_stroke_edges.csv")
    write.csv(stroke_specific_edges, edges_output, row.names=FALSE)
    cat(sprintf("Stroke特异性新边已保存到: %s\n", edges_output))
    flush.console()
    
    # 统计：stroke特异性边数
    cat(sprintf("Stroke特异性边数: %d\n", nrow(stroke_specific_edges)))
    flush.console()
    
    # 涉及的关键基因
    key_genes_in_edges <- unique(c(stroke_specific_edges$From, stroke_specific_edges$To))
    # 检查是否涉及铜死亡核心基因
    cuproptosis_core_genes <- gene_modules$cuproptosis_core
    key_genes_in_edges <- intersect(key_genes_in_edges, cuproptosis_core_genes)
    if (length(key_genes_in_edges) > 0) {
      cat("涉及的铜死亡核心基因: ", paste(key_genes_in_edges, collapse=", "), "\n")
      flush.console()
    }
  } else {
    cat("未发现Stroke特异性新边\n")
    flush.console()
    # 创建空文件
    edges_output <- file.path(output_dir, "novel_stroke_edges.csv")
    write.csv(data.frame(From = character(), To = character(), Type = character(), weight = numeric()), 
              edges_output, row.names=FALSE)
  }
} else {
  cat("至少一组网络构建失败，无法进行差异分析\n")
  flush.console()
  # 创建空文件
  edges_output <- file.path(output_dir, "novel_stroke_edges.csv")
  write.csv(data.frame(From = character(), To = character(), Type = character(), weight = numeric()), 
            edges_output, row.names=FALSE)
}

cat(">>>> 差异网络分析完成\n\n")

# 保存网络统计信息
if (nrow(network_statistics) > 0) {
  cat("保存网络统计信息...\n")
  flush.console()
  write.csv(network_statistics, file.path(output_dir, "network_statistics.csv"), row.names=FALSE)
  cat("网络统计信息已保存到: network_statistics.csv\n")
  flush.console()
}

# 验证所有输出文件是否生成
cat(">>>> 验证输出文件\n")
output_files <- c(
  file.path(output_dir, "gene_mapping_93.csv"),
  file.path(output_dir, "umap_overview.pdf"),
  file.path(output_dir, "network_control_edges.csv"),
  file.path(output_dir, "network_stroke_edges.csv"),
  file.path(output_dir, "network_control_plot.pdf"),
  file.path(output_dir, "network_stroke_plot.pdf"),
  file.path(output_dir, "network_statistics.csv"),
  file.path(output_dir, "novel_stroke_edges.csv")
)

for (file in output_files) {
  if (file.exists(file)) {
    cat(sprintf("✓ %s 已生成\n", basename(file)))
  } else {
    cat(sprintf("✗ %s 未生成\n", basename(file)))
  }
  flush.console()
}

cat("\n分析完成！\n")
flush.console()
