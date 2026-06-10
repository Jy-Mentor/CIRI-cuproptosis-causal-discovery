#!/usr/bin/env Rscript

# 设置全局选项
options(stringsAsFactors = FALSE, encoding = "UTF-8")
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

# 不再使用缓存，每次从步骤1开始运行
# 确保数据一致性和可重复性
goto_step5 <- FALSE
cat("Running all steps from beginning to ensure data consistency...\n")

# 93个目标基因（人类Symbol）
target_genes <- c("ACTA2","ADORA1","AIF1","ALDH1A1","ALDH9A1","AOC3","ATF4","BRD4","C3","CASP8","CCL2","CCND1","CCR5","CDC42","CNDP2","CNR2","COL1A1","CP","CPT1A","CPT2","CTSB","CTSD","CTSS","CXCR3","DDIT3","EGR1","F3","FABP3","FABP5","FAS","FASN","GAD1","GFAP","GPT","HMGCR","HMOX1","HSPA5","HTR2A","ICAM1","IGF1R","IL6","IRF1","JAK1","MAOB","MAPK9","MDM2","MGLL","NFE2L2","NFKB1","NOTCH1","NR1H3","PARP1","PLA2G4A","PRKCQ","PTGES","PTGS1","PTGS2","PTPN6","PTPRC","RELA","S100A6","S1PR1","SAT1","SOD2","SREBF1","STAT1","STAT3","STAT5A","TIMP1","TGFB1","TSPO","XDH","FDX1","DLAT","DLD","LIPT1","PDHX","PDHB","SLC31A1","ATP7B","ATP7A","ATOX1","COMMD1","MT2A","TLR4")

# 铜死亡核心基因和铜稳态调控基因（全局常量，小写）
CORE_CUPROPTOSIS <- c("fdx1", "lias", "dlat", "slc31a1", "dld", "lipt1", "pdhx", "pdhb")
CORE_HOMEOSTASIS <- c("atp7b", "atp7a", "atox1", "commd1", "mt2a")

# 步骤1：物种转换（人类→小鼠）
if (!goto_step5) {
  cat(">>>> Step 1: Species conversion (Human → Mouse)\n")

# 读取映射文件
mapping_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt"
if (!file.exists(mapping_file)) {
  stop("Species mapping file path error: ", mapping_file)
}

# 跳过前36行注释，读取数据（使用UTF-8编码避免中文路径乱码）
cat("Reading mapping file...\n")
mapping_data <- read.table(file(mapping_file, "r", encoding="UTF-8"), sep="\t", header=TRUE, skip=36, stringsAsFactors=FALSE)

# 检查必要的列是否存在
if (!all(c("HUMAN_ORTHOLOG_SYMBOL", "MOUSE_ORTHOLOG_SYMBOL") %in% colnames(mapping_data))) {
  stop("Mapping file missing required columns")
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
cat(sprintf("Mapping success: %d, Mapping failure: %d\n", success_count, failure_count))
if (failure_count > 0) {
  cat("Unmapped genes: ", paste(unmapped_genes, collapse=", "), "\n")
}

# 保存映射表
mapping_output <- file.path(output_dir, "gene_mapping_93.csv")
write.csv(gene_mapping, mapping_output, row.names=FALSE)
cat(sprintf("Mapping table saved to: %s\n", mapping_output))

cat(">>>> Species conversion completed | Key stats: Success =", success_count, "Failure =", failure_count, "\n\n")

# 检查GSE目录是否存在
cat("Checking GSE directory...\n")
gse_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/GSE"
if (!dir.exists(gse_dir)) {
  cat("GSE directory not found, skipping subsequent steps\n")
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
    cat(sprintf("Sample %s directory exists\n", sample_name))
    # 检查必要文件
    barcodes <- file.exists(file.path(sample_path, "barcodes.tsv"))
    features <- file.exists(file.path(sample_path, "features.tsv"))
    genes <- file.exists(file.path(sample_path, "genes.tsv"))
    matrix <- file.exists(file.path(sample_path, "matrix.mtx"))
    cat(sprintf("  barcodes.tsv: %s\n", ifelse(barcodes, "exists", "missing")))
    cat(sprintf("  features.tsv: %s\n", ifelse(features, "exists", "missing")))
    cat(sprintf("  genes.tsv: %s\n", ifelse(genes, "exists", "missing")))
    cat(sprintf("  matrix.mtx: %s\n", ifelse(matrix, "exists", "missing")))
  } else {
    cat(sprintf("Sample %s directory not exists\n", sample_name))
  }
}

# 尝试加载Seurat和ggplot2
cat("\nTrying to load required packages...\n")
if (require("Seurat", character.only = TRUE) && require("ggplot2", character.only = TRUE)) {
  cat("Seurat and ggplot2 packages loaded successfully\n")
  library(data.table)
} else {
  cat("Required packages not installed, trying to install...\n")
  install.packages(c("Seurat", "ggplot2"), dependencies = TRUE)
  if (require("Seurat", character.only = TRUE) && require("ggplot2", character.only = TRUE)) {
    cat("Packages installed and loaded successfully\n")
    library(data.table)
  } else {
    cat("Package installation failed, skipping subsequent steps\n")
    quit()
  }
}

# 步骤2：读取单细胞数据
cat("\n>>>> Step 2: Reading single-cell data\n")

seurat_objects <- list()

# 创建临时英文路径
 temp_dir <- "C:/Users/Jy-Mentor-7/Desktop/test/GSE"
 if (!dir.exists(temp_dir)) {
   dir.create(temp_dir, recursive = TRUE)
 }

for (sample_name in names(samples)) {
  sample_info <- samples[[sample_name]]
  sample_path <- sample_info$path
  
  cat(sprintf("Processing sample: %s (path: %s)\n", sample_name, sample_path))
  
  # 检查文件存在性
  barcodes_path <- file.path(sample_path, "barcodes.tsv")
  features_path <- file.path(sample_path, "features.tsv")
  genes_path <- file.path(sample_path, "genes.tsv")
  matrix_path <- file.path(sample_path, "matrix.mtx")
  
  # 检查barcodes.tsv
  if (!file.exists(barcodes_path)) {
    stop(sprintf("Sample %s missing file: %s", sample_name, barcodes_path))
  }
  
  # 检查matrix.mtx
  if (!file.exists(matrix_path)) {
    stop(sprintf("Sample %s missing file: %s", sample_name, matrix_path))
  }
  
  # 检查features.tsv或genes.tsv
  if (!file.exists(features_path) && !file.exists(genes_path)) {
    stop(sprintf("Sample %s missing features.tsv and genes.tsv", sample_name))
  }
  
  # 创建Seurat对象
  cat(sprintf("Creating Seurat object: %s\n", sample_name))
  
  # 尝试手动读取和处理文件
  cat("Trying to read and process files manually...\n")
  
  # 读取barcodes.tsv（去除空行）
  barcodes_path <- file.path(sample_path, "barcodes.tsv")
  genes_path <- file.path(sample_path, "genes.tsv")
  matrix_path <- file.path(sample_path, "matrix.mtx")
  
  cat("Reading barcodes.tsv...\n")
  barcodes <- readLines(barcodes_path, encoding = "UTF-8")
  barcodes <- barcodes[barcodes != ""]  # 去除空行
  cat(sprintf("Read %d barcodes\n", length(barcodes)))
  
  cat("Reading genes.tsv...\n")
  features <- read.delim(genes_path, header = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("Read %d genes\n", nrow(features)))
  
  cat("Reading matrix.mtx...\n")
  library(Matrix)
  counts <- readMM(matrix_path)
  cat(sprintf("Matrix dimensions: %d x %d\n", nrow(counts), ncol(counts)))
  
  # 确保barcodes数量与矩阵列数匹配
  if (length(barcodes) != ncol(counts)) {
    cat(sprintf("Warning: Barcodes count (%d) doesn't match matrix columns (%d)\n", length(barcodes), ncol(counts)))
    # 取最小长度
    min_length <- min(length(barcodes), ncol(counts))
    barcodes <- barcodes[1:min_length]
    counts <- counts[, 1:min_length]
    cat(sprintf("Adjusted to %d cells\n", min_length))
  }
  
  # 设置行名和列名
  # 关键修复：使用第2列Gene Symbol作为基因名（而非第1列Ensembl ID）
  if (ncol(features) >= 2) {
    gene_symbols <- features[, 2]  # 第二列是Symbol
  } else {
    gene_symbols <- features[, 1]
  }
  
  # 处理空值
  gene_symbols[gene_symbols == "" | is.na(gene_symbols)] <- "unknown"
  
  # 确保基因名唯一性，强制转为小写
  gene_symbols <- tolower(make.names(gene_symbols, unique = TRUE))
  
  # 再次检查并处理重复基因名
  if (any(duplicated(gene_symbols))) {
    cat(sprintf("Warning: Found %d duplicate gene symbols, adding suffixes\n", sum(duplicated(gene_symbols))))
    # 为重复的基因名添加后缀
    gene_symbols <- make.unique(gene_symbols, sep = ".")
  }
  
  # 后续所有基因匹配必须使用小写（包括target_genes_mouse转小写）
  rownames(counts) <- gene_symbols
  colnames(counts) <- barcodes
  
  # 创建Seurat对象
  cat("Creating Seurat object...\n")
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    min.cells = 3,
    min.features = 200,
    project = sample_name
  )
  cat("Seurat object created successfully\n")
  
  # 添加metadata
  seurat_obj$group <- sample_info$group
  seurat_obj$sample <- sample_info$sample_name
  
  seurat_objects[[sample_name]] <- seurat_obj
}

# 清理临时目录
unlink(temp_dir, recursive = TRUE)

# 合并6个样本
cat("Merging samples...\n")
sc_total <- merge(seurat_objects[[1]], seurat_objects[2:length(seurat_objects)])

# 输出合并结果
cat(sprintf("Merge completed, total cells: %d\n", ncol(sc_total)))
cat(">>>> Reading single-cell data completed | Key stats: Total cells =", ncol(sc_total), "\n\n")

# 步骤3：Seurat标准化流程
cat(">>>> Step 3: Seurat normalization workflow\n")

# 质控：计算 percent.mt 和 percent.rb
cat("Calculating quality control metrics...\n")
sc_total <- PercentageFeatureSet(sc_total, pattern = "^mt-", col.name = "percent.mt")
sc_total <- PercentageFeatureSet(sc_total, pattern = "^Rp[sl]", col.name = "percent.rb")

# 过滤
cat("Filtering cells...\n")
sc_total <- subset(sc_total, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 15)
cat(sprintf("Cells after filtering: %d\n", ncol(sc_total)))

# 标准化
cat("Normalizing data...\n")
sc_total <- NormalizeData(sc_total)

# 找可变基因
cat("Finding variable genes...\n")
sc_total <- FindVariableFeatures(sc_total, selection.method = "vst", nfeatures = 3000)

# 缩放数据
cat("Scaling data...\n")
sc_total <- ScaleData(sc_total)

# 降维
cat("Performing PCA...\n")
sc_total <- RunPCA(sc_total, dims = 1:30)

# 找邻居
cat("Finding neighbors...\n")
sc_total <- FindNeighbors(sc_total, dims = 1:30)

# 聚类
cat("Clustering...\n")
sc_total <- FindClusters(sc_total, resolution = 0.8)

# UMAP
cat("Performing UMAP...\n")
sc_total <- RunUMAP(sc_total, dims = 1:30)

# 输出UMAP可视化
umap_output <- file.path(output_dir, "umap_overview.pdf")
pdf(umap_output, width = 10, height = 8)
DimPlot(sc_total, group.by = "group", cols = c("control" = "#2E86AB", "stroke" = "#F24236")) + 
  ggtitle("UMAP Overview (Control vs Stroke)")
dev.off()

cat(sprintf("UMAP visualization saved to: %s\n", umap_output))
cat(">>>> Seurat normalization workflow completed | Key stats: Cells after filtering =", ncol(sc_total), "\n\n")

# 步骤4：细胞类型注释与神经元提取
cat(">>>> Step 4: Cell type annotation and neuron extraction\n")

# 检查基因名称格式
cat("Checking gene name format...\n")
head_genes <- head(rownames(sc_total), 10)
cat("First 10 gene names:\n")
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
    cat(sprintf("%s: Found %d marker genes\n", cell_type, length(present_markers)))
    print(present_markers)
    # 使用存在的标记基因
    markers[[cell_type]] <- present_markers
  } else {
    cat(sprintf("%s: No marker genes found\n", cell_type))
    # 移除没有标记基因的细胞类型
    markers[[cell_type]] <- NULL
  }
}

# 过滤掉没有标记基因的细胞类型
markers <- markers[!sapply(markers, is.null)]

# 添加ModuleScore
for (cell_type in names(markers)) {
  if (length(markers[[cell_type]]) >= 1) {
    cat(sprintf("Calculating ModuleScore for %s...\n", cell_type))
    sc_total <- AddModuleScore(
      sc_total,
      features = list(markers[[cell_type]]),
      name = paste0(cell_type, "_score")
    )
  } else {
    cat(sprintf("Skipping %s, insufficient marker genes\n", cell_type))
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

cat(sprintf("Neuron count: %d, Ratio: %.2f%%\n", neuron_count, neuron_ratio))

# control vs stroke 细胞数对比
control_neurons <- subset(neurons, subset = group == "control")
stroke_neurons <- subset(neurons, subset = group == "stroke")

cat(sprintf("Control group neurons: %d\n", ncol(control_neurons)))
cat(sprintf("Stroke group neurons: %d\n", ncol(stroke_neurons)))

# 样本量前置检查
if(ncol(control_neurons) < 30 | ncol(stroke_neurons) < 30) {
  cat(sprintf("Error: Insufficient neurons: Control=%d, Stroke=%d (min 30 required)\n",
               ncol(control_neurons), ncol(stroke_neurons)))
  flush.console()
  stop("Insufficient neuron count for network analysis")
}

# 检查神经元数量是否足够
if (neuron_count < 50) {
  cat("Warning: Insufficient neurons (<50), cannot perform subsequent analysis\n")
  flush.console()
  quit()
}

cat(">>>> Cell type annotation and neuron extraction completed | Key stats: Neuron count =", neuron_count, "Ratio =", round(neuron_ratio, 2), "%\n\n")

  # 释放内存
  cat("Releasing memory...\n")
  flush.console()
  rm(sc_total, seurat_objects)
  gc()
  
  # 不再使用磁盘缓存，确保每次运行数据一致性
  cat("Skipping cache save to ensure data consistency...\n")
  flush.console()
}

# 步骤5：93基因提取（修复版）
cat(">>>> Step 5: Gene extraction with format auto-detection\n")

# 首先诊断：显示Seurat对象中的实际基因名样本
cat("Sample gene names in Seurat object (first 20):\n")
print(head(rownames(neurons), 20))

# 铜死亡核心基因检查
cat("Copper death core gene check:\n")
core_candidates <- c("fdx1", "lias", "slc31a1", "dlat", "Fdx1", "Lias", "Slc31a1", "Dlat", 
                     "FDX1", "LIAS", "SLC31A1", "DLAT")
for(g in core_candidates) {
  if(g %in% rownames(neurons)) {
    cat(sprintf("  ✓ Found: %s\n", g))
  }
}

# 基因名诊断：在Step 5开头打印表达矩阵维度和名称，确认维度方向正确（应为细胞×基因）
cat("\nGene name diagnosis:\n")
if (exists("exp_matrix")) {
  cat(sprintf("exp_matrix dimensions: %d cells × %d genes\n", nrow(exp_matrix), ncol(exp_matrix)))
  cat("First 10 cell IDs (rownames):\n")
  print(head(rownames(exp_matrix), 10))
  cat("First 10 gene names (colnames):\n")
  print(head(colnames(exp_matrix), 10))
} else {
  cat("exp_matrix not yet created\n")
}

# 加载基因映射表，获取小鼠基因符号
cat("Loading gene mapping table...\n")
gene_mapping_file <- file.path(output_dir, "gene_mapping_93.csv")

# 检查文件存在性
if (!file.exists(gene_mapping_file)) {
  cat("Error: Gene mapping file not exists, please complete step 1 first\n")
  flush.console()
  quit()
}

gene_mapping <- read.csv(gene_mapping_file, stringsAsFactors = FALSE)

# 智能匹配：找出实际格式
target_genes_mouse <- gene_mapping$mouse_symbol[!is.na(gene_mapping$mouse_symbol)]
actual_genes <- rownames(neurons)

cat("Sample gene names in neurons:", paste(head(actual_genes, 5), collapse=", "), "\n")

# 智能匹配：找出实际格式
matched_genes <- c()
for(tg in target_genes_mouse) {
  # 忽略大小写匹配
  idx <- which(tolower(actual_genes) == tolower(tg))
  if(length(idx) > 0) {
    matched_genes <- c(matched_genes, actual_genes[idx[1]])
  }
}
matched_genes <- unique(matched_genes)

cat(sprintf("Matched %d genes\n", length(matched_genes)))

# 强制确保铜死亡核心基因被包含（即使表达低）
core_priority <- c("fdx1", "lias", "slc31a1", "dlat")  # 小写格式
core_found <- c()
for(core in core_priority) {
  # 忽略大小写匹配
  idx <- which(tolower(actual_genes) == tolower(core))
  if(length(idx) > 0) {
    core_found <- c(core_found, actual_genes[idx[1]])
  }
}

if(length(core_found) < 4) {
  cat(sprintf("Warning: Only found %d/4 copper death core genes, trying to find other formats...\n", length(core_found)))
  # 尝试查找任何包含这些基因名片段的基因
  for(core in c("fdx", "lias", "slc31a", "dlat")) {
    partial_match <- grep(core, actual_genes, value=TRUE, ignore.case=TRUE)
    if(length(partial_match) > 0) {
      cat(sprintf("  Found partial match '%s': %s\n", core, paste(partial_match, collapse=", ")))
      matched_genes <- unique(c(matched_genes, partial_match))
    }
  }
}

# 2. 使用Seurat v5安全提取（防御性代码）
cat("Extracting expression matrix...\n")

# 初始化expr_sparse为NULL
expr_sparse <- NULL

# 尝试使用GetAssayData获取完整的表达矩阵（合并数据）
tryCatch({
  cat("Trying GetAssayData for complete expression matrix...\n")
  # Seurat v5标准方法：获取合并后的data层（包含所有细胞）
  expr_sparse <- GetAssayData(neurons, assay = "RNA", layer = "data")
  
  # 验证维度：必须包含所有细胞
  if (ncol(expr_sparse) != ncol(neurons)) {
    stop(sprintf("Dimension mismatch: GetAssayData returned %d cells, but neurons has %d cells", 
                 ncol(expr_sparse), ncol(neurons)))
  }
  
  # 只保留匹配的基因（安全子集：取交集避免报错）
  available_genes <- rownames(expr_sparse)
  genes_to_extract <- intersect(matched_genes, available_genes)
  if (length(genes_to_extract) < length(matched_genes)) {
    cat(sprintf("Warning: Only %d/%d matched genes found in expression matrix\n", 
                length(genes_to_extract), length(matched_genes)))
  }
  expr_sparse <- expr_sparse[genes_to_extract, , drop=FALSE]
  
  cat(sprintf("GetAssayData extraction successful: %d genes × %d cells\n", 
              nrow(expr_sparse), ncol(expr_sparse)))
  
}, error = function(e) {
  cat(sprintf("GetAssayData extraction failed: %s\n", e$message))
  cat("Falling back to JoinLayers + LayerData...\n")
  
  # 降级方案：先合并所有层，再提取
  tryCatch({
    # 合并所有样本层（关键修复）
    neurons <<- JoinLayers(neurons, assay = "RNA")
    
    # 提取默认层（现在已包含所有细胞）
    expr_sparse <<- LayerData(neurons, assay="RNA", features=matched_genes)
    
    cat(sprintf("LayerData extraction successful after JoinLayers: %d genes × %d cells\n", 
                nrow(expr_sparse), ncol(expr_sparse)))
  }, error = function(e2) {
    stop(sprintf("All extraction methods failed: %s", e2$message))
  })
})

# 再次验证细胞数（关键质控）
if (is.null(expr_sparse) || ncol(expr_sparse) == 0) {
  stop("Expression matrix extraction failed: empty result")
}
if (ncol(expr_sparse) != ncol(neurons)) {
  stop(sprintf("Critical error: Extracted %d cells, but neurons object has %d cells", 
               ncol(expr_sparse), ncol(neurons)))
}

# 转置为细胞×基因格式（PC算法标准输入）
exp_matrix <- t(expr_sparse)

# 关键验证：转置后行数必须等于neurons细胞数
cat(sprintf("\nValidation: exp_matrix has %d cells (neurons has %d) - %s\n", 
            nrow(exp_matrix), ncol(neurons), 
            ifelse(nrow(exp_matrix)==ncol(neurons), "MATCH", "MISMATCH")))

# 基因名诊断
cat("\nGene name diagnosis after extraction:\n")
cat(sprintf("exp_matrix dimensions: %d cells × %d genes\n", 
            nrow(exp_matrix), ncol(exp_matrix)))
cat("First 10 cell IDs (rownames):\n")
print(head(rownames(exp_matrix), 10))
cat("First 10 gene names (colnames):\n")
print(head(colnames(exp_matrix), 10))

# 3. 智能过滤：保留高表达基因 + 强制保留核心基因（即使低表达）
gene_counts <- colSums(exp_matrix > 0) 
high_genes <- names(gene_counts[gene_counts >= 5]) 

# 确保核心基因常量存在（防御缓存缺失）
if(!exists("CORE_CUPROPTOSIS")) {
  CORE_CUPROPTOSIS <- c("fdx1", "lias", "dlat", "slc31a1", "dld", "lipt1", "pdhx", "pdhb")
}

# 核心基因：改为硬编码小写列表+intersect，确保8个核心基因全在
# 强制保留铜死亡核心基因（硬编码小写列表）
must_keep <- intersect(CORE_CUPROPTOSIS, colnames(exp_matrix))
cat(sprintf("Forcing to keep %d copper death core genes\n", length(must_keep)))

# 合并（使用实际检测到的名称，而非硬编码） 
keep_genes <- unique(c(high_genes, must_keep)) 
exp_matrix <- exp_matrix[, keep_genes, drop=FALSE] 

cat(sprintf("Final matrix: %d cells x %d genes (core genes: %d)\n", 
             nrow(exp_matrix), ncol(exp_matrix), length(must_keep)))

# 检查铜死亡核心基因是否存在
missing_core <- setdiff(CORE_CUPROPTOSIS, colnames(exp_matrix))
if(length(missing_core) > 0) {
  cat("Warning: Missing copper death core genes: ", paste(missing_core, collapse=", "), "\n")
  flush.console()
} else {
  cat("All copper death core genes found!\n")
  flush.console()
}

# 保存key_genes变量，避免后续代码出错
key_genes <- c()

cat(">>>> Gene extraction and preprocessing completed | Key stats: Final genes =", ncol(exp_matrix), "\n\n")

# 步骤6：分块PC网络（修复版）
cat(">>>> Step 6: Modular PC network analysis\n")

# 安装并加载必要的包
cat("Installing and loading pcalg and dependencies...\n")
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
  cat("pcalg package loaded successfully\n")
  flush.console()
}, error = function(e) {
  cat(sprintf("Package installation failed: %s\n", e$message))
  flush.console()
  stop("Cannot install necessary packages, please install pcalg and its dependencies manually")
})

# 创建与exp_matrix行数匹配的逻辑向量
exp_cells <- rownames(exp_matrix)
cat("Matching cell indices...\n")
cat(sprintf("exp_matrix cells: %d\n", length(exp_cells)))
cat(sprintf("neurons cells: %d\n", ncol(neurons)))

# 细胞对齐：使用细胞条形码直接索引，零NA风险
# 确保exp_cells中的所有细胞都在neurons中
common_cells <- intersect(exp_cells, colnames(neurons))
if (length(common_cells) < length(exp_cells)) {
  cat(sprintf("Warning: %d cells not found in neurons object, filtering\n", length(exp_cells) - length(common_cells)))
  exp_matrix <- exp_matrix[common_cells, , drop=FALSE]
  exp_cells <- common_cells
}

# 方法：命名向量索引确保顺序严格对应exp_cells
group_vec <- setNames(neurons$group, colnames(neurons))
cell_groups <- group_vec[exp_cells]  # 按exp_cells顺序提取，无NA风险
# 验证无NA（应全为control/stroke）
if(any(is.na(cell_groups))) stop("Group assignment failed: NA detected")

# 维度验证：使用stopifnot强制验证cells×genes格式
if(nrow(exp_matrix) != length(cell_groups)) stop("Cell groups length does not match exp_matrix rows")
if(nrow(exp_matrix) != length(exp_cells)) stop("exp_cells length does not match exp_matrix rows")

# 检查cell_groups的分布
cat("Cell groups distribution:\n")
print(table(cell_groups))

# 动态调整模块定义（匹配实际基因名格式，增加容错）
adjust_to_actual <- function(gene_list, available_genes) {
  # 输入输出均为小写
  result <- c()
  for(g in gene_list) {
    # 完全匹配
    idx <- which(available_genes == g)
    if(length(idx) == 0) {
      # 模糊匹配（针对基因名变体如mt-nd1）
      idx <- grep(paste0("^", g, "$"), available_genes, ignore.case=TRUE)
    }
    if(length(idx) > 0) result <- c(result, available_genes[idx[1]])
  }
  return(unique(result))
}

available <- colnames(exp_matrix)

# 基于网络药理学结果（KEGG通路+PPI Hub）的6模块精细划分
# 解决原others模块过大（71基因）导致PC卡死的问题

# 先定义核心基因列表（确保完整保留）
core_age_rage <- c("nfkb1", "rela", "stat1", "stat3", "ccl2", "icam1", "tgfb1", "pparg", "tlr4", "il6")
core_cuproptosis <- c("fdx1", "lias", "dlat", "dld", "lipt1", "pdhx", "pdhb")
core_homeostasis <- c("slc31a1", "atp7b", "atp7a", "atox1", "commd1", "mt2a", "cp")
core_inflammation <- c("ptgs2", "parp1", "hmox1", "fabp3", "fabp5", "cpt1a", "sod2")
core_death <- c("casp8", "fas", "f3", "mapk9", "gad1", "gfap")
core_transcription <- c("egr1", "atf4", "stat5a", "brd4", "notch1", "nfe2l2", "srebf1")

# 动态匹配实际存在的基因
gene_modules <- list(
  # 模块1：AGE-RAGE信号轴（Hub基因+通路基因，基于申报书2.2.2 KEGG结果）
  age_rage_axis = adjust_to_actual(core_age_rage, available),
  
  # 模块2：铜死亡执行核心（必须完整保留，申报书2.3.1）
  cuproptosis_execution = adjust_to_actual(core_cuproptosis, available),
  
  # 模块3：铜稳态调控（上游转运）
  copper_homeostasis = adjust_to_actual(core_homeostasis, available),
  
  # 模块4：炎症-代谢交叉（剩余Hub基因和对接靶点，申报书表2）
  inflammation_metabolism = adjust_to_actual(core_inflammation, available),
  
  # 模块5：细胞死亡调控（凋亡/焦亡交叉）
  cell_death_regulation = adjust_to_actual(core_death, available),
  
  # 模块6：转录调控（GO富集RNA聚合酶II相关，申报书2.2.1）
  transcription_regulation = adjust_to_actual(core_transcription, available)
)

# 打印模块统计（增加警告提示）
cat("Module gene format adjustment completed:\n")
for(n in names(gene_modules)) {
  cat(sprintf("  Module %s: %d genes", n, length(gene_modules[[n]])))
  if (length(gene_modules[[n]]) > 15) {
    cat(sprintf(" [WARNING: >15 genes, may be slow]"))
  } else if (length(gene_modules[[n]]) < 3) {
    cat(sprintf(" [WARNING: <3 genes, will skip PC]"))
  }
  cat("\n")
}
flush.console()

# 2. 分块构建网络（每块基因数<25，计算极快）
build_module_network <- function(module_genes, module_name, data) {
  start_time <- Sys.time()
  
  tryCatch({
    # 增加基因数硬性限制（保险机制）
    if (length(module_genes) > 20) {
      cat(sprintf("Module %s: Too many genes (%d), limiting to top 20 by variance\n", 
                  module_name, length(module_genes)))
      # 按方差排序取前20（保留变化最大的基因）
      gene_vars <- apply(as.matrix(data[, module_genes]), 2, var, na.rm=TRUE)
      module_genes <- names(sort(gene_vars, decreasing = TRUE))[1:20]
    }
    
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
      cat(sprintf("Module %s: Insufficient genes (%d), skipping\n", module_name, length(present_genes)))
      flush.console()
      return(NULL)
    }
    
    module_data <- data[, present_genes, drop=FALSE]
    
    # 根据样本量动态调整alpha
    n_cells <- nrow(module_data)
    alpha <- ifelse(n_cells < 100, 0.1, ifelse(n_cells > 500, 0.01, ALPHA_BASE))
    cat(sprintf("Module %s: %d cells, using alpha=%.3f\n", module_name, n_cells, alpha))
    flush.console()
    
    # 临时转为dense矩阵用于计算，若细胞数过多则采样
    if (n_cells > 10000) {
      cat("Too many cells, sampling 10000 cells for correlation calculation...\n")
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
    
    # 在return前增加运行时间检查
    elapsed <- difftime(Sys.time(), start_time, units="secs")
    if (elapsed > 60) {
      cat(sprintf("Module %s: Warning - took %.1f seconds\n", module_name, elapsed))
    }
    
    cat(sprintf("Module %s: %d genes, %d edges\n", module_name, length(present_genes), length(edges)))
    flush.console()
    return(list(adj=adj, genes=present_genes))
  }, error = function(e) {
    cat(sprintf("Module %s calculation failed: %s\n", module_name, e$message))
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
cat("\n==== Stroke group diagnosis ====\n")
stroke_idx <- cell_groups == "stroke"
cat(sprintf("Stroke group cell count: %d\n", sum(stroke_idx)))
if(sum(stroke_idx) > 0) {
    stroke_data <- exp_matrix[stroke_idx, , drop=FALSE]
    cat(sprintf("Stroke group expression matrix dimensions: %d cells × %d genes\n", nrow(stroke_data), ncol(stroke_data)))
    cat("Core gene expression in Stroke group:\n")
    for(g in must_keep) {  # 使用must_keep而非core_actual_names
      if(g %in% colnames(stroke_data)) {
        expr_cells <- sum(stroke_data[, g] > 0)
        cat(sprintf("  %s: %d/%d cells expressed (%.1f%%)\n", 
                    g, expr_cells, nrow(stroke_data), 100*expr_cells/nrow(stroke_data)))
      }
    }
    if(sum(stroke_idx) < 30) {
      cat("Warning: Stroke group cell count < 30, PC algorithm cannot run (insufficient sample size)\n")
    }
  }
cat("====================\n\n")

for (group in c("control", "stroke")) {
  tryCatch({
    idx <- cell_groups == group
    group_data <- exp_matrix[idx, , drop=FALSE]
    
    cat(sprintf("\nBuilding %s group network (%d cells):\n", group, nrow(group_data)))
    flush.console()
    
    # 单核顺序计算（4个模块，稳定性优先）
    cat("Processing 4 modules sequentially...\n")
    flush.console()
    
    networks <- list()
    for (mod_name in names(gene_modules)) {
      cat(sprintf("  → Module %s: ", mod_name))
      flush.console()
      
      networks[[mod_name]] <- build_module_network(
        gene_modules[[mod_name]],
        mod_name,
        group_data
      )
      
      # 如果模块成功，输出边数
      if (!is.null(networks[[mod_name]])) {
        adj <- networks[[mod_name]]$adj
        n_edges <- sum(adj[upper.tri(adj)] == 1)
        cat(sprintf("%d edges\n", n_edges))
      } else {
        cat("skipped\n")
      }
      flush.console()
    }
    
    # 合并边（原有逻辑保持不变）
    cat("Merging module edges...\n")
    flush.console()
    
    all_edges <- data.frame(from=character(), to=character(), weight=numeric())
    
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
    
    # 去重（原有逻辑）
    if (nrow(all_edges) > 0) {
      all_edges$edge_id <- apply(all_edges[, c("from", "to")], 1, 
                               function(x) paste(sort(x), collapse="_"))
      unique_indices <- !duplicated(all_edges$edge_id)
      all_edges <- all_edges[unique_indices, c("from", "to", "weight")]
    }
    
    # 保存（原有逻辑）
    write.csv(all_edges, 
              file.path(output_dir, sprintf("network_%s_edges.csv", group)), 
              row.names=FALSE)
    all_edges_list[[group]] <- all_edges
    
    cat(sprintf("%s group: %d unique edges saved\n", group, nrow(all_edges)))
    flush.console()
    
  }, error = function(e) {
    cat(sprintf("%s group network construction failed: %s\n", group, e$message))
    flush.console()
    # 保存空文件
    all_edges <- data.frame(from=character(), to=character(), weight=numeric())
    write.csv(all_edges, 
              file.path(output_dir, sprintf("network_%s_edges.csv", group)), 
              row.names=FALSE)
    all_edges_list[[group]] <- all_edges
  })
}

cat(">>>> PC causal network analysis completed\n\n")

# 步骤7：差异网络分析
cat(">>>> Step 7: Differential network analysis\n")

# 检查两组网络是否都构建成功
if (!is.null(all_edges_list$control) && !is.null(all_edges_list$stroke)) {
  # 提取两组的边
  control_edges <- all_edges_list$control
  stroke_edges <- all_edges_list$stroke
  
  # 创建边的唯一标识符（添加tolower确保FDX1-LIAS与fdx1-lias被视为同一边）
  control_edges$edge_id <- apply(control_edges[, c("from", "to")], 1, function(x) paste(sort(tolower(x)), collapse="_"))
  stroke_edges$edge_id <- apply(stroke_edges[, c("from", "to")], 1, function(x) paste(sort(tolower(x)), collapse="_"))
  
  # 识别Stroke特异性新边（在stroke中存在，在control中不存在）
  stroke_specific_edges <- stroke_edges[!stroke_edges$edge_id %in% control_edges$edge_id, ]
  
  if (nrow(stroke_specific_edges) > 0) {
    # 移除edge_id列，重命名列名
    stroke_specific_edges <- stroke_specific_edges[, c("from", "to", "weight")]
    colnames(stroke_specific_edges) <- c("From", "To", "weight")
    stroke_specific_edges$Type <- "Stroke_Specific"
    
    # 输出novel_stroke_edges.csv
    edges_output <- file.path(output_dir, "novel_stroke_edges.csv")
    write.csv(stroke_specific_edges, edges_output, row.names=FALSE)
    cat(sprintf("Stroke-specific edges saved to: %s\n", edges_output))
    flush.console()
    
    # 统计：stroke特异性边数
    cat(sprintf("Stroke-specific edge count: %d\n", nrow(stroke_specific_edges)))
    flush.console()
  } else {
    cat("No Stroke-specific edges found\n")
    flush.console()
    # 创建空文件
    edges_output <- file.path(output_dir, "novel_stroke_edges.csv")
    write.csv(data.frame(From = character(), To = character(), Type = character(), weight = numeric()), 
              edges_output, row.names=FALSE)
  }
} else {
  cat("At least one group network construction failed, cannot perform differential analysis\n")
  flush.console()
  # 创建空文件
  edges_output <- file.path(output_dir, "novel_stroke_edges.csv")
  write.csv(data.frame(From = character(), To = character(), Type = character(), weight = numeric()), 
            edges_output, row.names=FALSE)
}

cat(">>>> Differential network analysis completed\n\n")

# 验证所有输出文件是否生成
cat(">>>> Verifying output files\n")
output_files <- c(
  file.path(output_dir, "gene_mapping_93.csv"),
  file.path(output_dir, "umap_overview.pdf"),
  file.path(output_dir, "network_control_edges.csv"),
  file.path(output_dir, "network_stroke_edges.csv"),
  file.path(output_dir, "novel_stroke_edges.csv")
)

for (file in output_files) {
  if (file.exists(file)) {
    cat(sprintf("✓ %s generated\n", basename(file)))
  } else {
    cat(sprintf("✗ %s not generated\n", basename(file)))
  }
  flush.console()
}

cat("\nAnalysis completed!\n")
flush.console()