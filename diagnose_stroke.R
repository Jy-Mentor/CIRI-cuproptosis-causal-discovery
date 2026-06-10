#!/usr/bin/env Rscript

# 加载缓存数据进行诊断
cache_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/analysis_cache.rds"

# 确保Seurat包已加载
if (!require("Seurat", character.only = TRUE)) {
  install.packages("Seurat", dependencies = TRUE)
  library(Seurat)
}

# 加载缓存数据
cat("Loading cache data...\n")
cached_data <- readRDS(cache_file)
neurons <- cached_data$neurons

# 检查步骤4后的神经元对象
cat("\n=== Checking neuron group distribution ===\n")
group_cell_type_table <- table(neurons$group, neurons$cell_type)
print(group_cell_type_table)

# 计算各组细胞数
cat("\n=== Cell count per group ===\n")
group_counts <- table(neurons$group)
print(group_counts)

# 检查是否有stroke组细胞
if ("stroke" %in% names(group_counts)) {
  stroke_count <- group_counts["stroke"]
  cat(sprintf("Stroke group cell count: %d\n", stroke_count))
  if (stroke_count == 0) {
    cat("WARNING: No cells in stroke group!\n")
  }
} else {
  cat("WARNING: No stroke group found!\n")
}

# 检查细胞类型分布
cat("\n=== Cell type distribution ===\n")
cell_type_counts <- table(neurons$cell_type)
print(cell_type_counts)

# 检查原始sc_total对象（如果存在）
if (file.exists("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/sc_total.rds")) {
  cat("\n=== Checking original sc_total object ===\n")
  sc_total <- readRDS("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/sc_total.rds")
  sc_total_table <- table(sc_total$group, sc_total$cell_type)
  print(sc_total_table)
} else {
  cat("\n=== sc_total.rds not found, skipping ===\n")
  # 尝试从缓存中加载sc_total（如果存在）
  if ("sc_total" %in% names(cached_data)) {
    cat("Loading sc_total from cache...\n")
    sc_total <- cached_data$sc_total
    sc_total_table <- table(sc_total$group, sc_total$cell_type)
    print(sc_total_table)
  }
}

cat("\nDiagnosis completed!\n")