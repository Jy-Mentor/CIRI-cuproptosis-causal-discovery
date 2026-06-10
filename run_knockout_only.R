# 仅运行scTenifoldKnk虚拟敲除 (绕过QC)
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))

library(Seurat)
library(Matrix)
library(scTenifoldNet)
library(scTenifoldKnk)
library(dplyr)

result_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/analysis_results"

cat("=== scTenifoldKnk虚拟敲除 (跳过QC) ===\n\n")

# 1. 读取之前保存的矩阵
cat("1. 读取count矩阵...\n")
count_file <- file.path(result_dir, "sc_count_matrix.rds")

if(file.exists(count_file)) {
  count_matrix <- readRDS(count_file)
  cat(sprintf("  矩阵: %d 基因 x %d 细胞\n", nrow(count_matrix), ncol(count_matrix)))
} else {
  # 如果没有，从Seurat对象创建
  cat("  从Seurat对象提取...\n")
  sc_file <- file.path(result_dir, "sc_annotated.rds")
  if(!file.exists(sc_file)) {
    cat("错误: 找不到数据文件\n")
    quit(status=1)
  }
  sc_obj <- readRDS(sc_file)
  count_matrix <- as.matrix(GetAssayData(sc_obj, layer="counts"))
  keep_genes <- rowSums(count_matrix > 0) >= 50
  count_matrix <- count_matrix[keep_genes, ]
  saveRDS(count_matrix, count_file)
  cat(sprintf("  矩阵: %d 基因 x %d 细胞\n", nrow(count_matrix), ncol(count_matrix)))
}

# 2. 目标基因
target_genes <- c("Nfkb1", "Fdx1", "Tlr4")

# 3. 执行虚拟敲除 (qc=FALSE)
cat("\n2. 执行虚拟敲除...\n")

for(gene in target_genes) {
  if(gene %in% rownames(count_matrix)) {
    cat(sprintf("\n  虚拟敲除 %s...\n", gene))
    cat("    (这可能需要15-40分钟，请耐心等待)...\n")

    # 手动QC - 过滤低质量细胞
    cat("    预处理...\n")
    min_lib_size <- 1000
    lib_sizes <- colSums(count_matrix)
    valid_cells <- lib_sizes >= min_lib_size
    cm <- count_matrix[, valid_cells]
    cat(sprintf("    质控后: %d 细胞\n", ncol(cm)))

    # 过滤低表达基因
    gene_exp <- rowMeans(cm > 0)
    cm <- cm[gene_exp >= 0.05, ]
    cat(sprintf("    过滤基因后: %d 基因 x %d 细胞\n", nrow(cm), ncol(cm)))

    set.seed(666)
    knk_res <- scTenifoldKnk(
      countMatrix = cm,
      gKO = gene,
      qc = FALSE,  # 禁用内部QC
      nc_lambda = 0,
      nc_nNet = 10,
      nc_nCells = 500,
      nc_nComp = 3,
      nc_q = 0.9,
      td_K = 3,
      ma_nDim = 2,
      nCores = 2
    )

    # 处理结果
    cat("    保存结果...\n")
    degs <- knk_res$diffRegulation %>%
      filter(p.value < 0.05) %>%
      mutate(logFC = log2(FC)) %>%
      arrange(p.value)

    write.csv(degs, file.path(result_dir, paste0("knockout_", gene, "_DEGs.csv")), row.names=FALSE)
    cat(sprintf("    ✅ 差异基因数: %d\n", nrow(degs)))

    # 显示top10
    cat("    Top 10差异基因:\n")
    for(i in 1:min(10, nrow(degs))) {
      direction <- ifelse(degs$logFC[i] > 0, "↑", "↓")
      cat(sprintf("      %s: %s%.3f, p=%.2e\n", degs$gene[i], direction, abs(degs$logFC[i]), degs$p.value[i]))
    }

    saveRDS(knk_res, file.path(result_dir, paste0("knockout_", gene, "_result.rds")))
    cat("    ✅ 完成\n")

  } else {
    cat(sprintf("  [跳过] %s 不在数据中\n", gene))
  }
}

cat("\n=== ✅ 完成 ===\n")
cat("结果保存到:", result_dir, "\n")
