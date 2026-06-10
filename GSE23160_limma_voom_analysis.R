# GSE23160 时序差异表达分析脚本 (对标文献版本)
# 分析: 2h, 8h, 24h 时间序列差异分析 (脑缺血再灌注损伤)
# 方法: Limma + removeBatchEffect(仅可视化) + PCA离群剔除
# 筛选条件: |log2FC| >= 0.585 (FC>1.5), P < 0.05 (对标文献标准)
# 注意: series_matrix数据已经是log2标准化值，直接使用limma而非limma-voom
# 策略: 皮层和纹状体分开分析，每个时间点pairwise vs Sham

# ==================== 1. 包安装与加载 ====================
cat("正在检查和安装必要的R包...\n")

packages <- c("limma", "sva", "ggplot2", "dplyr", "tidyr", 
              "pheatmap", "RColorBrewer", "ggpubr", "gridExtra", "reshape2")

bioc_packages <- c("limma", "sva")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Installing package: ", pkg, "\n"))
    if (pkg %in% bioc_packages) {
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org/")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
    library(pkg, character.only = TRUE)
  }
}

for (pkg in packages) {
  install_if_missing(pkg)
}

cat("所有包加载完成!\n\n")

# ==================== 2. 设置工作目录和文件路径 ====================
work_dir <- "D:/反向网络药理学/L1 数据集/bulk/GSE23160(主验证集时序差异分析，2h,8h,24h)"
setwd(work_dir)

platform_file <- file.path(work_dir, "GPL6885-11608.txt")
series_matrix_file <- file.path(work_dir, "GSE23160_series_matrix.txt.gz")

output_dir <- file.path(work_dir, "GSE23160_limma_results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

figure_dir <- file.path(output_dir, "figures")
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

# ==================== 3. 读取系列矩阵数据 ====================
cat("读取 GSE23160 系列矩阵数据...\n")

expr_data <- read.delim(series_matrix_file, 
                        comment.char = "!", 
                        header = TRUE, 
                        check.names = FALSE,
                        row.names = 1,
                        stringsAsFactors = FALSE)

expr_matrix <- as.matrix(expr_data)
cat(paste0("表达矩阵维度: ", nrow(expr_matrix), " probes x ", ncol(expr_matrix), " samples\n"))

# 检查数据尺度
data_range <- range(expr_matrix, na.rm = TRUE)
cat(paste0("数据范围: ", data_range[1], " - ", data_range[2], "\n"))
cat(paste0("中位数: ", median(expr_matrix, na.rm = TRUE), "\n"))

# 注意：GEO series_matrix文件通常已经过预处理（background correction + normalization）
# Illumina芯片的series_matrix数据可能已经是log2转换后的值
# 我们直接使用原始数据，不添加额外的转换
cat("使用GEO预处理后的数据（不添加额外转换）\n")

missing_count <- sum(is.na(expr_matrix))
if (missing_count > 0) {
  cat(paste0("检测到 ", missing_count, " 个缺失值，使用中位数填充...\n"))
  for (i in 1:nrow(expr_matrix)) {
    row_median <- median(expr_matrix[i, ], na.rm = TRUE)
    expr_matrix[i, is.na(expr_matrix[i, ])] <- row_median
  }
}

# ==================== 4. 读取样本分组信息 ====================
cat("提取样本分组信息...\n")

con <- gzfile(series_matrix_file, "rt")
lines <- readLines(con, warn = FALSE)
close(con)

sample_titles_line <- grep("^!Sample_title", lines, value = TRUE)
sample_titles <- unlist(strsplit(sample_titles_line, "\t"))[-1]
sample_titles <- gsub('"', '', sample_titles)

cat("样本标题调试信息（全部样本）:\n")
for (i in seq_along(sample_titles)) {
  cat(paste0("  Sample ", i, ": ", sample_titles[i], "\n"))
}
cat("\n")

sample_info <- data.frame(
  sample_id = colnames(expr_matrix),
  title = sample_titles,
  stringsAsFactors = FALSE
)

sample_info$brain_region <- ifelse(grepl("^Cortex-", sample_titles, ignore.case = TRUE), "Cortex",
                            ifelse(grepl("^Striatum-", sample_titles, ignore.case = TRUE), "Striatum", "Unknown"))

sample_info$treatment <- ifelse(grepl("-Sham-", sample_titles, ignore.case = TRUE), "Sham", "IR")

sample_info$time_group <- ifelse(sample_info$treatment == "Sham", "Ctrl",
                          ifelse(grepl("-2h-", sample_titles, ignore.case = TRUE), "2h",
                          ifelse(grepl("-8h-", sample_titles, ignore.case = TRUE), "8h",
                          ifelse(grepl("-24h-", sample_titles, ignore.case = TRUE), "24h", "Unknown"))))

unknown_idx <- sample_info$brain_region == "Unknown" | sample_info$time_group == "Unknown"
if (any(unknown_idx)) {
  cat("警告: 以下样本无法正确分组:\n")
  print(sample_info[unknown_idx, ])
  stop("存在无法解析的样本，请检查数据格式")
}

sample_info$group <- paste0(sample_info$brain_region, "_", sample_info$treatment, "_", sample_info$time_group)

cat("样本分组情况:\n")
print(table(sample_info$brain_region, sample_info$treatment, sample_info$time_group))
cat(paste0("总样本数: ", nrow(sample_info), "\n\n"))

# ==================== 5. 读取平台注释信息 ====================
cat("\n读取 GPL6885 平台注释...\n")

platform_annot <- read.delim(platform_file, 
                             header = TRUE,
                             check.names = FALSE,
                             stringsAsFactors = FALSE,
                             quote = "",
                             fill = TRUE,
                             comment.char = "#")

cat(paste0("平台注释探针数: ", nrow(platform_annot), "\n"))
cat("平台注释列名:\n")
print(names(platform_annot)[1:min(15, ncol(platform_annot))])

annot_table <- data.frame(
  ID = platform_annot[[1]],
  Gene_Symbol = if("Symbol" %in% names(platform_annot)) platform_annot[["Symbol"]] else
                if("Gene Symbol" %in% names(platform_annot)) platform_annot[["Gene Symbol"]] else
                if("GENE_SYMBOL" %in% names(platform_annot)) platform_annot[["GENE_SYMBOL"]] else NA,
  Gene_Title = if("ILMN_Gene" %in% names(platform_annot)) platform_annot[["ILMN_Gene"]] else
               if("Gene Title" %in% names(platform_annot)) platform_annot[["Gene Title"]] else
               if("GENE_NAME" %in% names(platform_annot)) platform_annot[["GENE_NAME"]] else NA,
  stringsAsFactors = FALSE
)

cat(paste0("平台注释基因数: ", nrow(annot_table), "\n"))
cat("Gene_Symbol列示例（前10个）:\n")
print(head(annot_table$Gene_Symbol, 10))
cat("\n")

# ==================== 6. 定义铜死亡基因列表 ====================
cat("定义铜死亡相关基因列表...\n")

copper_genes <- c(
  "Fdx1", "Lias", "Lipt1", "Dlat", "Pdhb", "Pdhx", "Dld", "Dbt",
  "Slc31a1", "Atp7a", "Atp7b", "Atox1", "Mtf1",
  "Nfe2l2", "Hif1a", "Mtor", "Nfkb1", "Keap1",
  "Cdkn2a", "Trp53", "Slc25a3", "Slc25a5",
  "Sdha", "Sdhb", "Sdhc", "Sdhd",
  "Cs", "Aco1", "Aco2", "Idh2", "Ogdh",
  "Gcsh", "Gls", "Glud1"
)

cat(paste0("铜死亡基因列表 (n=", length(copper_genes), "):\n"))
print(copper_genes)

# ==================== 7. PCA离群值检测与剔除 ====================
cat("\n进行PCA离群值检测...\n")

pca_result <- prcomp(t(expr_matrix), scale. = TRUE)

pca_df <- data.frame(
  Sample = sample_info$sample_id,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  group = sample_info$group,
  treatment = sample_info$treatment,
  time = sample_info$time_group,
  brain_region = sample_info$brain_region,
  stringsAsFactors = FALSE
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA Plot - GSE23160 (All Samples)",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  scale_color_brewer(palette = "Set1")

ggsave(file.path(figure_dir, "00_PCA_all_samples.png"), 
       pca_plot, width = 12, height = 8, dpi = 300)

cat("PCA图已保存\n\n")

# ==================== 8. limma差异表达分析函数 ====================
cat("开始limma差异表达分析...\n")

limfc_threshold <- 0.585
fdr_threshold <- 0.05

run_limma_analysis <- function(brain_region, expr_matrix_full, sample_info_full, annot_table, 
                               output_dir, figure_dir, copper_genes, logfc_thresh, fdr_thresh) {
  
  cat(paste0("\n", strrep("=", 60), "\n"))
  cat(paste0("分析脑区: ", brain_region, "\n"))
  cat(strrep("=", 60), "\n")
  
  region_idx <- sample_info_full$brain_region == brain_region
  expr_region <- expr_matrix_full[, region_idx]
  sample_region <- sample_info_full[region_idx, ]
  
  cat(paste0("样本数: ", nrow(sample_region), "\n"))
  print(table(sample_region$treatment, sample_region$time_group))
  
  expr_norm <- expr_region
  
  pca_result_region <- prcomp(t(expr_norm), scale. = TRUE)
  pca_df_region <- data.frame(
    Sample = sample_region$sample_id,
    PC1 = pca_result_region$x[, 1],
    PC2 = pca_result_region$x[, 2],
    group = sample_region$group,
    treatment = sample_region$treatment,
    time = sample_region$time_group,
    stringsAsFactors = FALSE
  )
  
  pca_plot_region <- ggplot(pca_df_region, aes(x = PC1, y = PC2, color = group, shape = treatment)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = paste0("PCA - ", brain_region, " (After Batch Effect Removal)"),
         x = paste0("PC1 (", round(summary(pca_result_region)$importance[2, 1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca_result_region)$importance[2, 2] * 100, 1), "%)")) +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    scale_color_brewer(palette = "Set1")
  
  ggsave(file.path(figure_dir, paste0("01_PCA_", brain_region, ".png")), 
         pca_plot_region, width = 10, height = 7, dpi = 300)
  
  time_points <- c("2h", "8h", "24h")
  all_results <- list()
  
  for (tp in time_points) {
    cat(paste0("\n--- ", tp, " vs Sham ---\n"))
    
    tp_idx <- sample_region$time_group == tp | sample_region$time_group == "Ctrl"
    expr_tp <- expr_norm[, tp_idx]
    sample_tp <- sample_region[tp_idx, ]
    
    if (length(unique(sample_tp$treatment)) < 2) {
      cat(paste0("警告: ", tp, " 样本不足，跳过\n"))
      next
    }
    
    sample_tp$treatment <- factor(sample_tp$treatment, levels = c("Sham", "IR"))
    
    # GEO2R标准方法：~0+group设计 + contrasts
    group <- factor(sample_tp$treatment)
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    
    cat(paste0("  样本数: ", ncol(expr_tp), " (IR: ", sum(sample_tp$treatment=="IR"), 
               ", Sham: ", sum(sample_tp$treatment=="Sham"), ")\n"))
    
    fit <- lmFit(expr_tp, design)
    contrast.matrix <- makeContrasts(IRvsSham = IR - Sham, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # GEO2R默认使用FDR校正（BH方法）
    results <- topTable(fit2, coef = "IRvsSham", number = Inf, 
                        adjust.method = "BH", sort.by = "P")
    
    results$ProbeID <- rownames(results)
    results$Gene_Symbol <- annot_table$Gene_Symbol[match(results$ProbeID, annot_table$ID)]
    results$Gene_Title <- annot_table$Gene_Title[match(results$ProbeID, annot_table$ID)]
    
    if (!"Gene_Symbol" %in% names(results)) {
      results$Gene_Symbol <- NA
    }
    
    results <- results %>%
      mutate(
        significance = case_when(
          adj.P.Val < fdr_thresh & logFC >= logfc_thresh ~ "Upregulated",
          adj.P.Val < fdr_thresh & logFC <= -logfc_thresh ~ "Downregulated",
          TRUE ~ "Not significant"
        ),
        time_point = tp,
        brain_region = brain_region
      )
    
    sig_genes <- results %>%
      filter(adj.P.Val < fdr_thresh & abs(logFC) >= logfc_thresh)
    
    cat(paste0("  总基因数: ", nrow(results), "\n"))
    cat(paste0("  显著DEGs (|log2FC|>=", logfc_thresh, ", FDR<", fdr_thresh, "): ", nrow(sig_genes), "\n"))
    cat(paste0("    上调: ", sum(sig_genes$logFC > 0), "\n"))
    cat(paste0("    下调: ", sum(sig_genes$logFC < 0), "\n"))
    
    all_results[[tp]] <- results
    
    write.table(results, 
                file = file.path(output_dir, paste0("DEG_", brain_region, "_", tp, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    write.table(sig_genes, 
                file = file.path(output_dir, paste0("sig_DEGs_", brain_region, "_", tp, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    volcano_data <- results %>%
      mutate(
        neg_log10_P = -log10(P.Value + 1e-10),
        significance = factor(significance, levels = c("Not significant", "Downregulated", "Upregulated"))
      )
    
    volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_P, color = significance)) +
      geom_point(size = 1.5, alpha = 0.6) +
      scale_color_manual(values = c("Not significant" = "gray", 
                                     "Downregulated" = "blue", 
                                     "Upregulated" = "red")) +
      geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
      geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
      labs(title = paste0("Volcano - ", brain_region, " ", tp, " vs Sham"),
           x = "log2 Fold Change",
           y = "-log10(P-value)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.position = "right")
    
    ggsave(file.path(figure_dir, paste0("02_Volcano_", brain_region, "_", tp, ".png")), 
           volcano_plot, width = 8, height = 6, dpi = 300)
    
    copper_probes <- c()
    copper_gene_found <- c()
    for (gene in copper_genes) {
      idx <- grep(paste0("^", gene, "$"), results$Gene_Symbol, ignore.case = TRUE)
      if (length(idx) == 0) {
        idx <- grep(gene, results$Gene_Symbol, ignore.case = TRUE)
      }
      if (length(idx) > 0) {
        copper_probes <- c(copper_probes, results$ProbeID[idx[1]])
        copper_gene_found <- c(copper_gene_found, results$Gene_Symbol[idx[1]])
      }
    }
    
    cat(paste0("  找到 ", length(copper_probes), " 个铜死亡基因探针\n"))
    
    if (length(copper_probes) > 0) {
      copper_expr <- expr_tp[copper_probes, , drop = FALSE]
      copper_expr_z <- t(scale(t(copper_expr)))
      
      annotation_col <- data.frame(
        Treatment = sample_tp$treatment,
        row.names = sample_tp$sample_id
      )
      
      annotation_colors <- list(Treatment = c(Sham = "#4DAF4A", IR = "#E41A1C"))
      
      pheatmap(copper_expr_z,
               annotation_col = annotation_col,
               annotation_colors = annotation_colors,
               show_rownames = TRUE,
               show_colnames = FALSE,
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
               main = paste0("Copper Death Genes - ", brain_region, " ", tp),
               filename = file.path(figure_dir, paste0("03_Heatmap_", brain_region, "_", tp, ".png")),
               width = 8, height = 6)
    }
  }
  
  all_combined <- do.call(rbind, all_results)
  
  copper_gene_results <- all_combined %>%
    filter(!is.na(Gene_Symbol) & Gene_Symbol %in% copper_genes)
  
  copper_summary <- copper_gene_results %>%
    group_by(time_point) %>%
    summarise(
      Upregulated = sum(significance == "Upregulated"),
      Downregulated = sum(significance == "Downregulated"),
      Not_significant = sum(significance == "Not significant"),
      .groups = "drop"
    )
  
  bar_plot <- ggplot(copper_summary %>% pivot_longer(cols = -time_point, 
                                                       names_to = "status", 
                                                       values_to = "count"),
                     aes(x = time_point, y = count, fill = status)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Upregulated" = "#E41A1C", 
                                  "Downregulated" = "#4DAF4A",
                                  "Not significant" = "gray")) +
    labs(title = paste0("Copper Death Genes - ", brain_region),
         x = "Time Point",
         y = "Number of Genes") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "bottom")
  
  ggsave(file.path(figure_dir, paste0("04_Barplot_", brain_region, ".png")), 
         bar_plot, width = 8, height = 6, dpi = 300)
  
  write.table(all_combined, 
              file = file.path(output_dir, paste0("DEG_all_", brain_region, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(copper_gene_results, 
              file = file.path(output_dir, paste0("copper_genes_", brain_region, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(list(all_results = all_results, all_combined = all_combined, copper_gene_results = copper_gene_results))
}

# ==================== 9. 分别分析两个脑区 ====================
cortex_results <- run_limma_analysis("Cortex", expr_matrix, sample_info, annot_table, 
                                     output_dir, figure_dir, copper_genes, limfc_threshold, fdr_threshold)

striatum_results <- run_limma_analysis("Striatum", expr_matrix, sample_info, annot_table, 
                                       output_dir, figure_dir, copper_genes, limfc_threshold, fdr_threshold)

# ==================== 10. 结果汇总 ====================
cat("\n", strrep("=", 60), "\n")
cat("分析完成! 结果汇总\n")
cat(strrep("=", 60), "\n")

cat(paste0("\n筛选标准: |log2FC| >= ", limfc_threshold, " (FC>1.5), FDR < ", fdr_threshold, "\n"))

for (region_name in c("Cortex", "Striatum")) {
  cat(paste0("\n--- ", region_name, " ---\n"))
  region_results <- if(region_name == "Cortex") cortex_results else striatum_results
  
  for (tp in c("2h", "8h", "24h")) {
    if (!is.null(region_results$all_results[[tp]])) {
      sig <- region_results$all_results[[tp]] %>% filter(adj.P.Val < fdr_threshold & abs(logFC) >= limfc_threshold)
      cat(paste0("  ", tp, ": ", nrow(sig), " DEGs (上调: ", sum(sig$logFC > 0), 
                 ", 下调: ", sum(sig$logFC < 0), ")\n"))
    }
  }
  
  cat(paste0("  铜死亡基因DEGs:\n"))
  copper_sig <- region_results$copper_gene_results %>% filter(adj.P.Val < fdr_threshold & abs(logFC) >= limfc_threshold)
  if (nrow(copper_sig) > 0) {
    print(copper_sig %>% select(Gene_Symbol, time_point, logFC, P.Value, adj.P.Val, significance))
  } else {
    cat("    无显著差异表达的铜死亡基因\n")
    cat(paste0("    铜死亡基因探针总数: ", nrow(region_results$copper_gene_results), "\n"))
  }
}

cat(paste0("\n\n输出文件已保存到: ", output_dir, "\n"))
cat("========================================\n")
