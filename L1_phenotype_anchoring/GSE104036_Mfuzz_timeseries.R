# ==================== L1 表型锚定：GSE104036 多时序 Mfuzz 趋势聚类 ====================
# 数据来源: GSE104036 - MCAO小鼠脑组织bulk RNA-seq 多时序 (0/3/6/12/24h)
# 方法: Mfuzz 模糊C均值聚类 (Shen-Orr et al., Nature/Cell同款方法)
# 输出: 损伤侧 & 对侧 基因表达时序趋势图谱

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(42)

# ==================== 0. 包管理 ====================
required_packages <- c("Mfuzz", "data.table", "RColorBrewer", "ggplot2", "pheatmap", "clusterProfiler", "org.Mm.eg.db", "enrichplot")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ==================== 1. 配置参数 ====================
OUTPUT_DIR <- file.path(getwd(), "results", "L1_phenotype_anchoring", "GSE104036_Mfuzz")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

DATA_DIR <- "D:/反向网络药理学/L1 数据集/bulk/GSE104036（多时序）"

CLUSTER_NUM <- 8          # Mfuzz聚类数（参考Cell文章Fig S5用8类；也可设为6类）
EXPRESSION_FILTER <- 1    # 基因在所有时间点的平均CPM至少>1
FC_FILTER <- 1.5          # 至少一个时间点 vs Sham fold change > 1.5（筛选有变化的基因）

# 铜死亡核心基因（人）
CUPROPTOSIS_HUMAN <- c('FDX1', 'LIAS', 'LIPT1', 'DLD', 'DLAT', 'PDHA1', 'PDHB', 'MTF1', 'GLS', 'CDKN2A',
                       'SLC31A1', 'ATP7B', 'ATP7A', 'ATOX1', 'CCS', 'COX17', 'COMMD1', 'MT2A', 'SOD1', 'GCSH')

# 鼠对应同源基因（大写一致，但有些鼠名可能不同，Gene symbol基本一致）
CUPROPTOSIS_MOUSE <- CUPROPTOSIS_HUMAN

cat("============================================================\n")
cat("GSE104036 Mfuzz 时序趋势聚类分析\n")
cat("============================================================\n\n")

# ==================== 2. 数据加载 ====================
cat("=== 步骤2: 加载GSE104036计数数据 ===\n")

exp_raw <- fread(file.path(DATA_DIR, "GSE104036_TC-RNAseq_counts.txt.gz"), 
                 data.table = FALSE, header = TRUE)

gene_names <- exp_raw[, 1]
exp_raw <- exp_raw[, -1]

# 处理重复基因名：合并为平均值
dup_genes <- gene_names[duplicated(gene_names)]
if (length(dup_genes) > 0) {
  cat(sprintf("发现 %d 个重复基因名（如 %s），合并为平均值\n", 
              length(dup_genes), paste(head(unique(dup_genes), 3), collapse = ", ")))
  exp_raw$gene_symbol <- gene_names
  exp_agg <- aggregate(. ~ gene_symbol, data = exp_raw, FUN = mean)
  gene_names <- exp_agg$gene_symbol
  exp_raw <- exp_agg[, -1, drop = FALSE]
} else {
  exp_raw$gene_symbol <- gene_names
}

rownames(exp_raw) <- gene_names

cat(sprintf("原始数据: %d 基因 × %d 样本\n", nrow(exp_raw), ncol(exp_raw)))
cat(sprintf("列名: %s\n", paste(colnames(exp_raw)[1:min(5, ncol(exp_raw))], collapse = ", ")))

# ==================== 3. 提取分组数据 ====================
# 损伤侧 (Ipsilateral): Sham(0h) → I_3h → I_6h → I_12h → I_24h
# 对侧 (Contralateral): Sham(0h) → C_3h → C_6h → C_12h → C_24h

extract_time_series <- function(exp_matrix, side = c("ipsilateral", "contralateral")) {
  side <- match.arg(side)
  
  sham_cols <- c("S1", "S2", "S3")
  
  if (side == "ipsilateral") {
    time_cols <- list(
      "I_3h" = c("I1_3hr", "I2_3hr", "I3_3hr"),
      "I_6h" = c("I1_6hr", "I2_6hr", "I3_6hr"),
      "I_12h" = c("I1_12hr", "I2_12hr", "I3_12hr"),
      "I_24h" = c("I1_24hr", "I2_24hr", "I3_24hr")
    )
  } else {
    time_cols <- list(
      "C_3h" = c("C1_3hr", "C2_3hr", "C3_3hr"),
      "C_6h" = c("C1_6hr", "C2_6hr", "C3_6hr"),
      "C_12h" = c("C1_12hr", "C2_12hr", "C3_12hr"),
      "C_24h" = c("C1_24hr", "C2_24hr", "C3_24hr")
    )
  }
  
  # 提取Sham
  sham_data <- exp_matrix[, sham_cols, drop = FALSE]
  
  # 提取各时间点
  time_data_list <- lapply(time_cols, function(cols) {
    exp_matrix[, cols, drop = FALSE]
  })
  
  list(
    sham = sham_data,
    timepoints = time_data_list,
    time_labels = c("0h", "3h", "6h", "12h", "24h"),
    all_cols = c(sham_cols, unlist(time_cols)),
    data = exp_matrix[, c(sham_cols, unlist(time_cols)), drop = FALSE]
  )
}

ipsi_data <- extract_time_series(exp_raw, "ipsilateral")
contra_data <- extract_time_series(exp_raw, "contralateral")

# ==================== 4. CPM标准化 ====================
cat("=== 步骤4: CPM标准化 ===\n")

cpm_normalize <- function(count_matrix) {
  total_counts <- colSums(count_matrix)
  cpm <- t(t(count_matrix) / total_counts * 1e6)
  return(cpm)
}

ipsi_cpm <- cpm_normalize(ipsi_data$data)
contra_cpm <- cpm_normalize(contra_data$data)

# ==================== 5. 过滤低表达基因 ====================
cat("=== 步骤5: 过滤低表达基因 ===\n")

filter_genes <- function(cpm_matrix, group_labels, min_cpm = EXPRESSION_FILTER) {
  time_labels <- c("0h", "3h", "6h", "12h", "24h")
  n_per_group <- 3
  
  # 计算每个时间点的平均CPM
  avg_by_time <- matrix(NA, nrow = nrow(cpm_matrix), ncol = length(time_labels))
  colnames(avg_by_time) <- time_labels
  rownames(avg_by_time) <- rownames(cpm_matrix)
  
  for (i in seq_along(time_labels)) {
    cols <- ((i - 1) * n_per_group + 1):(i * n_per_group)
    avg_by_time[, i] <- rowMeans(cpm_matrix[, cols, drop = FALSE])
  }
  
  # 至少一个时间点平均CPM >= min_cpm
  keep <- rowSums(avg_by_time >= min_cpm) >= 1
  
  cat(sprintf("过滤前: %d 基因, 过滤后: %d 基因 (至少一个时间点CPM>=%.0f)\n", 
              nrow(cpm_matrix), sum(keep), min_cpm))
  
  return(list(
    filtered = cpm_matrix[keep, ],
    avg_by_time = avg_by_time[keep, ],
    keep_idx = keep
  ))
}

ipsi_filt <- filter_genes(ipsi_cpm, ipsi_data$time_labels)
contra_filt <- filter_genes(contra_cpm, contra_data$time_labels)

# ==================== 6. 筛选有显著时序变化的基因 ====================
cat("=== 步骤6: 筛选有显著时序变化的基因（Fold Change过滤）===\n")

filter_by_foldchange <- function(avg_by_time, fc_threshold = FC_FILTER) {
  sham_avg <- avg_by_time[, "0h"]
  max_fc <- apply(avg_by_time, 1, function(x) {
    max(abs(x[-1] / (sham_avg[match(names(x)[1], rownames(avg_by_time))] + 0.01)), na.rm = TRUE)
  })
  
  # 重新计算max FC
  time_cols <- c("3h", "6h", "12h", "24h")
  max_fc <- rep(1, nrow(avg_by_time))
  for (i in seq_len(nrow(avg_by_time))) {
    sham_val <- avg_by_time[i, "0h"] + 0.01
    fcs <- avg_by_time[i, time_cols] / sham_val
    max_fc[i] <- max(abs(log2(fcs)), na.rm = TRUE)
  }
  
  keep <- max_fc >= log2(fc_threshold)
  
  cat(sprintf("FC过滤后保留: %d / %d 基因 (|log2FC|>=%.2f)\n", 
              sum(keep), nrow(avg_by_time), log2(fc_threshold)))
  
  return(list(
    avg_by_time = avg_by_time[keep, ],
    keep_idx = which(keep)
  ))
}

ipsi_fc <- filter_by_foldchange(ipsi_filt$avg_by_time)
contra_fc <- filter_by_foldchange(contra_filt$avg_by_time)

# ==================== 7. Mfuzz 模糊C均值聚类 ====================
cat("\n=== 步骤7: Mfuzz 模糊C均值聚类 ===\n")

run_mfuzz <- function(avg_matrix, cluster_num = CLUSTER_NUM, side_label = "") {
  
  # 构建 ExpressionSet
  mfuzz_eset <- new("ExpressionSet", exprs = as.matrix(avg_matrix))
  
  # 过滤缺失值
  mfuzz_eset <- filter.NA(mfuzz_eset, thres = 0.25)
  mfuzz_eset <- fill.NA(mfuzz_eset, mode = "mean")
  
  # 过滤无变化基因
  mfuzz_eset <- filter.std(mfuzz_eset, min.std = 0)
  
  # 标准化
  mfuzz_eset <- standardise(mfuzz_eset)
  
  cat(sprintf("\n[%s] Mfuzz分析基因数: %d\n", side_label, nrow(exprs(mfuzz_eset))))
  
  n_genes <- nrow(exprs(mfuzz_eset))
  actual_clusters <- min(cluster_num, max(2, floor(n_genes / 10)))
  
  cat(sprintf("[%s] 使用聚类数: %d\n", side_label, actual_clusters))
  
  # 计算模糊度m
  m_value <- mestimate(mfuzz_eset)
  cat(sprintf("[%s] 估计模糊度 m = %.3f\n", side_label, m_value))
  
  # 执行Mfuzz聚类
  set.seed(42)
  mfuzz_result <- mfuzz(mfuzz_eset, c = actual_clusters, m = m_value)
  
  # 提取聚类成员分数
  membership <- mfuzz_result$membership
  colnames(membership) <- paste0("Cluster_", 1:actual_clusters)
  cluster_assignment <- apply(membership, 1, which.max)
  max_score <- apply(membership, 1, max)
  
  result_df <- data.frame(
    Gene = rownames(avg_matrix),
    Cluster = cluster_assignment,
    Membership_Score = max_score,
    stringsAsFactors = FALSE
  )
  
  for (k in 1:actual_clusters) {
    result_df[[paste0("Membership_C", k)]] <- membership[, k]
  }
  
  list(
    eset = mfuzz_eset,
    mfuzz_cluster = mfuzz_result,
    membership = membership,
    result_df = result_df,
    cluster_num = actual_clusters,
    m_value = m_value
  )
}

ipsi_mfuzz <- run_mfuzz(ipsi_fc$avg_by_time, CLUSTER_NUM, "损伤侧")
contra_mfuzz <- run_mfuzz(contra_fc$avg_by_time, CLUSTER_NUM, "对侧")

# ==================== 8. 可视化 ====================
cat("\n=== 步骤8: 时序趋势可视化 ===\n")

# 配色方案 - 模仿Cell/Nature风格
color_palette <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)

plot_mfuzz <- function(mfuzz_obj, side_label, cluster_num, time_labels) {
  
  ncols <- min(3, cluster_num)
  nrows <- ceiling(cluster_num / ncols)
  
  # 主图 - mfuzz.plot2 (使用默认配色)
  pdf(file.path(OUTPUT_DIR, sprintf("Mfuzz_Clusters_%s.pdf", side_label)), 
      width = ncols * 4, height = nrows * 3.5)
  
  mfuzz.plot2(mfuzz_obj$eset, 
              cl = mfuzz_obj$mfuzz_cluster, 
              mfrow = c(nrows, ncols),
              time.labels = time_labels,
              xlab = "Time after MCAO",
              ylab = "Standardized Expression")
  
  title(main = sprintf("Mfuzz Time-Series Clusters - %s (GSE104036)", side_label), 
        outer = TRUE, line = -1.5, cex.main = 1.2)
  
  dev.off()
  cat(sprintf("  → 保存: Mfuzz_Clusters_%s.pdf\n", side_label))
  
  # 每个聚类单独绘制（大图）
  pdf(file.path(OUTPUT_DIR, sprintf("Mfuzz_Individual_%s.pdf", side_label)), 
      width = 8, height = 6)
  
  for (k in 1:cluster_num) {
    genes_in_cluster <- rownames(mfuzz_obj$membership)[max.col(mfuzz_obj$membership) == k]
    n_genes <- length(genes_in_cluster)
    
    if (n_genes == 0) next
    
    # 提取该聚类的表达数据
    cluster_expr <- exprs(mfuzz_obj$eset)[genes_in_cluster, , drop = FALSE]
    
    # 绘制
    par(mar = c(4, 4, 3, 1))
    if (n_genes > 1) {
      matplot(t(cluster_expr), type = "l", lty = 1, col = adjustcolor("grey60", 0.3),
              xaxt = "n", xlab = "Time after MCAO", ylab = "Standardized Expression",
              main = sprintf("Cluster %d (%d genes)", k, n_genes))
      
      # 添加均值曲线
      mean_trend <- colMeans(cluster_expr)
      lines(mean_trend, col = "#E41A1C", lwd = 2.5)
    } else {
      plot(cluster_expr[1, ], type = "o", col = "#E41A1C", lwd = 2,
           xaxt = "n", xlab = "Time after MCAO", ylab = "Standardized Expression",
           main = sprintf("Cluster %d (%d gene)", k, n_genes))
    }
    axis(1, at = 1:5, labels = time_labels)
    grid(col = "grey90")
  }
  
  dev.off()
  cat(sprintf("  → 保存: Mfuzz_Individual_%s.pdf\n", side_label))
}

plot_mfuzz(ipsi_mfuzz, "Ipsilateral", ipsi_mfuzz$cluster_num, ipsi_data$time_labels)
plot_mfuzz(contra_mfuzz, "Contralateral", contra_mfuzz$cluster_num, contra_data$time_labels)

# ==================== 9. 铜死亡基因趋势提取 ====================
cat("\n=== 步骤9: 铜死亡基因时序趋势提取 ===\n")

extract_cuproptosis_trend <- function(cpm_matrix, avg_matrix, cuproptosis_genes, time_labels, side_label) {
  
  available_genes <- intersect(cuproptosis_genes, rownames(avg_matrix))
  cat(sprintf("[%s] 铜死亡基因可用: %d/%d\n", side_label, length(available_genes), length(cuproptosis_genes)))
  cat(sprintf("  可用基因: %s\n", paste(available_genes, collapse = ", ")))
  
  if (length(available_genes) == 0) return(NULL)
  
  trend_data <- avg_matrix[available_genes, , drop = FALSE]
  
  # 绘制铜死亡基因趋势
  pdf(file.path(OUTPUT_DIR, sprintf("Cuproptosis_Trend_%s.pdf", side_label)), 
      width = 10, height = 8)
  
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 8))
  
  colors <- RColorBrewer::brewer.pal(min(8, length(available_genes)), "Set1")
  if (length(available_genes) > 8) {
    colors <- colorRampPalette(colors)(length(available_genes))
  }
  
  plot(1:5, trend_data[1, ], type = "n", 
       xaxt = "n", xlim = c(1, 5), 
       ylim = range(trend_data, na.rm = TRUE) * 1.1,
       xlab = "Time after MCAO", ylab = "Average CPM",
       main = sprintf("Cuproptosis Gene Expression Trends - %s", side_label),
       cex.main = 1.1, cex.lab = 1.0)
  
  for (i in seq_len(nrow(trend_data))) {
    lines(1:5, trend_data[i, ], col = colors[i], lwd = 2.5, type = "o", pch = 19, cex = 0.8)
  }
  
  axis(1, at = 1:5, labels = time_labels)
  grid(col = "grey90")
  
  legend("topright", inset = c(-0.35, 0), 
         legend = rownames(trend_data),
         col = colors, lwd = 2, pch = 19, cex = 0.7,
         xpd = TRUE, bty = "n")
  
  dev.off()
  cat(sprintf("  → 保存: Cuproptosis_Trend_%s.pdf\n", side_label))
  
  # 保存趋势数据
  trend_df <- data.frame(
    Gene = rownames(trend_data),
    trend_data,
    stringsAsFactors = FALSE
  )
  write.csv(trend_df, file.path(OUTPUT_DIR, sprintf("Cuproptosis_Trend_Data_%s.csv", side_label)), 
            row.names = FALSE)
  
  return(trend_data)
}

ipsi_cupro <- extract_cuproptosis_trend(ipsi_cpm, ipsi_fc$avg_by_time, 
                                         CUPROPTOSIS_MOUSE, ipsi_data$time_labels, "Ipsilateral")
contra_cupro <- extract_cuproptosis_trend(contra_cpm, contra_fc$avg_by_time, 
                                           CUPROPTOSIS_MOUSE, contra_data$time_labels, "Contralateral")

# ==================== 10. 热图可视化每个聚类的代表性趋势 ====================
cat("\n=== 步骤10: 聚类热图 ===\n")

plot_cluster_heatmap <- function(mfuzz_obj, avg_matrix, side_label) {
  
  cluster_labels <- mfuzz_obj$result_df$Cluster
  names(cluster_labels) <- mfuzz_obj$result_df$Gene
  
  # 每个聚类取membership score最高的top基因
  cluster_order <- order(cluster_labels)
  sorted_avg <- avg_matrix[names(sort(cluster_labels)), , drop = FALSE]
  
  # 采样以防基因太多
  max_genes_per_cluster <- 50
  sampled_genes <- character(0)
  annotation_row <- character(0)
  
  for (k in sort(unique(cluster_labels))) {
    genes_k <- names(cluster_labels)[cluster_labels == k]
    if (length(genes_k) > max_genes_per_cluster) {
      genes_k <- sample(genes_k, max_genes_per_cluster)
    }
    sampled_genes <- c(sampled_genes, genes_k)
    annotation_row <- c(annotation_row, rep(paste0("C", k), length(genes_k)))
  }
  
  heat_data <- avg_matrix[sampled_genes, , drop = FALSE]
  
  # z-score标准化用于可视化
  heat_data_z <- t(scale(t(heat_data)))
  heat_data_z[heat_data_z > 2] <- 2
  heat_data_z[heat_data_z < -2] <- -2
  
  annot_row <- data.frame(Cluster = annotation_row, row.names = sampled_genes)
  
  cluster_colors <- RColorBrewer::brewer.pal(min(8, mfuzz_obj$cluster_num), "Set2")
  names(cluster_colors) <- paste0("C", 1:mfuzz_obj$cluster_num)
  annot_colors <- list(Cluster = cluster_colors)
  
  pdf(file.path(OUTPUT_DIR, sprintf("Cluster_Heatmap_%s.pdf", side_label)), 
      width = 8, height = 12)
  
  pheatmap(heat_data_z, 
           cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = FALSE,
           annotation_row = annot_row,
           annotation_colors = annot_colors,
           color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
           main = sprintf("Mfuzz Clusters Heatmap - %s (GSE104036)", side_label),
           fontsize = 10,
           angle_col = 45)
  
  dev.off()
  cat(sprintf("  → 保存: Cluster_Heatmap_%s.pdf\n", side_label))
}

plot_cluster_heatmap(ipsi_mfuzz, ipsi_fc$avg_by_time, "Ipsilateral")
plot_cluster_heatmap(contra_mfuzz, contra_fc$avg_by_time, "Contralateral")

# ==================== 11. GO/KEGG富集分析（各聚类功能注释） ====================
cat("\n=== 步骤11: 各聚类GO/KEGG功能注释 ===\n")

run_cluster_enrichment <- function(mfuzz_obj, side_label, top_n_pathways = 5) {
  
  cluster_labels <- mfuzz_obj$result_df$Cluster
  names(cluster_labels) <- mfuzz_obj$result_df$Gene
  
  all_enrich_results <- list()
  
  for (k in sort(unique(cluster_labels))) {
    genes_k <- names(cluster_labels)[cluster_labels == k]
    n_genes <- length(genes_k)
    
    if (n_genes < 5) {
      cat(sprintf("[%s] Cluster %d: 基因数过少(%d), 跳过富集分析\n", side_label, k, n_genes))
      next
    }
    
    cat(sprintf("[%s] Cluster %d: %d 基因 → GO/KEGG富集\n", side_label, k, n_genes))
    
    # GO BP富集
    tryCatch({
      ego <- enrichGO(gene = genes_k,
                      OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
      
      if (!is.null(ego) && nrow(ego@result) > 0) {
        all_enrich_results[[paste0("Cluster_", k)]] <- ego@result
        
        # 点图
        if (nrow(ego@result) >= 3) {
          pdf(file.path(OUTPUT_DIR, sprintf("GO_Cluster%d_%s.pdf", k, side_label)), 
              width = 10, height = 6)
          print(dotplot(ego, showCategory = min(top_n_pathways, nrow(ego@result)), 
                        title = sprintf("GO BP - Cluster %d (%d genes) - %s", k, n_genes, side_label)))
          dev.off()
        }
      }
    }, error = function(e) {
      cat(sprintf("[%s] Cluster %d GO富集失败: %s\n", side_label, k, e$message))
    })
    
    # KEGG富集
    tryCatch({
      entrez_ids <- bitr(genes_k, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
      
      if (nrow(entrez_ids) > 5) {
        kk <- enrichKEGG(gene = entrez_ids$ENTREZID,
                         organism = "mmu",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2)
        
        if (!is.null(kk) && nrow(kk@result) > 0) {
          all_enrich_results[[paste0("Cluster_", k, "_KEGG")]] <- kk@result
          
          if (nrow(kk@result) >= 3) {
            pdf(file.path(OUTPUT_DIR, sprintf("KEGG_Cluster%d_%s.pdf", k, side_label)), 
                width = 10, height = 6)
            print(dotplot(kk, showCategory = min(top_n_pathways, nrow(kk@result)),
                          title = sprintf("KEGG - Cluster %d (%d genes) - %s", k, n_genes, side_label)))
            dev.off()
          }
        }
      }
    }, error = function(e) {
      cat(sprintf("[%s] Cluster %d KEGG富集失败: %s\n", side_label, k, e$message))
    })
  }
  
  # 保存所有富集结果
  if (length(all_enrich_results) > 0) {
    saveRDS(all_enrich_results, file.path(OUTPUT_DIR, sprintf("Enrichment_Results_%s.rds", side_label)))
    
    for (name in names(all_enrich_results)) {
      write.csv(all_enrich_results[[name]], 
                file.path(OUTPUT_DIR, sprintf("Enrichment_%s_%s.csv", name, side_label)),
                row.names = FALSE)
    }
  }
  
  return(all_enrich_results)
}

ipsi_enrich <- run_cluster_enrichment(ipsi_mfuzz, "Ipsilateral")
contra_enrich <- run_cluster_enrichment(contra_mfuzz, "Contralateral")

# ==================== 12. 输出汇总结果 ====================
cat("\n=== 步骤12: 输出汇总结果 ===\n")

# 聚类统计
cluster_summary <- function(mfuzz_obj, side_label) {
  tab <- table(mfuzz_obj$result_df$Cluster)
  cat(sprintf("\n[%s] 聚类分布:\n", side_label))
  for (k in names(tab)) {
    cat(sprintf("  Cluster %s: %d genes\n", k, tab[k]))
  }
  
  df <- data.frame(
    Cluster = names(tab),
    Gene_Count = as.integer(tab),
    stringsAsFactors = FALSE
  )
  write.csv(df, file.path(OUTPUT_DIR, sprintf("Cluster_Summary_%s.csv", side_label)), 
            row.names = FALSE)
}

cluster_summary(ipsi_mfuzz, "Ipsilateral")
cluster_summary(contra_mfuzz, "Contralateral")

# 完整聚类成员表
write.csv(ipsi_mfuzz$result_df, 
          file.path(OUTPUT_DIR, "Mfuzz_Cluster_Membership_Ipsilateral.csv"),
          row.names = FALSE)
write.csv(contra_mfuzz$result_df, 
          file.path(OUTPUT_DIR, "Mfuzz_Cluster_Membership_Contralateral.csv"),
          row.names = FALSE)

# 各时间点平均表达数据
write.csv(as.data.frame(ipsi_fc$avg_by_time), 
          file.path(OUTPUT_DIR, "AvgExpression_By_Time_Ipsilateral.csv"),
          row.names = TRUE)
write.csv(as.data.frame(contra_fc$avg_by_time), 
          file.path(OUTPUT_DIR, "AvgExpression_By_Time_Contralateral.csv"),
          row.names = TRUE)

# ==================== 13. 损伤侧vs对侧趋势一致性分析 ====================
cat("\n=== 步骤13: 损伤侧vs对侧趋势对比 ===\n")

compare_trends <- function(ipsi_avg, contra_avg) {
  common_genes <- intersect(rownames(ipsi_avg), rownames(contra_avg))
  cat(sprintf("损伤侧与对侧共有显著变化基因: %d\n", length(common_genes)))
  
  if (length(common_genes) < 10) return(NULL)
  
  # 计算Spearman相关性
  cors <- numeric(length(common_genes))
  names(cors) <- common_genes
  
  for (gene in common_genes) {
    ct <- cor(as.numeric(ipsi_avg[gene, ]), as.numeric(contra_avg[gene, ]), 
              method = "spearman")
    cors[gene] <- ct
  }
  
  pdf(file.path(OUTPUT_DIR, "Ipsi_vs_Contra_Trend_Correlation.pdf"), 
      width = 8, height = 6)
  
  hist(cors, breaks = 50, col = "steelblue", border = "white",
       main = "Ipsilateral vs Contralateral Trend Correlation",
       xlab = "Spearman Correlation", ylab = "Gene Count")
  abline(v = 0, lty = 2, col = "red", lwd = 2)
  abline(v = median(cors, na.rm = TRUE), lty = 2, col = "darkgreen", lwd = 2)
  legend("topleft", 
         legend = c(sprintf("Median r = %.3f", median(cors, na.rm = TRUE)),
                    sprintf("Consistent (r>0.7): %d genes", sum(cors > 0.7, na.rm = TRUE)),
                    sprintf("Opposite (r< -0.3): %d genes", sum(cors < -0.3, na.rm = TRUE))),
         bty = "n", cex = 0.8)
  
  dev.off()
  
  # 保存相关性数据
  cor_df <- data.frame(Gene = common_genes, Spearman_R = cors, stringsAsFactors = FALSE)
  write.csv(cor_df, file.path(OUTPUT_DIR, "Ipsi_vs_Contra_Correlation.csv"), row.names = FALSE)
  
  cat(sprintf("趋势一致基因 (r>0.7): %d\n", sum(cors > 0.7, na.rm = TRUE)))
  cat(sprintf("趋势相反基因 (r< -0.3): %d\n", sum(cors < -0.3, na.rm = TRUE)))
  
  return(list(cors = cors, common_genes = common_genes))
}

compare_result <- compare_trends(ipsi_fc$avg_by_time, contra_fc$avg_by_time)

# ==================== 14. 损伤侧vs对侧聚类交叉表 ====================
cat("\n=== 步骤14: 损伤侧vs对侧交叉验证 ===\n")

cross_side_analysis <- function(ipsi_mfuzz_obj, contra_mfuzz_obj) {
  
  ipsi_cluster <- ipsi_mfuzz_obj$result_df
  contra_cluster <- contra_mfuzz_obj$result_df
  
  common_genes <- intersect(ipsi_cluster$Gene, contra_cluster$Gene)
  cat(sprintf("两侧共有聚类基因: %d\n", length(common_genes)))
  
  if (length(common_genes) < 10) return(NULL)
  
  ipsi_sub <- ipsi_cluster[ipsi_cluster$Gene %in% common_genes, ]
  contra_sub <- contra_cluster[contra_cluster$Gene %in% common_genes, ]
  
  ipsi_sub <- ipsi_sub[match(common_genes, ipsi_sub$Gene), ]
  contra_sub <- contra_sub[match(common_genes, contra_sub$Gene), ]
  
  cross_tab <- table(Ipsilateral = ipsi_sub$Cluster, 
                     Contralateral = contra_sub$Cluster)
  
  write.csv(as.data.frame.matrix(cross_tab), 
            file.path(OUTPUT_DIR, "CrossTab_Ipsi_vs_Contra_Cluster.csv"))
  
  cat("交叉聚类表:\n")
  print(cross_tab)
  
  return(list(cross_tab = cross_tab, ipsi_sub = ipsi_sub, contra_sub = contra_sub))
}

cross_result <- cross_side_analysis(ipsi_mfuzz, contra_mfuzz)

# ==================== 15. 总结 ====================
cat("\n============================================================\n")
cat("GSE104036 Mfuzz 时序趋势聚类分析 完成！\n")
cat("============================================================\n")
cat(sprintf("输出目录: %s\n", OUTPUT_DIR))
cat("\n生成文件列表:\n")
output_files <- list.files(OUTPUT_DIR, full.names = FALSE)
for (f in output_files) {
  cat(sprintf("  • %s\n", f))
}
cat("\n============================================================\n")

sessionInfo()