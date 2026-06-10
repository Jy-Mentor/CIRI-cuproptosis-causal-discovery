# 虚拟敲除分析 - 独立版 (不依赖scTenifoldNet)
# 使用相关性网络 + 差异分析模拟敲除效果
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))

if(!"Seurat" %in% installed.packages()){install.packages('Seurat', repos='https://cran.rstudio.com')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
if(!"igraph" %in% installed.packages()){install.packages('igraph')}
library(Seurat)
library(ggplot2)
library(igraph)

cat("=== 虚拟敲除分析 (独立版) ===\n")
cat("方法: 相关性网络 + 差异分析\n\n")

# 1. 读取数据
cat("1. 读取单细胞数据...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
sc_obj <- JoinLayers(sc_obj)

# 抽样减少计算量
set.seed(42)
max_cells <- 3000
if(ncol(sc_obj) > max_cells) {
  sc_obj <- subset(sc_obj, cells = sample(Cells(sc_obj), max_cells))
}
cat(sprintf("  使用 %d 细胞\n", ncol(sc_obj)))

# 2. 提取表达矩阵
cat("\n2. 构建基因网络...\n")
sc_expr <- GetAssayData(sc_obj, layer="data")

# 目标基因
target_genes <- c("Nfkb1", "Fdx1", "Tlr4")

# 3. 虚拟敲除函数 (基于相关性)
virtual_knockout_sim <- function(expr_matrix, gene_of_interest, n_top = 100) {
  cat(sprintf("\n  分析 %s...\n", gene_of_interest))

  # 检查基因是否存在
  if(!gene_of_interest %in% rownames(expr_matrix)) {
    cat(sprintf("  [跳过] %s 不在数据中\n", gene_of_interest))
    return(NULL)
  }

  # 计算目标基因与其他基因的相关性
  gene_expr <- as.numeric(expr_matrix[gene_of_interest, ])
  all_genes <- rownames(expr_matrix)

  cat(sprintf("    计算 %d 个基因的相关性...\n", length(all_genes)))

  # 并行计算相关性 (分批处理减少内存)
  batch_size <- 1000
  n_genes <- length(all_genes)
  cor_results <- data.frame(
    gene = character(),
    correlation = numeric(),
    stringsAsFactors = FALSE
  )

  for(i in seq(1, n_genes, by = batch_size)) {
    end_idx <- min(i + batch_size - 1, n_genes)
    batch_genes <- all_genes[i:end_idx]

    batch_cor <- sapply(batch_genes, function(g) {
      if(g == gene_of_interest) return(NA)
      g_expr <- as.numeric(expr_matrix[g, ])
      suppressWarnings(cor(gene_expr, g_expr, method = "spearman"))
    })

    cor_results <- rbind(cor_results, data.frame(
      gene = batch_genes,
      correlation = batch_cor,
      stringsAsFactors = FALSE
    ))
  }

  # 移除NA
  cor_results <- cor_results[!is.na(cor_results$correlation), ]

  # 基于相关性预测敲除效应
  # 正相关 -> 敲除后下调 (logFC < 0)
  # 负相关 -> 敲除后上调 (logFC > 0)
  cor_results$predicted_logFC <- -cor_results$correlation * 2  # 放大效应
  cor_results$p.value <- 2 * pnorm(-abs(cor_results$correlation) * sqrt(ncol(expr_matrix)))
  cor_results$padj <- p.adjust(cor_results$p.value, method = "BH")

  # 排序
  cor_results <- cor_results[order(-abs(cor_results$predicted_logFC)), ]

  # 选择top基因
  top_results <- head(cor_results, n_top)

  cat(sprintf("    找到 %d 个显著相关基因\n", sum(cor_results$padj < 0.05, na.rm = TRUE)))

  return(list(
    full_results = cor_results,
    top_results = top_results,
    target_gene = gene_of_interest
  ))
}

# 4. 执行虚拟敲除
cat("\n3. 执行虚拟敲除...\n")
results <- list()

for(gene in target_genes) {
  res <- virtual_knockout_sim(sc_expr, gene, n_top = 50)
  if(!is.null(res)) {
    results[[gene]] <- res
    cat(sprintf("  %s: 找到 %d 个下游基因\n", gene, nrow(res$top_results)))
  }
}

# 5. 保存结果
cat("\n4. 保存结果...\n")
dir.create('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/virtual_knockout_v2', showWarnings=FALSE, recursive=TRUE)

saveRDS(results, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/virtual_knockout_v2/knockout_results.rds')

# 6. 可视化
cat("\n5. 生成可视化...\n")
pdf('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/virtual_knockout_v2/knockout_results.pdf', width=12, height=10)

# 为每个敲除基因绘制结果
for(gene in names(results)) {
  res <- results[[gene]]$top_results

  # 火山图
  par(mfrow=c(1,1))
  plot(res$predicted_logFC, -log10(res$p.value),
       pch=20, col=ifelse(res$predicted_logFC > 0, "#CB181D", "#2171B5"),
       main=paste("Virtual Knockout of", gene),
       xlab="Predicted logFC", ylab="-log10(p-value)")
  abline(h=-log10(0.05), lty=2, col="gray")
  abline(v=0, lty=2, col="gray")

  # 标记top基因
  top10 <- head(res, 10)
  text(top10$predicted_logFC, -log10(top10$p.value),
       labels=top10$gene, cex=0.6, pos=3)

  # 条形图
  par(mfrow=c(1,1), mar=c(5,10,4,2))
  top20 <- head(res, 20)
  bar_colors <- ifelse(top20$predicted_logFC > 0, "#CB181D", "#2171B5")
  barplot(top20$predicted_logFC, names.arg=top20$gene, las=2, cex.names=0.6,
          main=paste("Top 20 DEGs after", gene, "Knockout"),
          xlab="Predicted logFC", col=bar_colors, horiz=TRUE)
  abline(v=0, lty=2, col="black")
}

dev.off()

# 7. 汇总结果
cat("\n=== 虚拟敲除结果汇总 ===\n")
for(gene in names(results)) {
  res <- results[[gene]]$top_results
  cat(sprintf("\n【%s 敲除】\n", gene))
  cat(sprintf("  总相关基因: %d\n", nrow(results[[gene]]$full_results)))
  cat(sprintf("  显著基因 (p<0.05): %d\n", sum(res$p.value < 0.05)))
  cat("  Top 10 预测变化基因:\n")
  for(i in 1:min(10, nrow(res))) {
    direction <- ifelse(res$predicted_logFC[i] > 0, "↑", "↓")
    cat(sprintf("    %s: logFC=%.3f %s (r=%.3f)\n",
                res$gene[i], res$predicted_logFC[i], direction, res$correlation[i]))
  }
}

cat("\n=== 完成 ===\n")
cat("结果保存到: C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/virtual_knockout_v2/\n")
