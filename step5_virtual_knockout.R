# 虚拟敲除分析 - 简化版 (优化计算效率)
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))

if(!"Seurat" %in% installed.packages()){install.packages('Seurat', repos='https://cran.rstudio.com')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
library(Seurat)
library(ggplot2)

cat("=== 虚拟敲除分析 (优化版) ===\n")

# 1. 读取单细胞数据
cat("1. 读取单细胞数据...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')

# 随机抽样减少计算量
set.seed(42)
max_cells <- 2000
if(ncol(sc_obj) > max_cells) {
  sc_obj <- subset(sc_obj, cells = sample(Cells(sc_obj), max_cells))
  cat(sprintf("  随机抽样至: %d 细胞\n", max_cells))
}

sc_expr <- GetAssayData(sc_obj, layer="data")
cat(sprintf("  表达矩阵: %d 基因 x %d 细胞\n", nrow(sc_expr), ncol(sc_expr)))

# 2. PC网络推断的因果下游
pc_downstream <- list(
  Nfkb1 = c("Stat1", "Tgfbr1", "Pdha1"),
  Fdx1 = c("Atox1", "Slc31a1", "Sod1"),
  Tlr4 = c("Nfe2l2", "Atox1"),
  Stat1 = c("Atox1"),
  Nfe2l2 = c("Tlr4", "Tgfbr1", "Sod1")
)

# 3. 快速虚拟敲除函数 (只计算关键基因的相关性)
virtual_knockout_fast <- function(expr_matrix, gene_of_interest, downstream_genes) {
  gene_idx <- which(rownames(expr_matrix) == gene_of_interest)
  if(length(gene_idx) == 0) return(NULL)

  gene_expr <- as.numeric(expr_matrix[gene_idx, ])

  # 只计算下游基因的相关性
  cor_vals <- sapply(downstream_genes, function(g) {
    g_idx <- which(rownames(expr_matrix) == g)
    if(length(g_idx) == 0) return(NA)
    g_expr <- as.numeric(expr_matrix[g_idx, ])
    suppressWarnings(cor(gene_expr, g_expr, method = "spearman"))
  })

  return(cor_vals)
}

# 4. 对每个目标基因进行虚拟敲除
cat("\n2. 执行虚拟敲除...\n")
knockout_targets <- c("Nfkb1", "Fdx1", "Tlr4")
results <- list()

for(target in knockout_targets) {
  cat(sprintf("\n  虚拟敲除 %s:\n", target))
  downstream <- pc_downstream[[target]]
  cor_vals <- virtual_knockout_fast(sc_expr, target, downstream)
  results[[target]] <- cor_vals

  for(g in downstream) {
    if(!is.na(cor_vals[g])) {
      direction <- ifelse(cor_vals[g] > 0, "正相关(敲除后↓)", "负相关(敲除后↑)")
      cat(sprintf("    %s → %s: r=%.3f (%s)\n", target, g, cor_vals[g], direction))
    }
  }
}

# 5. 保存结果
cat("\n3. 保存结果...\n")
dir.create('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/virtual_knockout', showWarnings=FALSE, recursive=TRUE)
saveRDS(results, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/virtual_knockout/knockout_results.rds')

# 6. 可视化
cat("\n4. 生成可视化...\n")
pdf('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/virtual_knockout/knockout_results.pdf', width=12, height=6)

par(mfrow=c(1, 3), mar=c(5, 4, 3, 2))
colors <- c("#2171B5", "#CB181D")

for(target in knockout_targets) {
  cor_vals <- results[[target]]
  cor_vals <- cor_vals[!is.na(cor_vals)]

  if(length(cor_vals) > 0) {
    bar_colors <- ifelse(cor_vals > 0, colors[1], colors[2])
    barplot(cor_vals, main = sprintf("敲除 %s", target),
            horiz = TRUE, las = 1, col = bar_colors,
            xlim = c(-1, 1), xlab = "Spearman r")
    abline(v = 0, lty = 2, col = "gray")
  }
}

dev.off()

# 7. 总结
cat("\n=== 虚拟敲除结果总结 ===\n")
for(target in knockout_targets) {
  cat(sprintf("\n【%s 敲除】\n", target))
  downstream <- pc_downstream[[target]]
  cor_vals <- results[[target]]
  for(g in downstream) {
    if(!is.na(cor_vals[g])) {
      direction <- ifelse(cor_vals[g] > 0, "↓(正相关)", "↑(负相关)")
      cat(sprintf("  预测 %s: r=%.3f %s\n", g, cor_vals[g], direction))
    }
  }
}

cat("\n=== 完成 ===\n")
cat("结果保存到: C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/virtual_knockout/\n")