# ============================================================================
# 第3步：虚拟敲除干预验证 - 快速版（Top 1000细胞）
# ============================================================================

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

cat("=== 第3步：虚拟敲除干预验证（快速版 - Top 2000细胞）===\n\n")

result_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/causal_analysis_results"
knockout_dir <- file.path(result_dir, "virtual_knockout_results")
dir.create(knockout_dir, showWarnings=FALSE, recursive=TRUE)

# ============================================================================
# 1. 读取单细胞数据
# ============================================================================
cat("【3.1】读取单细胞数据...\n")

sc_file <- file.path(result_dir, "SCISSOR_results", "sc_final_with_scores.rds")
sc_obj <- readRDS(sc_file)
cat(sprintf("  原始细胞数: %d\n", ncol(sc_obj)))

# ============================================================================
# 2. 选择PC Hub基因表达最显著的Top 1000细胞
# ============================================================================
cat("\n【3.2】选择PC Hub基因表达最显著的Top 1000细胞...\n")

# 使用PC_Hub_Score排序
hub_scores <- sc_obj$PC_Hub_Score
names(hub_scores) <- colnames(sc_obj)

# 选择top 2000细胞
top_cells <- names(sort(hub_scores, decreasing=TRUE))[1:2000]
sc_obj_top <- subset(sc_obj, cells=top_cells)

cat(sprintf("  选择Top 2000细胞（PC_Hub_Score范围: %.3f - %.3f）\n",
            min(hub_scores[top_cells]), max(hub_scores[top_cells])))

# 检查分组分布
cat("  Top 2000细胞分组分布:\n")
print(table(sc_obj_top$group))

# ============================================================================
# 3. 提取count矩阵
# ============================================================================
cat("\n【3.3】准备count矩阵...\n")

count_matrix <- as.matrix(GetAssayData(sc_obj_top, layer="counts"))

# 过滤低表达基因
keep_genes <- rowSums(count_matrix > 0) >= 10
count_matrix <- count_matrix[keep_genes, ]

cat(sprintf("  Count矩阵: %d 基因 x 2000 细胞\n", nrow(count_matrix), ncol(count_matrix)))

# ============================================================================
# 4. 确定敲除目标基因
# ============================================================================
cat("\n【3.4】确定敲除目标基因...\n")

bcp_genes <- c("Nfkb1", "Stat3", "Fdx1", "Tlr4")
knockout_targets <- bcp_genes[bcp_genes %in% rownames(count_matrix)]

if(length(knockout_targets) < length(bcp_genes)) {
  missing <- setdiff(bcp_genes, knockout_targets)
  cat(sprintf("  [警告] 以下基因不在数据中: %s\n", paste(missing, collapse=", ")))
}

cat(sprintf("  敲除目标基因 (%d 个): %s\n", 
            length(knockout_targets), paste(knockout_targets, collapse=", ")))

# ============================================================================
# 5. 快速虚拟敲除（基于相关性）
# ============================================================================
cat("\n【3.5】执行快速虚拟敲除（基于相关性）...\n")

knockout_results <- list()

for(target_gene in knockout_targets) {
  cat(sprintf("\n  虚拟敲除 %s...\n", target_gene))
  
  # 获取目标基因表达
  target_expr <- as.numeric(count_matrix[target_gene, ])
  
  # 计算与所有基因的相关性
  cat("    计算基因相关性...\n")
  all_genes <- setdiff(rownames(count_matrix), target_gene)
  
  # 分批计算避免内存问题
  batch_size <- 500
  n_genes <- length(all_genes)
  cor_results <- data.frame(
    gene = character(),
    correlation = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(i in seq(1, n_genes, by=batch_size)) {
    end_idx <- min(i + batch_size - 1, n_genes)
    batch_genes <- all_genes[i:end_idx]
    
    batch_cor <- sapply(batch_genes, function(g) {
      g_expr <- as.numeric(count_matrix[g, ])
      suppressWarnings(cor(target_expr, g_expr, method = "spearman"))
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
  cor_results$predicted_logFC <- -cor_results$correlation * 2
  cor_results$p.value <- 2 * pnorm(-abs(cor_results$correlation) * sqrt(ncol(count_matrix)))
  cor_results$padj <- p.adjust(cor_results$p.value, method = "BH")
  
  # 筛选显著差异基因
  degs <- cor_results %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
  cat(sprintf("    ✅ 差异基因数: %d\n", nrow(degs)))
  
  # 显示top10差异基因
  cat("    Top 10差异基因:\n")
  for(i in 1:min(10, nrow(degs))) {
    direction <- ifelse(degs$predicted_logFC[i] > 0, "↑", "↓")
    cat(sprintf("      %2d. %-15s %s%.3f (p=%.2e)\n", 
                i, degs$gene[i], direction, abs(degs$predicted_logFC[i]), degs$padj[i]))
  }
  
  # 保存结果
  write.csv(degs, file.path(knockout_dir, paste0("knockout_", target_gene, "_DEGs_fast.csv")), row.names=FALSE)
  
  knockout_results[[target_gene]] <- list(
    degs = degs,
    n_degs = nrow(degs),
    cor_results = cor_results
  )
}

# ============================================================================
# 6. 汇总分析
# ============================================================================
cat("\n【3.6】敲除结果汇总...\n")

cat("\n  各基因敲除后的差异基因数:\n")
for(gene in names(knockout_results)) {
  cat(sprintf("    %s: %d 个差异基因\n", gene, knockout_results[[gene]]$n_degs))
}

# 寻找共同的下游基因
cat("\n  寻找共同的下游调控基因...\n")
all_degs <- lapply(knockout_results, function(x) x$degs$gene)

if(length(all_degs) >= 2) {
  # 两两比较找共同基因
  for(i in 1:(length(all_degs)-1)) {
    for(j in (i+1):length(all_degs)) {
      gene1 <- names(all_degs)[i]
      gene2 <- names(all_degs)[j]
      common <- intersect(all_degs[[i]], all_degs[[j]])
      cat(sprintf("    %s ∩ %s: %d 个共同基因\n", gene1, gene2, length(common)))
    }
  }
  
  # 所有基因的共同下游
  common_all <- Reduce(intersect, all_degs)
  cat(sprintf("\n    所有敲除共同下游基因: %d 个\n", length(common_all)))
  
  if(length(common_all) > 0) {
    cat("    共同基因列表:", paste(head(common_all, 10), collapse=", "), "...\n")
    write.csv(data.frame(gene=common_all), 
              file.path(knockout_dir, "common_downstream_genes_fast.csv"), row.names=FALSE)
  }
}

# ============================================================================
# 7. 可视化
# ============================================================================
cat("\n【3.7】生成可视化...\n")

pdf(file.path(knockout_dir, "knockout_summary_fast.pdf"), width=14, height=10)

# 差异基因数柱状图
n_degs_list <- sapply(knockout_results, function(x) x$n_degs)
par(mar=c(5,4,4,2))
barplot(n_degs_list, main="Number of DEGs after Virtual Knockout (Top 1000 cells)",
        ylab="Number of DEGs", xlab="Knockout Gene",
        col="#2171B5", las=2)

# 每个敲除基因的top差异基因柱状图
for(gene in names(knockout_results)) {
  degs <- knockout_results[[gene]]$degs
  
  if(nrow(degs) > 0) {
    par(mfrow=c(1,1), mar=c(5,10,4,2))
    
    top20 <- head(degs, 20)
    bar_colors <- ifelse(top20$predicted_logFC > 0, "#CB181D", "#2171B5")
    
    barplot(top20$predicted_logFC, names.arg=top20$gene, las=2, cex.names=0.6,
            main=paste("Top 20 DEGs after", gene, "Knockout"),
            xlab="Predicted log2 Fold Change", col=bar_colors, horiz=TRUE)
    abline(v=0, lty=2, col="black")
  }
}

dev.off()
cat("  ✅ 可视化已保存\n")

# 保存汇总结果
saveRDS(knockout_results, file.path(knockout_dir, "all_knockout_results_fast.rds"))

cat("\n"); cat(rep("=", 60), sep=""); cat("\n")
cat("第3步（快速版）完成！\n")
cat("结果目录:", knockout_dir, "\n")
cat("\n关键发现:\n")
for(gene in names(knockout_results)) {
  cat(sprintf("  - %s敲除: %d 个差异基因\n", gene, knockout_results[[gene]]$n_degs))
}
if(exists("common_all") && length(common_all) > 0) {
  cat(sprintf("  - 发现 %d 个共同下游调控基因\n", length(common_all)))
}
cat("  - 基于Top 2000 PC Hub基因高表达细胞\n")
cat(rep("=", 60), sep=""); cat("\n")
