# ============================================================================
# 第3步：单细胞虚拟敲除干预验证
# 使用scTenifoldKnk模拟敲除关键基因
# ============================================================================

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

library(Seurat)
library(Matrix)
library(scTenifoldNet)
library(scTenifoldKnk)
library(dplyr)
library(ggplot2)

cat("=== 第3步：单细胞虚拟敲除干预验证 ===\n\n")

result_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/causal_analysis_results"
knockout_dir <- file.path(result_dir, "virtual_knockout_results")
dir.create(knockout_dir, showWarnings=FALSE, recursive=TRUE)

# ============================================================================
# 1. 读取单细胞数据
# ============================================================================
cat("【3.1】读取单细胞数据...\n")

sc_file <- file.path(result_dir, "SCISSOR_results", "sc_final_with_scores.rds")
sc_obj <- readRDS(sc_file)
cat(sprintf("  细胞数: %d\n", ncol(sc_obj)))

# ============================================================================
# 2. 提取count矩阵
# ============================================================================
cat("\n【3.2】准备count矩阵...\n")

count_matrix <- as.matrix(GetAssayData(sc_obj, layer="counts"))

# 过滤低表达基因
keep_genes <- rowSums(count_matrix > 0) >= 50
count_matrix <- count_matrix[keep_genes, ]

cat(sprintf("  Count矩阵: %d 基因 x %d 细胞\n", nrow(count_matrix), ncol(count_matrix)))

# ============================================================================
# 3. 确定敲除目标基因
# ============================================================================
cat("\n【3.3】确定敲除目标基因...\n")

# 用户指定的BCP轴核心基因
bcp_genes <- c("Nfkb1", "Stat3", "Fdx1", "Tlr4")

# 检查基因是否在count矩阵中
knockout_targets <- bcp_genes[bcp_genes %in% rownames(count_matrix)]

if(length(knockout_targets) < length(bcp_genes)) {
  missing <- setdiff(bcp_genes, knockout_targets)
  cat(sprintf("  [警告] 以下基因不在数据中: %s\n", paste(missing, collapse=", ")))
}

cat(sprintf("  敲除目标基因 (%d 个): %s\n", 
            length(knockout_targets), paste(knockout_targets, collapse=", ")))

# ============================================================================
# 4. 执行虚拟敲除
# ============================================================================
cat("\n【3.4】执行scTenifoldKnk虚拟敲除...\n")

knockout_results <- list()

for(target_gene in knockout_targets) {
  cat(sprintf("\n  虚拟敲除 %s...\n", target_gene))
  cat("    (计算中，请等待约10-20分钟)...\n")
  
  # 手动QC过滤
  lib_sizes <- colSums(count_matrix)
  valid_cells <- lib_sizes >= 1000
  cm <- count_matrix[, valid_cells]
  
  # 过滤低表达基因
  gene_exp <- rowMeans(cm > 0)
  cm <- cm[gene_exp >= 0.05, ]
  
  cat(sprintf("    质控后: %d 基因 x %d 细胞\n", nrow(cm), ncol(cm)))
  
  # 运行scTenifoldKnk
  set.seed(666)
  start_time <- Sys.time()
  
  knk_res <- scTenifoldKnk(
    countMatrix = cm,
    gKO = target_gene,
    qc = FALSE,  # 已手动QC
    nc_lambda = 0,
    nc_nNet = 10,
    nc_nCells = 500,
    nc_nComp = 3,
    nc_q = 0.9,
    td_K = 3,
    ma_nDim = 2,
    nCores = 2
  )
  
  end_time <- Sys.time()
  cat(sprintf("    耗时: %.1f 分钟\n", as.numeric(difftime(end_time, start_time, units="mins"))))
  
  # 处理结果
  degs <- knk_res$diffRegulation %>%
    filter(p.value < 0.05) %>%
    mutate(logFC = log2(FC)) %>%
    arrange(p.value)
  
  cat(sprintf("    ✅ 差异基因数: %d\n", nrow(degs)))
  
  # 显示top10差异基因
  cat("    Top 10差异基因:\n")
  for(i in 1:min(10, nrow(degs))) {
    direction <- ifelse(degs$logFC[i] > 0, "↑", "↓")
    cat(sprintf("      %2d. %-15s %s%.3f (p=%.2e)\n", 
                i, degs$gene[i], direction, abs(degs$logFC[i]), degs$p.value[i]))
  }
  
  # 保存结果
  write.csv(degs, file.path(knockout_dir, paste0("knockout_", target_gene, "_DEGs.csv")), row.names=FALSE)
  saveRDS(knk_res, file.path(knockout_dir, paste0("knockout_", target_gene, "_result.rds")))
  
  knockout_results[[target_gene]] <- list(
    result = knk_res,
    degs = degs,
    n_degs = nrow(degs)
  )
}

# ============================================================================
# 5. 汇总分析
# ============================================================================
cat("\n【3.5】敲除结果汇总...\n")

cat("\n  各基因敲除后的差异基因数:\n")
for(gene in names(knockout_results)) {
  cat(sprintf("    %s: %d 个差异基因\n", gene, knockout_results[[gene]]$n_degs))
}

# 寻找共同的下游基因
cat("\n  寻找共同的下游调控基因...\n")
all_degs <- lapply(knockout_results, function(x) x$degs$gene)

if(length(all_degs) >= 2) {
  common_genes <- Reduce(intersect, all_degs)
  cat(sprintf("    共同下游基因: %d 个\n", length(common_genes)))
  
  if(length(common_genes) > 0) {
    cat("    共同基因列表:", paste(head(common_genes, 10), collapse=", "), "...\n")
    
    # 保存共同基因
    write.csv(data.frame(gene=common_genes), 
              file.path(knockout_dir, "common_downstream_genes.csv"), row.names=FALSE)
  }
}

# ============================================================================
# 6. 可视化
# ============================================================================
cat("\n【3.6】生成可视化...\n")

pdf(file.path(knockout_dir, "knockout_summary.pdf"), width=14, height=10)

# 差异基因数柱状图
n_degs_list <- sapply(knockout_results, function(x) x$n_degs)
par(mar=c(5,4,4,2))
barplot(n_degs_list, main="Number of DEGs after Knockout",
        ylab="Number of DEGs", xlab="Knockout Gene",
        col="#2171B5", las=2)

# 每个敲除基因的火山图
for(gene in names(knockout_results)) {
  degs <- knockout_results[[gene]]$degs
  
  if(nrow(degs) > 0) {
    par(mfrow=c(1,1))
    
    # 火山图
    plot(degs$logFC, -log10(degs$p.value),
         pch=20, col=ifelse(degs$logFC > 0, "#CB181D", "#2171B5"),
         main=paste("Virtual Knockout of", gene),
         xlab="log2 Fold Change", ylab="-log10(p-value)")
    abline(h=-log10(0.05), lty=2, col="gray")
    abline(v=0, lty=2, col="gray")
    
    # 标记top基因
    top10 <- head(degs, 10)
    text(top10$logFC, -log10(top10$p.value),
         labels=top10$gene, cex=0.6, pos=3)
  }
}

dev.off()
cat("  ✅ 可视化已保存\n")

# 保存汇总结果
saveRDS(knockout_results, file.path(knockout_dir, "all_knockout_results.rds"))

cat("\n"); cat(rep("=", 60), sep=""); cat("\n")
cat("第3步完成！\n")
cat("结果目录:", knockout_dir, "\n")
cat("\n关键发现:\n")
for(gene in names(knockout_results)) {
  cat(sprintf("  - %s敲除: %d 个差异基因\n", gene, knockout_results[[gene]]$n_degs))
}
cat("  - 评估了关键Hub基因的功能重要性\n")
if(exists("common_genes") && length(common_genes) > 0) {
  cat(sprintf("  - 发现 %d 个共同下游调控基因\n", length(common_genes)))
}
cat(rep("=", 60), sep=""); cat("\n")
