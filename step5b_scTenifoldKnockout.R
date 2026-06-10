# scTenifoldKnk虚拟敲除分析 - 专业版
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

# 加载scTenifoldKnk
library(scTenifoldKnk)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(parallel)

cat("=== scTenifoldKnk虚拟敲除分析 ===\n")
cat("目标: 验证RAGE-NFKB1-FDX1因果轴\n\n")

# 1. 读取单细胞count矩阵
cat("1. 读取单细胞数据...\n")
sc_Matrix <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_Matrix.rds')
cat(sprintf("  Count矩阵: %d 基因 x %d 细胞\n", nrow(sc_Matrix), ncol(sc_Matrix)))

# 2. 定义敲除基因 (基于PC网络)
knk_genes <- c("Nfkb1", "Fdx1", "Tlr4")

# 3. 检查基因是否存在
cat("\n2. 检查目标基因...\n")
for(g in knk_genes) {
  if(g %in% rownames(sc_Matrix)){
    cat(sprintf("  [OK] %s 在数据中\n", g))
  } else {
    cat(sprintf("  [WARN] %s 不在数据中\n", g))
  }
}

# 4. 执行虚拟敲除
cat("\n3. 执行scTenifoldKnk虚拟敲除...\n")
dir.create('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scTenifoldKnk_results', showWarnings=FALSE, recursive=TRUE)

for(knk_gene in knk_genes) {
  cat(sprintf("\n  正在敲除: %s\n", knk_gene))

  if(!knk_gene %in% rownames(sc_Matrix)) {
    cat(sprintf("  [跳过] %s不在数据中\n", knk_gene))
    next
  }

  set.seed(666)
  knk_res <- scTenifoldKnk(
    countMatrix = sc_Matrix,
    gKO = knk_gene,
    qc = TRUE,
    qc_mtThreshold = 0.1,
    qc_minLSize = 1000,
    nc_lambda = 0,
    nc_nNet = 10,
    nc_nCells = 500,
    nc_nComp = 3,
    nc_scaleScores = TRUE,
    nc_symmetric = FALSE,
    nc_q = 0.9,
    td_K = 3,
    td_maxIter = 1000,
    td_maxError = 1e-05,
    td_nDecimal = 3,
    ma_nDim = 2,
    nCores = parallel::detectCores()
  )

  # 5. 提取差异结果
  knk_df <- knk_res$diffRegulation %>%
    mutate(logFC = log2(.data$FC)) %>%
    filter(gene != knk_gene) %>%
    filter(p.value < 0.05)
  knk_df[,2:7] <- sapply(knk_df[,2:7], as.numeric)

  # 6. 保存结果
  write.xlsx(knk_df, paste0('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scTenifoldKnk_results/', knk_gene, '.xlsx'))

  cat(sprintf("  差异基因数: %d\n", nrow(knk_df)))

  # 7. 可视化 - 柱状图
  p1 <- ggplot(knk_df, aes(x = reorder(gene, logFC), y = logFC)) +
    geom_bar(stat = 'identity', fill='orange') +
    coord_flip() +
    labs(title = paste("Differentially Regulated Genes after", knk_gene, "Knockout"),
         x = "", y = "logFC") +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12, hjust = .5))

  ggsave(plot = p1,
         device = 'png',
         width = 8, height = 10,
         dpi = 300,
         filename = paste0('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scTenifoldKnk_results/', knk_gene, '_bar.png'))

  # 8. 火山图
  df <- knk_res$diffRegulation %>%
    mutate(label_genes = ifelse(gene %in% knk_df$gene, gene, NA))

  p2 <- ggplot(df, aes(x = log2(FC), y = -log10(p.value))) +
    geom_point(alpha = ifelse(df$gene %in% df$label_genes, 0.9, 0.4),
               fill = ifelse(df$gene %in% df$label_genes, '#ffd6a5', 'gray'),
               shape = 21, color = 'white',
               aes(size = -log10(p.value))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text_repel(data = df, aes(label = label_genes),
                    size = 2, parse = F, fontface = 'bold.italic',
                    color = '#c7522a', max.overlaps = 20) +
    labs(title = paste('Differentially Regulated Genes:', kng_gene),
         x = "logFC", y = "-log10(Pvalue)") +
    theme_bw() +
    theme(plot.title = element_text(size = 12, hjust = .5)) +
    guides(size = 'none')

  ggsave(plot = p2,
         device = 'png',
         width = 8, height = 8,
         dpi = 300,
         filename = paste0('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scTenifoldKnk_results/', knk_gene, '_volcano.png'))

  cat(sprintf("  结果已保存: %s\n", knk_gene))
}

cat("\n=== 完成 ===\n")
cat("结果保存到: C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scTenifoldKnk_results/\n")