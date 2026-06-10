# ============================================================================
# 第1步：PC因果发现算法分析
# 数据集：GSE163614（大鼠MCAO模型Bulk RNA-seq）
# 目标：构建基因调控网络，识别MCAO关键驱动基因
# ============================================================================

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"bnlearn" %in% installed.packages()){install.packages('bnlearn')}
if(!"limma" %in% installed.packages()){BiocManager::install('limma')}
if(!"readxl" %in% installed.packages()){install.packages('readxl')}
library(bnlearn)
library(limma)
library(readxl)
library(ggplot2)
library(igraph)

cat("=== 第1步：PC因果发现算法（GSE163614）===\n\n")

# 设置路径
result_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/causal_analysis_results"
dir.create(result_dir, showWarnings=FALSE, recursive=TRUE)

# ============================================================================
# 1. 读取GSE163614 Bulk数据
# ============================================================================
cat("【1.1】读取GSE163614数据...\n")

bulk_file <- "C:/Users/Jy-Mentor-7/Downloads/GSE163614_mRNA_Expression_Profiling.xlsx"
bulk_df <- read_excel(bulk_file, skip=5)

cat(sprintf("  原始数据: %d 行 x %d 列\n", nrow(bulk_df), ncol(bulk_df)))

# 提取样本列 (假设第7-12列是样本)
colnames(bulk_df)[1] <- "gene_symbol"
colnames(bulk_df)[7:12] <- c("MCAO1","MCAO2","MCAO3","Sham1","Sham2","Sham3")

# 构建表达矩阵
bulk_expr <- as.matrix(bulk_df[,7:12])
rownames(bulk_expr) <- bulk_df$gene_symbol
bulk_expr <- apply(bulk_expr, 2, as.numeric)
rownames(bulk_expr) <- bulk_df$gene_symbol

cat(sprintf("  表达矩阵: %d 基因 x %d 样本\n", nrow(bulk_expr), ncol(bulk_expr)))

# 定义分组
bulk_pheno <- c(rep(1,3), rep(0,3))  # MCAO=1, Sham=0
names(bulk_pheno) <- colnames(bulk_expr)

cat(sprintf("  分组: MCAO=%d, Sham=%d\n", sum(bulk_pheno==1), sum(bulk_pheno==0)))

# ============================================================================
# 2. 差异分析（识别MCAO相关基因）
# ============================================================================
cat("\n【1.2】limma差异分析...\n")

design <- model.matrix(~0 + factor(bulk_pheno))
colnames(design) <- c("Sham", "MCAO")
contrast.matrix <- makeContrasts(MCAO - Sham, levels=design)

fit <- lmFit(bulk_expr, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg_res <- topTable(fit2, adjust="BH", number=Inf)

# 筛选显著差异基因 (FDR<0.05, |log2FC|>0.5)
sig_genes <- rownames(deg_res)[deg_res$adj.P.Val < 0.05 & abs(deg_res$logFC) > 0.5]
cat(sprintf("  显著差异基因: %d\n", length(sig_genes)))

# 选择top基因用于PC网络 (太多基因会导致计算不可行)
# 选择标准：FDR最小的前200个基因
top_genes <- head(sig_genes, 200)
cat(sprintf("  用于PC网络的基因: %d\n", length(top_genes)))

# ============================================================================
# 3. PC因果网络构建
# ============================================================================
cat("\n【1.3】PC因果网络构建...\n")

# 提取top基因的表达数据
pc_data <- t(bulk_expr[top_genes, ])
pc_data <- as.data.frame(pc_data)

# 确保所有列都是数值型
for(col in colnames(pc_data)) {
  pc_data[[col]] <- as.numeric(pc_data[[col]])
}

cat(sprintf("  PC输入矩阵: %d 样本 x %d 基因\n", nrow(pc_data), ncol(pc_data)))

# 使用pc.stable进行稳定的PC算法（带bootstrap）
set.seed(42)
cat("  运行PC-stable算法（带bootstrap）...\n")

tryCatch({
  # PC算法 (alpha=0.05)
  pc_result <- pc.stable(pc_data, alpha=0.05, test="mi-g-sh", B=100)

  # 提取边
  edges <- pc_result$arcs
  cat(sprintf("  发现 %d 条因果边\n", nrow(edges)))

  if(nrow(edges) > 0) {
    # 转换为数据框
    edges_df <- data.frame(
      from = edges[,1],
      to = edges[,2],
      stringsAsFactors = FALSE
    )

    # 计算边的强度 (基于偏相关)
    cat("  计算边强度...\n")
    edge_strengths <- sapply(1:nrow(edges_df), function(i) {
      from_expr <- pc_data[[edges_df$from[i]]]
      to_expr <- pc_data[[edges_df$to[i]]]
      abs(cor(from_expr, to_expr, method="spearman"))
    })
    edges_df$strength <- edge_strengths

    # 保存结果
    write.csv(edges_df, file.path(result_dir, "PC_causal_network_edges.csv"), row.names=FALSE)

    # 识别hub基因（连接数最多的基因）
    node_degrees <- table(c(edges_df$from, edges_df$to))
    hub_genes <- names(sort(node_degrees, decreasing=TRUE))[1:10]
    cat("\n  Top 10 Hub基因（潜在关键驱动基因）:\n")
    for(i in 1:10) {
      cat(sprintf("    %d. %s (连接数: %d)\n", i, hub_genes[i], node_degrees[hub_genes[i]]))
    }

    # 保存hub基因
    hub_df <- data.frame(
      gene = hub_genes,
      degree = as.numeric(node_degrees[hub_genes]),
      stringsAsFactors = FALSE
    )
    write.csv(hub_df, file.path(result_dir, "PC_hub_genes.csv"), row.names=FALSE)

    # 绘制网络图
    cat("\n  绘制因果网络图...\n")
    pdf(file.path(result_dir, "PC_causal_network.pdf"), width=14, height=12)

    g <- graph_from_data_frame(edges_df[,1:2], directed=TRUE)

    # 设置节点大小基于degree
    V(g)$size <- sqrt(node_degrees[V(g)$name]) * 5
    V(g)$color <- ifelse(V(g)$name %in% hub_genes[1:5], "#CB181D", "#2171B5")
    E(g)$width <- edges_df$strength * 2
    E(g)$arrow.size <- 0.5

    layout <- layout_with_fr(g)
    plot(g, layout=layout,
         vertex.label.cex=0.8,
         vertex.label.color="black",
         vertex.label.font=2,
         main="PC Causal Network - GSE163614 (MCAO)")

    dev.off()
    cat("  网络图已保存\n")

    # 保存完整结果
    saveRDS(list(
      pc_result = pc_result,
      edges = edges_df,
      hub_genes = hub_genes,
      top_genes = top_genes,
      deg_res = deg_res
    ), file.path(result_dir, "PC_analysis_result.rds"))

  } else {
    cat("  [警告] 未发现显著因果边\n")
  }

}, error = function(e) {
  cat("  [错误] PC算法运行失败:", conditionMessage(e), "\n")
  cat("  尝试简化版本...\n")

  # 简化：仅基于相关性构建网络
  cat("  使用相关性网络作为替代...\n")
  cor_matrix <- cor(pc_data, method="spearman")

  # 筛选强相关 (|r| > 0.8)
  strong_cor <- which(abs(cor_matrix) > 0.8 & abs(cor_matrix) < 1, arr.ind=TRUE)
  if(nrow(strong_cor) > 0) {
    edges_df <- data.frame(
      from = colnames(pc_data)[strong_cor[,1]],
      to = colnames(pc_data)[strong_cor[,2]],
      strength = cor_matrix[strong_cor],
      stringsAsFactors = FALSE
    )
    write.csv(edges_df, file.path(result_dir, "correlation_network_edges.csv"), row.names=FALSE)
    cat(sprintf("  相关性网络: %d 条边\n", nrow(edges_df)))
  }
})

# ============================================================================
# 4. BCP轴基因分析
# ============================================================================
cat("\n【1.4】BCP轴基因分析...\n")

bcp_genes <- c("Ager","Nfkb1","Fdx1","Tlr4","Stat1","Stat3","Tgfbr1","Nfe2l2")
bcp_found <- bcp_genes[bcp_genes %in% rownames(deg_res)]

if(length(bcp_found) > 0) {
  bcp_df <- data.frame(
    gene = bcp_found,
    log2FC = deg_res[bcp_found, "logFC"],
    FDR = deg_res[bcp_found, "adj.P.Val"],
    stringsAsFactors = FALSE
  )
  bcp_df <- bcp_df[order(-abs(bcp_df$log2FC)), ]

  cat("  BCP轴基因差异表达:\n")
  print(bcp_df)

  write.csv(bcp_df, file.path(result_dir, "BCP_axis_genes_DEG.csv"), row.names=FALSE)
}

cat("\n=== 第1步完成 ===\n")
cat("结果保存到:", result_dir, "\n")
cat("\n关键发现：\n")
cat("- PC因果网络已构建\n")
cat("- Hub基因已识别（潜在关键驱动基因）\n")
cat("- 可用于下一步SCISSOR-like分析\n")
