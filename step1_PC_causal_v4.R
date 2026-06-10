# ============================================================================
# 第1步：PC因果发现算法分析 - v4 (修复数据读取)
# ============================================================================

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

library(bnlearn)
library(limma)
library(readxl)
library(ggplot2)
library(igraph)

cat("=== 第1步：PC因果发现算法（GSE163614）===\n\n")

result_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/causal_analysis_results"
dir.create(result_dir, showWarnings=FALSE, recursive=TRUE)

# ============================================================================
# 1. 读取GSE163614 Bulk数据
# ============================================================================
cat("【1.1】读取GSE163614数据...\n")

bulk_file <- "C:/Users/Jy-Mentor-7/Downloads/GSE163614_mRNA_Expression_Profiling.xlsx"

# 跳过前5行，第6行是列名
bulk_df <- read_excel(bulk_file, skip=5)

cat(sprintf("  数据维度: %d 行 x %d 列\n", nrow(bulk_df), ncol(bulk_df)))
cat("  列名:\n")
print(head(colnames(bulk_df), 15))

# 提取基因名和表达值
# 假设第1列是基因名，第7-12列是表达值（MCAO1-3, Sham1-3）
gene_col <- 1
expr_cols <- 7:12  # FPKM列

# 提取数据
bulk_expr <- as.data.frame(bulk_df[, expr_cols])
rownames(bulk_expr) <- make.names(bulk_df[[gene_col]], unique=TRUE)
colnames(bulk_expr) <- c("MCAO1","MCAO2","MCAO3","Sham1","Sham2","Sham3")

# 转换为数值
for(col in colnames(bulk_expr)) {
  bulk_expr[[col]] <- as.numeric(as.character(bulk_expr[[col]]))
}

# 移除NA和零方差基因
bulk_expr <- bulk_expr[!apply(is.na(bulk_expr), 1, any), ]
bulk_expr <- bulk_expr[apply(bulk_expr, 1, var) > 0, ]

cat(sprintf("  过滤后: %d 基因 x %d 样本\n", nrow(bulk_expr), ncol(bulk_expr)))
cat(sprintf("  示例基因: %s\n", paste(head(rownames(bulk_expr), 5), collapse=", ")))

# 定义分组
bulk_pheno <- c(rep(1,3), rep(0,3))
names(bulk_pheno) <- colnames(bulk_expr)
cat(sprintf("  分组: MCAO=%d, Sham=%d\n", sum(bulk_pheno==1), sum(bulk_pheno==0)))

# ============================================================================
# 2. 差异分析
# ============================================================================
cat("\n【1.2】limma差异分析...\n")

bulk_expr_matrix <- as.matrix(bulk_expr)

design <- model.matrix(~0 + factor(bulk_pheno))
colnames(design) <- c("Sham", "MCAO")
contrast.matrix <- makeContrasts(MCAO - Sham, levels=design)

fit <- lmFit(bulk_expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg_res <- topTable(fit2, adjust="BH", number=Inf)

cat(sprintf("  差异分析完成，共 %d 个基因\n", nrow(deg_res)))

# 筛选显著差异基因
sig_idx <- deg_res$adj.P.Val < 0.05 & abs(deg_res$logFC) > 0.3
sig_genes <- rownames(deg_res)[sig_idx]
cat(sprintf("  显著差异基因 (FDR<0.05, |log2FC|>0.3): %d\n", length(sig_genes)))

# 确保基因在矩阵中
sig_genes <- sig_genes[sig_genes %in% rownames(bulk_expr)]
cat(sprintf("  在表达矩阵中的基因: %d\n", length(sig_genes)))

# 选择前80个基因用于PC网络（样本只有6个，不能太多）
n_pc_genes <- min(80, length(sig_genes))
if(n_pc_genes < 10) {
  cat("  [警告] 显著基因太少，使用top 100变异最大基因\n")
  gene_vars <- apply(bulk_expr, 1, var)
  sig_genes <- names(sort(gene_vars, decreasing=TRUE))[1:100]
  n_pc_genes <- 80
}

top_genes <- head(sig_genes, n_pc_genes)
cat(sprintf("  用于PC网络的基因: %d\n", length(top_genes)))

# ============================================================================
# 3. PC因果网络构建
# ============================================================================
cat("\n【1.3】PC因果网络构建...\n")

# 提取top基因的表达数据
pc_data <- t(bulk_expr[top_genes, , drop=FALSE])
pc_data <- as.data.frame(pc_data)

cat(sprintf("  PC输入矩阵: %d 样本 x %d 基因\n", nrow(pc_data), ncol(pc_data)))

# 使用Hill-Climbing算法
cat("  运行Hill-Climbing算法...\n")
set.seed(42)

tryCatch({
  hc_result <- hc(pc_data, score="bic-g", max.iter=50)
  
  # 提取边
  edges <- hc_result$arcs
  cat(sprintf("  发现 %d 条边\n", nrow(edges)))
  
  if(nrow(edges) > 0) {
    # 转换为数据框
    edges_df <- data.frame(
      from = edges[,1],
      to = edges[,2],
      stringsAsFactors = FALSE
    )
    
    # 计算边的强度
    edge_strengths <- sapply(1:nrow(edges_df), function(i) {
      from_expr <- pc_data[[edges_df$from[i]]]
      to_expr <- pc_data[[edges_df$to[i]]]
      cor(from_expr, to_expr, method="spearman")
    })
    edges_df$correlation <- edge_strengths
    edges_df$strength <- abs(edge_strengths)
    
    # 保存结果
    write.csv(edges_df, file.path(result_dir, "PC_causal_network_edges.csv"), row.names=FALSE)
    cat(sprintf("  保存了 %d 条边\n", nrow(edges_df)))
    
    # 识别hub基因
    node_degrees <- table(c(edges_df$from, edges_df$to))
    sorted_nodes <- sort(node_degrees, decreasing=TRUE)
    n_hub <- min(15, length(sorted_nodes))
    hub_genes <- names(sorted_nodes)[1:n_hub]
    
    cat("\n  Top Hub基因（关键驱动基因）:\n")
    for(i in 1:n_hub) {
      gene_name <- hub_genes[i]
      degree <- sorted_nodes[i]
      
      # 获取差异表达信息
      if(gene_name %in% rownames(deg_res)) {
        logfc <- deg_res[gene_name, "logFC"]
        fdr <- deg_res[gene_name, "adj.P.Val"]
        direction <- ifelse(logfc > 0, "↑", "↓")
        cat(sprintf("    %2d. %-15s (连接: %2d, %s%.2f, FDR: %.2e)\n", 
                    i, gene_name, degree, direction, abs(logfc), fdr))
      } else {
        cat(sprintf("    %2d. %-15s (连接: %2d)\n", i, gene_name, degree))
      }
    }
    
    # 保存hub基因
    hub_df <- data.frame(
      gene = hub_genes,
      degree = as.numeric(sorted_nodes[1:n_hub]),
      stringsAsFactors = FALSE
    )
    hub_df$log2FC <- deg_res[hub_genes, "logFC"]
    hub_df$FDR <- deg_res[hub_genes, "adj.P.Val"]
    write.csv(hub_df, file.path(result_dir, "PC_hub_genes.csv"), row.names=FALSE)
    
    # 绘制网络图
    cat("\n  绘制因果网络图...\n")
    pdf(file.path(result_dir, "PC_causal_network.pdf"), width=14, height=12)
    
    g <- graph_from_data_frame(edges_df[,1:2], directed=TRUE)
    
    # 设置节点属性
    node_colors <- ifelse(V(g)$name %in% hub_genes[1:5], "#CB181D", "#2171B5")
    V(g)$size <- sqrt(node_degrees[V(g)$name]) * 6
    V(g)$color <- node_colors
    V(g)$label.cex <- 0.8
    
    # 设置边属性
    E(g)$width <- edges_df$strength * 3
    edge_colors <- ifelse(edges_df$correlation > 0, "#2171B5", "#CB181D")
    E(g)$color <- edge_colors
    E(g)$arrow.size <- 0.4
    
    layout <- layout_with_fr(g)
    plot(g, layout=layout,
         vertex.label.color="black",
         vertex.label.font=2,
         main="PC Causal Network - GSE163614 (MCAO vs Sham)")
    
    dev.off()
    cat("  ✅ 网络图已保存\n")
    
    # 保存完整结果
    saveRDS(list(
      network = hc_result,
      edges = edges_df,
      hub_genes = hub_genes,
      top_genes = top_genes,
      deg_res = deg_res,
      node_degrees = node_degrees
    ), file.path(result_dir, "PC_analysis_result.rds"))
    
    cat("\n✅ PC网络构建成功!\n")
    
  } else {
    cat("  [警告] 未发现边\n")
  }
  
}, error = function(e) {
  cat("  [错误]", conditionMessage(e), "\n")
  cat("  尝试使用correlation网络...\n")
  
  # 备用方案：相关性网络
  cor_matrix <- cor(pc_data, method="spearman")
  cor_threshold <- 0.7
  strong_cor <- which(abs(cor_matrix) > cor_threshold & abs(cor_matrix) < 1, arr.ind=TRUE)
  
  if(nrow(strong_cor) > 0) {
    edges_df <- data.frame(
      from = colnames(pc_data)[strong_cor[,1]],
      to = colnames(pc_data)[strong_cor[,2]],
      correlation = cor_matrix[strong_cor],
      strength = abs(cor_matrix[strong_cor]),
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

bcp_genes <- c("Ager","Nfkb1","Fdx1","Tlr4","Stat1","Stat3","Tgfbr1","Nfe2l2","Sod1","Cat")
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
  
  # 检查是否有hub基因是BCP轴基因
  bcp_hub <- intersect(bcp_found, hub_genes)
  if(length(bcp_hub) > 0) {
    cat("\n  ✅ BCP轴基因中的Hub基因:", paste(bcp_hub, collapse=", "), "\n")
  }
}

cat("\n"); cat(rep("=", 60), sep=""); cat("\n")
cat("第1步完成！\n")
cat("结果目录:", result_dir, "\n")
cat("关键文件:\n")
cat("  - PC_causal_network_edges.csv (因果网络边)\n")
cat("  - PC_hub_genes.csv (关键驱动基因)\n")
cat("  - PC_causal_network.pdf (网络图)\n")
cat(rep("=", 60), sep=""); cat("\n")
