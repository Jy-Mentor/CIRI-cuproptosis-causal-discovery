# ============================================================================
# 第1步：PC因果发现算法分析 - v3 (修复版)
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
bulk_df <- read_excel(bulk_file, skip=5)

cat(sprintf("  原始数据: %d 行 x %d 列\n", nrow(bulk_df), ncol(bulk_df)))

# 提取基因符号和表达列
colnames(bulk_df)[1] <- "gene_symbol"

# 找到表达列（数值列）
expr_cols <- 7:12
colnames(bulk_df)[expr_cols] <- c("MCAO1","MCAO2","MCAO3","Sham1","Sham2","Sham3")

# 构建表达矩阵 - 使用正确的基因名列
bulk_expr <- as.matrix(bulk_df[, expr_cols])
rownames(bulk_expr) <- make.names(bulk_df$gene_symbol, unique=TRUE)

# 转换为数值并处理NA
bulk_expr <- apply(bulk_expr, 2, function(x) as.numeric(as.character(x)))
bulk_expr <- bulk_expr[!is.na(rownames(bulk_expr)), ]
bulk_expr <- bulk_expr[!apply(is.na(bulk_expr), 1, any), ]

cat(sprintf("  表达矩阵: %d 基因 x %d 样本\n", nrow(bulk_expr), ncol(bulk_expr)))

# 定义分组
bulk_pheno <- c(rep(1,3), rep(0,3))
names(bulk_pheno) <- colnames(bulk_expr)

cat(sprintf("  分组: MCAO=%d, Sham=%d\n", sum(bulk_pheno==1), sum(bulk_pheno==0)))

# ============================================================================
# 2. 差异分析
# ============================================================================
cat("\n【1.2】limma差异分析...\n")

design <- model.matrix(~0 + factor(bulk_pheno))
colnames(design) <- c("Sham", "MCAO")
contrast.matrix <- makeContrasts(MCAO - Sham, levels=design)

fit <- lmFit(bulk_expr, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg_res <- topTable(fit2, adjust="BH", number=Inf)

# 筛选显著差异基因 (FDR<0.05, |log2FC|>0.3)
sig_idx <- deg_res$adj.P.Val < 0.05 & abs(deg_res$logFC) > 0.3
sig_genes <- rownames(deg_res)[sig_idx]
cat(sprintf("  显著差异基因: %d\n", length(sig_genes)))

# 确保这些基因在表达矩阵中
sig_genes <- sig_genes[sig_genes %in% rownames(bulk_expr)]
cat(sprintf("  在矩阵中的基因: %d\n", length(sig_genes)))

# 选择前100个基因用于PC网络（样本只有6个，不能太多）
n_pc_genes <- min(100, length(sig_genes))
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
    
    # 识别hub基因
    node_degrees <- table(c(edges_df$from, edges_df$to))
    sorted_nodes <- sort(node_degrees, decreasing=TRUE)
    n_hub <- min(10, length(sorted_nodes))
    hub_genes <- names(sorted_nodes)[1:n_hub]
    
    cat("\n  Top Hub基因（关键驱动基因）:\n")
    for(i in 1:n_hub) {
      gene_name <- hub_genes[i]
      degree <- sorted_nodes[i]
      
      # 获取差异表达信息
      if(gene_name %in% rownames(deg_res)) {
        logfc <- deg_res[gene_name, "logFC"]
        fdr <- deg_res[gene_name, "adj.P.Val"]
        cat(sprintf("    %d. %s (连接数: %d, log2FC: %.2f, FDR: %.2e)\n", 
                    i, gene_name, degree, logfc, fdr))
      } else {
        cat(sprintf("    %d. %s (连接数: %d)\n", i, gene_name, degree))
      }
    }
    
    # 保存hub基因
    hub_df <- data.frame(
      gene = hub_genes,
      degree = as.numeric(sorted_nodes[1:n_hub]),
      stringsAsFactors = FALSE
    )
    
    # 添加差异表达信息
    hub_df$log2FC <- deg_res[hub_genes, "logFC"]
    hub_df$FDR <- deg_res[hub_genes, "adj.P.Val"]
    
    write.csv(hub_df, file.path(result_dir, "PC_hub_genes.csv"), row.names=FALSE)
    
    # 绘制网络图
    cat("\n  绘制因果网络图...\n")
    pdf(file.path(result_dir, "PC_causal_network.pdf"), width=14, height=12)
    
    g <- graph_from_data_frame(edges_df[,1:2], directed=TRUE)
    
    # 设置节点属性
    node_colors <- ifelse(V(g)$name %in% hub_genes[1:5], "#CB181D", "#2171B5")
    V(g)$size <- sqrt(node_degrees[V(g)$name]) * 5
    V(g)$color <- node_colors
    V(g)$label.cex <- 0.7
    
    # 设置边属性
    E(g)$width <- edges_df$strength * 3
    edge_colors <- ifelse(edges_df$correlation > 0, "#2171B5", "#CB181D")
    E(g)$color <- edge_colors
    E(g)$arrow.size <- 0.4
    
    layout <- layout_with_fr(g)
    plot(g, layout=layout,
         vertex.label.color="black",
         vertex.label.font=2,
         main="PC Causal Network - GSE163614 (MCAO)")
    
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

cat("\n=== 第1步完成 ===\n")
cat("结果保存到:", result_dir, "\n")
cat("\n下一步：使用PC识别的Hub基因进行SCISSOR-like分析\n")
