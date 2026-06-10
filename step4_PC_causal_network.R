# PC因果网络分析 - RAGE→NFKB1→FDX1轴
# 使用bnlearn包实现PC算法
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"bnlearn" %in% installed.packages()){install.packages('bnlearn')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
if(!"igraph" %in% installed.packages()){install.packages('igraph')}
if(!"gplots" %in% installed.packages()){install.packages('gplots')}
if(!"Seurat" %in% installed.packages()){install.packages('Seurat', repos='https://cran.rstudio.com')}
library(bnlearn)
library(ggplot2)
library(igraph)
library(gplots)
library(Seurat)

cat("=== PC因果网络分析 ===\n")
cat("目标: 验证 RAGE→NFKB1→FDX1 因果轴\n\n")

# 1. 读取单细胞数据
cat("1. 读取单细胞数据...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
cat(sprintf("  单细胞: %d 细胞 x %d 基因\n", ncol(sc_obj), nrow(sc_obj)))

# 2. 读取DEG结果获取基因列表
deg_list <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/scissor_plots_v10/DEGs.rds')
lfc <- deg_list[['log2fc']]
fdr <- deg_list[['fdr']]

# 3. 选择目标基因构建网络
cat("\n2. 选择目标基因...\n")
target_genes <- c("Ager","Nfkb1","Fdx1","Tlr4","Stat1","Stat3","Tgfbr1","Nfe2l2",
                  "Slc31a1","Atp7a","Atox1","Lias","Dlat","Pdha1","Pdhb","Dld",
                  "Cat","Sod1","Lipt1","Gcsh")

cat("  构建网络的基因:", paste(target_genes, collapse=", "), "\n")

# 4. 提取单细胞表达矩阵
cat("\n3. 提取表达矩阵...\n")
sc_expr <- GetAssayData(sc_obj, layer="data")
common_genes <- intersect(target_genes, rownames(sc_expr))
sc_expr <- sc_expr[common_genes, ]
cat(sprintf("  表达矩阵: %d 基因 x %d 细胞\n", nrow(sc_expr), ncol(sc_expr)))

# 随机抽样控制计算量
set.seed(42)
max_cells <- 3000
if(ncol(sc_expr) > max_cells) {
  cell_idx <- sample(ncol(sc_expr), max_cells)
  sc_expr <- sc_expr[, cell_idx]
  cat(sprintf("  随机抽样至: %d 细胞\n", max_cells))
}

# 5. 转置为样本x基因矩阵
expr_df <- as.data.frame(t(as.matrix(sc_expr)))
colnames(expr_df) <- make.names(colnames(expr_df))
expr_df <- expr_df[, colSums(is.na(expr_df)) == 0]
cat(sprintf("  最终矩阵: %d 样本 x %d 基因\n", nrow(expr_df), ncol(expr_df)))

# 6. 运行PC算法
cat("\n4. 运行PC算法...\n")

pc_result <- pc.stable(
  x = expr_df,
  alpha = 0.05,
  test = "cor",
  max.sx = 3
)
cat("  PC算法完成\n")

# 7. 构建网络
cat("\n5. 构建因果网络...\n")

edges_df <- as.data.frame(arcs(pc_result))
if(nrow(edges_df) > 0) {
  colnames(edges_df) <- c("from", "to")
  edges_df$strength <- 0.8
}

# 8. 保存网络结果
cat("\n6. 保存结果...\n")
dir.create('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/pc_network', showWarnings=FALSE, recursive=TRUE)

saveRDS(list(
  pc_result = pc_result,
  edges = edges_df
), 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/pc_network/network_result.rds')

# 9. 绘制网络图
cat("\n7. 绘制因果网络图...\n")

pdf('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/pc_network/causal_network.pdf', width=12, height=10)

if(nrow(edges_df) > 0) {
  g <- graph_from_data_frame(edges_df, directed=TRUE)
  layout <- layout_with_fr(g)

  node_colors <- rep("#A8D5E5", vcount(g))
  V(g)$color <- node_colors
  V(g)$size <- 30
  E(g)$width <- 2
  E(g)$color <- "#2171B5"

  plot(g, layout=layout,
       vertex.label.color="black",
       vertex.label.font=2,
       edge.arrow.size=0.5,
       main="PC Causal Network: BCP-Cuproptosis Axis")
} else {
  plot.new()
  text(0.5, 0.5, "No significant causal edges found\n(Try adjusting alpha or increasing sample size)")
}

dev.off()

# 10. 输出关键因果关系
cat("\n=== 关键因果关系 ===\n")
if(nrow(edges_df) > 0) {
  cat("\n所有因果边:\n")
  print(edges_df)

  # 关注RAGE→NFKB1→FDX1轴
  axis_edges <- edges_df[grep("Ager|Nfkb1|Fdx1", edges_df$from) |
                          grep("Ager|Nfkb1|Fdx1", edges_df$to), ]
  if(nrow(axis_edges) > 0) {
    cat("\nBCP轴相关因果边:\n")
    print(axis_edges)
  }
} else {
  cat("  未发现显著因果边\n")
}

if(nrow(edges_df) > 0) {
  write.csv(edges_df, 'C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/pc_network/causal_edges.csv', row.names=FALSE)
}

cat("\n=== 完成 ===\n")
cat("结果保存到: C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/pc_network/\n")