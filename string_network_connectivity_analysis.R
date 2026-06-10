# ============================================
# String PPI网络连通性诊断分析
# ============================================

cat("正在加载必要的R包...\n")

if (!require("igraph", quietly = TRUE)) {
  install.packages("igraph", repos = "https://cloud.r-project.org/")
  library(igraph)
}

# 设置文件路径
work_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
edges_file <- "C:/Users/Jy-Mentor-7/Downloads/string_interactions (4).tsv"
nodes_file <- "C:/Users/Jy-Mentor-7/Downloads/string_node_degrees (2).tsv"
coords_file <- "C:/Users/Jy-Mentor-7/Downloads/string_network_coordinates (2).tsv"

output_dir <- file.path(work_dir, "String_Network_Connectivity_Results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 读取边列表
cat("\n读取String网络边列表...\n")
edges <- read.delim(edges_file, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
# 设置列名
if (ncol(edges) >= 13) {
  names(edges) <- c("node1", "node2", "node1_string_id", "node2_string_id", 
                    "neighborhood_on_chromosome", "gene_fusion", "phylogenetic_cooccurrence",
                    "homology", "coexpression", "experimentally_determined_interaction",
                    "database_annotated", "automated_textmining", "combined_score")
} else {
  names(edges) <- c("node1", "node2", paste0("V", 3:ncol(edges)))
}
cat(paste0("边列表列名: ", paste(names(edges), collapse = ", "), "\n"))
cat(paste0("总边数: ", nrow(edges), "\n"))

# 读取节点列表
cat("\n读取String网络节点列表...\n")
nodes <- read.delim(nodes_file, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
# 设置列名
if (ncol(nodes) >= 3) {
  names(nodes) <- c("node", "identifier", "node_degree")
} else {
  names(nodes) <- c("node", paste0("V", 2:ncol(nodes)))
}
cat(paste0("节点列表列名: ", paste(names(nodes), collapse = ", "), "\n"))
cat(paste0("总节点数: ", nrow(nodes), "\n"))

# 构建igraph图对象
cat("\n构建igraph图对象...\n")

# 创建边数据框（只需要前两列：node1和node2）
edge_df <- edges[, c("node1", "node2")]

# 创建节点数据框（包含属性）
node_df <- nodes

# 构建图
g <- graph_from_data_frame(edge_df, vertices = node_df, directed = FALSE)

cat(paste0("图构建完成\n"))

# ==================== 连通性诊断 ====================
cat("\n========================================\n")
cat("         String网络连通性诊断报告\n")
cat("========================================\n")

# 1. 基本统计
cat("\n【基本网络统计】\n")
cat(paste0("总节点数: ", vcount(g), "\n"))
cat(paste0("总边数: ", ecount(g), "\n"))

# 2. 连通分量分析
components <- components(g)
cat("\n【连通分量分析】\n")
cat(paste0("连通分量数: ", components$no, "\n"))
cat(paste0("最大连通分量节点数: ", max(components$csize), "\n"))
cat(paste0("最大连通分量占比: ", round(max(components$csize)/vcount(g)*100, 2), "%\n"))
cat(paste0("第二大连通分量节点数: ", sort(components$csize, decreasing = TRUE)[2], "\n"))

# 3. 网络密度
density <- edge_density(g, loops = FALSE)
cat("\n【网络密度】\n")
cat(paste0("网络密度: ", round(density, 6), "\n"))
cat(paste0("可能的最大边数: ", vcount(g) * (vcount(g) - 1) / 2, "\n"))

# 4. 平均路径长度（仅对最大连通分量）
cat("\n【路径长度分析（最大连通分量）】\n")
largest_component_id <- which.max(components$csize)
largest_component_nodes <- names(which(components$membership == largest_component_id))
sub_g <- induced_subgraph(g, largest_component_nodes)

avg_path_length <- average.path.length(sub_g, directed = FALSE)
diameter <- diameter(sub_g, directed = FALSE)

cat(paste0("最大连通分量节点数: ", vcount(sub_g), "\n"))
cat(paste0("最大连通分量边数: ", ecount(sub_g), "\n"))
cat(paste0("平均最短路径长度: ", round(avg_path_length, 3), "\n"))
cat(paste0("网络直径: ", diameter, "\n"))

# 5. 聚类系数
cat("\n【聚类系数】\n")
transitivity_global <- transitivity(g, type = "global")
transitivity_avg <- transitivity(g, type = "average")
cat(paste0("全局聚类系数: ", round(transitivity_global, 4), "\n"))
cat(paste0("平均局部聚类系数: ", round(transitivity_avg, 4), "\n"))

# 6. 度分布统计
cat("\n【度分布统计】\n")
degree_values <- degree(g)
cat(paste0("平均度: ", round(mean(degree_values), 2), "\n"))
cat(paste0("度中位数: ", median(degree_values), "\n"))
cat(paste0("最大度: ", max(degree_values), "\n"))
cat(paste0("最小度: ", min(degree_values), "\n"))
cat(paste0("度标准差: ", round(sd(degree_values), 2), "\n"))

# 7. 连通分量详细分布
cat("\n【连通分量大小分布】\n")
component_sizes <- table(components$csize)
component_sizes_df <- as.data.frame(component_sizes)
names(component_sizes_df) <- c("Component_Size", "Count")
cat("分量大小分布:\n")
print(component_sizes_df)

# 8. 孤立节点
cat("\n【孤立节点分析】\n")
isolated_nodes <- names(which(degree_values == 0))
cat(paste0("孤立节点数: ", length(isolated_nodes), "\n"))
if (length(isolated_nodes) > 0 && length(isolated_nodes) <= 20) {
  cat("孤立节点列表:\n")
  print(isolated_nodes)
} else if (length(isolated_nodes) > 0) {
  cat(paste0("孤立节点示例 (前20个): ", paste(head(isolated_nodes, 20), collapse = ", "), "\n"))
}

# 9. 中心性指标（Top 10）
cat("\n【网络中心性指标 (Top 10)】\n")

# 度中心性
degree_centrality <- degree(g)
top_degree <- head(sort(degree_centrality, decreasing = TRUE), 10)
cat("\n度中心性 Top 10:\n")
print(top_degree)

# 介数中心性（只对最大连通分量计算，因为其他分量可能有问题）
cat("\n计算介数中心性（最大连通分量）...\n")
betweenness_centrality <- betweenness(sub_g, directed = FALSE, normalized = TRUE)
top_betweenness <- head(sort(betweenness_centrality, decreasing = TRUE), 10)
cat("介数中心性 Top 10:\n")
print(top_betweenness)

# 接近中心性
cat("\n计算接近中心性（最大连通分量）...\n")
closeness_centrality <- closeness(sub_g, normalized = TRUE)
top_closeness <- head(sort(closeness_centrality, decreasing = TRUE), 10)
cat("接近中心性 Top 10:\n")
print(top_closeness)

# 特征向量中心性
cat("\n计算特征向量中心性（最大连通分量）...\n")
eigenvector_centrality <- eigen_centrality(sub_g, directed = FALSE)$vector
top_eigen <- head(sort(eigenvector_centrality, decreasing = TRUE), 10)
cat("特征向量中心性 Top 10:\n")
print(top_eigen)

# 10. 网络鲁棒性分析
cat("\n【网络鲁棒性分析】\n")

# 计算边介数（边的关键性）
edge_betweenness <- edge.betweenness(g)
cat(paste0("平均边介数: ", round(mean(edge_betweenness), 2), "\n"))
cat(paste0("最大边介数: ", max(edge_betweenness), "\n"))

# 桥接边（高介数的边可能连接不同模块）
high_betweenness_edges <- which(edge_betweenness > quantile(edge_betweenness, 0.95))
cat(paste0("高介数边数 (>95%分位数): ", length(high_betweenness_edges), "\n"))

# ==================== 保存结果 ====================
cat("\n保存分析结果...\n")

# 1. 连通分量成员
membership_df <- data.frame(
  Node = names(components$membership),
  Component_ID = components$membership,
  stringsAsFactors = FALSE
)
write.table(membership_df, file = file.path(output_dir, "component_membership.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 2. 连通分量统计
component_stats <- data.frame(
  Component_ID = 1:components$no,
  Size = components$csize,
  stringsAsFactors = FALSE
)
write.table(component_stats, file = file.path(output_dir, "component_statistics.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 3. 中心性指标汇总
centrality_df <- data.frame(
  Node = names(degree_centrality),
  Degree = as.numeric(degree_centrality),
  stringsAsFactors = FALSE
)
# 添加其他中心性（对于最大连通分量的节点）
centrality_df$Betweenness <- NA
centrality_df$Closeness <- NA
centrality_df$Eigenvector <- NA

for (node in names(betweenness_centrality)) {
  centrality_df$Betweenness[centrality_df$Node == node] <- betweenness_centrality[node]
  centrality_df$Closeness[centrality_df$Node == node] <- closeness_centrality[node]
  centrality_df$Eigenvector[centrality_df$Node == node] <- eigenvector_centrality[node]
}

write.table(centrality_df, file = file.path(output_dir, "centrality_measures.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 4. 孤立节点
if (length(isolated_nodes) > 0) {
  write.table(data.frame(Isolated_Node = isolated_nodes), 
              file = file.path(output_dir, "isolated_nodes.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# 5. 综合报告
report <- data.frame(
  Metric = c("总节点数", "总边数", "连通分量数", "最大连通分量节点数", 
             "最大连通分量占比(%)", "网络密度", "平均最短路径长度", 
             "网络直径", "全局聚类系数", "平均度", "最大度", 
             "孤立节点数"),
  Value = c(vcount(g), ecount(g), components$no, max(components$csize),
            round(max(components$csize)/vcount(g)*100, 2), round(density, 6),
            round(avg_path_length, 3), diameter, round(transitivity_global, 4),
            round(mean(degree_values), 2), max(degree_values), length(isolated_nodes))
)
write.table(report, file = file.path(output_dir, "connectivity_summary.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n========================================\n")
cat("         分析结果已保存\n")
cat(paste0("输出目录: ", output_dir, "\n"))
cat("生成文件:\n")
cat("  - component_membership.txt (节点所属连通分量)\n")
cat("  - component_statistics.txt (连通分量统计)\n")
cat("  - centrality_measures.txt (中心性指标)\n")
cat("  - isolated_nodes.txt (孤立节点列表)\n")
cat("  - connectivity_summary.txt (综合报告)\n")
cat("========================================\n")
