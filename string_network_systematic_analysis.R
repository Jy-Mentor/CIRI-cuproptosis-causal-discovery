# ============================================
# String PPI网络系统性生物信息学分析
# 包含：K-core分解、中心性指标计算、RRA综合评分
# ============================================

cat("正在加载必要的R包...\n")

packages <- c("igraph", "dplyr", "tidyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

# 设置文件路径
work_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
edges_file <- "C:/Users/Jy-Mentor-7/Downloads/string_interactions (4).tsv"
nodes_file <- "C:/Users/Jy-Mentor-7/Downloads/string_node_degrees (2).tsv"

output_dir <- file.path(work_dir, "String_Network_Systematic_Analysis")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 读取边列表
cat("\n读取String网络边列表...\n")
edges <- read.delim(edges_file, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
if (ncol(edges) >= 13) {
  names(edges) <- c("node1", "node2", "node1_string_id", "node2_string_id", 
                    "neighborhood_on_chromosome", "gene_fusion", "phylogenetic_cooccurrence",
                    "homology", "coexpression", "experimentally_determined_interaction",
                    "database_annotated", "automated_textmining", "combined_score")
}
cat(paste0("总边数: ", nrow(edges), "\n"))

# 读取节点列表
cat("\n读取String网络节点列表...\n")
nodes <- read.delim(nodes_file, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
if (ncol(nodes) >= 3) {
  names(nodes) <- c("node", "identifier", "node_degree")
}
cat(paste0("总节点数: ", nrow(nodes), "\n"))

# 构建igraph图对象
cat("\n构建igraph图对象...\n")
edge_df <- edges[, c("node1", "node2")]
g <- graph_from_data_frame(edge_df, vertices = nodes, directed = FALSE)
cat(paste0("原始图 - 节点数: ", vcount(g), ", 边数: ", ecount(g), "\n"))

# ==================== 步骤1：删除孤立节点 ====================
cat("\n========================================\n")
cat("步骤1: 删除孤立节点（度为0的节点）\n")
cat("========================================\n")

# 计算度
degree_values <- degree(g)
isolated_nodes <- names(which(degree_values == 0))
cat(paste0("孤立节点数: ", length(isolated_nodes), "\n"))
if (length(isolated_nodes) > 0) {
  cat("孤立节点列表:\n")
  print(isolated_nodes)
}

# 删除孤立节点
g_clean <- delete_vertices(g, isolated_nodes)
cat(paste0("\n预处理后图 - 节点数: ", vcount(g_clean), ", 边数: ", ecount(g_clean), "\n"))

# 保存孤立节点列表
if (length(isolated_nodes) > 0) {
  write.table(data.frame(Isolated_Node = isolated_nodes), 
              file = file.path(output_dir, "01_isolated_nodes_removed.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ==================== 步骤2：K-core分解分析 ====================
cat("\n========================================\n")
cat("步骤2: K-core分解分析\n")
cat("========================================\n")

# 执行K-core分解
cores <- coreness(g_clean)
cat("K-core分布:\n")
core_table <- table(cores)
print(core_table)

# 保存K-core结果
kcore_df <- data.frame(
  Node = names(cores),
  K_core = as.numeric(cores),
  stringsAsFactors = FALSE
)
write.table(kcore_df, file = file.path(output_dir, "02_kcore_decomposition.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 提取K=3的核心子图（g_k3）
g_k3 <- induced_subgraph(g_clean, names(cores)[cores >= 3])
cat(paste0("\nK>=3核心子图 - 节点数: ", vcount(g_k3), ", 边数: ", ecount(g_k3), "\n"))

# 保存K>=3子图的节点和边
write.table(data.frame(Node = V(g_k3)$name), 
            file = file.path(output_dir, "02_kcore_k3_nodes.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 保存g_k3边列表
k3_edges <- as.data.frame(get.edgelist(g_k3))
write.table(k3_edges, file = file.path(output_dir, "02_kcore_k3_edges.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("source", "target"))

# 保存igraph对象
saveRDS(g_k3, file = file.path(output_dir, "02_g_k3_network_object.rds"))
cat("g_k3网络对象已保存\n")

# ==================== 步骤3：网络中心性指标计算 ====================
cat("\n========================================\n")
cat("步骤3: 网络中心性指标计算\n")
cat("========================================\n")

# 使用原始图（保留所有节点）进行中心性计算
# 1. 度中心性 (DC)
cat("计算度中心性 (DC)...\n")
dc <- degree(g_clean)

# 2. 介数中心性 (BC)
cat("计算介数中心性 (BC)...\n")
bc <- betweenness(g_clean, directed = FALSE, normalized = TRUE)

# 3. 接近中心性 (CC)
cat("计算接近中心性 (CC)...\n")
cc <- closeness(g_clean, normalized = TRUE)

# 4. 特征向量中心性 (EC)
cat("计算特征向量中心性 (EC)...\n")
ec <- eigen_centrality(g_clean, directed = FALSE)$vector

# 合并所有中心性指标
centrality_df <- data.frame(
  Node = names(dc),
  DC = as.numeric(dc),
  BC = as.numeric(bc[names(dc)]),
  CC = as.numeric(cc[names(dc)]),
  EC = as.numeric(ec[names(dc)]),
  stringsAsFactors = FALSE
)

# 添加K-core信息
centrality_df$K_core <- cores[centrality_df$Node]

# 保存中心性结果
write.table(centrality_df, file = file.path(output_dir, "03_centrality_measures.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("中心性指标已保存\n")

# ==================== 步骤4：RRA综合评分分析 ====================
cat("\n========================================\n")
cat("步骤4: Robust Rank Aggregation (RRA) 综合评分\n")
cat("========================================\n")

# 为每种中心性指标创建排序列表（降序）
# 注意：RRA需要升序排名（排名1 = 最好），所以我们需要反转

# 创建排名（升序，1=最高值）
get_rank <- function(x) {
  return(rank(-x, ties.method = "min"))
}

centrality_df$DC_rank <- get_rank(centrality_df$DC)
centrality_df$BC_rank <- get_rank(centrality_df$BC)
centrality_df$CC_rank <- get_rank(centrality_df$CC)
centrality_df$EC_rank <- get_rank(centrality_df$EC)

# 计算平均排名（简单RRA近似）
centrality_df$Avg_Rank <- rowMeans(centrality_df[, c("DC_rank", "BC_rank", "CC_rank", "EC_rank")])

# 计算RRA综合评分（使用几何平均的变体）
# RRA_score = (DC_rank * BC_rank * CC_rank * EC_rank)^(1/4)
centrality_df$RRA_Score <- apply(centrality_df[, c("DC_rank", "BC_rank", "CC_rank", "EC_rank")], 1, 
                                 function(x) prod(x)^(1/4))

# 按RRA评分排序（越小越好）
centrality_df <- centrality_df[order(centrality_df$RRA_Score), ]

# 保存完整结果
write.table(centrality_df, file = file.path(output_dir, "04_rra_centrality_all.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 提取Top 20 Hub节点
top20 <- head(centrality_df, 20)
write.table(top20, file = file.path(output_dir, "04_top20_hub_nodes_rra.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nTop 20 Hub节点 (RRA综合评分):\n")
print(top20[, c("Node", "DC", "BC", "CC", "EC", "RRA_Score")])

# ==================== 步骤5：准备Cytoscape导入文件 ====================
cat("\n========================================\n")
cat("步骤5: 准备Cytoscape导入文件\n")
cat("========================================\n")

# 1. 节点属性文件（包含中心性指标）
node_attributes <- centrality_df[, c("Node", "DC", "BC", "CC", "EC", "K_core", "RRA_Score")]
write.table(node_attributes, file = file.path(output_dir, "05_cytoscape_node_attributes.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 2. 边列表文件
edge_list <- as.data.frame(get.edgelist(g_clean))
write.table(edge_list, file = file.path(output_dir, "05_cytoscape_edge_list.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("source", "target"))

# 3. 用于MCODE的节点列表（仅包含节点名）
write.table(data.frame(Node = V(g_clean)$name), 
            file = file.path(output_dir, "05_mcode_input_nodes.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 4. 用于cytoHubba的输入文件
cytohubba_df <- centrality_df[, c("Node", "DC", "BC", "CC", "EC")]
write.table(cytohubba_df, file = file.path(output_dir, "05_cytohubba_input.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# ==================== 步骤6：生成Cytoscape操作指南 ====================
cat("\n========================================\n")
cat("步骤6: Cytoscape操作指南\n")
cat("========================================\n")

guide <- c(
  "========================================",
  "    Cytoscape网络分析操作指南",
  "========================================",
  "",
  "【MCODE聚类分析】",
  "1. 打开Cytoscape软件",
  "2. 导入网络文件: 05_cytoscape_edge_list.txt",
  "   - File -> Import -> Network from File",
  "   - 选择05_cytoscape_edge_list.txt",
  "   - Source Column: source, Target Column: target",
  "",
  "3. 导入节点属性: 05_cytoscape_node_attributes.txt",
  "   - File -> Import -> Table from File",
  "   - 选择05_cytoscape_node_attributes.txt",
  "   - Key Column for Network: 选择Node列",
  "",
  "4. 运行MCODE插件:",
  "   - Apps -> MCODE",
  "   - 参数设置:",
  "     * Network Scoring: Degree Cutoff = 2, K-Core = 2",
  "     * Cluster Finding: Node Score Cutoff = 0.2, Haircut = TRUE",
  "   - 点击'Analyze Cluster'",
  "",
  "5. 导出MCODE结果:",
  "   - MCODE -> View/Export Clusters",
  "   - 保存为CSV或Cytoscape会话文件",
  "",
  "【cytoHubba中心性分析】",
  "1. 确保网络已导入",
  "2. 运行cytoHubba:",
  "   - Apps -> cytoHubba",
  "   - 选择要计算的指标:",
  "     * SC (Subgraph Centrality)",
  "     * NC (Network Centrality)",
  "     * LAC (Local Average Connectivity)",
  "     * MCC (Maximum Clique Centrality)",
  "   - 选择网络 -> 选择节点 -> 点击'Calculate'",
  "",
  "3. 导出cytoHubba结果:",
  "   - cytoHubba -> Export Results",
  "   - 保存为TXT或CSV格式",
  "",
  "【综合分析建议】",
  "1. 将cytoHubba结果与R计算的中心性指标合并",
  "2. 使用RRA方法整合所有8种中心性指标:",
  "   - DC (igraph)",
  "   - BC (igraph)",
  "   - CC (igraph)",
  "   - EC (igraph)",
  "   - SC (cytoHubba)",
  "   - NC (cytoHubba)",
  "   - LAC (cytoHubba)",
  "   - MCC (cytoHubba)",
  "3. 最终确定Top 20 Hub基因",
  "",
  "========================================"
)

writeLines(guide, file.path(output_dir, "06_cytoscape_guide.txt"))
cat("Cytoscape操作指南已保存\n")

# ==================== 最终总结 ====================
cat("\n========================================\n")
cat("         系统性分析完成\n")
cat("========================================\n")
cat(paste0("输出目录: ", output_dir, "\n\n"))
cat("生成文件列表:\n")
cat("【预处理】\n")
cat("  01_isolated_nodes_removed.txt - 删除的孤立节点\n\n")
cat("【K-core分解】\n")
cat("  02_kcore_decomposition.txt - K-core分解结果\n")
cat("  02_kcore_k3_nodes.txt - K>=3核心节点\n")
cat("  02_kcore_k3_edges.txt - K>=3核心边列表\n")
cat("  02_g_k3_network_object.rds - g_k3网络对象\n\n")
cat("【中心性分析】\n")
cat("  03_centrality_measures.txt - 4种中心性指标\n\n")
cat("【RRA综合评分】\n")
cat("  04_rra_centrality_all.txt - 所有节点的RRA评分\n")
cat("  04_top20_hub_nodes_rra.txt - Top 20 Hub节点\n\n")
cat("【Cytoscape输入文件】\n")
cat("  05_cytoscape_node_attributes.txt - 节点属性\n")
cat("  05_cytoscape_edge_list.txt - 边列表\n")
cat("  05_mcode_input_nodes.txt - MCODE输入节点\n")
cat("  05_cytohubba_input.txt - cytoHubba输入\n\n")
cat("【操作指南】\n")
cat("  06_cytoscape_guide.txt - Cytoscape操作指南\n")
cat("========================================\n")
