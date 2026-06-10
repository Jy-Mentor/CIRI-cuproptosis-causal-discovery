#!/usr/bin/env Rscript

# ============================================================================
# PPI网络拓扑分析脚本 - 修正版本
# 图8: 66个核心交集基因的网络拓扑图 (string_interactions_short 6)
# 图9: 核心交集基因∪铜死亡核心基因的扩展网络拓扑图 (string_interactions_short 5)
# ============================================================================

# 一、环境设置与包加载 ----
set.seed(42)

# 智能包安装与加载
packages <- c("igraph", "dplyr", "ggplot2", "ggraph", "tidygraph", 
              "RColorBrewer", "scales", "gridExtra", "patchwork")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# 二、数据读取函数 ----

read_string_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("文件不存在:", file_path))
  }
  
  # 读取文件内容的第一行来获取列名
  con <- file(file_path, "r")
  header_line <- readLines(con, n = 1)
  close(con)
  
  # 处理列名：移除开头的#号并按制表符分割
  col_names <- strsplit(sub("^#", "", header_line), "\t")[[1]]
  
  # 读取数据，跳过第一行（因为那是注释/列名）
  df <- read.delim(file_path, stringsAsFactors = FALSE, comment.char = "#", header = FALSE, skip = 1)
  names(df) <- col_names
  
  return(df)
}

# 三、读取两个网络数据 ----

# 图8: 66个核心基因网络 (文件6)
edges_core <- read_string_file("C:/Users/Jy-Mentor-7/Downloads/string_interactions_short (6).tsv")

# 图9: 扩展网络 (文件5) - 并入铜死亡基因
edges_extended <- read_string_file("C:/Users/Jy-Mentor-7/Downloads/string_interactions_short (5).tsv")

# 去重处理
edges_core <- edges_core %>%
  mutate(pair_id = paste(pmin(node1, node2), pmax(node1, node2), sep = "_")) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(-pair_id)

edges_extended <- edges_extended %>%
  mutate(pair_id = paste(pmin(node1, node2), pmax(node1, node2), sep = "_")) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(-pair_id)

# 四、定义基因集合 ----

# 66个核心基因（从图8网络提取）
core_genes <- unique(c(edges_core$node1, edges_core$node2))
cat("核心交集基因数:", length(core_genes), "\n")

# 15个铜死亡核心基因
cuproptosis_genes <- c("FDX1", "SLC31A1", "NFKB1", "RELA", "TNF", "IL6", 
                       "STAT3", "JAK1", "STAT1", "HIF1A", "NFE2L2", "HMOX1",
                       "PPARG", "MTOR", "EGFR")

# 检查铜死亡基因在两个网络中的存在情况
cuproptosis_in_core <- cuproptosis_genes[cuproptosis_genes %in% core_genes]
cuproptosis_not_in_core <- cuproptosis_genes[!(cuproptosis_genes %in% core_genes)]

cat("核心网络中已存在的铜死亡基因:", length(cuproptosis_in_core), "\n")
cat("核心网络中缺失的铜死亡基因:", paste(cuproptosis_not_in_core, collapse = ", "), "\n")

# 扩展网络基因
extended_genes <- unique(c(edges_extended$node1, edges_extended$node2))
cat("扩展网络基因数:", length(extended_genes), "\n")

# 五、构建图8网络（66个核心交集基因） ----

g_core <- graph_from_data_frame(
  d = edges_core %>% select(node1, node2, combined_score),
  directed = FALSE,
  vertices = data.frame(name = core_genes)
)

# 移除孤立节点
g_core <- delete_vertices(g_core, which(degree(g_core) == 0))

# 计算图8网络拓扑指标
cat("\n========== 图8: 66个核心交集基因网络拓扑特征 ==========\n")
cat("网络节点数:", vcount(g_core), "\n")
cat("网络边数:", ecount(g_core), "\n")
cat("网络密度:", round(edge_density(g_core), 4), "\n")
cat("平均节点度:", round(mean(degree(g_core)), 2), "\n")
cat("平均路径长度:", round(mean_distance(g_core), 2), "\n")
cat("网络直径:", diameter(g_core), "\n")
cat("聚类系数:", round(transitivity(g_core, type = "global"), 4), "\n")

# 计算各节点拓扑指标
core_node_metrics <- data.frame(
  gene = V(g_core)$name,
  degree = degree(g_core),
  betweenness = round(betweenness(g_core, normalized = TRUE), 4),
  closeness = round(closeness(g_core), 4),
  eigenvector = round(eigen_centrality(g_core)$vector, 4),
  stringsAsFactors = FALSE
) %>% arrange(desc(degree))

cat("\n图8网络中心性排名前10的节点:\n")
print(head(core_node_metrics, 10))

# 六、构建图9扩展网络 ----

g_extended <- graph_from_data_frame(
  d = edges_extended %>% select(node1, node2, combined_score),
  directed = FALSE,
  vertices = data.frame(name = extended_genes)
)

# 移除孤立节点
g_extended <- delete_vertices(g_extended, which(degree(g_extended) == 0))

# 计算图9扩展网络拓扑指标
cat("\n========== 图9: 扩展网络拓扑特征 ==========\n")
cat("扩展网络节点数:", vcount(g_extended), "\n")
cat("扩展网络边数:", ecount(g_extended), "\n")
cat("扩展网络密度:", round(edge_density(g_extended), 4), "\n")
cat("扩展网络平均节点度:", round(mean(degree(g_extended)), 2), "\n")
cat("扩展网络平均路径长度:", round(mean_distance(g_extended), 2), "\n")
cat("扩展网络直径:", diameter(g_extended), "\n")
cat("扩展网络聚类系数:", round(transitivity(g_extended, type = "global"), 4), "\n")

# 七、铜死亡基因与核心基因的拓扑交叉分析 ----

# 识别扩展网络中的铜死亡基因
cuproptosis_in_extended <- cuproptosis_genes[cuproptosis_genes %in% V(g_extended)$name]
cat("\n扩展网络中的铜死亡基因数:", length(cuproptosis_in_extended), "\n")

# 计算铜死亡基因与核心基因之间的直接互作边
# 注意：这里核心基因是66个，不是扩展网络中的所有基因
cuproptosis_edges <- edges_extended %>%
  filter(
    (node1 %in% cuproptosis_in_extended & node2 %in% core_genes) |
    (node2 %in% cuproptosis_in_extended & node1 %in% core_genes)
  )

cat("铜死亡基因与66个核心基因之间的直接互作边数:", nrow(cuproptosis_edges), "\n")
cat("占扩展网络总边数比例:", round(nrow(cuproptosis_edges) / ecount(g_extended) * 100, 2), "%\n")

# 计算扩展网络各节点拓扑指标
extended_node_metrics <- data.frame(
  gene = V(g_extended)$name,
  degree = degree(g_extended),
  betweenness = round(betweenness(g_extended, normalized = TRUE), 4),
  closeness = round(closeness(g_extended), 4),
  eigenvector = round(eigen_centrality(g_extended)$vector, 4),
  is_cuproptosis = V(g_extended)$name %in% cuproptosis_genes,
  is_core = V(g_extended)$name %in% core_genes,
  stringsAsFactors = FALSE
) %>% arrange(desc(degree))

cat("\n扩展网络中心性排名前15的节点:\n")
print(head(extended_node_metrics, 15))

# 八、识别关键枢纽基因 ----

# 找出连接度最高的桥接基因
top_hub <- extended_node_metrics$gene[1]
top_hub_degree <- extended_node_metrics$degree[1]

cat("\n========== 关键枢纽基因分析 ==========\n")
cat("最高连接度枢纽基因:", top_hub, "(Degree =", top_hub_degree, ")\n")

# 分析铜死亡核心执行蛋白FDX1
if ("FDX1" %in% extended_node_metrics$gene) {
  fdx1_info <- extended_node_metrics %>% filter(gene == "FDX1")
  cat("\nFDX1 (铜死亡核心执行蛋白) 连接度:", fdx1_info$degree, "\n")
  fdx1_neighbors <- neighbors(g_extended, "FDX1")$name
  cat("与FDX1直接互作的基因:", paste(fdx1_neighbors, collapse = ", "), "\n")
}

# 分析NFKB1
if ("NFKB1" %in% extended_node_metrics$gene) {
  nfkb1_info <- extended_node_metrics %>% filter(gene == "NFKB1")
  cat("\nNFKB1 连接度:", nfkb1_info$degree, "\n")
}

# 九、可视化 ----

output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/PPI_network_plots_corrected"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 图8: 核心交集基因网络拓扑图
png(file.path(output_dir, "Figure8_core66_PPI_network.png"), 
    width = 1600, height = 1400, res = 150)

set.seed(42)
layout_core <- layout_with_fr(g_core, niter = 2000)

node_degrees_core <- degree(g_core)
node_size_core <- rescale(node_degrees_core, to = c(5, 20))
degree_colors_core <- colorRampPalette(c("#FADBD8", "#E74C3C", "#922B21"))(max(node_degrees_core))
node_color_core <- degree_colors_core[node_degrees_core]

par(mar = c(0, 0, 2, 0))
plot(g_core,
     layout = layout_core,
     vertex.size = node_size_core,
     vertex.color = node_color_core,
     vertex.frame.color = "#641E16",
     vertex.label.color = "black",
     vertex.label.cex = 0.7,
     vertex.label.dist = 1.2,
     vertex.label.font = 2,
     edge.color = alpha("#7F8C8D", 0.4),
     edge.width = 0.6,
     main = "图8 66个核心交集基因PPI网络拓扑图")

dev.off()

# 图9: 扩展网络拓扑图
png(file.path(output_dir, "Figure9_extended_PPI_network.png"), 
    width = 1800, height = 1600, res = 150)

set.seed(42)
layout_extended <- layout_with_fr(g_extended, niter = 2000)

node_degrees_extended <- degree(g_extended)
node_size_extended <- rescale(node_degrees_extended, to = c(5, 22))

# 设置节点颜色：铜死亡基因为金色，核心基因为红色，其他为灰色
node_colors_extended <- ifelse(V(g_extended)$name %in% cuproptosis_genes, "#F39C12",
                               ifelse(V(g_extended)$name %in% core_genes, "#E74C3C", "#95A5A6"))
node_frame_colors_extended <- ifelse(V(g_extended)$name %in% cuproptosis_genes, "#D68910",
                                      ifelse(V(g_extended)$name %in% core_genes, "#C0392B", "#7F8C8D"))

par(mar = c(0, 0, 2, 0))
plot(g_extended,
     layout = layout_extended,
     vertex.size = node_size_extended,
     vertex.color = node_colors_extended,
     vertex.frame.color = node_frame_colors_extended,
     vertex.label.color = "black",
     vertex.label.cex = 0.65,
     vertex.label.dist = 1.2,
     vertex.label.font = 2,
     edge.color = alpha("#7F8C8D", 0.35),
     edge.width = 0.5,
     main = "图9 核心交集基因∪铜死亡核心基因PPI网络拓扑图")

legend("bottomleft", 
       legend = c("核心交集基因", "铜死亡核心基因", "新增基因"),
       col = c("#E74C3C", "#F39C12", "#95A5A6"),
       pch = 21,
       pt.bg = c("#E74C3C", "#F39C12", "#95A5A6"),
       pt.cex = c(1.5, 1.5, 1.5),
       cex = 0.9,
       title = "基因类型")

dev.off()

# 十、保存结果 ----

write.csv(core_node_metrics, 
          file.path(output_dir, "Figure8_core66_node_metrics.csv"),
          row.names = FALSE)

write.csv(extended_node_metrics,
          file.path(output_dir, "Figure9_extended_node_metrics.csv"),
          row.names = FALSE)

write.csv(cuproptosis_edges,
          file.path(output_dir, "cuproptosis_interaction_edges.csv"),
          row.names = FALSE)

# 保存网络统计摘要
sink(file.path(output_dir, "network_summary.txt"))
cat("========== PPI网络拓扑分析摘要 ==========\n\n")
cat("【图8】66个核心交集基因网络 (n=", vcount(g_core), "):\n")
cat("  - 有效互作边数:", ecount(g_core), "\n")
cat("  - 平均节点度:", round(mean(degree(g_core)), 2), "\n")
cat("  - 网络密度:", round(edge_density(g_core), 4), "\n")
cat("  - 平均路径长度:", round(mean_distance(g_core), 2), "\n")
cat("  - 网络直径:", diameter(g_core), "\n")
cat("  - 聚类系数:", round(transitivity(g_core, type = "global"), 4), "\n\n")

cat("【图9】扩展网络 (n=", vcount(g_extended), "):\n")
cat("  - 有效互作边数:", ecount(g_extended), "\n")
cat("  - 平均节点度:", round(mean(degree(g_extended)), 2), "\n")
cat("  - 网络密度:", round(edge_density(g_extended), 4), "\n")
cat("  - 平均路径长度:", round(mean_distance(g_extended), 2), "\n")
cat("  - 网络直径:", diameter(g_extended), "\n")
cat("  - 聚类系数:", round(transitivity(g_extended, type = "global"), 4), "\n")
cat("  - 铜死亡基因与66个核心基因直接互作边数:", nrow(cuproptosis_edges), "\n")
cat("  - 占扩展网络总边数比例:", round(nrow(cuproptosis_edges) / ecount(g_extended) * 100, 2), "%\n\n")

cat("【关键枢纽基因】\n")
cat("  - 最高连接度基因:", top_hub, "(Degree =", top_hub_degree, ")\n")
if ("FDX1" %in% extended_node_metrics$gene) {
  fdx1_deg <- extended_node_metrics$degree[extended_node_metrics$gene == "FDX1"]
  cat("  - FDX1 (铜死亡核心执行蛋白) Degree =", fdx1_deg, "\n")
}
if ("NFKB1" %in% extended_node_metrics$gene) {
  nfkb1_deg <- extended_node_metrics$degree[extended_node_metrics$gene == "NFKB1"]
  cat("  - NFKB1 (炎症信号关键基因) Degree =", nfkb1_deg, "\n")
}
sink()

cat("\n所有结果已保存到:", output_dir, "\n")
cat("分析完成!\n")
