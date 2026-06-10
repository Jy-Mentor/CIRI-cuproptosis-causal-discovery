#!/usr/bin/env Rscript

# ============================================================================
# PPI网络拓扑分析脚本
# 图8: 66个核心交集基因的网络拓扑图
# 图9: 核心交集基因∪铜死亡核心基因的扩展网络拓扑图
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

# 二、数据读取与预处理 ----

# 读取所有STRING交互文件
file_paths <- c(
  "C:/Users/Jy-Mentor-7/Downloads/string_interactions (7).tsv",
  "C:/Users/Jy-Mentor-7/Downloads/string_interactions (6).tsv",
  "C:/Users/Jy-Mentor-7/Downloads/string_interactions_short (4).tsv",
  "C:/Users/Jy-Mentor-7/Downloads/string_interactions (5).tsv",
  "C:/Users/Jy-Mentor-7/Downloads/string_interactions_short (3).tsv"
)

# 读取并合并所有交互数据
all_edges <- lapply(file_paths, function(f) {
  if (file.exists(f)) {
    # 读取文件内容的第一行来获取列名
    con <- file(f, "r")
    header_line <- readLines(con, n = 1)
    close(con)
    
    # 处理列名：移除开头的#号并按制表符分割
    col_names <- strsplit(sub("^#", "", header_line), "\t")[[1]]
    
    # 读取数据，跳过第一行（因为那是注释/列名）
    df <- read.delim(f, stringsAsFactors = FALSE, comment.char = "#", header = FALSE, skip = 1)
    names(df) <- col_names
    
    return(df)
  } else {
    warning(paste("文件不存在:", f))
    return(NULL)
  }
}) %>% bind_rows()

# 去重：对于无向网络，确保node1-node2和node2-node1只保留一条
all_edges <- all_edges %>%
  mutate(
    pair_id = paste(pmin(node1, node2), pmax(node1, node2), sep = "_")
  ) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(-pair_id)

# 验证数据
stopifnot(nrow(all_edges) > 0)
cat("总交互边数（去重后）:", nrow(all_edges), "\n")

# 提取所有唯一基因
all_genes <- unique(c(all_edges$node1, all_edges$node2))
cat("网络中总基因数:", length(all_genes), "\n")

# 三、定义基因集合 ----

# 66个核心交集基因（从网络数据中提取的所有基因）
core_genes <- all_genes
cat("核心交集基因数:", length(core_genes), "\n")

# 15个铜死亡核心基因（从135个铜死亡相关基因中遴选）
cuproptosis_genes <- c("FDX1", "SLC31A1", "NFKB1", "RELA", "TNF", "IL6", 
                       "STAT3", "JAK1", "STAT1", "HIF1A", "NFE2L2", "HMOX1",
                       "PPARG", "MTOR", "EGFR")

# 检查铜死亡基因在网络中的存在情况
cuproptosis_in_network <- cuproptosis_genes[cuproptosis_genes %in% all_genes]
cuproptosis_not_in_network <- cuproptosis_genes[!(cuproptosis_genes %in% all_genes)]

cat("网络中已存在的铜死亡基因:", length(cuproptosis_in_network), "\n")
cat("网络中缺失的铜死亡基因:", paste(cuproptosis_not_in_network, collapse = ", "), "\n")

# 扩展基因集合：核心基因 ∪ 铜死亡基因
extended_genes <- unique(c(core_genes, cuproptosis_genes))
cat("扩展网络基因数:", length(extended_genes), "\n")

# 四、构建图8网络（核心交集基因） ----

# 筛选核心基因之间的交互
edges_core <- all_edges %>%
  filter(node1 %in% core_genes, node2 %in% core_genes)

# 创建igraph对象
g_core <- graph_from_data_frame(
  d = edges_core %>% select(node1, node2, combined_score),
  directed = FALSE,
  vertices = data.frame(name = core_genes)
)

# 移除孤立节点（如果有）
g_core <- delete_vertices(g_core, which(degree(g_core) == 0))

# 计算图8网络拓扑指标
cat("\n========== 图8: 核心交集基因网络拓扑特征 ==========\n")
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

# 五、构建图9扩展网络（核心基因 ∪ 铜死亡基因） ----

# 对于缺失的铜死亡基因，我们需要从原始数据中获取其交互
# 由于缺失基因不在当前网络中，我们标记它们并分析其与现有网络的连接

# 筛选扩展基因集合之间的交互（包括现有边）
edges_extended <- all_edges %>%
  filter(node1 %in% extended_genes, node2 %in% extended_genes)

# 创建扩展网络
g_extended <- graph_from_data_frame(
  d = edges_extended %>% select(node1, node2, combined_score),
  directed = FALSE,
  vertices = data.frame(name = extended_genes[extended_genes %in% unique(c(edges_extended$node1, edges_extended$node2))])
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

# 六、铜死亡基因与交集基因的拓扑交叉分析 ----

# 识别扩展网络中的铜死亡基因
cuproptosis_in_extended <- cuproptosis_genes[cuproptosis_genes %in% V(g_extended)$name]
cat("\n扩展网络中的铜死亡基因数:", length(cuproptosis_in_extended), "\n")
cat("铜死亡基因列表:", paste(cuproptosis_in_extended, collapse = ", "), "\n")

# 计算铜死亡基因与核心基因之间的直接互作边
cuproptosis_edges <- edges_extended %>%
  filter(
    (node1 %in% cuproptosis_in_extended & node2 %in% core_genes) |
    (node2 %in% cuproptosis_in_extended & node1 %in% core_genes)
  )

cat("\n铜死亡基因与核心基因之间的直接互作边数:", nrow(cuproptosis_edges), "\n")
cat("占总边数比例:", round(nrow(cuproptosis_edges) / ecount(g_extended) * 100, 2), "%\n")

# 计算扩展网络各节点拓扑指标
extended_node_metrics <- data.frame(
  gene = V(g_extended)$name,
  degree = degree(g_extended),
  betweenness = round(betweenness(g_extended, normalized = TRUE), 4),
  closeness = round(closeness(g_extended), 4),
  eigenvector = round(eigen_centrality(g_extended)$vector, 4),
  is_cuproptosis = V(g_extended)$name %in% cuproptosis_genes,
  stringsAsFactors = FALSE
) %>% arrange(desc(degree))

cat("\n扩展网络中心性排名前15的节点:\n")
print(head(extended_node_metrics, 15))

# 七、识别关键枢纽基因 ----

# 找出连接度最高的桥接基因
top_hub <- extended_node_metrics$gene[1]
top_hub_degree <- extended_node_metrics$degree[1]

cat("\n========== 关键枢纽基因分析 ==========\n")
cat("最高连接度枢纽基因:", top_hub, "(Degree =", top_hub_degree, ")\n")

# 分析铜死亡核心执行蛋白（如FDX1）的连接情况
if ("FDX1" %in% extended_node_metrics$gene) {
  fdx1_info <- extended_node_metrics %>% filter(gene == "FDX1")
  cat("\nFDX1 (铜死亡核心执行蛋白) 连接度:", fdx1_info$degree, "\n")
  
  # 找出与FDX1直接互作的基因
  fdx1_neighbors <- neighbors(g_extended, "FDX1")$name
  cat("与FDX1直接互作的基因:", paste(fdx1_neighbors, collapse = ", "), "\n")
}

# 分析NFKB1的连接情况（炎症信号关键基因）
if ("NFKB1" %in% extended_node_metrics$gene) {
  nfkb1_info <- extended_node_metrics %>% filter(gene == "NFKB1")
  cat("\nNFKB1 连接度:", nfkb1_info$degree, "\n")
  
  nfkb1_neighbors <- neighbors(g_extended, "NFKB1")$name
  cat("与NFKB1直接互作的基因:", paste(nfkb1_neighbors, collapse = ", "), "\n")
}

# 八、可视化 ----

# 设置输出目录
output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/PPI_network_plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 图8: 核心交集基因网络拓扑图
png(file.path(output_dir, "Figure8_core_PPI_network.png"), 
    width = 1200, height = 1000, res = 150)

set.seed(42)
layout_core <- layout_with_fr(g_core, niter = 1000)

# 根据度设置节点大小
node_size_core <- degree(g_core)
node_size_core <- rescale(node_size_core, to = c(3, 15))

# 绘制网络
plot(g_core,
     layout = layout_core,
     vertex.size = node_size_core,
     vertex.color = "#E74C3C",
     vertex.frame.color = "#C0392B",
     vertex.label.color = "black",
     vertex.label.cex = 0.6,
     vertex.label.dist = 1,
     edge.color = alpha("#7F8C8D", 0.5),
     edge.width = 0.5,
     main = "图8 核心交集基因PPI网络拓扑图")

dev.off()

# 图9: 扩展网络拓扑图
png(file.path(output_dir, "Figure9_extended_PPI_network.png"), 
    width = 1400, height = 1200, res = 150)

set.seed(42)
layout_extended <- layout_with_fr(g_extended, niter = 1000)

# 设置节点颜色：铜死亡基因为金色，核心基因为红色
node_colors <- ifelse(V(g_extended)$name %in% cuproptosis_genes, 
                      "#F39C12", "#E74C3C")
node_frame_colors <- ifelse(V(g_extended)$name %in% cuproptosis_genes,
                            "#D68910", "#C0392B")

node_size_extended <- degree(g_extended)
node_size_extended <- rescale(node_size_extended, to = c(3, 18))

plot(g_extended,
     layout = layout_extended,
     vertex.size = node_size_extended,
     vertex.color = node_colors,
     vertex.frame.color = node_frame_colors,
     vertex.label.color = "black",
     vertex.label.cex = 0.55,
     vertex.label.dist = 1,
     edge.color = alpha("#7F8C8D", 0.4),
     edge.width = 0.4,
     main = "图9 核心交集基因∪铜死亡核心基因PPI网络拓扑图")

# 添加图例
legend("topright", 
       legend = c("核心交集基因", "铜死亡核心基因"),
       col = c("#E74C3C", "#F39C12"),
       pch = 21,
       pt.bg = c("#E74C3C", "#F39C12"),
       cex = 0.8)

dev.off()

# 九、生成网络统计摘要 ----

cat("\n========== 网络统计摘要 ==========\n")

cat("\n【图8】核心交集基因网络:\n")
cat("  - 有效互作边数:", ecount(g_core), "\n")
cat("  - 平均节点度:", round(mean(degree(g_core)), 2), "\n")
cat("  - 网络密度:", round(edge_density(g_core), 4), "\n")

cat("\n【图9】扩展网络:\n")
cat("  - 有效互作边数:", ecount(g_extended), "\n")
cat("  - 平均节点度:", round(mean(degree(g_extended)), 2), "\n")
cat("  - 网络密度:", round(edge_density(g_extended), 4), "\n")
cat("  - 铜死亡基因与交集基因直接互作边数:", nrow(cuproptosis_edges), "\n")
cat("  - 占总边数比例:", round(nrow(cuproptosis_edges) / ecount(g_extended) * 100, 2), "%\n")

cat("\n【关键发现】\n")
cat("  - 桥接枢纽基因:", top_hub, "(Degree =", top_hub_degree, ")\n")

# 十、保存结果 ----

# 保存节点拓扑指标
write.csv(core_node_metrics, 
          file.path(output_dir, "Figure8_core_node_metrics.csv"),
          row.names = FALSE)

write.csv(extended_node_metrics,
          file.path(output_dir, "Figure9_extended_node_metrics.csv"),
          row.names = FALSE)

# 保存铜死亡基因互作边
write.csv(cuproptosis_edges,
          file.path(output_dir, "cuproptosis_interaction_edges.csv"),
          row.names = FALSE)

cat("\n所有结果已保存到:", output_dir, "\n")
cat("分析完成!\n")
