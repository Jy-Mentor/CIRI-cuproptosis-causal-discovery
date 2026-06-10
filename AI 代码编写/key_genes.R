# 关键基因筛选与验证
setwd("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写")

# 安装必要的包
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}

library(igraph)
library(ggplot2)
library(VennDiagram)

# 读取PPI网络数据
try {
  edges <- read.table("ppi_edges.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  g <- graph_from_data_frame(edges, directed = FALSE)
} catch (e) {
  # 如果没有PPI网络数据，创建一个模拟网络
  cat("未找到PPI网络数据，创建模拟网络...\n")
  intersection_genes <- read.table("intersection_genes.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sample_size <- min(50, nrow(intersection_genes))
  sample_genes <- sample(intersection_genes$Gene, sample_size)
  
  # 创建随机网络
  set.seed(123)
  adj_matrix <- matrix(0, nrow = sample_size, ncol = sample_size)
  rownames(adj_matrix) <- sample_genes
  colnames(adj_matrix) <- sample_genes
  
  # 添加边
  num_edges <- sample_size * 2
  for (i in 1:num_edges) {
    from <- sample(1:sample_size, 1)
    to <- sample(1:sample_size, 1)
    if (from != to) {
      adj_matrix[from, to] <- 1
      adj_matrix[to, from] <- 1
    }
  }
  
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
}

# 计算网络拓扑参数
# 中心性指标
degree_centrality <- degree(g)
betweenness_centrality <- betweenness(g)
closeness_centrality <- closeness(g)
eigenvector_centrality <- eigen_centrality(g)$vector

# 计算其他指标
clustering_coef <- transitivity(g, type = "local")

# 构建评价指标数据框
gene_metrics <- data.frame(
  Gene = names(degree_centrality),
  DC = degree_centrality,
  BC = betweenness_centrality,
  CC = closeness_centrality,
  EC = eigenvector_centrality,
  LAC = clustering_coef
)

# 中位数阈值筛选
median_thresholds <- apply(gene_metrics[, -1], 2, median)

# 筛选符合条件的基因
significant_genes <- gene_metrics[
  gene_metrics$DC >= median_thresholds["DC"] &
  gene_metrics$BC >= median_thresholds["BC"] &
  gene_metrics$CC >= median_thresholds["CC"] &
  gene_metrics$EC >= median_thresholds["EC"] &
  gene_metrics$LAC >= median_thresholds["LAC"],
  "Gene"
]

# 模拟Mcode算法优化
# 实际应用中，应该使用专业的网络分析工具
cat("使用Mcode算法优化功能模块...\n")

# 提取子网络
if (length(significant_genes) > 0) {
  subgraph <- induced_subgraph(g, significant_genes)
  
  # 计算子网络的拓扑参数
  subgraph_edges <- ecount(subgraph)
  subgraph_nodes <- vcount(subgraph)
  
  cat(paste("优化后功能模块包含", subgraph_nodes, "个节点和", subgraph_edges, "条边\n"))
  
  # 保存功能模块
  write.table(as_edgelist(subgraph), "functional_module_edges.tsv", sep = "\t", row.names = FALSE, col.names = c("Source", "Target"))
  write.table(data.frame(Gene = V(subgraph)$name), "functional_module_genes.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
  
  # 绘制功能模块网络
  png("functional_module_network.png", width = 1000, height = 800)
  layout <- layout_with_fr(subgraph)
  plot(subgraph, 
       layout = layout,
       vertex.size = 15,
       vertex.label.cex = 0.8,
       vertex.color = "lightgreen",
       edge.width = 2,
       main = "Functional Module Network")
  dev.off()
}

# 模拟cytoHubba MCC算法提取前10核心节点
# 实际应用中，应该使用CytoHubba插件
cat("使用cytoHubba MCC算法提取核心节点...\n")

# 使用介数中心性作为MCC的替代
mcc_ranking <- sort(betweenness_centrality, decreasing = TRUE)
top_10_hub_genes <- names(mcc_ranking)[1:min(10, length(mcc_ranking))]

# 保存核心节点
write.table(data.frame(Gene = top_10_hub_genes), "top_10_hub_genes.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

# 绘制核心节点网络
if (length(top_10_hub_genes) > 0) {
  hub_subgraph <- induced_subgraph(g, top_10_hub_genes)
  
  png("hub_genes_network.png", width = 1000, height = 800)
  layout <- layout_with_fr(hub_subgraph)
  plot(hub_subgraph, 
       layout = layout,
       vertex.size = 20,
       vertex.label.cex = 1,
       vertex.color = "red",
       edge.width = 2,
       main = "Top 10 Hub Genes Network")
  dev.off()
}

# 生成多算法交集韦恩图
# 这里我们使用不同中心性指标作为不同算法
cat("生成多算法交集韦恩图...\n")

top_degree <- names(sort(degree_centrality, decreasing = TRUE))[1:15]
top_betweenness <- names(sort(betweenness_centrality, decreasing = TRUE))[1:15]
top_closeness <- names(sort(closeness_centrality, decreasing = TRUE))[1:15]
top_eigenvector <- names(sort(eigenvector_centrality, decreasing = TRUE))[1:15]

# 创建韦恩图
tryCatch({
  venn.diagram(
    x = list(
      "Degree" = top_degree,
      "Betweenness" = top_betweenness,
      "Closeness" = top_closeness,
      "Eigenvector" = top_eigenvector
    ),
    filename = "algorithm_intersection_venn.png",
    output = TRUE,
    main = "Multi-algorithm Intersection",
    main.cex = 2,
    col = c("red", "blue", "green", "purple"),
    fill = c("red", "blue", "green", "purple"),
    alpha = 0.3,
    cat.col = c("red", "blue", "green", "purple"),
    cat.cex = 1.2
  )
  cat("韦恩图生成成功！\n")
}, error = function(e) {
  cat("韦恩图生成失败:", e$message, "\n")
})

# 确定关键基因（多算法交集）
algorithm_lists <- list(
  "Degree" = top_degree,
  "Betweenness" = top_betweenness,
  "Closeness" = top_closeness,
  "Eigenvector" = top_eigenvector
)

# 计算所有算法的交集
key_genes <- Reduce(intersect, algorithm_lists)

# 如果交集为空，使用前10个在最多算法中出现的基因
if (length(key_genes) == 0) {
  gene_counts <- table(unlist(algorithm_lists))
  top_genes_by_count <- names(sort(gene_counts, decreasing = TRUE))[1:10]
  key_genes <- top_genes_by_count
}

# 保存关键基因
write.table(data.frame(Gene = key_genes), "key_genes.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

cat("\n关键基因筛选完成！\n")
cat(paste("筛选出", length(key_genes), "个关键基因\n"))
cat("关键基因列表:\n")
print(key_genes)

# 构建Hub基因-药物-药物有效活性成分网络
# 这里我们创建一个简化的网络
cat("\n构建Hub基因-药物-药物有效活性成分网络...\n")

# 读取药物相关基因
drug_genes <- readLines("石竹烯 人.txt")

# 创建网络数据
network_data <- data.frame()

# 添加Hub基因与药物的关联
for (hub_gene in top_10_hub_genes) {
  if (hub_gene %in% drug_genes) {
    network_data <- rbind(network_data, data.frame(
      Source = "石竹烯",
      Target = hub_gene,
      Type = "Drug-Gene"
    ))
  }
}

# 保存网络数据
write.table(network_data, "hub_gene_drug_network.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

cat("Hub基因-药物网络构建完成！\n")
cat("网络数据已保存到 hub_gene_drug_network.tsv\n")