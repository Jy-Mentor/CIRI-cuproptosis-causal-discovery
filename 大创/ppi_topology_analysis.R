# PPI 网络拓扑分析脚本
# 分两个阶段分析：基础网络（78基因）和扩展网络（78+15基因）

# 设置随机种子保证可重复性
set.seed(123)

# 智能检查并安装包
install.packages.if.needed <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# 安装必要的包
install.packages.if.needed(c("igraph"))

# 加载包
library(igraph)

# 定义文件路径
base_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/78ppi.tsv"
extended_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/78+15 ppi.tsv"

# 定义铜死亡基因列表
cuproptosis_genes <- c("FDX1", "SLC31A1", "NFKB1", "DLAT", "DLST", "PDHA1", "PDHB", "LIAS", "LIPT1", "LIPT2", "GLS", "CDKN2A", "GLS2", "DLD", "ATP7B")

# 数据处理函数
process_ppi_data <- function(file_path) {
  # 读取数据，确保正确解析表头
  # 使用 read.delim 并设置 comment.char = "" 以避免 # 被视为注释
  data <- read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
  
  # 打印列名以调试
  cat("文件列名:", colnames(data), "\n")
  
  # 自动识别节点列
  if ("#node1" %in% colnames(data)) {
    colnames(data)[colnames(data) == "#node1"] <- "node1"
  }
  
  # 确保 combined_score 是数值型
  if ("combined_score" %in% colnames(data)) {
    # 尝试转换为数值型
    data$combined_score <- as.numeric(as.character(data$combined_score))
    
    # 打印转换结果
    cat("combined_score 列类型:", class(data$combined_score), "\n")
    cat("非NA值数量:", sum(!is.na(data$combined_score)), "\n")
    cat("大于等于0.4的值数量:", sum(data$combined_score >= 0.4, na.rm = TRUE), "\n")
    
    # 过滤 combined_score < 0.4
    data <- data[!is.na(data$combined_score) & data$combined_score >= 0.4, ]
  } else {
    cat("警告: 文件中没有 combined_score 列!\n")
  }
  
  # 打印过滤后的数据量
  cat("过滤后的数据行数:", nrow(data), "\n")
  
  # 统一转为大写 Symbol
  if ("node1" %in% colnames(data) && "node2" %in% colnames(data)) {
    data$node1 <- toupper(data$node1)
    data$node2 <- toupper(data$node2)
    
    # 去除自环
    data <- data[data$node1 != data$node2, ]
    
    # 打印去除自环后的数据量
    cat("去除自环后的数据行数:", nrow(data), "\n")
    
    # 去除重复边（无向图）
    # 先排序每行的节点，然后去重
    if (nrow(data) > 0) {
      edge_pairs <- t(apply(data[, c("node1", "node2")], 1, sort))
      data <- data[!duplicated(edge_pairs), ]
      
      # 打印去除重复边后的数据量
      cat("去除重复边后的数据行数:", nrow(data), "\n")
    }
  }
  
  return(data)
}

# 计算网络拓扑参数
calculate_topology <- function(graph) {
  # 计算有效互作边总数
  edges <- ecount(graph)
  
  # 计算平均节点度
  avg_degree <- mean(degree(graph))
  
  # 计算网络密度
  density <- edge_density(graph)
  
  # 计算平均路径长度（若图不连通则计算最大连通子图）
  if (is_connected(graph)) {
    apl <- mean_distance(graph)
  } else {
    comps <- components(graph)
    largest_comp_id <- which.max(comps$csize)
    largest_comp <- induced_subgraph(graph, which(comps$membership == largest_comp_id))
    apl <- mean_distance(largest_comp)
  }
  
  # 计算平均聚类系数
  clustering <- transitivity(graph, type = "average")
  
  return(list(
    edges = edges,
    avg_degree = avg_degree,
    density = density,
    apl = apl,
    clustering = clustering
  ))
}

# 识别 Hub 基因
identify_hub_genes <- function(graph, top_n = 10) {
  # 计算三种中心性
  degree_centrality <- degree(graph)
  betweenness_centrality <- betweenness(graph)
  closeness_centrality <- closeness(graph)
  
  # 按 Degree 降序排序
  sorted_genes <- names(sort(degree_centrality, decreasing = TRUE))[1:top_n]
  
  # 构建结果数据框
  hub_df <- data.frame(
    Gene = sorted_genes,
    Degree = degree_centrality[sorted_genes],
    Betweenness = betweenness_centrality[sorted_genes],
    Closeness = closeness_centrality[sorted_genes]
  )
  
  return(hub_df)
}

# 阶段一：基础网络（78 基因）拓扑解析
cat("阶段一：基础网络拓扑解析\n")

# 处理基础网络数据
base_data <- process_ppi_data(base_file)

# 构建基础网络
base_graph <- graph_from_data_frame(base_data, directed = FALSE)

# 计算基础网络拓扑参数
base_topology <- calculate_topology(base_graph)

# 识别基础网络 Hub 基因
base_hub_genes <- identify_hub_genes(base_graph)

# 阶段二：扩展网络（78+15 基因）交叉拓扑解析
cat("\n阶段二：扩展网络拓扑解析\n")

# 处理扩展网络数据
extended_data <- process_ppi_data(extended_file)

# 构建扩展网络
extended_graph <- graph_from_data_frame(extended_data, directed = FALSE)

# 计算扩展网络拓扑参数
extended_topology <- calculate_topology(extended_graph)

# 识别扩展网络 Hub 基因
extended_hub_genes <- identify_hub_genes(extended_graph)

# 铜死亡基因专项统计
cat("\n铜死亡基因专项统计\n")

# 提取铜死亡基因子图
cuproptosis_subgraph <- induced_subgraph(extended_graph, cuproptosis_genes[cuproptosis_genes %in% V(extended_graph)$name])

# 计算铜死亡基因在全局网络中的 Degree 分布
cuproptosis_degree <- degree(extended_graph, v = cuproptosis_genes[cuproptosis_genes %in% V(extended_graph)$name])

# 铜死亡基因 Degree Top 3
cuproptosis_top3 <- sort(cuproptosis_degree, decreasing = TRUE)[1:3]

# 跨组互作量化：计算铜死亡基因与 78 个交集基因之间的直接连接边数
# 假设 78 个交集基因为基础网络的节点
intersection_genes <- V(base_graph)$name
cross_edges <- 0

# 获取扩展网络中的所有节点
extended_nodes <- V(extended_graph)$name

for (gene1 in cuproptosis_genes) {
  # 只检查存在于扩展网络中的铜死亡基因
  if (gene1 %in% extended_nodes) {
    for (gene2 in intersection_genes) {
      # 只检查存在于扩展网络中的交集基因
      if (gene2 %in% extended_nodes) {
        if (are_adjacent(extended_graph, gene1, gene2)) {
          cross_edges <- cross_edges + 1
        }
      }
    }
  }
}

# 统计铜死亡基因之间的内部连接边数
internal_edges <- ecount(cuproptosis_subgraph)

# 计算铜死亡基因与 78 交集基因的连接数
cuproptosis_connected_count <- sapply(cuproptosis_genes, function(gene) {
  if (gene %in% V(extended_graph)$name) {
    neighbors <- neighbors(extended_graph, gene)
    sum(neighbors$name %in% intersection_genes)
  } else {
    0
  }
})

# 生成铜死亡基因中心性数据
cuproptosis_centrality <- data.frame(
  Gene = cuproptosis_genes,
  Degree = ifelse(cuproptosis_genes %in% names(degree(extended_graph)), 
                  degree(extended_graph)[cuproptosis_genes], 0),
  Betweenness = ifelse(cuproptosis_genes %in% names(betweenness(extended_graph)), 
                       betweenness(extended_graph)[cuproptosis_genes], 0),
  Closeness = ifelse(cuproptosis_genes %in% names(closeness(extended_graph)), 
                     closeness(extended_graph)[cuproptosis_genes], 0),
  Connected_78Genes_Count = cuproptosis_connected_count
)

# 生成拓扑对比数据
topology_comparison <- data.frame(
  Network_Type = c("基础网络", "扩展网络"),
  Nodes = c(vcount(base_graph), vcount(extended_graph)),
  Edges = c(base_topology$edges, extended_topology$edges),
  Avg_Degree = c(round(base_topology$avg_degree, 3), round(extended_topology$avg_degree, 3)),
  Density = c(round(base_topology$density, 3), round(extended_topology$density, 3)),
  APL = c(round(base_topology$apl, 3), round(extended_topology$apl, 3)),
  Clustering_Coefficient = c(round(base_topology$clustering, 3), round(extended_topology$clustering, 3))
)

# 生成 Hub 基因对比数据
base_hub_df <- data.frame(
  Network = "基础网络",
  Rank = 1:nrow(base_hub_genes),
  Gene = base_hub_genes$Gene,
  Degree = base_hub_genes$Degree,
  Betweenness = round(base_hub_genes$Betweenness, 3)
)

extended_hub_df <- data.frame(
  Network = "扩展网络",
  Rank = 1:nrow(extended_hub_genes),
  Gene = extended_hub_genes$Gene,
  Degree = extended_hub_genes$Degree,
  Betweenness = round(extended_hub_genes$Betweenness, 3)
)

hub_genes_comparison <- rbind(base_hub_df, extended_hub_df)

# 生成统计摘要
edge_increase <- extended_topology$edges - base_topology$edges
cross_edges_percent <- round((cross_edges / extended_topology$edges) * 100, 1)

# 构建统计摘要文本
cross_talk_text <- paste0(
  "基础网络包含 ", base_topology$edges, " 条有效互作边，平均节点度 ", round(base_topology$avg_degree, 3), "，网络密度 ", round(base_topology$density, 3), "。\n",
  "扩展网络包含 ", extended_topology$edges, " 条有效互作边（较基础网络增加 ", edge_increase, " 条），平均节点度 ", round(extended_topology$avg_degree, 3), "，网络密度 ", round(extended_topology$density, 3), "；其中铜死亡基因与交集基因之间存在 ", cross_edges, " 条直接互作边，占总边数的 ", cross_edges_percent, "%。\n",
  "铜死亡基因 Degree Top 3 为 ", names(cuproptosis_top3)[1], "（", cuproptosis_top3[1], "）、", names(cuproptosis_top3)[2], "（", cuproptosis_top3[2], "）、", names(cuproptosis_top3)[3], "（", cuproptosis_top3[3], "）。"
)

# 输出文件
output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 输出拓扑对比表
write.csv(topology_comparison, file = paste0(output_dir, "/topology_comparison.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 输出铜死亡基因中心性
write.csv(cuproptosis_centrality, file = paste0(output_dir, "/cuproptosis_centrality.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 输出 Hub 基因对比
write.csv(hub_genes_comparison, file = paste0(output_dir, "/hub_genes_comparison.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 输出统计摘要
writeLines(cross_talk_text, con = paste0(output_dir, "/cross_talk_statistics.txt"))

# 打印结果摘要
cat("\n分析完成！\n")
cat("输出文件：\n")
cat(paste0(output_dir, "/topology_comparison.csv\n"))
cat(paste0(output_dir, "/cuproptosis_centrality.csv\n"))
cat(paste0(output_dir, "/hub_genes_comparison.csv\n"))
cat(paste0(output_dir, "/cross_talk_statistics.txt\n"))

cat("\n统计摘要：\n")
cat(cross_talk_text)
