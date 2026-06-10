# 网络药理学分析脚本 - 替代Cytoscape插件流程
# 作者：生物信息学分析师
# 日期：2026-03-15
# R版本：4.5.2

# 设置随机种子以确保可重复性
set.seed(123)
options(stringsAsFactors = FALSE)

# 检查并安装必要的包
check_and_install_packages <- function() {
  required_packages <- c("igraph", "tidyverse", "VennDiagram", "ggVennDiagram")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
  
  # 加载包
  library(igraph)
  library(tidyverse)
  library(VennDiagram)
  library(ggVennDiagram)
}

# 执行包检查和安装
check_and_install_packages()

# 1. 数据读入与网络构建
cat("\nStep 1: 数据读入与网络构建...\n")

# 读取PPI数据
ppi_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/78ppi.tsv"
ppi_data <- read.delim(ppi_file, header = TRUE, sep = "\t")

# 检查列名并标准化
if ("source" %in% colnames(ppi_data) && "target" %in% colnames(ppi_data)) {
  colnames(ppi_data) <- c("node1", "node2")
}

# 构建无向图
g <- graph_from_data_frame(ppi_data, directed = FALSE)

# 去除自环和重复边
g <- simplify(g)

cat("网络构建完成：")
cat("节点数 =", vcount(g), "")
cat("边数 =", ecount(g), "\n")

# 2. 替代cytoNCA：8指标拓扑筛选（中位数阈值）
cat("\nStep 2: 替代cytoNCA - 8指标拓扑筛选...\n")

# 计算8种中心性指标
calculate_centrality <- function(g) {
  # BC (Betweenness)
  bc <- betweenness(g)
  
  # CC (Closeness)
  cc <- closeness(g)
  
  # DC (Degree)
  dc <- degree(g)
  
  # EC (Eigenvector)
  ec <- eigen_centrality(g)$vector
  
  # LAC (Local Average Connectivity) - 节点邻居的平均度
  lac <- sapply(V(g), function(v) {
    neighbors <- neighbors(g, v)
    if (length(neighbors) == 0) return(0)
    mean(degree(g, neighbors))
  })
  
  # NC (Network) - 中心度
  nc <- centr_degree(g)$res
  
  # SC (Subgraph) - 用subgraph_centrality()
  sc <- tryCatch(
    subgraph_centrality(g),
    error = function(e) {
      # 如果不支持，用alpha_centrality()替代
      alpha_centrality(g)
    }
  )
  
  # IC (Information) - 用harmonic_centrality()作为近似
  ic <- harmonic_centrality(g)
  
  # 整理结果
  centrality_df <- data.frame(
    gene = V(g)$name,
    BC = bc,
    CC = cc,
    DC = dc,
    EC = ec,
    LAC = lac,
    NC = nc,
    SC = sc,
    IC = ic
  )
  
  return(centrality_df)
}

centrality_df <- calculate_centrality(g)

# 计算各指标中位数并筛选
screen_centrality <- function(centrality_df) {
  thresholds <- apply(centrality_df[, -1], 2, median)
  
  # 筛选满足≥中位数的节点（至少满足6个指标）
  selected <- apply(centrality_df[, -1], 1, function(row) {
    sum(row >= thresholds) >= 6
  })
  
  return(list(
    selected_genes = centrality_df$gene[selected],
    centrality_df = centrality_df,
    thresholds = thresholds
  ))
}

cytoNCA_result <- screen_centrality(centrality_df)
list_cytoNCA <- cytoNCA_result$selected_genes

cat("cytoNCA 筛选后节点数:", length(list_cytoNCA), "\n")

# 3. 替代MCODE：密度聚类精炼
cat("\nStep 3: 替代MCODE - 密度聚类精炼...\n")

# K-core 分解，筛选 coreness >= 2 的节点
coreness_values <- coreness(g)
k_core_nodes <- V(g)$name[coreness_values >= 2]

# 密度筛选：计算局部聚类系数，保留系数 ≥ 0.2 的节点
transitivity_values <- transitivity(g, type = "local")
density_nodes <- V(g)$name[!is.na(transitivity_values) & transitivity_values >= 0.2]

# 取交集
filtered_nodes <- intersect(k_core_nodes, density_nodes)

# 在筛选后的子图上运行模块识别
if (length(filtered_nodes) > 0) {
  subgraph <- induced_subgraph(g, filtered_nodes)
  
  # 尝试使用cluster_louvain()
  tryCatch({
    clusters <- cluster_louvain(subgraph)
  }, error = function(e) {
    # 如果失败，使用cluster_fast_greedy()
    clusters <- cluster_fast_greedy(subgraph)
  })
  
  # 提取模块内节点数 ≥ 3 的簇
  module_sizes <- sizes(clusters)
  valid_modules <- which(module_sizes >= 3)
  
  if (length(valid_modules) > 0) {
    # 收集所有有效模块的节点
    module_nodes <- c()
    module_assignments <- data.frame()
    
    for (module_id in valid_modules) {
      module_genes <- V(subgraph)$name[clusters$membership == module_id]
      module_nodes <- c(module_nodes, module_genes)
      
      # 收集模块信息
      module_info <- data.frame(
        gene = module_genes,
        module = module_id,
        local_density = transitivity_values[module_genes]
      )
      module_assignments <- rbind(module_assignments, module_info)
    }
    
    # 取最大模块或前3个模块
    if (length(valid_modules) > 3) {
      # 按模块大小排序，取前3个
      sorted_modules <- sort(module_sizes, decreasing = TRUE)[1:3]
      top_module_nodes <- c()
      for (module_id in names(sorted_modules)) {
        top_module_nodes <- c(top_module_nodes, V(subgraph)$name[clusters$membership == as.numeric(module_id)])
      }
      list_MCODE <- unique(top_module_nodes)
    } else {
      list_MCODE <- unique(module_nodes)
    }
  } else {
    # 如果没有有效模块，使用所有筛选后的节点
    list_MCODE <- filtered_nodes
    module_assignments <- data.frame(
      gene = filtered_nodes,
      module = 1,
      local_density = transitivity_values[filtered_nodes]
    )
  }
} else {
  # 如果没有筛选后的节点，使用原始节点
  list_MCODE <- V(g)$name
  module_assignments <- data.frame(
    gene = V(g)$name,
    module = 1,
    local_density = transitivity_values[V(g)$name]
  )
}

cat("MCODE 筛选后节点数:", length(list_MCODE), "\n")

# 4. 替代cytoHubba MCC：最大团中心性
cat("\nStep 4: 替代cytoHubba - 最大团中心性(MCC)...\n")

# 计算MCC (Maximal Clique Centrality)
calculate_mcc <- function(g) {
  # 找到所有极大团
  cliques <- max_cliques(g)
  
  # 统计每个节点参与的团数
  mcc_scores <- rep(0, vcount(g))
  names(mcc_scores) <- V(g)$name
  
  for (clique in cliques) {
    clique_nodes <- V(g)$name[clique]
    mcc_scores[clique_nodes] <- mcc_scores[clique_nodes] + 1
  }
  
  # 归一化
  mcc_scores <- mcc_scores / max(mcc_scores)
  
  return(mcc_scores)
}

mcc_scores <- calculate_mcc(g)

# 提取Top 10节点
list_cytoHubba <- names(sort(mcc_scores, decreasing = TRUE))[1:10]

cat("cytoHubba Top 10 节点:", length(list_cytoHubba), "\n")

# 5. 三算法交集与可视化
cat("\nStep 5: 三算法交集与可视化...\n")

# 计算三算法交集
key_genes <- intersect(intersect(list_cytoNCA, list_MCODE), list_cytoHubba)

cat("三算法交集关键基因数量:", length(key_genes), "\n")
cat("三算法交集关键基因:", paste(key_genes, collapse = ", "), "\n")

# 6. 输出结果
cat("\nStep 6: 结果输出...\n")

# 创建output目录
output_dir <- "output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 保存cytoNCA结果
cytoNCA_output <- centrality_df %>%
  filter(gene %in% list_cytoNCA)
write.csv(cytoNCA_output, file.path(output_dir, "cytoNCA_genes.csv"), row.names = FALSE)
cat("cytoNCA结果已保存: output/cytoNCA_genes.csv\n")

# 保存MCODE结果
write.csv(module_assignments, file.path(output_dir, "MCODE_module_genes.csv"), row.names = FALSE)
cat("MCODE结果已保存: output/MCODE_module_genes.csv\n")

# 保存cytoHubba结果
cytoHubba_output <- data.frame(
  gene = list_cytoHubba,
  mcc_score = mcc_scores[list_cytoHubba],
  rank = 1:10
)
write.csv(cytoHubba_output, file.path(output_dir, "cytoHubba_Top10.csv"), row.names = FALSE)
cat("cytoHubba结果已保存: output/cytoHubba_Top10.csv\n")

# 保存关键基因交集
if (length(key_genes) > 0) {
  key_genes_output <- data.frame(gene = key_genes)
  write.csv(key_genes_output, file.path(output_dir, "key_genes_intersection.csv"), row.names = FALSE)
  cat("关键基因交集已保存: output/key_genes_intersection.csv\n")
}

# 7. 可视化
cat("\nStep 7: 可视化...\n")

# 绘制三集合韦恩图
venn_list <- list(
  cytoNCA = list_cytoNCA,
  MCODE = list_MCODE,
  cytoHubba = list_cytoHubba
)

# 使用ggVennDiagram绘制韦恩图
venn_plot <- ggVennDiagram(venn_list, label_alpha = 0.7)
# 添加标题
venn_plot <- venn_plot + 
  ggtitle("三算法筛选结果韦恩图") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

# 保存韦恩图
ggsave(file.path(output_dir, "venn_diagram.pdf"), venn_plot, width = 10, height = 8, dpi = 300)
cat("韦恩图已保存: output/venn_diagram.pdf\n")

# 绘制网络图，用不同颜色标注三算法筛选出的节点
# 为每个节点分配颜色属性
node_colors <- rep("gray", vcount(g))
names(node_colors) <- V(g)$name

# cytoNCA=蓝色, MCODE=绿色, cytoHubba=红色，交集=紫色
node_colors[list_cytoNCA] <- "blue"
node_colors[list_MCODE] <- "green"
node_colors[list_cytoHubba] <- "red"
node_colors[key_genes] <- "purple"

# 设置节点大小（基于度）
node_sizes <- degree(g) * 3

# 绘制网络图
pdf(file.path(output_dir, "network_topology_plot.pdf"), width = 12, height = 10)
tryCatch({
  plot(g,
       vertex.color = node_colors,
       vertex.size = node_sizes,
       vertex.label = NA,  # 不显示标签以避免混乱
       edge.width = 1,
       edge.color = "gray80",
       layout = layout_with_kk,
       main = "网络拓扑图（三算法筛选结果）"
  )
  
  # 添加图例
  legend("bottomright",
         legend = c("cytoNCA", "MCODE", "cytoHubba", "交集"),
         col = c("blue", "green", "red", "purple"),
         pch = 19,
         cex = 0.8
  )
}, finally = {
  dev.off()
})

cat("网络图已保存: output/network_topology_plot.pdf\n")

# 绘制条形图展示各算法的基因数量及交集大小
algorithm_counts <- data.frame(
  algorithm = c("cytoNCA", "MCODE", "cytoHubba", "三算法交集"),
  count = c(length(list_cytoNCA), length(list_MCODE), length(list_cytoHubba), length(key_genes))
)

bar_plot <- ggplot(algorithm_counts, aes(x = algorithm, y = count, fill = algorithm)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = -0.5, size = 4) +
  labs(title = "各算法筛选基因数量", x = "算法", y = "基因数量") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

ggsave(file.path(output_dir, "algorithm_counts_barplot.pdf"), bar_plot, width = 8, height = 6, dpi = 300)
cat("算法计数条形图已保存: output/algorithm_counts_barplot.pdf\n")

# 获取 cytoNCA 子网的边数
g_nca <- induced_subgraph(g, list_cytoNCA)
nca_edges <- ecount(g_nca)  # 这就是第一个【】中的边数

# 获取 MCODE 最大模块的边数（假设取 module 4）
mcode_genes <- module_assignments %>% filter(module == 4) %>% pull(gene)
if (length(mcode_genes) == 0) {
  # 如果没有 module 4，取最大的模块
  module_sizes <- table(module_assignments$module)
  max_module <- as.numeric(names(module_sizes)[which.max(module_sizes)])
  mcode_genes <- module_assignments %>% filter(module == max_module) %>% pull(gene)
}
g_mcode <- induced_subgraph(g, mcode_genes)
mcode_edges <- ecount(g_mcode)  # 这就是第二个【】中的边数

# 控制台输出汇总信息
cat("\n=== 分析结果汇总 ===\n")
cat("cytoNCA 筛选后节点数:", length(list_cytoNCA), "\n")
cat("MCODE 筛选后节点数:", length(list_MCODE), "\n")
cat("cytoHubba Top 10 节点数:", length(list_cytoHubba), "\n")
cat("三算法交集关键基因数量:", length(key_genes), "\n")
if (length(key_genes) > 0) {
  cat("三算法交集关键基因:", paste(key_genes, collapse = ", "), "\n")
}
cat("cytoNCA 子网边数:", nca_edges, "\n")
cat("MCODE 最大模块边数:", mcode_edges, "\n")

# 分析完成
cat("\n分析完成！所有结果文件已输出至output目录。\n")
