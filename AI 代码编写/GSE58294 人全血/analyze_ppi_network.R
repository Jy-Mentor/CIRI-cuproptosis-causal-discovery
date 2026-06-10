# 蛋白质互作网络分析

# 设置CRAN镜像
options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

# 安装和加载必要的包
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("igraph")) install.packages("igraph")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("RColorBrewer")) install.packages("RColorBrewer")

library(tidyverse)
library(igraph)
library(ggplot2)
library(RColorBrewer)

# 设置工作目录
setwd("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写\\GSE58294 人全血")

# 1. 读取STRING互作网络数据
read_string_data <- function(file_path) {
  data <- read_tsv(file_path)
  # 查看列名
  cat("列名：", paste(colnames(data), collapse = ", "), "\n")
  # 重命名列
  data_renamed <- data %>%
    rename(
      from = `#node1`,
      to = node2
    )
  # 不使用阈值筛选，保留所有互作
  data_filtered <- data_renamed
  return(data_filtered)
}

# 文件路径
ppi_files <- list(
  "3H" = "PPI STRING\\string_interactions (3h).tsv",
  "5H" = "PPI STRING\\string_interactions (5h).tsv",
  "24H" = "PPI STRING\\string_interactions (24h).tsv"
)

# 读取并处理数据
ppi_data_list <- list()
for (time_point in names(ppi_files)) {
  cat("读取", time_point, "时间点的PPI数据...\n")
  ppi_data_list[[time_point]] <- read_string_data(ppi_files[[time_point]])
  cat("筛选后互作数：", nrow(ppi_data_list[[time_point]]), "\n")
}

# 2. 构建蛋白质互作网络
build_network <- function(ppi_data) {
  # 创建边列表
  edges <- ppi_data %>%
    select(from, to, combined_score)
  
  # 构建图
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # 添加节点属性
  V(g)$degree <- degree(g)
  V(g)$betweenness <- betweenness(g)
  V(g)$closeness <- closeness(g)
  V(g)$eigenvector <- eigen_centrality(g)$vector
  
  return(g)
}

# 构建所有时间点的网络
network_list <- list()
for (time_point in names(ppi_data_list)) {
  cat("构建", time_point, "时间点的网络...\n")
  network_list[[time_point]] <- build_network(ppi_data_list[[time_point]])
  cat("节点数：", vcount(network_list[[time_point]]), "，边数：", ecount(network_list[[time_point]]), "\n")
}

# 3. 识别Hub基因
identify_hub_genes <- function(g, top_n = 10) {
  # 计算MCC (Maximum Clique Centrality) - 这里使用度中心性作为替代
  # 在实际应用中，应使用Cytoscape的cytoHubba插件计算MCC
  
  # 综合多种中心性指标
  hub_scores <- data.frame(
    gene = V(g)$name,
    degree = V(g)$degree,
    betweenness = V(g)$betweenness,
    closeness = V(g)$closeness,
    eigenvector = V(g)$eigenvector
  )
  
  # 标准化各项指标
  hub_scores$degree_norm <- scale(hub_scores$degree)[, 1]
  hub_scores$betweenness_norm <- scale(hub_scores$betweenness)[, 1]
  hub_scores$closeness_norm <- scale(hub_scores$closeness)[, 1]
  hub_scores$eigenvector_norm <- scale(hub_scores$eigenvector)[, 1]
  
  # 计算综合得分
  hub_scores$综合得分 <- rowMeans(hub_scores[, c("degree_norm", "betweenness_norm", "closeness_norm", "eigenvector_norm")], na.rm = TRUE)
  
  # 排序并选择top_n个Hub基因
  top_hubs <- hub_scores %>%
    arrange(desc(综合得分)) %>%
    head(top_n)
  
  return(list(
    hub_scores = hub_scores,
    top_hubs = top_hubs
  ))
}

# 识别所有时间点的Hub基因
hub_analysis_list <- list()
for (time_point in names(network_list)) {
  cat("识别", time_point, "时间点的Hub基因...\n")
  hub_analysis_list[[time_point]] <- identify_hub_genes(network_list[[time_point]], top_n = 10)
  cat("Top 10 Hub基因：", paste(hub_analysis_list[[time_point]]$top_hubs$gene, collapse = ", "), "\n")
}

# 4. 可视化网络
output_dir <- "output/ppi_analysis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 绘制网络可视化
plot_network <- function(g, time_point, top_hubs) {
  # 设置布局
  layout <- layout_with_fr(g)
  
  # 设置节点大小和颜色
  V(g)$size <- V(g)$degree * 2 + 3
  V(g)$color <- ifelse(V(g)$name %in% top_hubs$gene, "red", "skyblue")
  
  # 绘制网络
  png(file.path(output_dir, paste0("network_visualization_", time_point, ".png")), width = 1200, height = 1000)
  plot(g, 
       layout = layout,
       vertex.label = ifelse(V(g)$name %in% top_hubs$gene, V(g)$name, NA),
       vertex.label.cex = 1.2,
       vertex.label.color = "black",
       edge.width = E(g)$combined_score / 200,
       edge.color = "gray80",
       main = paste("蛋白质互作网络（", time_point, "）", sep = ""))
  
  # 添加图例
  legend("bottomright", 
         legend = c("Hub基因", "非Hub基因"),
         col = c("red", "skyblue"),
         pch = 19,
         pt.cex = 2,
         cex = 1.2,
         bty = "n")
  
dev.off()
}

# 绘制所有时间点的网络
for (time_point in names(network_list)) {
  cat("绘制", time_point, "时间点的网络可视化...\n")
  plot_network(network_list[[time_point]], time_point, hub_analysis_list[[time_point]]$top_hubs)
}

# 5. 绘制Hub基因排名图
plot_hub_ranking <- function(hub_scores, time_point, top_n = 10) {
  # 排序并选择top_n
  top_hubs <- hub_scores %>%
    arrange(desc(综合得分)) %>%
    head(top_n)
  
  # 绘制条形图
  png(file.path(output_dir, paste0("hub_ranking_", time_point, ".png")), width = 1000, height = 600)
  ggplot(top_hubs, aes(x = reorder(gene, 综合得分), y = 综合得分)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste("Hub基因排名（", time_point, "）", sep = ""),
         x = "基因",
         y = "综合得分") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 10))
  dev.off()
}

# 绘制所有时间点的Hub基因排名
for (time_point in names(hub_analysis_list)) {
  cat("绘制", time_point, "时间点的Hub基因排名...\n")
  plot_hub_ranking(hub_analysis_list[[time_point]]$hub_scores, time_point)
}

# 6. 分析Hub基因的时间变化
# 提取所有时间点的Hub基因
hub_genes_all <- list()
for (time_point in names(hub_analysis_list)) {
  hub_genes_all[[time_point]] <- hub_analysis_list[[time_point]]$top_hubs$gene
}

# 找出在多个时间点出现的Hub基因
common_hubs <- Reduce(intersect, hub_genes_all)
cat("\n在所有时间点均为Hub基因的基因：", paste(common_hubs, collapse = ", "), "\n")

# 7. 保存分析结果
saveRDS(list(
  ppi_data_list = ppi_data_list,
  network_list = network_list,
  hub_analysis_list = hub_analysis_list,
  common_hubs = common_hubs
), file.path(output_dir, "ppi_analysis_results.rds"))

# 8. 生成分析报告
report_file <- file.path(output_dir, "ppi_analysis_report.txt")
sink(report_file)

cat("蛋白质互作网络分析报告\n")
cat("======================\n\n")

# 报告网络基本信息
cat("1. 网络基本信息\n")
cat("================\n")
for (time_point in names(network_list)) {
  cat(paste0("\n时间点：", time_point, "\n"))
  cat("- 节点数：", vcount(network_list[[time_point]]), "\n")
  cat("- 边数：", ecount(network_list[[time_point]]), "\n")
  cat("- 平均度：", mean(degree(network_list[[time_point]])), "\n")
  cat("- 网络密度：", graph.density(network_list[[time_point]]), "\n")
}
cat("\n")

# 报告Hub基因
cat("2. Hub基因分析\n")
cat("===============\n")
for (time_point in names(hub_analysis_list)) {
  cat(paste0("\n时间点：", time_point, "\n"))
  cat("Top 10 Hub基因：\n")
  top_hubs <- hub_analysis_list[[time_point]]$top_hubs
  for (i in 1:nrow(top_hubs)) {
    cat(sprintf("%d. %s (综合得分: %.4f, 度: %d)\n", 
                i, top_hubs$gene[i], top_hubs$综合得分[i], top_hubs$degree[i]))
  }
}
cat("\n")

# 报告共同Hub基因
cat("3. 共同Hub基因\n")
cat("===============\n")
if (length(common_hubs) > 0) {
  cat("在所有时间点均为Hub基因的基因：\n")
  for (gene in common_hubs) {
    cat(paste0("- ", gene, "\n"))
  }
} else {
  cat("无共同Hub基因\n")
}
cat("\n")

# 报告总结
cat("4. 分析总结\n")
cat("============\n")
cat("- 分析了3个时间点的蛋白质互作网络\n")
cat("- 识别了每个时间点的Top 10 Hub基因\n")
cat("- 可视化了网络结构和Hub基因\n")
cat("- 结果保存在：", output_dir, "\n")
cat("\n可视化结果包括：\n")
for (time_point in names(network_list)) {
  cat(paste0("- network_visualization_", time_point, ".png: ", time_point, "时间点的网络可视化\n"))
  cat(paste0("- hub_ranking_", time_point, ".png: ", time_point, "时间点的Hub基因排名\n"))
}

sink()

print("PPI网络分析完成！结果保存在output/ppi_analysis目录中。")
