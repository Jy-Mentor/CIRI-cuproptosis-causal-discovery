# PPI 网络分析脚本

# 智能检查并安装包
install.packages.if.needed <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# 安装必要的包
install.packages.if.needed(c("igraph", "ggraph", "tidyverse"))

# 加载包
library(igraph)
library(ggraph)
library(tidyverse)

# 读取数据
file_path <- "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创\\78ppi.tsv"
data <- read_tsv(file_path)

# 过滤 combined_score > 0.4
data_filtered <- data %>% filter(combined_score > 0.4)

# 创建图对象
g <- graph_from_data_frame(data_filtered, directed = FALSE)

# 计算 Degree 和 Betweenness
degree_values <- degree(g)
betweenness_values <- betweenness(g)

# 准备节点属性
g <- set_vertex_attr(g, "degree", value = degree_values)
g <- set_vertex_attr(g, "betweenness", value = betweenness_values)

# 基因类型分类（这里需要根据实际情况调整，假设需要手动指定或从其他文件读取）
# 这里暂时创建一个示例分类，实际使用时需要修改
copper_death_genes <- c("TP53", "SOD1", "ATP7A", "ATP7B", "MTF1")  # 示例铜死亡基因
bcp_target_genes <- c("AKT1", "EGFR", "VEGFA", "MAPK1", "TP53")  # 示例BCP靶点

# 为节点添加颜色属性
V(g)$color <- ifelse(V(g)$name %in% copper_death_genes, "#7b3294",  # 铜死亡紫色
                     ifelse(V(g)$name %in% bcp_target_genes, "#1a9850",  # BCP靶点绿色
                            "#999999"))  # 其他基因灰色

# 计算 Degree 前 10% 的 Hub 基因
top_10_percent <- round(vcount(g) * 0.1)
top_genes <- names(sort(degree_values, decreasing = TRUE))[1:top_10_percent]

# 创建标签属性（仅标注前10%的Hub基因）
V(g)$label <- ifelse(V(g)$name %in% top_genes, V(g)$name, "")

# 绘制网络
p <- ggraph(g, layout = "fr") +  # Fruchterman-Reingold 布局
  geom_edge_link(alpha = 0.3) +
  geom_node_point(aes(size = degree, color = color)) +
  scale_size_continuous(range = c(2, 10)) +
  scale_color_identity() +
  geom_node_text(aes(label = label), size = 3, repel = TRUE) +
  theme_void() +
  labs(title = "PPI Network Analysis")

# 输出 PDF
pdf("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创\\ppi_network.pdf", width = 12, height = 10)
print(p)
dev.off()

# 输出 300dpi PNG
png("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创\\ppi_network.png", width = 1200, height = 1000, res = 300)
print(p)
dev.off()

# 打印摘要信息
cat("网络分析完成！\n")
cat(paste("节点数量:", vcount(g), "\n"))
cat(paste("边数量:", ecount(g), "\n"))
cat("前10% Hub基因:", paste(top_genes, collapse = ", "), "\n")
