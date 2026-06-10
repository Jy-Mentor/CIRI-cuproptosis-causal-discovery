#!/usr/bin/env Rscript
# Hub基因筛选分析脚本
# 严格遵循AI代码编写规则

set.seed(123)

# ==============================================================================
# 0. 环境准备：自动安装并加载所需包
# ==============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quiet = TRUE)

required_packages <- c("STRINGdb", "igraph", "ggplot2", "VennDiagram", "dplyr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg == "STRINGdb") {
      BiocManager::install(pkg, ask = FALSE, quiet = TRUE)
    } else {
      install.packages(pkg, quiet = TRUE, repos = "https://cloud.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# 1. 基因交集分析
# ==============================================================================
cat("=== 步骤1：基因交集分析 ===\n")

gene_list1 <- c("ACTA2", "ADORA1", "AIF1", "ALDH1A1", "ALDH9A1", "AOC3", "ATF4", "BRD4", "C3", "CASP8", "CCL2", "CCND1", "CCR5", "CDC42", "CNDP2", "CNR2", "COL1A1", "CP", "CPT1A", "CPT2", "CTSB", "CTSD", "CTSS", "CXCR3", "DDIT3", "EGR1", "F3", "FABP3", "FABP5", "FAS", "FASN", "GAD1", "GFAP", "GPT", "HMGCR", "HMOX1", "HSPA5", "HTR2A", "ICAM1", "IGF1R", "IL6", "IRF1", "JAK1", "MAOB", "MAPK9", "MDM2", "MGLL", "NFE2L2", "NFKB1", "NOTCH1", "NR1H3", "PARP1", "PLA2G4A", "PRKCQ", "PTGES", "PTGS1", "PTGS2", "PTPN6", "PTPRC", "RELA", "S100A6", "S1PR1", "SAT1", "SOD2", "SREBF1", "STAT1", "STAT3", "STAT5A", "TIMP1", "TGFB1", "TSPO", "XDH")

gene_list2 <- c("FDX1", "DLAT", "DLD", "LIPT1", "PDHX", "PDHB", "SLC31A1", "ATP7B", "ATP7A", "ATOX1", "COMMD1", "MT2A", "NFKB1", "ATF4", "TLR4")

direct_intersect <- intersect(gene_list1, gene_list2)
cat("直接交集基因:", paste(direct_intersect, collapse = ", "), "\n")

expand_genes <- union(gene_list1, gene_list2)
cat("合并去重后基因数量:", length(expand_genes), "\n\n")

# ==============================================================================
# 2. PPI数据获取 (STRING数据库)
# ==============================================================================
cat("=== 步骤2：获取PPI数据 ===\n")

string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 700,
  input_directory = tempdir()
)

gene_mapping <- string_db$map(
  data.frame(Gene = expand_genes),
  "Gene",
  removeUnmappedRows = TRUE
)

if (nrow(gene_mapping) == 0) {
  stop("错误：没有基因能成功映射到STRING数据库")
}

cat("成功映射基因数:", nrow(gene_mapping), "\n")

ppi_data <- string_db$get_interactions(gene_mapping$STRING_id)

if (nrow(ppi_data) == 0) {
  stop("错误：未找到基因间的互作关系")
}

cat("获得互作关系数:", nrow(ppi_data), "\n\n")

# ==============================================================================
# 3. 网络构建与初筛
# ==============================================================================
cat("=== 步骤3：网络构建与初筛 ===\n")

ppi_network <- graph_from_data_frame(
  d = ppi_data[, c("from", "to")],
  directed = FALSE
)

id2gene <- setNames(gene_mapping$Gene, gene_mapping$STRING_id)
V(ppi_network)$gene_name <- id2gene[V(ppi_network)$name]

cat("正在计算8种中心性指标...\n")
cat("  1. Degree (度中心性)\n")
cat("  2. Betweenness (介数中心性)\n")
cat("  3. Closeness (接近中心性)\n")
cat("  4. Eigenvector (特征向量中心性)\n")
cat("  5. PageRank\n")
cat("  6. Subgraph (子图中心性)\n")
cat("  7. Authority (权威性)\n")
cat("  8. HubScore\n")

centrality_df <- data.frame(
  Gene = V(ppi_network)$gene_name,
  Degree = degree(ppi_network),
  Betweenness = betweenness(ppi_network, normalized = TRUE),
  Closeness = closeness(ppi_network, normalized = TRUE),
  Eigenvector = eigen_centrality(ppi_network)$vector,
  PageRank = page_rank(ppi_network)$vector,
  Subgraph = subgraph_centrality(ppi_network),
  Authority = authority_score(ppi_network)$vector,
  HubScore = hub_score(ppi_network)$vector,
  stringsAsFactors = FALSE
)

centrality_df <- na.omit(centrality_df)

medians <- apply(centrality_df[, -1], 2, median)
filter_pass <- apply(sweep(centrality_df[, -1], 2, medians, ">="), 1, all)
selected_genes <- centrality_df$Gene[filter_pass]

if (length(selected_genes) == 0) {
  stop("错误：初筛后没有剩余节点！")
}

subnet_initial <- induced_subgraph(
  ppi_network,
  V(ppi_network)$gene_name %in% selected_genes
)

cat("\n初筛后节点数:", vcount(subnet_initial), "边数:", ecount(subnet_initial), "\n\n")

# ==============================================================================
# 4. 网络密度验证
# ==============================================================================
cat("=== 步骤4：网络密度验证 ===\n")

real_density <- edge_density(subnet_initial)
cat("真实子网密度:", round(real_density, 4), "\n")

num_nodes <- vcount(subnet_initial)
num_edges <- ecount(subnet_initial)

random_densities <- replicate(1000, {
  g <- erdos.renyi.game(n = num_nodes, m = num_edges, type = "gnm", directed = FALSE)
  edge_density(g)
})

p_value <- t.test(random_densities, mu = real_density, alternative = "less")$p.value
cat("随机网络密度均值:", round(mean(random_densities), 4), "\n")
cat("P值 (真实 > 随机):", format.pval(p_value, digits = 3), "\n")

if (p_value < 0.05) {
  cat("结论：子网密度显著高于随机网络 (P < 0.05)\n\n")
} else {
  cat("结论：子网密度与随机网络无显著差异\n\n")
}

# ==============================================================================
# 5. 子网精炼 (K-core + MCODE)
# ==============================================================================
cat("=== 步骤5：子网精炼 ===\n")

core_levels <- coreness(subnet_initial)
subnet_kcore <- induced_subgraph(subnet_initial, core_levels >= 2)

if (vcount(subnet_kcore) == 0) {
  stop("错误：K-core分解后没有剩余节点！")
}

cat("K-core (K=2) 后节点数:", vcount(subnet_kcore), "\n")
cat("正在进行MCODE聚类...\n")

mcode_result <- mcode(
  subnet_kcore,
  vwp = 0.5,
  haircut = TRUE
)

module_sizes <- sapply(mcode_result$modules, length)
if (length(module_sizes) == 0) {
  subnet_refined <- subnet_kcore
} else {
  largest_module <- mcode_result$modules[[which.max(module_sizes)]]
  subnet_refined <- induced_subgraph(subnet_kcore, largest_module)
}

cat("精炼核心网络信息:\n")
cat("  节点数:", vcount(subnet_refined), "\n")
cat("  边数:", ecount(subnet_refined), "\n")
cat("  平均节点度:", round(mean(degree(subnet_refined)), 2), "\n")
cat("  网络密度:", round(edge_density(subnet_refined), 3), "\n\n")

# ==============================================================================
# 6. Hub基因筛选 (MCC + 韦恩图)
# ==============================================================================
cat("=== 步骤6：Hub基因筛选 ===\n")

calculate_mcc <- function(graph) {
  nodes <- V(graph)$gene_name
  mcc_scores <- numeric(length(nodes))
  cliques_info <- max_cliques(graph, min = 3)

  for (i in seq_along(nodes)) {
    node <- nodes[i]
    count <- sum(sapply(cliques_info, function(clq) node %in% as.numeric(clq)))
    mcc_scores[i] <- count
  }
  return(mcc_scores)
}

mcc_scores <- calculate_mcc(subnet_refined)
hub_candidates <- data.frame(
  Gene = V(subnet_refined)$gene_name,
  MCC = mcc_scores,
  stringsAsFactors = FALSE
)
hub_candidates <- hub_candidates[order(-hub_candidates$MCC), ]

top8_hub <- head(hub_candidates$Gene, 8)
cat("最终Top 8 Hub基因 (按MCC排序):\n")
print(hub_candidates[1:8, ])

refined_centrality <- data.frame(
  Gene = V(subnet_refined)$gene_name,
  Degree = degree(subnet_refined),
  Betweenness = betweenness(subnet_refined, normalized = TRUE),
  Closeness = closeness(subnet_refined, normalized = TRUE),
  MCC = mcc_scores,
  stringsAsFactors = FALSE
)

get_top <- function(x, n = 8) {
  head(refined_centrality$Gene[order(-x)], n)
}

venn_list <- list(
  Degree = get_top(refined_centrality$Degree),
  Betweenness = get_top(refined_centrality$Betweenness),
  Closeness = get_top(refined_centrality$Closeness),
  MCC = top8_hub
)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1,
  main = "多算法Top 8基因交集"
)

# ==============================================================================
# 7. 可视化输出
# ==============================================================================
cat("\n=== 步骤7：生成可视化 ===\n")

pdf("1_Initial_Subnetwork.pdf", width = 10, height = 10)
set.seed(123)
plot(
  subnet_initial,
  vertex.label = V(subnet_initial)$gene_name,
  vertex.size = 8,
  vertex.color = "lightblue",
  vertex.label.cex = 0.7,
  edge.width = 0.8,
  main = "初筛功能子网"
)
dev.off()

node_colors <- ifelse(V(subnet_refined)$gene_name %in% top8_hub, "red", "gold")
node_sizes <- ifelse(V(subnet_refined)$gene_name %in% top8_hub, 15, 10)

pdf("2_Refined_Core_Network.pdf", width = 10, height = 10)
set.seed(123)
plot(
  subnet_refined,
  vertex.label = V(subnet_refined)$gene_name,
  vertex.size = node_sizes,
  vertex.color = node_colors,
  vertex.label.cex = 0.8,
  vertex.label.color = "black",
  edge.width = 1,
  main = "精炼核心网络 (红色为Hub基因)"
)
legend(
  "topleft",
  legend = c("Hub基因", "其他基因"),
  pch = 21,
  pt.cex = 2,
  pt.bg = c("red", "gold")
)
dev.off()

pdf("3_Venn_Diagram.pdf", width = 8, height = 8)
grid.draw(venn.plot)
dev.off()

cat("分析完成！已生成3个PDF文件:\n")
cat("  1. 1_Initial_Subnetwork.pdf (初筛子网)\n")
cat("  2. 2_Refined_Core_Network.pdf (核心网络)\n")
cat("  3. 3_Venn_Diagram.pdf (韦恩图)\n")
cat("\n最终Hub基因列表:", paste(top8_hub, collapse = ", "), "\n")