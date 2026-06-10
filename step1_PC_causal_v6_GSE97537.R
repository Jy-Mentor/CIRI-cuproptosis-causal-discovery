# ============================================================================
# Step 1: PC因果网络分析 - 使用GSE97537 DEGs
# 研究问题: BCP是否通过RAGE→NFKB1→FDX1轴发挥作用
# 数据集: GSE163614 (大鼠MCAO Bulk RNA-seq) + GSE97537 DEGs
# 注意: GSE97537是大鼠数据，与GSE163614同物种，无需跨物种映射
# ============================================================================

# 清除环境
rm(list = ls())
gc()

# 设置工作目录
setwd("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙")

# 加载必要的包
library(readxl)
library(pcalg)
library(igraph)

# 设置随机种子确保可复现
set.seed(42)

# ============================================================================
# 第一部分：读取Bulk数据和GSE97537 DEG结果
# ============================================================================

cat("=== 第一步：读取GSE163614 Bulk RNA-seq数据 ===\n")

# 读取GSE163614 Bulk数据
bulk_file <- "C:/Users/Jy-Mentor-7/Downloads/GSE163614_mRNA_Expression_Profiling.xlsx"
if (!file.exists(bulk_file)) {
  stop("错误: 找不到GSE163614数据文件!")
}

# 读取数据（跳过前5行表头）
bulk_df <- read_excel(bulk_file, skip = 5)

# 设置列名
colnames(bulk_df)[1] <- "gene_symbol"
colnames(bulk_df)[7:12] <- c("MCAO1", "MCAO2", "MCAO3", "Sham1", "Sham2", "Sham3")

# 提取表达矩阵
bulk_expr <- as.matrix(bulk_df[, 7:12])
rownames(bulk_expr) <- make.names(bulk_df$gene_symbol, unique = TRUE)
bulk_expr <- apply(bulk_expr, 2, as.numeric)
rownames(bulk_expr) <- make.names(bulk_df$gene_symbol, unique = TRUE)

cat(sprintf("Bulk数据: %d 基因 x %d 样本\n", nrow(bulk_expr), ncol(bulk_expr)))
cat("基因名示例:", head(rownames(bulk_expr), 5), "\n")

# ============================================================================
# 第二部分：读取并筛选GSE97537 DEGs
# ============================================================================

cat("\n=== 第二步：读取GSE97537 DEG结果并应用分层筛选 ===\n")

# 读取GSE97537 GEO2R结果（大鼠DEGs，无需跨物种映射）
deg_file <- "C:/Users/Jy-Mentor-7/Downloads/GSE97537.top.table (1).tsv"
if (!file.exists(deg_file)) {
  stop("错误: 找不到GSE97537 DEG文件!")
}

deg_df <- read.delim(deg_file, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("GSE97537总基因数: %d\n", nrow(deg_df)))

# 标准化列名
colnames(deg_df)[colnames(deg_df) == "logFC"] <- "log2FC"
colnames(deg_df)[colnames(deg_df) == "adj.P.Val"] <- "adj.P"

# 去除空Gene.symbol的记录
deg_df <- deg_df[deg_df$Gene.symbol != "" & !is.na(deg_df$Gene.symbol), ]
cat(sprintf("有基因符号的记录: %d\n", nrow(deg_df)))

# ============================================================================
# 分层筛选策略
# ============================================================================

# 第一层：全基因组DEGs
# 标准：|log₂FC| > 0.585 (FC>1.5)，adj.P < 0.05
layer1_degs <- deg_df[abs(deg_df$log2FC) > 0.585 & deg_df$adj.P < 0.05, ]
layer1_genes <- unique(layer1_degs$Gene.symbol)
cat(sprintf("第一层（全基因组DEGs）: %d 个基因 (|log2FC|>0.585, adj.P<0.05)\n", length(layer1_genes)))

# 第二层：铜死亡相关基因（核心机制层）
cat("\n=== 铜死亡基因分层筛选 ===\n")

# 铜死亡基因（使用大鼠基因名格式）
cuproptosis_genes <- c("Fdx1", "Lias", "Lipt1", "Dld", "Dlat", "Pdha1", "Pdhb", 
                       "Mtf1", "Gls", "Cdkn2a", "Slc31a1", "Atp7a", "Atp7b", 
                       "Sod1", "Ccs", "Cox17", "Sco1", "Sco2", "Cycs")

# 1. 核心调控基因（Fdx1, Lias, Slc31a1）
core_regulators <- c("Fdx1", "Lias", "Slc31a1")
core_regulator_degs <- deg_df[deg_df$Gene.symbol %in% core_regulators & 
                                (abs(deg_df$log2FC) > 0.3 | deg_df$adj.P < 0.1), ]
cat(sprintf("  核心调控基因: %d 个\n", nrow(core_regulator_degs)))

# 2. 铜伴侣蛋白
copper_chaperones <- c("Atox1", "Ccs", "Cox17")
chaperone_degs <- deg_df[deg_df$Gene.symbol %in% copper_chaperones & 
                           (abs(deg_df$log2FC) > 0.3 | deg_df$adj.P < 0.1), ]
cat(sprintf("  铜伴侣蛋白: %d 个\n", nrow(chaperone_degs)))

# 3. 铜离子转运体（不限制FC，确保捕获）
copper_transporters <- c("Atp7a", "Atp7b", "Slc31a1")
transporter_degs <- deg_df[deg_df$Gene.symbol %in% copper_transporters, ]
cat(sprintf("  铜离子转运体: %d 个\n", nrow(transporter_degs)))

layer2_genes <- unique(c(core_regulator_degs$Gene.symbol, 
                         chaperone_degs$Gene.symbol, 
                         transporter_degs$Gene.symbol))
cat(sprintf("第二层（铜死亡基因）: %d 个基因\n", length(layer2_genes)))

# 第三层：石竹烯（BCP）靶点基因
cat("\n=== 石竹烯（BCP）靶点基因 ===\n")

bcp_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/石竹烯 人.txt"
bcp_targets <- readLines(bcp_file)
bcp_targets <- trimws(bcp_targets[bcp_targets != ""])

# BCP靶点筛选：放宽标准以获取更多基因
bcp_in_deg <- deg_df[deg_df$Gene.symbol %in% bcp_targets & 
                       (abs(deg_df$log2FC) > 0.3 | deg_df$adj.P < 0.1), ]
layer3_genes <- unique(bcp_in_deg$Gene.symbol)
cat(sprintf("第三层（BCP靶点）: %d 个基因\n", length(layer3_genes)))

# BCP核心轴基因（使用大鼠基因名格式，强制包含）
bcp_core_axis_human <- c("AGER", "NFKB1", "FDX1", "STAT3", "TGFB1", "IL6", "PTGS2", "TLR4")
bcp_core_axis <- c("Ager", "Nfkb1", "Fdx1", "Stat3", "Tgfb1", "Il6", "Ptgs2", "Tlr4")
bcp_core_in_deg <- bcp_core_axis[bcp_core_axis %in% deg_df$Gene.symbol]
cat(sprintf("BCP核心轴基因在GSE97537中: %s\n", paste(bcp_core_in_deg, collapse = ", ")))

# ============================================================================
# 整合基因并筛选在Bulk数据中存在的基因
# ============================================================================

cat("\n=== 整合基因列表用于PC网络构建 ===\n")

# 合并所有层级基因
all_pc_genes <- unique(c(layer1_genes, layer2_genes, layer3_genes, bcp_core_axis))
cat(sprintf("整合后PC网络基因总数: %d\n", length(all_pc_genes)))

# 筛选在Bulk数据中存在的基因（大鼠数据，直接匹配）
available_genes <- all_pc_genes[all_pc_genes %in% rownames(bulk_expr)]
cat(sprintf("在GSE163614 Bulk数据中可映射的基因: %d\n", length(available_genes)))

# 检查BCP核心基因是否在Bulk中
bcp_core_in_bulk <- bcp_core_axis[bcp_core_axis %in% rownames(bulk_expr)]
cat(sprintf("BCP核心基因在Bulk中: %s\n", paste(bcp_core_in_bulk, collapse = ", ")))

# 强制添加BCP核心基因（即使不在DEG列表中，也要检查它们在Bulk中的表达）
bcp_core_to_add <- bcp_core_axis[bcp_core_axis %in% rownames(bulk_expr)]
if (length(bcp_core_to_add) > 0) {
  cat(sprintf("将强制添加%d个BCP核心基因到PC网络\n", length(bcp_core_to_add)))
}

# 如果基因太少，从第一层扩展
if (length(available_genes) < 20) {
  cat("基因数量太少，从第一层DEGs扩展...\n")
  # 按|log2FC|排序，取前50个
  layer1_sorted <- layer1_degs[order(abs(layer1_degs$log2FC), decreasing = TRUE), ]
  additional_genes <- head(layer1_sorted$Gene.symbol, 50)
  available_genes <- unique(c(available_genes, additional_genes))
  available_genes <- available_genes[available_genes %in% rownames(bulk_expr)]
  cat(sprintf("扩展后可用的基因: %d\n", length(available_genes)))
}

# 确保BCP核心基因被包含（强制添加）
available_genes <- unique(c(bcp_core_to_add, available_genes))
cat(sprintf("强制添加BCP核心基因后: %d 个基因\n", length(available_genes)))

# 限制基因数量（PC算法适合50-200个基因）
if (length(available_genes) > 150) {
  cat("基因数量超过150，按优先级筛选...\n")
  # 优先级：BCP核心 > 铜死亡 > BCP其他 > 第一层DEGs
  priority_order <- c(bcp_core_axis, cuproptosis_genes, layer1_genes)
  priority_genes <- priority_order[priority_order %in% rownames(bulk_expr)]
  available_genes <- unique(c(priority_genes, available_genes))[1:150]
  cat(sprintf("筛选后PC网络基因: %d\n", length(available_genes)))
}

# 检查是否有可用基因
if (length(available_genes) == 0) {
  stop("错误: 没有可用的基因进行PC分析！")
}

# 提取表达矩阵
pc_data <- t(bulk_expr[available_genes, ])
rownames(pc_data) <- colnames(bulk_expr)

cat(sprintf("\nPC分析数据矩阵: %d 样本 x %d 基因\n", nrow(pc_data), ncol(pc_data)))

# 保存基因列表
gene_info <- data.frame(
  Gene = available_genes,
  In_Layer1 = available_genes %in% layer1_genes,
  In_Layer2 = available_genes %in% layer2_genes,
  In_Layer3 = available_genes %in% layer3_genes,
  Is_BCP_Core = available_genes %in% bcp_core_axis,
  Is_Cuproptosis = available_genes %in% cuproptosis_genes,
  stringsAsFactors = FALSE
)

# 添加DEG统计
gene_info$log2FC <- deg_df$log2FC[match(gene_info$Gene, deg_df$Gene.symbol)]
gene_info$adj_P <- deg_df$adj.P[match(gene_info$Gene, deg_df$Gene.symbol)]

write.table(gene_info, "PC_network_genes_GSE97537.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# ============================================================================
# 第三部分：PC算法因果网络推断
# ============================================================================

cat("\n=== 第三步：PC算法因果网络推断 ===\n")

# 数据标准化
pc_data_scaled <- scale(pc_data)

# PC算法参数
alpha_level <- 0.1  # 提高alpha以获取更多边

suffStat <- list(C = cor(pc_data_scaled), n = nrow(pc_data_scaled))
pc_fit <- pc(suffStat,
             indepTest = gaussCItest,
             labels = colnames(pc_data_scaled),
             alpha = alpha_level,
             verbose = FALSE)

cat(sprintf("PC算法完成，alpha=%.2f\n", alpha_level))

# ============================================================================
# 第四部分：网络提取与可视化
# ============================================================================

cat("\n=== 第四步：提取因果网络边 ===\n")

# 提取网络边
net_edges <- as(pc_fit@graph, "matrix")
edge_list <- data.frame(
  from = character(),
  to = character(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(net_edges)) {
  for (j in 1:ncol(net_edges)) {
    if (net_edges[i, j] != 0 && i != j) {
      edge_list <- rbind(edge_list, data.frame(
        from = rownames(net_edges)[i],
        to = colnames(net_edges)[j],
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 去重
edge_list <- edge_list[!duplicated(t(apply(edge_list[, 1:2], 1, sort))), ]
cat(sprintf("网络边总数: %d\n", nrow(edge_list)))

write.table(edge_list, "PC_network_edges_GSE97537.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# ============================================================================
# 第五部分：RAGE→NFKB1→FDX1轴验证
# ============================================================================

cat("\n=== 第五步：验证RAGE→NFKB1→FDX1调控轴 ===\n")

# 检查基因在网络中（使用大鼠基因名格式）
axis_genes <- c("Ager", "Nfkb1", "Fdx1")
axis_status <- data.frame(
  Gene = axis_genes,
  In_Network = axis_genes %in% available_genes,
  Is_Source = axis_genes %in% edge_list$from,
  Is_Target = axis_genes %in% edge_list$to,
  stringsAsFactors = FALSE
)
print(axis_status)

# 查找直接边（使用大鼠基因名格式）
rage_to_nfkb <- edge_list[edge_list$from == "Ager" & edge_list$to == "Nfkb1", ]
nfkb_to_fdx1 <- edge_list[edge_list$from == "Nfkb1" & edge_list$to == "Fdx1", ]

cat("\n查找RAGE→NFKB1→FDX1路径:\n")
if (nrow(rage_to_nfkb) > 0) {
  cat("✓ 发现 Ager → Nfkb1 因果边\n")
} else {
  cat("✗ 未发现 Ager → Nfkb1 直接因果边\n")
}

if (nrow(nfkb_to_fdx1) > 0) {
  cat("✓ 发现 Nfkb1 → Fdx1 因果边\n")
} else {
  cat("✗ 未发现 Nfkb1 → Fdx1 直接因果边\n")
}

# 查找间接路径
if (nrow(edge_list) > 0) {
  cat("\n所有网络边:\n")
  for (i in 1:min(nrow(edge_list), 20)) {
    cat(sprintf("  %s → %s\n", edge_list$from[i], edge_list$to[i]))
  }
}

# ============================================================================
# 第六部分：网络可视化
# ============================================================================

cat("\n=== 第六步：生成网络可视化 ===\n")

if (nrow(edge_list) > 0) {
  g <- graph_from_data_frame(edge_list, directed = TRUE)
  
  # 设置节点颜色
  V(g)$color <- "lightblue"
  bcp_nodes <- which(V(g)$name %in% bcp_core_axis)
  V(g)$color[bcp_nodes] <- "red"
  cupro_nodes <- which(V(g)$name %in% cuproptosis_genes)
  V(g)$color[cupro_nodes] <- "gold"
  both_nodes <- intersect(bcp_nodes, cupro_nodes)
  V(g)$color[both_nodes] <- "orange"
  
  # 保存PDF
  pdf("PC_causal_network_GSE97537.pdf", width = 14, height = 10)
  set.seed(42)
  plot(g, 
       layout = layout_with_fr(g, niter = 1000),
       vertex.size = 8,
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       edge.arrow.size = 0.3,
       edge.curved = 0.2,
       main = "PC Causal Network (GSE97537 DEGs + GSE163614 Expression)\nRed: BCP Core, Gold: Cuproptosis")
  dev.off()
  
  cat("✓ 网络图已保存: PC_causal_network_GSE97537.pdf\n")
} else {
  cat("警告: 网络边为空，无法生成可视化\n")
}

# ============================================================================
# 第七部分：保存结果
# ============================================================================

cat("\n=== 分析完成，生成结果摘要 ===\n")

summary_stats <- data.frame(
  Metric = c(
    "Total Genes in PC Network",
    "Layer 1 (DEGs |log2FC|>0.585)",
    "Layer 2 (Cuproptosis)",
    "Layer 3 (BCP Targets)",
    "Total Network Edges",
    "Ager (RAGE) in Network",
    "Nfkb1 in Network",
    "Fdx1 in Network"
  ),
  Count = c(
    length(available_genes),
    sum(available_genes %in% layer1_genes),
    sum(available_genes %in% layer2_genes),
    sum(available_genes %in% layer3_genes),
    nrow(edge_list),
    ifelse("Ager" %in% available_genes, "Yes", "No"),
    ifelse("Nfkb1" %in% available_genes, "Yes", "No"),
    ifelse("Fdx1" %in% available_genes, "Yes", "No")
  ),
  stringsAsFactors = FALSE
)

print(summary_stats)
write.table(summary_stats, "PC_analysis_summary_GSE97537.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("\n✓ PC因果网络分析（GSE97537 DEGs）完成!\n")
cat("\n生成文件:\n")
cat("1. PC_network_genes_GSE97537.txt - 基因详细信息\n")
cat("2. PC_network_edges_GSE97537.txt - 网络边列表\n")
cat("3. PC_analysis_summary_GSE97537.txt - 分析摘要\n")
cat("4. PC_causal_network_GSE97537.pdf - 网络图\n")
