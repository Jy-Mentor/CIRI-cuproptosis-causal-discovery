# ============================================================================
# Step 1: PC因果网络分析 - 分层筛选策略
# 研究问题: BCP是否通过RAGE→NFKB1→FDX1轴发挥作用
# 数据集: GSE163614 (大鼠MCAO Bulk RNA-seq)
# ============================================================================

# 清除环境
rm(list = ls())
gc()

# 设置工作目录
setwd("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙")

# 设置CRAN镜像
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# 加载必要的包（假设已安装）
library(readxl)
library(pcalg)
library(igraph)

# 设置随机种子确保可复现
set.seed(42)

# ============================================================================
# 第一部分：数据读取与基因筛选（分层筛选策略）
# ============================================================================

cat("=== 第一步：读取Bulk RNA-seq数据 ===\n")

# 读取GSE163614 Bulk数据
bulk_file <- "C:/Users/Jy-Mentor-7/Downloads/GSE163614_mRNA_Expression_Profiling.xlsx"
if (!file.exists(bulk_file)) {
  stop("错误: 找不到GSE163614数据文件!")
}

# 读取数据（跳过前5行表头）
bulk_df <- read_excel(bulk_file, skip = 5)

# 设置列名
colnames(bulk_df)[1] <- "gene_symbol"
# 样本列：MCAO 1-3 和 Sham 1-3
colnames(bulk_df)[7:12] <- c("MCAO1", "MCAO2", "MCAO3", "Sham1", "Sham2", "Sham3")

# 提取表达矩阵
bulk_expr <- as.matrix(bulk_df[, 7:12])
rownames(bulk_expr) <- make.names(bulk_df$gene_symbol, unique = TRUE)

# 调试: 检查rownames
if (is.null(rownames(bulk_expr)) || length(rownames(bulk_expr)) == 0) {
  cat("警告: rownames为空，使用基因符号列直接设置\n")
  rownames(bulk_expr) <- bulk_df$gene_symbol
}

bulk_expr <- apply(bulk_expr, 2, as.numeric)

# 重新设置rownames（因为apply会移除rownames）
rownames(bulk_expr) <- make.names(bulk_df$gene_symbol, unique = TRUE)

cat(sprintf("Bulk数据: %d 基因 x %d 样本\n", nrow(bulk_expr), ncol(bulk_expr)))
cat("Bulk数据基因名示例:", head(rownames(bulk_expr), 10), "\n")

# 读取基因映射文件（人→大鼠）
cat("\n读取基因映射文件...\n")
mapping_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt"
if (file.exists(mapping_file)) {
  # 跳过前53行注释，读取数据
  mapping_df <- read.delim(mapping_file, sep = "\t", stringsAsFactors = FALSE, 
                           header = FALSE, skip = 53)
  # 设置列名
  colnames(mapping_df) <- c("RAT_GENE_SYMBOL", "RAT_GENE_RGD_ID", "RAT_GENE_NCBI_GENE_ID",
                            "HUMAN_ORTHOLOG_SYMBOL", "HUMAN_ORTHOLOG_RGD_ID", "HUMAN_ORTHOLOG_NCBI_GENE_ID",
                            "HUMAN_ORTHOLOG_SOURCE", "MOUSE_ORTHOLOG_SYMBOL", "MOUSE_ORTHOLOG_RGD_ID",
                            "MOUSE_ORTHOLOG_NCBI_GENE_ID", "MOUSE_ORTHOLOG_MGI_ID", 
                            "MOUSE_ORTHOLOG_SOURCE", "HUMAN_ORTHOLOG_HGNC_ID")
  cat(sprintf("映射库: %d 条记录\n", nrow(mapping_df)))
  
  # 创建人类→大鼠映射表（处理多个人类基因映射到同一个大鼠基因的情况）
  # 只保留一对一映射
  mapping_df_clean <- mapping_df[mapping_df$HUMAN_ORTHOLOG_SYMBOL != "" & 
                                   mapping_df$RAT_GENE_SYMBOL != "" &
                                   !grepl("\\|", mapping_df$HUMAN_ORTHOLOG_SYMBOL), ]
  human_to_rat <- setNames(mapping_df_clean$RAT_GENE_SYMBOL, mapping_df_clean$HUMAN_ORTHOLOG_SYMBOL)
  cat(sprintf("有效一对一映射: %d 条\n", length(human_to_rat)))
} else {
  cat("警告: 未找到映射文件，使用大小写转换规则\n")
  human_to_rat <- NULL
}

# 基因名转换函数：人类→大鼠
convert_human_to_rat <- function(human_genes) {
  if (!is.null(human_to_rat)) {
    # 使用映射表 - 逐个查找
    rat_genes <- c()
    for (g in human_genes) {
      if (g %in% names(human_to_rat)) {
        rat_genes <- c(rat_genes, human_to_rat[g])
      }
    }
    rat_genes <- unique(rat_genes)
  } else {
    # 使用大小写转换（首字母大写）
    rat_genes <- sapply(human_genes, function(g) {
      if (nchar(g) > 0) {
        paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))))
      } else {
        g
      }
    })
    names(rat_genes) <- NULL
  }
  return(rat_genes)
}

# ============================================================================
# 第一层筛选：全基因组DEGs（构建背景网络）
# 标准：|log₂FC| > 1，adj.P < 0.05
# ============================================================================

cat("\n=== 第二层筛选：读取GEO2R DEG结果并应用分层筛选 ===\n")

# 读取用户提供的GEO2R结果 (使用GSE97537，包含adj.P.Val)
deg_file <- "C:/Users/Jy-Mentor-7/Downloads/GSE97537.top.table (1).tsv"
if (!file.exists(deg_file)) {
  stop("错误: 找不到GEO2R结果文件!")
}

deg_df <- read.delim(deg_file, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("GEO2R结果: %d 个基因\n", nrow(deg_df)))

# 检查列名
cat("GEO2R列名:", colnames(deg_df), "\n")

# 标准化列名（支持不同GEO2R输出格式）
if ("logFC" %in% colnames(deg_df) && !("log2FC" %in% colnames(deg_df))) {
  deg_df$log2FC <- deg_df$logFC
}
if ("adj.P.Val" %in% colnames(deg_df) && !("adj.P" %in% colnames(deg_df))) {
  deg_df$adj.P <- deg_df$adj.P.Val
}

# 检查是否有P值列
has_pvalue <- "adj.P" %in% colnames(deg_df)
if (!has_pvalue) {
  cat("警告: DEG文件缺少adj.P列，将仅使用log2FC筛选\n")
}

# 第一层筛选：全基因组DEGs
# 标准调整为：|log₂FC| > 0.5，adj.P < 0.05（放宽以获取更多基因用于PC网络）
if (has_pvalue) {
  layer1_degs <- deg_df[abs(deg_df$log2FC) > 0.5 & deg_df$adj.P < 0.05, ]
} else {
  layer1_degs <- deg_df[abs(deg_df$log2FC) > 0.5, ]
}
cat(sprintf("第一层（全基因组DEGs）: %d 个基因 (|log2FC|>0.5, adj.P<0.05)\n", nrow(layer1_degs)))

# 提取第一层基因列表
layer1_genes <- unique(na.omit(layer1_degs$Gene.symbol))

# ============================================================================
# 第二层筛选：铜死亡相关基因（核心机制层）
# ============================================================================

cat("\n=== 第二层：铜死亡基因分层筛选 ===\n")

# 19个核心铜死亡基因
cuproptosis_genes <- c("FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB", 
                       "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7A", "ATP7B", 
                       "SOD1", "CCS", "COX17", "SCO1", "SCO2", "CYCS")

# 1. 核心调控基因（FDX1, LIAS, SLC31A1）
# 筛选标准放宽：|log2FC| > 0.3 或 adj.P < 0.1（捕获更多潜在基因）
core_regulators <- c("FDX1", "LIAS", "SLC31A1")
if (has_pvalue) {
  core_regulator_degs <- deg_df[deg_df$Gene.symbol %in% core_regulators & 
                                  (abs(deg_df$log2FC) > 0.3 | deg_df$adj.P < 0.1), ]
} else {
  core_regulator_degs <- deg_df[deg_df$Gene.symbol %in% core_regulators & 
                                  abs(deg_df$log2FC) > 0.3, ]
}
cat(sprintf("  核心调控基因: %d 个 (|log2FC|>0.3 或 adj.P<0.1)\n", nrow(core_regulator_degs)))

# 2. 铜伴侣蛋白（ATOX1, CCS, COX17）
# 筛选标准放宽：|log2FC| > 0.3 或 adj.P < 0.1
copper_chaperones <- c("ATOX1", "CCS", "COX17")
if (has_pvalue) {
  chaperone_degs <- deg_df[deg_df$Gene.symbol %in% copper_chaperones & 
                             (abs(deg_df$log2FC) > 0.3 | deg_df$adj.P < 0.1), ]
} else {
  chaperone_degs <- deg_df[deg_df$Gene.symbol %in% copper_chaperones & 
                             abs(deg_df$log2FC) > 0.3, ]
}
cat(sprintf("  铜伴侣蛋白: %d 个 (|log2FC|>0.3 或 adj.P<0.1)\n", nrow(chaperone_degs)))

# 3. 铜离子转运体（ATP7A, ATP7B, SLC31A1）
# 不设筛选限制，确保捕获
# 注意：ATP7B可能不存在于映射库，只列出存在的基因
copper_transporters <- c("ATP7A", "SLC31A1")
transporter_degs <- deg_df[deg_df$Gene.symbol %in% copper_transporters, ]
cat(sprintf("  铜离子转运体: %d 个 (捕获所有)\n", nrow(transporter_degs)))

# 合并所有铜死亡基因
layer2_genes <- unique(c(core_regulator_degs$Gene.symbol, 
                         chaperone_degs$Gene.symbol, 
                         transporter_degs$Gene.symbol))
cat(sprintf("第二层（铜死亡基因）: %d 个基因\n", length(layer2_genes)))

# ============================================================================
# 第三层：石竹烯（BCP）靶点基因
# ============================================================================

cat("\n=== 第三层：石竹烯（BCP）靶点基因 ===\n")

# 读取石竹烯靶点
bcp_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/石竹烯 人.txt"
if (!file.exists(bcp_file)) {
  stop("错误: 找不到石竹烯靶点文件!")
}

bcp_targets <- readLines(bcp_file)
bcp_targets <- trimws(bcp_targets[bcp_targets != ""])
cat(sprintf("石竹烯靶点: %d 个基因\n", length(bcp_targets)))

# 筛选在DEG列表中的BCP靶点
if (has_pvalue) {
  bcp_in_deg <- deg_df[deg_df$Gene.symbol %in% bcp_targets & 
                         abs(deg_df$log2FC) > 0.5 & 
                         deg_df$adj.P < 0.1, ]
} else {
  bcp_in_deg <- deg_df[deg_df$Gene.symbol %in% bcp_targets & 
                         abs(deg_df$log2FC) > 0.5, ]
}
layer3_genes <- unique(bcp_in_deg$Gene.symbol)
cat(sprintf("第三层（BCP靶点）: %d 个基因 (|log2FC|>0.5, adj.P<0.1)\n", length(layer3_genes)))

# 重点关注的核心BCP靶点（RAGE→NFKB1→FDX1轴）
bcp_core_axis <- c("AGER", "NFKB1", "FDX1", "STAT3", "TGFB1", "IL6", "PTGS2", "TLR4")
bcp_core_in_data <- bcp_core_axis[bcp_core_axis %in% deg_df$Gene.symbol]
cat(sprintf("BCP核心轴基因在数据中: %s\n", paste(bcp_core_in_data, collapse = ", ")))

# ============================================================================
# 整合所有基因用于PC网络构建
# ============================================================================

cat("\n=== 整合基因列表用于PC网络构建 ===\n")

# 合并所有层级基因（人类基因名）
all_pc_genes_human <- unique(c(layer1_genes, layer2_genes, layer3_genes, bcp_core_axis))
cat(sprintf("整合后PC网络基因总数（人类）: %d\n", length(all_pc_genes_human)))

# 转换为大鼠基因名
all_pc_genes_rat <- convert_human_to_rat(all_pc_genes_human)
cat(sprintf("转换为大鼠基因名: %d 个\n", length(all_pc_genes_rat)))
if (length(all_pc_genes_rat) > 0) {
  cat("转换后的基因示例:", head(all_pc_genes_rat, 10), "\n")
}

# 筛选在Bulk数据中存在的基因
available_genes <- all_pc_genes_rat[all_pc_genes_rat %in% rownames(bulk_expr)]
cat(sprintf("在Bulk数据中可映射的基因: %d\n", length(available_genes)))

# 如果没有找到基因，尝试大小写转换
if (length(available_genes) == 0) {
  cat("尝试大小写转换...\n")
  all_pc_genes_rat_simple <- sapply(all_pc_genes_human, function(g) {
    if (nchar(g) > 0) {
      paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))))
    } else { g }
  })
  available_genes <- all_pc_genes_rat_simple[all_pc_genes_rat_simple %in% rownames(bulk_expr)]
  available_genes <- unique(available_genes)
  cat(sprintf("大小写转换后在Bulk数据中可映射的基因: %d\n", length(available_genes)))
}

# 如果基因太多，优先选择核心基因
if (length(available_genes) > 500) {
  cat("基因数量超过500，优先选择核心基因...\n")
  # 优先级：BCP核心轴 > 铜死亡核心调控 > 第一层DEGs
  priority_genes <- c(bcp_core_axis[bcp_core_axis %in% rownames(bulk_expr)],
                      core_regulators[core_regulators %in% rownames(bulk_expr)],
                      layer1_genes[layer1_genes %in% rownames(bulk_expr)])
  available_genes <- unique(priority_genes)[1:min(500, length(unique(priority_genes)))]
  cat(sprintf("筛选后PC网络基因: %d\n", length(available_genes)))
}

# 检查是否有可用基因
if (length(available_genes) == 0) {
  stop("错误: 没有可用的基因进行PC分析！请检查DEG文件格式和基因名匹配。")
}

# 提取表达矩阵
pc_data <- t(bulk_expr[available_genes, ])
rownames(pc_data) <- colnames(bulk_expr)

cat(sprintf("\nPC分析数据矩阵: %d 样本 x %d 基因\n", nrow(pc_data), ncol(pc_data)))

# 保存基因列表用于后续分析
write.table(data.frame(
  Gene = available_genes,
  Layer = ifelse(available_genes %in% layer1_genes, "Layer1_DEG", 
                 ifelse(available_genes %in% layer2_genes, "Layer2_Cuproptosis",
                        ifelse(available_genes %in% bcp_targets, "Layer3_BCP", "Other")))
), "PC_network_genes_stratified.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# ============================================================================
# 第二部分：PC算法因果网络推断
# ============================================================================

cat("\n=== 第二步：PC算法因果网络推断 ===\n")

# 数据标准化（z-score）
pc_data_scaled <- scale(pc_data)

# 定义样本分组（用于条件独立检验）
# MCAO = 处理组, Sham = 对照组
groups <- factor(c("MCAO", "MCAO", "MCAO", "Sham", "Sham", "Sham"))

# 使用Fisher's Z检验进行条件独立检验
suffStat <- list(C = cor(pc_data_scaled), n = nrow(pc_data_scaled))

# PC算法参数
# alpha: 显著性水平（控制I类错误）
# indepTest: 独立性检验方法
# labels: 变量名
# verbose: 输出详细信息

# PC算法参数调整：提高alpha值以增加网络边
# alpha从0.05提高到0.1，平衡I类错误和网络连接性
alpha_level <- 0.1
pc_fit <- pc(suffStat,
             indepTest = gaussCItest,
             labels = colnames(pc_data_scaled),
             alpha = alpha_level,
             verbose = FALSE)

cat(sprintf("PC算法完成，alpha=%.2f (已优化以获取更多网络边)\n", alpha_level))

# ============================================================================
# 第三部分：网络提取与可视化
# ============================================================================

cat("\n=== 第三步：提取因果网络边 ===\n")

# 提取网络边
net_edges <- as(pc_fit@graph, "matrix")
edge_list <- data.frame(
  from = character(),
  to = character(),
  edge_type = character(),
  stringsAsFactors = FALSE
)

# 遍历邻接矩阵提取边
for (i in 1:nrow(net_edges)) {
  for (j in 1:ncol(net_edges)) {
    if (net_edges[i, j] != 0) {
      edge_list <- rbind(edge_list, data.frame(
        from = rownames(net_edges)[i],
        to = colnames(net_edges)[j],
        edge_type = ifelse(net_edges[i, j] == 1 && net_edges[j, i] == 1, "undirected", "directed"),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 去重（只保留一个方向的有向边）
edge_list <- edge_list[!duplicated(t(apply(edge_list[, 1:2], 1, sort))), ]

cat(sprintf("网络边总数: %d\n", nrow(edge_list)))

# 保存边列表
write.table(edge_list, "PC_network_edges_stratified.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# ============================================================================
# 第四部分：RAGE→NFKB1→FDX1轴验证
# ============================================================================

cat("\n=== 第四步：验证RAGE→NFKB1→FDX1调控轴 ===\n")

# 基因名映射（大鼠基因名可能与人类不同）
# 在Bulk数据中查找这些基因
check_axis_genes <- function(genes, edge_df) {
  results <- data.frame(
    Gene = genes,
    In_Network = genes %in% c(edge_df$from, edge_df$to),
    Is_Source = genes %in% edge_df$from,
    Is_Target = genes %in% edge_df$to,
    stringsAsFactors = FALSE
  )
  return(results)
}

axis_genes <- c("AGER", "NFKB1", "FDX1")
axis_status <- check_axis_genes(axis_genes, edge_list)
print(axis_status)

# 查找特定路径
cat("\n查找RAGE→NFKB1→FDX1路径:\n")

# 直接边
rage_to_nfkb <- edge_list[edge_list$from == "AGER" & edge_list$to == "NFKB1", ]
nfkb_to_fdx1 <- edge_list[edge_list$from == "NFKB1" & edge_list$to == "FDX1", ]

if (nrow(rage_to_nfkb) > 0) {
  cat("✓ 发现 AGER → NFKB1 因果边\n")
} else {
  cat("✗ 未发现 AGER → NFKB1 直接因果边\n")
}

if (nrow(nfkb_to_fdx1) > 0) {
  cat("✓ 发现 NFKB1 → FDX1 因果边\n")
} else {
  cat("✗ 未发现 NFKB1 → FDX1 直接因果边\n")
}

# 查找间接路径（两步）
cat("\n查找间接调控路径:\n")
for (gene in axis_genes) {
  # 查找以该基因为起点的边
  outgoing <- edge_list[edge_list$from == gene, ]
  if (nrow(outgoing) > 0) {
    cat(sprintf("  %s 调控: %s\n", gene, paste(outgoing$to, collapse = ", ")))
  }
  
  # 查找以该基因为终点的边
  incoming <- edge_list[edge_list$to == gene, ]
  if (nrow(incoming) > 0) {
    cat(sprintf("  %s 被调控于: %s\n", gene, paste(incoming$from, collapse = ", ")))
  }
}

# 保存轴验证结果
axis_validation <- list(
  axis_status = axis_status,
  direct_edges = list(rage_to_nfkb = rage_to_nfkb, nfkb_to_fdx1 = nfkb_to_fdx1),
  all_edges = edge_list
)

# ============================================================================
# 第五部分：网络可视化（PDF格式，文字可编辑）
# ============================================================================

cat("\n=== 第五步：生成网络可视化 ===\n")

# 使用igraph进行可视化
library(igraph)

# 创建图对象
if (nrow(edge_list) > 0) {
  g <- graph_from_data_frame(edge_list[, 1:2], directed = TRUE)
  
  # 设置节点颜色（按层级）
  node_colors <- rep("lightblue", vcount(g))
  V(g)$color <- node_colors
  
  # 标记BCP核心轴基因
  bcp_nodes <- which(V(g)$name %in% bcp_core_axis)
  V(g)$color[bcp_nodes] <- "red"
  
  # 标记铜死亡基因
  cupro_nodes <- which(V(g)$name %in% cuproptosis_genes)
  V(g)$color[cupro_nodes] <- "gold"
  
  # 标记同时是BCP靶点和铜死亡基因
  both_nodes <- intersect(bcp_nodes, cupro_nodes)
  V(g)$color[both_nodes] <- "orange"
  
  # 保存PDF（文字可编辑）
  pdf("PC_causal_network_stratified.pdf", width = 16, height = 12)
  
  # 使用fr布局
  set.seed(42)
  layout <- layout_with_fr(g, niter = 1000)
  
  # 绘制网络
  plot(g, 
       layout = layout,
       vertex.size = 8,
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       vertex.frame.color = "gray50",
       edge.arrow.size = 0.3,
       edge.curved = 0.2,
       main = "PC Causal Network - Stratified Screening\n(Red: BCP targets, Gold: Cuproptosis, Orange: Both)")
  
  dev.off()
  
  cat("✓ 网络图已保存: PC_causal_network_stratified.pdf\n")
  
  # 生成子网络（仅BCP核心轴相关基因）
  bcp_subgraph <- induced_subgraph(g, which(V(g)$name %in% c(bcp_core_axis, cuproptosis_genes)))
  
  if (vcount(bcp_subgraph) > 0) {
    pdf("PC_network_BCP_axis_stratified.pdf", width = 12, height = 10)
    
    set.seed(42)
    layout_sub <- layout_with_fr(bcp_subgraph, niter = 1000)
    
    plot(bcp_subgraph,
         layout = layout_sub,
         vertex.size = 12,
         vertex.label.cex = 0.8,
         vertex.label.color = "black",
         vertex.frame.color = "gray50",
         edge.arrow.size = 0.4,
         edge.curved = 0.2,
         main = "BCP Core Axis & Cuproptosis Subnetwork")
    
    dev.off()
    
    cat("✓ BCP轴子网络图已保存: PC_network_BCP_axis_stratified.pdf\n")
  }
} else {
  cat("警告: 网络边为空，无法生成可视化\n")
}

# ============================================================================
# 第六部分：保存结果与统计
# ============================================================================

cat("\n=== 分析完成，生成结果摘要 ===\n")

# 生成结果摘要
summary_stats <- data.frame(
  Metric = c(
    "Total Genes in PC Network",
    "Layer 1 (DEGs |log2FC|>1)",
    "Layer 2 (Cuproptosis)",
    "Layer 3 (BCP Targets)",
    "Total Network Edges",
    "AGER in Network",
    "NFKB1 in Network",
    "FDX1 in Network"
  ),
  Count = c(
    length(available_genes),
    sum(available_genes %in% layer1_genes),
    sum(available_genes %in% layer2_genes),
    sum(available_genes %in% bcp_targets),
    nrow(edge_list),
    ifelse("AGER" %in% c(edge_list$from, edge_list$to), "Yes", "No"),
    ifelse("NFKB1" %in% c(edge_list$from, edge_list$to), "Yes", "No"),
    ifelse("FDX1" %in% c(edge_list$from, edge_list$to), "Yes", "No")
  ),
  stringsAsFactors = FALSE
)

print(summary_stats)

# 保存摘要
write.table(summary_stats, "PC_analysis_summary_stratified.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# 保存各层级基因详情
gene_details <- data.frame(
  Gene = available_genes,
  In_Layer1 = available_genes %in% layer1_genes,
  In_Layer2 = available_genes %in% layer2_genes,
  In_Layer3 = available_genes %in% bcp_targets,
  Is_BCP_Core = available_genes %in% bcp_core_axis,
  Is_Cuproptosis = available_genes %in% cuproptosis_genes,
  stringsAsFactors = FALSE
)

# 添加DEG统计信息
gene_details$log2FC <- deg_df$log2FC[match(gene_details$Gene, deg_df$Gene.symbol)]
if (has_pvalue) {
  gene_details$adj_P <- deg_df$adj.P[match(gene_details$Gene, deg_df$Gene.symbol)]
} else {
  gene_details$adj_P <- NA
}

write.table(gene_details, "PC_network_gene_details_stratified.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("\n=== 所有结果文件 ===\n")
cat("1. PC_network_edges_stratified.txt - 网络边列表\n")
cat("2. PC_network_genes_stratified.txt - 网络基因分层信息\n")
cat("3. PC_network_gene_details_stratified.txt - 基因详细信息\n")
cat("4. PC_analysis_summary_stratified.txt - 分析摘要\n")
cat("5. PC_causal_network_stratified.pdf - 完整网络图\n")
cat("6. PC_network_BCP_axis_stratified.pdf - BCP轴子网络图\n")

cat("\n✓ PC因果网络分析（分层筛选策略）完成!\n")
