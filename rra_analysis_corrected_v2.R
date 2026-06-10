# ============================================
# 修正版RRA分析：使用RobustRankAggreg包
# 只整合4种中心性（DC, BC, CC, EC），K_core用于预筛选
# ============================================

cat("正在加载必要的R包...\n")

# 安装和加载RobustRankAggreg
if (!require("RobustRankAggreg", quietly = TRUE)) {
  install.packages("RobustRankAggreg", repos = "https://cloud.r-project.org/")
}

# 检查包中的函数
cat("检查RobustRankAggreg包中的函数...\n")
print(ls("package:RobustRankAggreg"))

cat("\n加载RobustRankAggreg包...\n")
library(RobustRankAggreg)

# 设置路径
work_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
result_dir <- file.path(work_dir, "String_Network_Systematic_Analysis")

# 读取中心性数据
cat("\n读取中心性数据...\n")
centrality_df <- read.delim(file.path(result_dir, "03_centrality_measures.txt"), 
                            stringsAsFactors = FALSE)
cat(paste0("读取了 ", nrow(centrality_df), " 个节点的中心性数据\n"))

# ==================== 修正版RRA：只整合4种中心性 ====================
cat("\n========================================\n")
cat("步骤: 修正版RRA分析（只整合4种中心性）\n")
cat("========================================\n")

# 创建4个排序列表（按中心性值降序排列）
dc_list <- centrality_df$Node[order(-centrality_df$DC)]
bc_list <- centrality_df$Node[order(-centrality_df$BC)]
cc_list <- centrality_df$Node[order(-centrality_df$CC)]
ec_list <- centrality_df$Node[order(-centrality_df$EC)]

cat("\n各中心性列表长度:\n")
cat(paste0("  DC列表: ", length(dc_list), "\n"))
cat(paste0("  BC列表: ", length(bc_list), "\n"))
cat(paste0("  CC列表: ", length(cc_list), "\n"))
cat(paste0("  EC列表: ", length(ec_list), "\n"))

# 创建列表集合
glist <- list(DC=dc_list, BC=bc_list, CC=cc_list, EC=ec_list)

# 运行RRA分析 - 使用aggregateRanks函数
cat("\n运行RRA分析...\n")

if (exists("aggregateRanks")) {
  rra_result <- aggregateRanks(glist, N=nrow(centrality_df), method="RRA")
} else if (exists("rankAggreg")) {
  rra_result <- rankAggreg(glist, N=nrow(centrality_df), method="RRA")
} else {
  cat("错误: 找不到RRA函数\n")
  cat("可用的函数:\n")
  cat(paste0("  ", paste(ls("package:RobustRankAggreg"), collapse=", "), "\n"))
  stop("RRA函数不可用")
}

# 查看RRA结果结构
cat(paste0("\nRRA结果共 ", nrow(rra_result), " 个节点\n"))
cat("RRA结果列名:\n")
print(names(rra_result))

# ==================== 查看FDX1排名 ====================
cat("\n========================================\n")
cat("FDX1基因RRA分析结果\n")
cat("========================================\n")

if ("FDX1" %in% rra_result$Name) {
  fdx1_row <- rra_result[rra_result$Name == "FDX1", ]
  cat("\nFDX1在RRA结果中的信息:\n")
  print(fdx1_row)
  
  fdx1_rank <- which(rra_result$Name == "FDX1")
  cat(paste0("\nFDX1排名: ", fdx1_rank, " / ", nrow(rra_result), "\n"))
  
  if ("Score" %in% names(fdx1_row)) {
    cat(paste0("FDX1得分 (Score): ", fdx1_row$Score, "\n"))
  }
  if ("Pvalue" %in% names(fdx1_row)) {
    cat(paste0("FDX1 P-value: ", fdx1_row$Pvalue, "\n"))
  }
  if ("pvalue" %in% names(fdx1_row)) {
    cat(paste0("FDX1 p-value: ", fdx1_row$pvalue, "\n"))
  }
} else {
  cat("\n警告: FDX1不在RRA结果中\n")
  cat("检查数据中是否有FDX1...\n")
  cat(paste0("数据中FDX1存在: ", "FDX1" %in% centrality_df$Node, "\n"))
}

# ==================== 查看铜死亡基因整体排名 ====================
cat("\n========================================\n")
cat("铜死亡基因RRA分析结果\n")
cat("========================================\n")

# 铜死亡基因列表
cup_genes <- c("FDX1", "LIAS", "SLC31A1", "DLAT", "PDHB", "PDHX", 
               "GPX4", "CP", "ATP7A", "ATOX1", "HIF1A", "NFKB1")

cat(paste0("\n铜死亡基因列表 (", length(cup_genes), " 个):\n"))
print(cup_genes)

# 筛选铜死亡基因的RRA结果
cup_in_rra <- rra_result[rra_result$Name %in% cup_genes, ]

# 添加排名信息
cup_in_rra$Rank <- sapply(cup_in_rra$Name, function(x) which(rra_result$Name == x))

# 按排名排序
cup_in_rra <- cup_in_rra[order(cup_in_rra$Rank), ]

cat(paste0("\n在RRA结果中找到 ", nrow(cup_in_rra), " 个铜死亡基因\n"))
cat("\n铜死亡基因RRA排名:\n")
print(cup_in_rra)

# 查看未找到的基因
not_found <- setdiff(cup_genes, rra_result$Name)
if (length(not_found) > 0) {
  cat(paste0("\n未在RRA结果中找到的铜死亡基因 (", length(not_found), " 个):\n"))
  print(not_found)
  
  # 检查是否在原始数据中
  cat("\n检查原始数据中的铜死亡基因:\n")
  for (gene in cup_genes) {
    in_data <- gene %in% centrality_df$Node
    cat(paste0("  ", gene, ": ", ifelse(in_data, "在数据中", "不在数据中"), "\n"))
  }
}

# ==================== 保存完整RRA结果 ====================
cat("\n保存RRA结果...\n")

# 添加排名列
rra_result$Rank <- 1:nrow(rra_result)

# 重新排列列顺序
cols_to_keep <- c("Name", "Rank", "Score")
if ("Pvalue" %in% names(rra_result)) cols_to_keep <- c(cols_to_keep, "Pvalue")
if ("pvalue" %in% names(rra_result)) cols_to_keep <- c(cols_to_keep, "pvalue")

rra_output <- rra_result[, cols_to_keep]

write.table(rra_output, 
            file = file.path(result_dir, "04_rra_corrected_all_nodes.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 保存铜死亡基因结果
if (nrow(cup_in_rra) > 0) {
  write.table(cup_in_rra, 
              file = file.path(result_dir, "04_rra_copper_genes.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# 保存Top 50 Hub节点
top50 <- head(rra_output, 50)
write.table(top50, 
            file = file.path(result_dir, "04_rra_top50_hub_nodes.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n========================================\n")
cat("         修正版RRA分析完成\n")
cat("========================================\n")
cat("\n生成文件:\n")
cat("  - 04_rra_corrected_all_nodes.txt (所有节点的RRA排名)\n")
cat("  - 04_rra_copper_genes.txt (铜死亡基因RRA排名)\n")
cat("  - 04_rra_top50_hub_nodes.txt (Top 50 Hub节点)\n")

# 输出Top 20 Hub节点
cat("\n========================================\n")
cat("Top 20 Hub节点 (修正版RRA)\n")
cat("========================================\n")
top20 <- head(rra_output, 20)
print(top20)
