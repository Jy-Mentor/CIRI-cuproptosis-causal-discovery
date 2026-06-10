# 展示虚拟敲除候选基因
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))

library(Seurat)
library(Matrix)

cat("=== 虚拟敲除候选基因展示 ===\n\n")

result_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/causal_analysis_results"

# 1. 读取PC Hub基因
cat("【PC因果网络识别的Hub基因】\n")
pc_result <- readRDS(file.path(result_dir, "PC_analysis_result.rds"))
hub_genes <- pc_result$hub_genes

cat("Top 15 Hub基因:\n")
hub_df <- data.frame(
  rank = 1:length(hub_genes),
  gene = hub_genes,
  degree = as.numeric(pc_result$node_degrees[hub_genes]),
  log2FC = pc_result$deg_res[hub_genes, "logFC"],
  FDR = pc_result$deg_res[hub_genes, "adj.P.Val"]
)
print(hub_df)

# 2. 读取SCISSOR结果中的表达差异
cat("\n【SCISSOR-like分析 - 单细胞中MCAO vs Sham差异】\n")
hub_expr <- read.csv(file.path(result_dir, "SCISSOR_results/hub_genes_expression_comparison.csv"))
hub_expr <- hub_expr[order(-abs(hub_expr$log2FC)), ]
cat("Hub基因在单细胞中的表达差异 (Top 10):\n")
print(head(hub_expr, 10))

# 3. 检查BCP轴基因
cat("\n【BCP轴基因在数据中的情况】\n")
bcp_genes <- c("Ager","Nfkb1","Fdx1","Tlr4","Stat1","Stat3","Tgfbr1","Nfe2l2","Sod1","Cat")

sc_file <- file.path(result_dir, "SCISSOR_results/sc_final_with_scores.rds")
sc_obj <- readRDS(sc_file)

# 检查哪些BCP基因存在
bcp_in_sc <- bcp_genes[bcp_genes %in% rownames(sc_obj)]
cat(sprintf("BCP轴基因在单细胞数据中: %d/%d\n", length(bcp_in_sc), length(bcp_genes)))
cat("存在的基因:", paste(bcp_in_sc, collapse=", "), "\n")

# 检查哪些在count矩阵中（可用于敲除）
count_matrix <- as.matrix(GetAssayData(sc_obj, layer="counts"))
bcp_in_count <- bcp_genes[bcp_genes %in% rownames(count_matrix)]
cat(sprintf("可用于敲除的BCP基因: %d\n", length(bcp_in_count)))
cat("可敲除:", paste(bcp_in_count, collapse=", "), "\n")

cat("\n"); cat(rep("=", 60), sep=""); cat("\n")
cat("建议敲除基因选择:\n")
cat("\n选项A: PC Hub基因（验证因果关系）\n")
cat("  - S100a4: 差异最大 (log2FC=1.57)\n")
cat("  - Bcl3: 转录因子，免疫调控\n")
cat("  - Aldoc: 星形胶质细胞marker\n")
cat("\n选项B: BCP轴核心基因（验证BCP机制）\n")
cat("  - Nfkb1: 核心炎症转录因子\n")
cat("  - Stat3: 显著上调 (log2FC=45.45)\n")
cat("  - Fdx1: 铁死亡相关\n")
cat("\n选项C: 组合敲除\n")
cat("  - 同时敲除Hub基因 + BCP轴基因\n")
cat(rep("=", 60), sep=""); cat("\n")
