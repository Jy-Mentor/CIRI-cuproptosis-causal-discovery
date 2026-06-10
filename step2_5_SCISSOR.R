# SCISSOR分析脚本
# 目的：识别与卒中表型相关的细胞亚群
# 安装: devtools::install_github('sunduanchen/Scissor')

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

if(!"devtools" %in% installed.packages()){install.packages('devtools')}
if(!"Scissor" %in% installed.packages()){devtools::install_github('sunduanchen/Scissor')}
if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"survival" %in% installed.packages()){install.packages('survival')}

options(Seurat.object.assay.version = "v3")
library(Scissor)
library(Seurat)
library(survival)

# 读取已注释的单细胞数据
sc_obj <- readRDS('./result/sc_annotated.rds')

# 检查数据
cat("=== 单细胞数据信息 ===\n")
cat(sprintf("细胞数: %d, 基因数: %d\n", ncol(sc_obj), nrow(sc_obj)))
cat("\n细胞类型分布:\n")
print(table(sc_obj$cell_type))

# 查看样本来源
cat("\n样本条件分布:\n")
print(table(sc_obj$condition))

# ============================================================
# SCISSOR需要三种数据:
# 1. 单细胞数据 (已准备好: sc_obj)
# 2. Bulk表达矩阵 (需要准备)
# 3. 表型数据 (需要准备)
# ============================================================

# 对于卒中研究，你需要准备:
# - Bulk表达矩阵: 卒中患者 vs 对照的基因表达
# - 表型数据: 疾病状态 (0=对照, 1=卒中) 或生存数据

# 示例: 创建模拟的表型数据用于演示
# 实际使用时替换为真实的bulk数据

cat("\n=== SCISSOR分析准备 ===\n")
cat("注意: SCISSOR需要匹配的bulk表达数据和表型数据\n")
cat("对于卒中研究，建议使用:\n")
cat("1. Blood transcriptome data (如from GEO: GSE58294, GSE22255)\n")
cat("2. 卒中GWAS summary statistics\n")
cat("3. 或使用PSEUCO等方法推断bulk profile\n")

# 如果你有bulk数据，按以下格式准备:
# bulk_data: genes x samples 矩阵
# phenotype: named vector (sample_id = phenotype_value)

# ============================================================
# SCISSOR主函数使用示例 (注释掉，需要真实数据)
# ============================================================

# # 加载bulk表达矩阵 (需要自行准备)
# load("bulk_expression.RData")  # bulk_dataset: genes x samples

# # 准备表型数据 (二元: 0=对照, 1=卒中)
# phenotype <- c(rep(0, 100), rep(1, 100))  # 示例
# names(phenotype) <- colnames(bulk_dataset)

# # 或使用生存数据 (需要time和status列)
# survival_data <- data.frame(
#   time = c(...),   # 生存时间
#   status = c(...)  # 0=存活, 1=死亡
# )

# # 运行SCISSOR
# Scissor_result <- Scissor(
#   sc_data = sc_obj,
#   bulk_data = bulk_dataset,
#   phenotype = phenotype,
#   alpha = 0.05,
#   family = "binomial"  # 或 "cox" for survival
# )

# # 获取结果
# sc_obj <- Scissor_result[[1]]
# head(metadata(sc_obj)$Scissor_result)

# # 可视化
# p1 <- DimPlot(sc_obj, reduction = "umap",
#               group.by = "Scissor",
#               cols = c("grey", "red", "blue"))
# p1

# ============================================================
# 替代方案: 使用差异基因进行表型关联
# ============================================================

cat("\n=== 替代方案: 基于DEG的表型关联 ===\n")
cat("如果你没有匹配的bulk数据，可以使用已完成的DEG结果进行表型关联分析\n")

# 读取DEG结果
deg_file <- './result/DEG_significant.xlsx'
if(file.exists(deg_file)) {
  library(openxlsx)
  deg_sig <- read.xlsx(deg_file)
  cat(sprintf("显著差异基因数: %d\n", nrow(deg_sig)))

  # BCP轴基因在DEG中的情况
  bcp_genes <- c("Ager", "Nfkb1", "Fdx1")
  cat("\nBCP轴基因差异状态:\n")
  for(g in bcp_genes) {
    g_deg <- deg_sig[deg_sig$gene == g, ]
    if(nrow(g_deg) > 0) {
      cat(sprintf("  %s: log2FC=%.3f, p=%.2e (%s)\n",
                  g, g_deg$avg_log2FC[1], g_deg$p_val_adj[1],
                  ifelse(g_deg$avg_log2FC[1] > 0, "上调", "下调")))
    } else {
      cat(sprintf("  %s: 无显著差异\n", g))
    }
  }
}

cat("\n=== SCISSOR分析说明 ===\n")
cat("运行SCISSOR需要:\n")
cat("1. 单细胞数据 (已完成注释)\n")
cat("2. 匹配的bulk表达矩阵 (卒中患者 vs 对照)\n")
cat("3. 表型数据 (二元/连续/生存)\n\n")
cat("准备好数据后，取消上面SCISSOR代码的注释并运行\n")
