# SCISSOR分析脚本
# 整合单细胞数据和Bulk数据，识别与卒中表型相关的细胞

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

cat("=== SCISSOR分析 ===\n")

# 安装Scissor包
if(!"devtools" %in% installed.packages()){install.packages('devtools')}
if(!"Scissor" %in% installed.packages()){
  cat("安装Scissor包...\n")
  devtools::install_github('sunduanchen/Scissor')
}
if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"survival" %in% installed.packages()){install.packages('survival')}

library(Scissor)
library(Seurat)
library(survival)

# 读取单细胞数据
cat("\n1. 读取单细胞数据...\n")
sc_obj <- readRDS('C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_annotated.rds')
sc_obj <- JoinLayers(sc_obj)
cat(sprintf("  单细胞: %d 细胞 x %d 基因\n", ncol(sc_obj), nrow(sc_obj)))
cat("  细胞类型分布:\n")
print(table(sc_obj$cell_type))

# 读取Bulk表达矩阵
cat("\n2. 读取Bulk表达矩阵...\n")
bulk_expr <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/expr_matrix_gene.rds')
bulk_pheno <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE58294_result/phenotype.rds')
cat(sprintf("  Bulk: %d 基因 x %d 样本\n", nrow(bulk_expr), ncol(bulk_expr)))
cat(sprintf("  表型: Control=%d, Stroke=%d\n", sum(bulk_pheno==0), sum(bulk_pheno==1)))

# 找出单细胞和Bulk数据中共有的基因
cat("\n3. 匹配单细胞和Bulk数据的基因...\n")
sc_genes <- rownames(sc_obj)
bulk_genes <- rownames(bulk_expr)
common_genes <- intersect(sc_genes, bulk_genes)
cat(sprintf("  单细胞基因数: %d\n", length(sc_genes)))
cat(sprintf("  Bulk基因数: %d\n", length(bulk_genes)))
cat(sprintf("  共有基因数: %d\n", length(common_genes)))

# 检查BCP轴基因是否在共有基因中
bcp_genes <- c("AGER", "RAGE", "NFKB1", "FDX1", "TLR4", "STAT1", "STAT3", "TGFB1", "NFE2L2")
bcp_in_common <- bcp_genes[bcp_genes %in% common_genes]
cat(sprintf("  BCP轴基因在共有基因中: %d/%d\n", length(bcp_in_common), length(bcp_genes)))

# 准备SCISSOR输入
cat("\n4. 准备SCISSOR输入...\n")
bulk_data <- as.matrix(bulk_expr[common_genes, ])
colnames(bulk_data) <- names(bulk_pheno)

# 确保样本顺序一致
pheno_numeric <- as.numeric(bulk_pheno)
names(pheno_numeric) <- names(bulk_pheno)

# SCISSOR需要特定的表型格式
phenotype <- pheno_numeric

cat(sprintf("  Bulk矩阵: %d x %d\n", nrow(bulk_data), ncol(bulk_data)))
cat(sprintf("  表型向量长度: %d\n", length(phenotype)))

# 运行SCISSOR
cat("\n5. 运行SCISSOR分析...\n")
set.seed(666)

tryCatch({
  Scissor_result <- Scissor(
    sc_data = sc_obj,
    bulk_data = bulk_data,
    phenotype = phenotype,
    alpha = 0.05,
    family = "binomial",
    save.file = "Scissor_result.RData"
  )

  # 获取结果
  sc_obj <- Scissor_result[[1]]

  # 查看Scissor标签分布
  cat("\n6. SCISSOR结果...\n")
  if("Scissor" %in% colnames(sc_obj@meta.data)) {
    cat("  Scissor标签分布:\n")
    print(table(sc_obj@meta.data$Scissor, useNA = "ifany"))

    # 按细胞类型分析Scissor分布
    cat("\n  各细胞类型的Scissor分布:\n")
    scissor_by_cell <- table(sc_obj$cell_type, sc_obj$Scissor, useNA = "ifany")
    print(scissor_by_cell)

    # 计算每个细胞类型的Scissor+比例
    cat("\n  各细胞类型的Scissor+比例:\n")
    for(ct in levels(sc_obj$cell_type)) {
      cells_ct <- sc_obj@meta.data[sc_obj@meta.data$cell_type == ct, ]
      n_total <- nrow(cells_ct)
      n_pos <- sum(cells_ct$Scissor == 2, na.rm = TRUE)  # 2 = Scissor+
      n_neg <- sum(cells_ct$Scissor == 1, na.rm = TRUE)  # 1 = Scissor-
      pct_pos <- n_pos / n_total * 100
      pct_neg <- n_neg / n_total * 100
      cat(sprintf("    %s: Scissor+ = %.1f%% (%d/%d), Scissor- = %.1f%% (%d/%d)\n",
                  ct, pct_pos, n_pos, n_total, pct_neg, n_neg, n_total))
    }

    # 保存结果
    saveRDS(sc_obj, file = "C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/sc_scissor.rds")

    # 保存Scissor信息
    scissor_info <- sc_obj@meta.data[, c("cell_type", "condition", "Scissor")]
    write.xlsx(scissor_info, "C:/Users/Jy-Mentor-7/Desktop/虚拟敲除/result/Scissor_cell_info.xlsx")

    cat("\n  结果已保存:\n")
    cat("    - sc_scissor.rds: 带有Scissor标签的Seurat对象\n")
    cat("    - Scissor_cell_info.xlsx: 每个细胞的Scissor信息\n")

  } else {
    cat("  警告: Scissor分析未成功生成标签\n")
  }

}, error = function(e) {
  cat(sprintf("  SCISSOR运行出错: %s\n", conditionMessage(e)))
  cat("\n  尝试简化分析...\n")
})

cat("\n=== SCISSOR分析完成 ===\n")
