# GSE58294芯片数据分析
# 处理Affymetrix CEL文件

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

cat("=== GSE58294 芯片数据处理 ===\n")

# 安装必要包
if(!"BiocManager" %in% installed.packages()){install.packages('BiocManager')}
if(!"affy" %in% installed.packages()){BiocManager::install('affy')}
if(!"affyio" %in% installed.packages()){BiocManager::install('affyio')}
if(!"Biobase" %in% installed.packages()){BiocManager::install('Biobase')}
if(!"openxlsx" %in% installed.packages()){install.packages('openxlsx')}

library(affy)
library(Biobase)
library(openxlsx)
library(dplyr)

# 设置数据目录
dir_data <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创"
dir_save <- paste0(dir_data, "/GSE58294_result")
if(!dir.exists(dir_save)){dir.create(dir_save, recursive = T)}

# 查找CEL文件
cel_files <- list.files(dir_data, pattern = "\\.CEL$", full.names = TRUE)
cat(sprintf("找到 %d 个CEL文件\n", length(cel_files)))

# 读取样本信息（用于分组的）
sample_info <- read.table(paste0(dir_data, "/GSE58294_series_matrix.txt"),
                          header = TRUE, sep = "\t", comment.char = "!",
                          stringsAsFactors = FALSE, nrows = 100)
cat("样本信息列名:\n")
print(colnames(sample_info)[1:10])

# 读取芯片注释
cat("\n读取GPL570注释...\n")
annot <- read.table(paste0(dir_data, "/GPL570-55999.txt"),
                    header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE,
                    quote = "", comment.char = "#")

cat(sprintf("注释探针数: %d\n", nrow(annot)))

# 提取探针→基因映射
probe2gene <- annot[, c("ID", "Gene.Symbol")]
probe2gene <- probe2gene[probe2gene$Gene.Symbol != "", ]
probe2gene <- probe2gene[!duplicated(probe2gene$ID), ]
probe2gene <- probe2gene[!is.na(probe2gene$Gene.Symbol), ]
cat(sprintf("有效探针数: %d\n", nrow(probe2gene)))

# 读取CEL文件
cat("\n读取CEL文件（这可能需要几分钟）...\n")
set.seed(666)
raw_data <- ReadAffy(filenames = cel_files)
cat(sprintf("CEL文件读取完成: %d 样本 x %d 探针\n",
            ncol(raw_data), nrow(raw_data)))

# 质量控制
cat("\n执行RMA归一化...\n")
eset <- rma(raw_data)

# 提取表达矩阵
expr_matrix <- exprs(eset)
cat(sprintf("表达矩阵: %d 探针 x %d 样本\n", nrow(expr_matrix), ncol(expr_matrix)))

# 探针ID转换为基因符号
cat("\n探针ID转换为基因符号...\n")
expr_df <- as.data.frame(expr_matrix)
expr_df$probe_id <- rownames(expr_df)
expr_df <- merge(expr_df, probe2gene, by.x = "probe_id", by.y = "ID", all.x = TRUE)

# 对于多个探针对应同一基因，取平均值
expr_df <- expr_df[!is.na(expr_df$Gene.Symbol), ]
expr_matrix_gene <- expr_df %>%
  dplyr::group_by(Gene.Symbol) %>%
  dplyr::summarise(across(dplyr::starts_with("GSM"), mean, na.rm = TRUE))
expr_matrix_gene <- as.data.frame(expr_matrix_gene)
rownames(expr_matrix_gene) <- expr_matrix_gene$Gene.Symbol
expr_matrix_gene <- expr_matrix_gene[, -1]

cat(sprintf("基因级别表达矩阵: %d 基因 x %d 样本\n",
            nrow(expr_matrix_gene), ncol(expr_matrix_gene)))

# 准备表型数据
cat("\n准备表型数据...\n")
sample_names <- colnames(expr_matrix_gene)

# 根据样本名判断条件（需要根据实际样本名调整）
# 通常Control样本vsStroke样本
phenotype <- ifelse(grepl("Control", sample_names, ignore.case = TRUE), 0, 1)
names(phenotype) <- sample_names

cat("表型分布:\n")
cat(sprintf("  Control (0): %d\n", sum(phenotype == 0)))
cat(sprintf("  Stroke (1): %d\n", sum(phenotype == 1)))

# 保存结果
saveRDS(expr_matrix_gene, file = paste0(dir_save, "/expr_matrix_gene.rds"))
saveRDS(phenotype, file = paste0(dir_save, "/phenotype.rds"))

# 保存为CSV
write.xlsx(expr_matrix_gene, paste0(dir_save, "/expr_matrix_gene.xlsx"))
write.xlsx(data.frame(Gene = names(phenotype), Phenotype = phenotype),
           paste0(dir_save, "/phenotype.xlsx"))

# BCP轴基因检查
cat("\n=== BCP轴基因表达检查 ===\n")
bcp_genes <- c("AGER", "RAGE", "NFKB1", "FDX1", "TLR4", "STAT1", "STAT3", "TGFB1", "NFE2L2")
bcp_genes_found <- rownames(expr_matrix_gene)[toupper(rownames(expr_matrix_gene)) %in% toupper(bcp_genes)]

cat("在Bulk数据中检测到的BCP轴基因:\n")
for(g in bcp_genes_found) {
  cat(sprintf("  %s\n", g))
}

# 计算BCP基因在Control vs Stroke中的差异
if(length(bcp_genes_found) > 0) {
  bcp_expr <- expr_matrix_gene[bcp_genes_found, ]
  bcp_diff <- data.frame(
    Gene = rownames(bcp_expr),
    Control_mean = rowMeans(bcp_expr[, phenotype == 0]),
    Stroke_mean = rowMeans(bcp_expr[, phenotype == 1]),
    log2FC = rowMeans(bcp_expr[, phenotype == 1]) - rowMeans(bcp_expr[, phenotype == 0])
  )
  bcp_diff <- bcp_diff[order(abs(bcp_diff$log2FC), decreasing = TRUE), ]
  write.xlsx(bcp_diff, paste0(dir_save, "/BCP_genes_expression.xlsx"))
  cat("\nBCP基因差异表达:\n")
  print(bcp_diff)
}

cat("\n=== 处理完成 ===\n")
cat(sprintf("结果保存在: %s\n", dir_save))
cat("1. expr_matrix_gene.rds - 基因表达矩阵\n")
cat("2. phenotype.rds - 表型数据\n")
cat("3. BCP_genes_expression.xlsx - BCP轴基因表达\n")
