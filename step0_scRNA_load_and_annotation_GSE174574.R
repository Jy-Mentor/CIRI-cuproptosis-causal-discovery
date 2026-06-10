# 国内镜像设置----
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

# 安装程序包----
if(!"BiocManager" %in% installed.packages()){install.packages('BiocManager')}
if(!"openxlsx" %in% installed.packages()){install.packages('openxlsx')}
if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"SeuratObject" %in% installed.packages()){BiocManager::install('SeuratObject')}
if(!"dplyr" %in% installed.packages()){install.packages('dplyr')}

# 加载程序包----
options(warn = -1)
library(Seurat)
library(SeuratObject)
library(openxlsx)
library(dplyr)

# 读取单细胞注释数据----
sc_obj <- readRDS('./result/sc_annotated.rds')

# 查看细胞类型和条件----
cat("=== 细胞类型分布 ===\n")
print(table(sc_obj$cell_type))

cat("\n=== 样本条件分布 ===\n")
print(table(sc_obj$condition))

# 选择关注的细胞类型（可修改）----
# 根据之前的分析，Nfkb1在Microglia中显著下调，选择Microglia
target_celltype <- "Microglia"
sc_obj_sub <- subset(sc_obj, cell_type == target_celltype)

cat(sprintf("\n=== 选择 %s 细胞 ===\n", target_celltype))
print(table(sc_obj_sub$condition))

# 选择基因：高变基因 + BCP轴基因----
# 1. 选择前1000个高变基因
hvg <- Seurat::VariableFeatures(sc_obj)[1:1000]

# 2. BCP轴基因
bcp_genes <- c("Ager", "Nfkb1", "Fdx1", "Tlr4", "Stat1", "Stat3", "Tgfb1", "Nfe2l2", "Jak1", "Ccl2", "Icam1", "Hmox1")

# 3. 合并并去重
selected_genes <- unique(c(hvg, bcp_genes[bcp_genes %in% rownames(sc_obj_sub)]))
sc_obj_sub <- subset(sc_obj_sub, features = selected_genes)

cat(sprintf("选择的基因数: %d\n", length(selected_genes)))

# 随机提取部分细胞----
set.seed(666)
sc_obj_select <- subset(sc_obj_sub, downsample = 1000)
cat(sprintf("下采样后细胞数: %d\n", ncol(sc_obj_select)))
print(table(sc_obj_select$condition))

# 提取表达矩阵----
sc_Matrix <- SeuratObject::LayerData(sc_obj_select, assay="RNA", layer='counts')

# 结果保存----
dir_save <- './result/knk_input'
if(!dir.exists(dir_save)){dir.create(dir_save, recursive = T)}

saveRDS(sc_Matrix, file = paste0(dir_save,'/sc_Matrix.rds'))
saveRDS(sc_obj_select, file = paste0(dir_save,'/sc_obj_select.rds'))

cat("\n=== 数据准备完成 ===\n")
cat(sprintf("表达矩阵保存至: %s/sc_Matrix.rds\n", dir_save))
cat(sprintf("细胞数: %d, 基因数: %d\n", ncol(sc_Matrix), nrow(sc_Matrix)))
