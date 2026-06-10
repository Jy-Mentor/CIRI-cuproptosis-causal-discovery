# 安装GEO研究必备包

# 安装Bioconductor核心包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 安装GEO研究必备包
BiocManager::install(c(
    "GEOquery",        # 从GEO数据库下载和处理数据
    "limma",           # 差异表达分析
    "Biobase",         # 基础生物信息学数据结构
    "edgeR",           # RNA-seq差异表达分析
    "DESeq2",          # RNA-seq差异表达分析
    "clusterProfiler", # 基因富集分析
    "org.Hs.eg.db",    # 人类基因注释
    "org.Mm.eg.db"     # 小鼠基因注释
), dependencies = TRUE)

# 安装CRAN包
install.packages(c(
    "ggplot2",         # 数据可视化
    "dplyr",           # 数据处理
    "tidyr",           # 数据整理
    "readxl",          # 读取Excel文件
    "stringr"          # 字符串处理
), dependencies = TRUE)

# 验证安装
cat("\n=== 安装验证 ===\n")
packages <- c(
    "GEOquery", "limma", "Biobase", "edgeR", "DESeq2", 
    "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db",
    "ggplot2", "dplyr", "tidyr", "readxl", "stringr"
)

for (pkg in packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("✓ %s 已安装\n", pkg))
    } else {
        cat(sprintf("✗ %s 安装失败\n", pkg))
    }
}

cat("\n安装完成！\n")
