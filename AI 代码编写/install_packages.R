# 安装必要的R包

# 安装BiocManager（如果尚未安装）
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# 安装CRAN包
cran_packages <- c(
    "limma",
    "ggplot2",
    "VennDiagram",
    "ComplexHeatmap",
    "reshape2",
    "ggpubr",
    "pheatmap",
    "stringr",
    "caret",
    "pROC",
    "tidyverse"
)

for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
    }
}

# 安装Bioconductor包
bioc_packages <- c(
    "org.Hs.eg.db",
    "enrichplot",
    "GSEABase",
    "GSVA",
    "clusterProfiler"
)

for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg)
    }
}

# 验证安装
print("安装完成，验证包是否安装成功:")
for (pkg in c(cran_packages, bioc_packages)) {
    status <- ifelse(requireNamespace(pkg, quietly = TRUE), "✓", "✗")
    print(paste(status, pkg))
}
