# 检查已安装的包

# 检查CRAN包
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

# 检查Bioconductor包
bioc_packages <- c(
    "org.Hs.eg.db",
    "enrichplot",
    "GSEABase",
    "GSVA",
    "clusterProfiler"
)

# 验证安装
print("已安装的包:")
for (pkg in c(cran_packages, bioc_packages)) {
    status <- ifelse(requireNamespace(pkg, quietly = TRUE), "✓", "✗")
    print(paste(status, pkg))
}
