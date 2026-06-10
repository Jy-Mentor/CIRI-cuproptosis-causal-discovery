# 检查neurons Seurat对象的基因
library(Seurat)

cache_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/analysis_cache.rds"
cache_data <- readRDS(cache_file)

neurons <- cache_data$neurons

cat("=== Seurat Object Info ===\n")
cat("Number of cells:", ncol(neurons), "\n")
cat("Number of genes:", nrow(neurons), "\n")

# 获取基因名
genes <- rownames(neurons)
cat("\nFirst 30 gene names:\n")
print(head(genes, 30))

# 检查是否有目标基因
target_genes <- c("IL6", "STAT3", "NFKB1", "PPARG", "TGFB1", "CCL2", "TLR4", "PTGS2",
                  "CCND1", "STAT1", "ICAM1", "PTPRC", "RELA", "CASP8", "CXCR4",
                  "NOTCH1", "MAPK1", "MDM2", "HSPA5", "PARP1", "JAK1", "CREBBP",
                  "MMP2", "SREBF1", "CDC42", "STAT5A", "NFE2L2", "IRF1", "HMOX1")

target_upper <- toupper(target_genes)
genes_upper <- toupper(genes)
matched <- target_genes[target_upper %in% genes_upper]
cat("\nMatched target genes:", length(matched), "\n")
print(matched)

# 显示未匹配的基因
unmatched <- target_genes[!target_upper %in% genes_upper]
cat("\nUnmatched genes:\n")
print(unmatched)
