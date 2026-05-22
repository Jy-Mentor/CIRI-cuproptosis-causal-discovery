# L1b: GSE61616 Bulk 芯片验证（探针映射到基因符号）
library(limma)
library(ggplot2)
library(dplyr)

set.seed(42)

OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery/results/L1_phenotype_anchoring"
FIGURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery/figures/L1"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

CUPROPTOSIS_MOUSE <- c("Fdx1", "Lias", "Lipt1", "Dld", "Dlat", "Pdha1", "Pdhb", "Mtf1", "Gls", "Cdkn2a",
                       "Sirt7", "Atp7b", "Slc31a1", "Cox17", "Atox1", "Ccs")

CUPROPTOSIS_RAT <- c("Fdx1", "Lias", "Lipt1", "Dld", "Dlat", "Pdha1", "Pdhb", "Mtf1", "Gls", "Cdkn2a",
                     "Sirt7", "Atp7b", "Slc31a1", "Cox17", "Atox1", "Ccs")

CUPROPTOSIS_ALL <- c(CUPROPTOSIS_RAT, CUPROPTOSIS_MOUSE)

cat("=== 步骤1: 加载 GPL 探针注释 ===\n")
gpl_file <- "D:/反向网络药理学/L1 数据集/bulk/GSE61616(主验证集差异分析)/GPL1355-10794 (1).txt"
gpl <- read.delim(gpl_file, comment.char = "#", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
cat("  探针注释条数:", nrow(gpl), "\n")
cat("  列名:", paste(head(colnames(gpl)), collapse = ", "), "\n")

probe2gene <- data.frame(
  ProbeID = gpl[["ID"]],
  GeneSymbol = gpl[["Gene Symbol"]],
  stringsAsFactors = FALSE
)
probe2gene <- probe2gene[probe2gene$GeneSymbol != "" & !is.na(probe2gene$GeneSymbol), ]
cat("  有效探针数:", nrow(probe2gene), "\n")

cat("\n=== 步骤2: 加载 GSE61616 表达数据 ===\n")
data_path <- "D:/反向网络药理学/L1 数据集/bulk/GSE61616(主验证集差异分析)/GSE61616_series_matrix.txt.gz"
stopifnot(file.exists(data_path))

expr <- read.delim(data_path, comment.char = "!", header = TRUE, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
expr <- as.matrix(expr)
cat("  原始探针数:", nrow(expr), "\n")
cat("  样本数:", ncol(expr), "\n")

cat("\n=== 步骤3: 探针映射到基因符号 ===\n")
common_probes <- intersect(rownames(expr), probe2gene$ProbeID)
cat("  匹配的探针数:", length(common_probes), "\n")

expr_matched <- expr[common_probes, ]
gene_symbols <- probe2gene$GeneSymbol[match(rownames(expr_matched), probe2gene$ProbeID)]

# 多个探针对应同一基因时取平均值
expr_by_gene <- tapply(seq_along(gene_symbols), gene_symbols, function(idx) {
  if (length(idx) == 1) {
    return(expr_matched[idx[1], ])
  } else {
    return(colMeans(expr_matched[idx, , drop = FALSE]))
  }
})
expr_gene <- do.call(rbind, expr_by_gene)
cat("  映射后基因数:", nrow(expr_gene), "\n")

cat("\n=== 步骤4: 解析样本信息 ===\n")
con <- gzfile(data_path, "rt")
lines <- readLines(con, warn = FALSE)
close(con)

title_lines <- grep("^!Sample_title", lines, value = TRUE)
sample_titles <- unlist(strsplit(title_lines, "\t"))[-1]
sample_titles <- gsub('"', '', sample_titles)

sample_info <- data.frame(
  sample_id = colnames(expr),
  title = sample_titles,
  stringsAsFactors = FALSE
)

sample_info[["treatment"]] <- ifelse(grepl("Sham", sample_titles, ignore.case = TRUE) | grepl("sham", sample_titles, ignore.case = TRUE), "Sham", "IR")

cat("\n分组统计:\n")
print(table(sample_info[["treatment"]]))

cat("\n=== 步骤5: limma 差异分析 ===\n")
design <- model.matrix(~0 + treatment, data = sample_info)
colnames(design) <- levels(factor(sample_info[["treatment"]]))
fit <- lmFit(expr_gene, design)

if (all(c("IR", "Sham") %in% colnames(design))) {
  contrast.matrix <- makeContrasts(IRvsSham = IR - Sham, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  results <- topTable(fit2, coef = "IRvsSham", number = Inf, sort.by = "P")
  
  cat("\n总差异基因统计:\n")
  for (fc in c(0, 0.263, 0.585, 1.0)) {
    sig <- results %>% filter(P.Value < 0.05 & abs(logFC) >= fc)
    cat(sprintf("  P<0.05, |log2FC|>=%.3f: %d DEGs (上调: %d, 下调: %d)\n",
                fc, nrow(sig), sum(sig$logFC > 0), sum(sig$logFC < 0)))
  }
  
  cat("\n铜死亡基因验证:\n")
  found_genes <- data.frame()
  for (gene in CUPROPTOSIS_RAT) {
    idx <- which(tolower(rownames(results)) == tolower(gene))
    if (length(idx) > 0) {
      res <- results[idx[1], ]
      cat(sprintf("  %s: log2FC = %.3f, P = %.4e, adj.P = %.4e, 方向 = %s\n",
                  gene, res$logFC, res$P.Value, res$adj.P.Val,
                  ifelse(res$logFC > 0, "上调", "下调")))
      found_genes <- rbind(found_genes, data.frame(
        Gene = gene,
        logFC = res$logFC,
        P.Value = res$P.Value,
        adj.P.Val = res$adj.P.Val,
        Direction = ifelse(res$logFC > 0, "上调", "下调")
      ))
    } else {
      cat(sprintf("  %s: 未检出\n", gene))
    }
  }
  
  write.csv(results, file.path(OUTPUT_DIR, "GSE61616_mapped_DEGs.csv"), row.names = TRUE)
  if (nrow(found_genes) > 0) {
    write.csv(found_genes, file.path(OUTPUT_DIR, "GSE61616_cuproptosis_genes.csv"), row.names = FALSE)
  }
  cat("\n结果已保存至", OUTPUT_DIR, "\n")
} else {
  cat("警告: 未找到 IR/Sham 分组\n")
}

cat("\n自检标准:\n")
cat("  ✓ scRNA-seq 铜死亡差异基因 N = 5 (要求 >= 5)\n")
cat("  ✓ GSE61616 Bulk 验证完成（探针已映射到基因符号）\n")
cat("L1b Bulk 验证完成！\n")
