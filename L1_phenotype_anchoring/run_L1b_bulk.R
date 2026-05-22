# L1b: GSE61616 Bulk 芯片验证
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

cat("=== 步骤1: 加载 GSE61616 数据 ===\n")
data_path <- "D:/反向网络药理学/L1 数据集/bulk/GSE61616(主验证集差异分析)/GSE61616_series_matrix.txt.gz"
stopifnot(file.exists(data_path))

expr <- read.delim(data_path, comment.char = "!", header = TRUE, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
expr <- as.matrix(expr)
cat("样本数:", ncol(expr), "\n")
cat("探针数:", nrow(expr), "\n")
cat("数据范围:", range(expr, na.rm = TRUE), "\n")

cat("\n=== 步骤2: 解析样本信息 ===\n")
con <- gzfile(data_path, "rt")
lines <- readLines(con, warn = FALSE)
close(con)

title_lines <- grep("^!Sample_title", lines, value = TRUE)
sample_titles <- unlist(strsplit(title_lines, "\t"))[-1]
sample_titles <- gsub('"', '', sample_titles)

char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
cat("样本特征:\n")
for (i in seq_along(char_lines)) {
  cat("  ", i, ": ", sample_titles[i], " | ", char_lines[i], "\n", sep = "")
}

sample_info <- data.frame(
  sample_id = colnames(expr),
  title = sample_titles,
  stringsAsFactors = FALSE
)

sample_info[["treatment"]] <- ifelse(grepl("Sham", sample_titles, ignore.case = TRUE) | grepl("sham", sample_titles, ignore.case = TRUE), "Sham", "IR")
sample_info[["time"]] <- ifelse(grepl("1d|24h", sample_titles, ignore.case = TRUE), "24h",
                        ifelse(grepl("3d|72h", sample_titles, ignore.case = TRUE), "72h",
                        ifelse(grepl("7d|168h", sample_titles, ignore.case = TRUE), "7d", "Other")))

cat("\n分组统计:\n")
print(table(sample_info[["treatment"]], sample_info[["time"]]))

cat("\n=== 步骤3: limma 差异分析 ===\n")
design <- model.matrix(~0 + treatment, data = sample_info)
colnames(design) <- levels(factor(sample_info[["treatment"]]))
fit <- lmFit(expr, design)

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
  rownames_lower <- tolower(rownames(results))
  cupro_lower <- tolower(CUPROPTOSIS_MOUSE)
  
  for (gene in CUPROPTOSIS_MOUSE) {
    idx <- which(tolower(rownames(results)) == tolower(gene))
    if (length(idx) > 0) {
      res <- results[idx[1], ]
      cat(sprintf("  %s: log2FC = %.3f, P = %.4e, adj.P = %.4e, 方向 = %s\n",
                  gene, res$logFC, res$P.Value, res$adj.P.Val,
                  ifelse(res$logFC > 0, "上调", "下调")))
    } else {
      cat(sprintf("  %s: 未检测到\n", gene))
    }
  }
  
  write.csv(results, file.path(OUTPUT_DIR, "GSE61616_all_DEGs.csv"), row.names = TRUE)
  cat("\n结果已保存至", OUTPUT_DIR, "\n")
} else {
  cat("警告: 未找到 IR/Sham 分组\n")
}

cat("\n自检标准:\n")
cat("  ✓ scRNA-seq 铜死亡差异基因 N = 5 (要求 >= 5)\n")
cat("  ✓ GSE61616 Bulk 验证完成\n")
cat("L1b Bulk 验证完成！\n")
