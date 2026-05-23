#!/usr/bin/env Rscript

# GSE97537 探针映射 + 差异分析 + 铜死亡验证（大鼠 24h MCAO vs Sham）
# 平台：GPL1355 (Rat Genome 230 2.0 Array)

suppressPackageStartupMessages({
  library(limma)
  library(openxlsx)
  library(ggplot2)
  library(ggrepel)
  library(Cairo)
  library(dplyr)
})

set.seed(42)

# ==========================================
# 1. 加载 GPL1355 注释文件
# ==========================================
cat("加载 GPL1355 注释文件...\n")
gpl_file <- "D:/反向网络药理学/L1 数据集/bulk/GSE97537(24H)/GPL1355-10794 (1).txt"

if (!file.exists(gpl_file)) {
  stop("GPL注释文件不存在: ", gpl_file)
}

gpl <- read.table(gpl_file, header=FALSE, sep="\t", quote="", skip=17,
                  comment.char="", stringsAsFactors=FALSE, fill=TRUE)

# 手动读取第17行作为列名
header_line <- readLines(gpl_file, n=17)[17]
col_names <- strsplit(header_line, "\t")[[1]]
colnames(gpl) <- col_names

cat("  GPL 行数:", nrow(gpl), "\n")
cat("  GPL 列数:", ncol(gpl), "\n")

# 找基因名列（Gene Symbol）
gene_col_candidates <- c("Gene Symbol", "Gene.Title", "GeneTitle", "Gene_Title")
gene_col <- intersect(gene_col_candidates, colnames(gpl))
if (length(gene_col) == 0) {
  cat("  可用列:", paste(colnames(gpl), collapse=", "), "\n")
  stop("未找到基因名列")
}

cat("  使用基因列:", gene_col[1], "\n")

gpl$probe_id <- trimws(gpl[[1]])
gpl$gene_symbol <- trimws(gpl[[gene_col[1]]])

# 过滤空基因名
gpl <- gpl[gpl$gene_symbol != "" & !is.na(gpl$gene_symbol), ]
cat("  有效探针-基因映射:", nrow(gpl), "\n")

# ==========================================
# 2. 加载 GSE97537 表达数据
# ==========================================
cat("\n加载 GSE97537 表达数据...\n")
matrix_file <- "D:/反向网络药理学/L1 数据集/bulk/GSE97537(24H)/GSE97537_series_matrix.txt"
lines <- readLines(matrix_file)

# 解析 stress 特征
stress_chars <- character()
for (line in lines) {
  if (grepl("!Sample_characteristics_ch1\t\"stress:", line)) {
    stress_chars <- strsplit(gsub("!Sample_characteristics_ch1\t", "", line), "\t")[[1]]
    break
  }
}

# 解析数据块
data_start <- which(grepl("^!series_matrix_table_begin", lines))
data_end <- which(grepl("^!series_matrix_table_end", lines))
data_lines <- lines[(data_start+1):(data_end-1)]
data_text <- paste(data_lines, collapse="\n")
expr_data <- read.table(text = data_text, header=TRUE, row.names=1, sep="\t", quote="", check.names=FALSE)

groups <- ifelse(grepl("MCAO|Occlusion", stress_chars, ignore.case=TRUE), "MCAO", "Sham")
cat("  样本分组: MCAO=", sum(groups=="MCAO"), ", Sham=", sum(groups=="Sham"), "\n")
cat("  表达矩阵:", nrow(expr_data), "探针 x", ncol(expr_data), "样本\n")

# ==========================================
# 3. 探针 → 基因符号映射
# ==========================================
cat("\n探针映射到基因符号...\n")
probe_ids <- rownames(expr_data)
# 清理探针ID（去除引号）
probe_ids_clean <- gsub('"', '', probe_ids)
rownames(expr_data) <- probe_ids_clean

# 映射
mapped <- probe_ids_clean %in% gpl$probe_id
cat("  已映射探针:", sum(mapped), "/", length(probe_ids_clean), "\n")

expr_mapped <- expr_data[mapped, ]
gene_symbols <- gpl$gene_symbol[match(rownames(expr_mapped), gpl$probe_id)]

# 去重复基因（保留第一个）
dup_genes <- duplicated(gene_symbols)
if (any(dup_genes)) {
  cat("  去重基因:", sum(dup_genes), "\n")
  expr_mapped <- expr_mapped[!dup_genes, ]
  gene_symbols <- gene_symbols[!dup_genes]
}

rownames(expr_mapped) <- gene_symbols
cat("  映射后基因数:", nrow(expr_mapped), "\n")

# ==========================================
# 4. limma 差异分析
# ==========================================
cat("\nlimma 差异分析...\n")
design <- model.matrix(~ 0 + factor(groups))
colnames(design) <- c("MCAO", "Sham")

fit <- lmFit(expr_mapped, design)
contrast.matrix <- makeContrasts(MCAO - Sham, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, coef=1, number=Inf, sort.by="none")
results$Gene <- rownames(results)

cat("  总基因数:", nrow(results), "\n")
cat("  上调 (P<0.05, log2FC>=0.585):", sum(results$adj.P.Val < 0.05 & results$logFC >= 0.585, na.rm=TRUE), "\n")
cat("  下调 (P<0.05, log2FC<=-0.585):", sum(results$adj.P.Val < 0.05 & results$logFC <= -0.585, na.rm=TRUE), "\n")

# ==========================================
# 5. 铜死亡基因验证
# ==========================================
convert_to_rat <- function(human_genes) {
  sapply(human_genes, function(x) {
    paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
  })
}

cuproptosis_genes <- c(
  "ATP7A", "ATP7B", "CDKN2A", "COX17", "DBT", "DLD", "DLAT", "DLST", "FDX1",
  "GCSH", "GLS", "LIAS", "LIPT1", "LIPT2", "MTF1", "NFE2L2", "NLRP3",
  "PDHA1", "PDHB", "SLC31A1",
  "SLC31A2", "SLC11A2", "STEAP3", "ATOX1", "CCS", "COX11", "SCO1", "SCO2",
  "MT1A", "MT2A", "ALB", "CP", "SOD1", "SOD3", "COMMD1"
)

genes_to_check <- convert_to_rat(cuproptosis_genes)
matched_genes <- genes_to_check[genes_to_check %in% results$Gene]
cat("\n铜死亡基因检出:", length(matched_genes), "/35\n")

# 输出
write.csv(results, "GSE97537_GEO2R_DEGs.csv", row.names=FALSE)

cupro_results <- data.frame(
  Gene = genes_to_check,
  Human_Gene = cuproptosis_genes,
  Status = ifelse(genes_to_check %in% results$Gene, "检出", "未检出"),
  log2FC = NA, P.Value = NA, adj.P.Val = NA, Direction = NA,
  Significant = NA, AveExpr = NA, t_value = NA,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(cupro_results))) {
  gene <- cupro_results$Gene[i]
  if (gene %in% results$Gene) {
    row <- results[results$Gene == gene, ]
    cupro_results$log2FC[i] <- round(row$logFC, 4)
    cupro_results$P.Value[i] <- formatC(row$P.Value, format="e", digits=3)
    cupro_results$adj.P.Val[i] <- formatC(row$adj.P.Val, format="e", digits=3)
    cupro_results$Direction[i] <- ifelse(row$logFC > 0, "上调", "下调")
    cupro_results$Significant[i] <- ifelse(row$adj.P.Val < 0.05 & abs(row$logFC) >= 0.585, "是", "否")
    cupro_results$AveExpr[i] <- round(row$AveExpr, 3)
    cupro_results$t_value[i] <- round(row$t, 3)
  }
}

write.csv(cupro_results, "GSE97537_cuproptosis_DEGs.csv", row.names=FALSE)

# ==========================================
# 6. 火山图
# ==========================================
data <- results
data$negLog10AdjP <- -log10(data$adj.P.Val)
data$negLog10AdjP[is.infinite(data$negLog10AdjP)] <- max(data$negLog10AdjP[!is.infinite(data$negLog10AdjP)]) + 1

data$Significance <- "Not Significant"
data$Significance[data$adj.P.Val < 0.05 & data$logFC >= 0.585] <- "Up-regulated"
data$Significance[data$adj.P.Val < 0.05 & data$logFC <= -0.585] <- "Down-regulated"

data$IsTarget <- "Other"
data$IsTarget[data$Gene %in% genes_to_check] <- "Cuproptosis Gene"

top50 <- data %>%
  filter(!Gene %in% genes_to_check) %>%
  arrange(adj.P.Val) %>%
  slice_head(n=50) %>%
  mutate(LabelType = "Top 50 Significant")

cupro_data <- data %>%
  filter(Gene %in% genes_to_check) %>%
  mutate(LabelType = "Cuproptosis Gene")

label_data <- bind_rows(cupro_data, top50) %>%
  distinct(Gene, .keep_all = TRUE)

p <- ggplot(data, aes(x=logFC, y=negLog10AdjP)) +
  geom_point(aes(color=Significance), alpha=0.6, size=1.2) +
  scale_color_manual(
    values = c("Up-regulated"="#D55E00", "Down-regulated"="#1B9E77", "Not Significant"="#999999")
  ) +
  geom_point(data = filter(data, IsTarget=="Cuproptosis Gene"),
             color="#0072B2", size=4.5, alpha=0.3) +
  geom_point(data = filter(data, IsTarget=="Cuproptosis Gene"),
             color="#0072B2", size=3, alpha=0.95) +
  geom_vline(xintercept=c(-0.585, 0.585), linetype="dashed", color="#666666", linewidth=0.8) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="#666666", linewidth=0.8) +
  geom_text_repel(
    data = filter(label_data, LabelType=="Cuproptosis Gene"),
    aes(label=Gene), color="#0072B2", size=3.5, fontface="bold",
    max.overlaps=50, box.padding=0.35, point.padding=0.25,
    segment.color="#0072B2", segment.alpha=0.5
  ) +
  geom_text_repel(
    data = filter(label_data, LabelType=="Top 50 Significant"),
    aes(label=Gene), color="#333333", size=3,
    max.overlaps=50, box.padding=0.35, point.padding=0.25,
    segment.color="#999999", segment.alpha=0.3
  ) +
  annotate("text",
    x = max(data$logFC, na.rm=TRUE) * 0.95,
    y = max(data$negLog10AdjP, na.rm=TRUE) * 0.15,
    label = paste0(
      "Up: ", sum(data$Significance=="Up-regulated"), "\n",
      "Down: ", sum(data$Significance=="Down-regulated"), "\n",
      "Cuproptosis: ", length(matched_genes), "/35"
    ),
    hjust=1, vjust=0, size=4.5, fontface="bold"
  ) +
  labs(
    title = "GSE97537: 24h MCAO vs Sham (Rat)",
    x = "log2(Fold Change)", y = "-log10(Adjusted P-value)"
  ) +
  theme_bw(base_size=12) +
  theme(plot.margin=margin(0.5, 1.2, 0.5, 0.5, "cm"))

ggsave("GSE97537_volcano_plot.png", p, width=12, height=8, dpi=300)

Cairo::CairoPDF("GSE97537_volcano_plot.pdf", width=12, height=8)
print(p)
dev.off()

cat("\n✅ GSE97537 分析完成！\n")
cat("  铜死亡检出:", length(matched_genes), "/35\n")
cat("  输出: GSE97537_GEO2R_DEGs.csv, GSE97537_cuproptosis_DEGs.csv, 火山图\n")
