#!/usr/bin/env Rscript

# L1 交叉验证 v3.1：修正 scRNA-seq 标准化流程
# 修正：使用 Wolf et al. 2018 标准流程 normalize_total(1e4) → log1p 计算 log2FC
# 新增：纳入 Wilcoxon p 值、表达百分比统计

suppressPackageStartupMessages({
  library(openxlsx)
})

set.seed(42)

base_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery"
results_dir <- file.path(base_dir, "results/L1_phenotype_anchoring")
l1_dir <- file.path(base_dir, "L1_phenotype_anchoring")

# 基因分类
gene_categories <- c(
  FDX1="铜死亡核心", LIAS="铜死亡核心", DLD="铜死亡核心", DLAT="铜死亡核心",
  DLST="铜死亡核心", PDHA1="铜死亡核心", PDHB="铜死亡核心", GLS="铜死亡核心",
  GCSH="铜死亡核心", LIPT1="铜死亡核心", LIPT2="铜死亡核心", CDKN2A="铜死亡核心",
  NFE2L2="铜死亡核心", NLRP3="铜死亡核心",
  SLC31A1="铜离子转运", SLC31A2="铜离子转运", SLC11A2="铜离子转运",
  STEAP3="铜离子转运", ATP7A="铜离子转运", ATP7B="铜离子转运",
  ATOX1="铜伴侣蛋白", CCS="铜伴侣蛋白", COX17="铜伴侣蛋白",
  COX11="铜伴侣蛋白", SCO1="铜伴侣蛋白", SCO2="铜伴侣蛋白",
  MT1A="铜储存缓冲", MT2A="铜储存缓冲", ALB="铜储存缓冲",
  CP="铜储存缓冲", SOD1="铜储存缓冲", SOD3="铜储存缓冲",
  COMMD1="铜代谢调控", MTF1="铜代谢调控"
)

# ==========================================
# 1. 加载数据
# ==========================================

# scRNA-seq（从 .raw 原始计数按 Wolf2018 标准流程重新计算）
scrna <- read.csv(file.path(results_dir, "scRNA_cuproptosis_all_genes.csv"),
  stringsAsFactors = FALSE)
scrna$gene <- toupper(scrna$gene)
scrna$direction <- ifelse(!is.na(scrna$log2FC), ifelse(scrna$log2FC > 0, "上调", "下调"), NA)
scrna$category <- gene_categories[scrna$gene]
scrna$category[is.na(scrna$category)] <- "其他"
scrna$found_logical <- scrna$found == "True"

# 添加 BH 校正 p 值（处理未找到基因的 NA）
scrna$p_adjust <- NA_real_
scrna$p_adjust[!is.na(scrna$p_value)] <- p.adjust(
  scrna$p_value[!is.na(scrna$p_value)], method = "BH")
scrna$significant <- ifelse(!is.na(scrna$p_adjust) & scrna$p_adjust < 0.05, TRUE, FALSE)

cat("scRNA-seq: 找到", sum(scrna$found_logical), "/", nrow(scrna), "基因\n")
cat("  - 显著差异 (BH adj.p<0.05):", sum(scrna$significant, na.rm=TRUE), "基因\n")

# GSE97537 (24h Rat)
gse97537 <- read.csv(file.path(l1_dir, "GSE97537_cuproptosis_DEGs.csv"), stringsAsFactors = FALSE)
gse97537$Human_Gene <- toupper(gse97537$Human_Gene)
gse97537$Direction <- ifelse(!is.na(gse97537$log2FC), ifelse(gse97537$log2FC > 0, "上调", "下调"), NA)

# GSE61616 (7d Mouse)
gse61616_excel <- file.path(results_dir, "L1_Bulk_GSE61616_Summary.xlsx")
gse61616_raw <- read.xlsx(gse61616_excel, sheet = "Cuproptosis_Genes", startRow = 3, colNames = TRUE)
gse61616 <- data.frame(
  gene = toupper(as.character(gse61616_raw[[2]])),
  log2FC = as.numeric(gse61616_raw[[4]]),
  direction = ifelse(!is.na(as.numeric(gse61616_raw[[4]])),
    ifelse(as.numeric(gse61616_raw[[4]]) > 0, "上调", "下调"), NA),
  category = gene_categories[toupper(as.character(gse61616_raw[[2]]))],
  stringsAsFactors = FALSE
)
gse61616 <- gse61616[!is.na(gse61616$gene) & gse61616$gene != "", ]
gse61616$category[is.na(gse61616$category)] <- "其他"

# 构建映射
scrna_map <- setNames(scrna$log2FC, scrna$gene)
scrna_dir_map <- setNames(scrna$direction, scrna$gene)
scrna_p_map <- setNames(scrna$p_adjust, scrna$gene)
scrna_pct_map <- setNames(paste0(scrna$pct_mcao, "%/", scrna$pct_sham, "%"), scrna$gene)
gse61616_map <- setNames(gse61616$log2FC, gse61616$gene)
gse61616_dir_map <- setNames(gse61616$direction, gse61616$gene)

# ==========================================
# 2. 横向对比：GSE97537 (24h) vs scRNA-seq (24h)
# ==========================================
cat("\n=== 横向对比：24h 急性期 (GSE97537 vs scRNA-seq) ===\n")

common_genes_h <- intersect(gse97537$Human_Gene[gse97537$Status == "检出"],
  scrna$gene[scrna$found_logical])
cat("  共同检出基因:", length(common_genes_h), "\n")

horizontal <- data.frame(
  Gene = common_genes_h,
  Category = unname(gene_categories[common_genes_h]),
  GSE97537_log2FC = gse97537$log2FC[match(common_genes_h, gse97537$Human_Gene)],
  GSE97537_Dir = gse97537$Direction[match(common_genes_h, gse97537$Human_Gene)],
  GSE97537_Sig = gse97537$Significant[match(common_genes_h, gse97537$Human_Gene)],
  scRNA_log2FC = unname(scrna_map[common_genes_h]),
  scRNA_Dir = unname(scrna_dir_map[common_genes_h]),
  scRNA_padj = unname(scrna_p_map[common_genes_h]),
  scRNA_pctExpr = unname(scrna_pct_map[common_genes_h]),
  stringsAsFactors = FALSE
)
horizontal$Category[is.na(horizontal$Category)] <- "其他"

# 计算一致性
horizontal$Dir_Consistent <- ifelse(
  !is.na(horizontal$GSE97537_Dir) & !is.na(horizontal$scRNA_Dir),
  ifelse(horizontal$GSE97537_Dir == horizontal$scRNA_Dir, "一致", "不一致"),
  "NA"
)

# Spearman相关（排除NA）
valid_h <- !is.na(horizontal$GSE97537_log2FC) & !is.na(horizontal$scRNA_log2FC)
if (sum(valid_h) >= 3) {
  cor_test <- cor.test(horizontal$GSE97537_log2FC[valid_h],
    horizontal$scRNA_log2FC[valid_h], method = "spearman")
  cat(sprintf("  Spearman rho = %.3f, p = %.4f\n", cor_test$estimate, cor_test$p.value))
} else {
  cat("  有效数据点不足，无法计算Spearman相关\n")
  cor_test <- list(estimate = NA, p.value = NA)
}

concordant <- sum(horizontal$Dir_Consistent == "一致", na.rm = TRUE)
total_valid <- sum(!is.na(horizontal$Dir_Consistent))
cat(sprintf("  方向一致性: %d/%d (%.1f%%)\n\n", concordant, total_valid, 100*concordant/total_valid))

for (i in seq_len(nrow(horizontal))) {
  cat(sprintf("  %-10s | 24hRat: %+7.4f (%s) | scRNA: %+7.4f (%s, p_adj=%.2e) | %s\n",
    horizontal$Gene[i],
    horizontal$GSE97537_log2FC[i],
    horizontal$GSE97537_Dir[i],
    horizontal$scRNA_log2FC[i],
    horizontal$scRNA_Dir[i],
    horizontal$scRNA_padj[i],
    horizontal$Dir_Consistent[i]))
}

# ==========================================
# 3. 纵向对比：GSE97537 (24h) vs GSE61616 (7d)
# ==========================================
cat("\n=== 纵向对比：时间动态 (24h vs 7d) ===\n")

common_genes_v <- intersect(gse97537$Human_Gene[gse97537$Status == "检出"], gse61616$gene)
cat("  双时间点检出:", length(common_genes_v), "\n")

vertical <- data.frame(
  Gene = common_genes_v,
  Category = unname(gene_categories[common_genes_v]),
  GSE97537_log2FC = gse97537$log2FC[match(common_genes_v, gse97537$Human_Gene)],
  GSE97537_Dir = gse97537$Direction[match(common_genes_v, gse97537$Human_Gene)],
  GSE97537_Sig = gse97537$Significant[match(common_genes_v, gse97537$Human_Gene)],
  GSE61616_log2FC = unname(gse61616_map[common_genes_v]),
  GSE61616_Dir = unname(gse61616_dir_map[common_genes_v]),
  stringsAsFactors = FALSE
)
vertical$Category[is.na(vertical$Category)] <- "其他"

vertical$Time_Dynamic <- ifelse(
  !is.na(vertical$GSE97537_Dir) & !is.na(vertical$GSE61616_Dir),
  ifelse(vertical$GSE97537_Dir == vertical$GSE61616_Dir, "持续响应", "方向反转"),
  ifelse(!is.na(vertical$GSE97537_Dir), "仅24h响应",
    ifelse(!is.na(vertical$GSE61616_Dir), "仅7d响应", "无响应")))

cat("  持续响应:", sum(vertical$Time_Dynamic == "持续响应"), "\n")
cat("  方向反转:", sum(vertical$Time_Dynamic == "方向反转"), "\n")
cat("  仅24h响应:", sum(vertical$Time_Dynamic == "仅24h响应"), "\n")
cat("  仅7d响应:", sum(vertical$Time_Dynamic == "仅7d响应"), "\n\n")

for (i in seq_len(nrow(vertical))) {
  cat(sprintf("  %-10s | 24h: %+7.4f (%s) | 7d: %+7.4f (%s) | %s\n",
    vertical$Gene[i],
    vertical$GSE97537_log2FC[i],
    vertical$GSE97537_Dir[i],
    vertical$GSE61616_log2FC[i],
    vertical$GSE61616_Dir[i],
    vertical$Time_Dynamic[i]))
}

# ==========================================
# 4. 输出 Excel
# ==========================================
out_file <- file.path(results_dir, "L1_Cross_Validation_v3.xlsx")
wb <- createWorkbook()

addWorksheet(wb, "Horizontal_24h")
writeData(wb, "Horizontal_24h", horizontal)

addWorksheet(wb, "Vertical_TimeCourse")
writeData(wb, "Vertical_TimeCourse", vertical)

addWorksheet(wb, "scRNA_Stats")
writeData(wb, "scRNA_Stats", scrna[, c("gene", "sham_mean_log1p", "mcao_mean_log1p",
  "log2FC", "p_value", "p_adjust", "significant", "pct_mcao", "pct_sham",
  "category")])

addWorksheet(wb, "Summary")
summary_stats <- data.frame(
  Metric = c(
    "数据版本",
    "scRNA-seq 标准化方法",
    "参考文献",
    "横向: 共同检出基因数",
    "横向: Spearman rho",
    "横向: Spearman p值",
    "横向: 方向一致性",
    "纵向: 双时间点检出基因数",
    "纵向: 持续响应",
    "纵向: 方向反转",
    "纵向: 仅24h响应",
    "纵向: 仅7d响应"
  ),
  Value = c(
    "v3.1 — 修正log2FC计算方法",
    "normalize_total(1e4) → log1p (Wolf et al. 2018)",
    "Wolf et al. 2018 Genome Biology; Luecken & Theis 2019 Mol Syst Biol",
    length(common_genes_h),
    ifelse(is.na(cor_test$estimate), "NA", sprintf("%.3f", cor_test$estimate)),
    ifelse(is.na(cor_test$p.value), "NA", sprintf("%.4f", cor_test$p.value)),
    sprintf("%d/%d (%.1f%%)", concordant, total_valid, 100*concordant/total_valid),
    length(common_genes_v),
    sum(vertical$Time_Dynamic == "持续响应"),
    sum(vertical$Time_Dynamic == "方向反转"),
    sum(vertical$Time_Dynamic == "仅24h响应"),
    sum(vertical$Time_Dynamic == "仅7d响应")
  ),
  stringsAsFactors = FALSE
)
writeData(wb, "Summary", summary_stats)

saveWorkbook(wb, out_file, overwrite = TRUE)

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("  交叉验证 v3.1 完成！\n")
cat("  输出:", out_file, "\n")
cat("  scRNA-seq 标准化: normalize_total(1e4) → log1p (Wolf et al. 2018)\n")
cat(paste(rep("=", 60), collapse=""), "\n")