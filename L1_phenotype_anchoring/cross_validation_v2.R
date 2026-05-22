#!/usr/bin/env Rscript

# L1 交叉验证 v2.0：优化版
# 参考：Nature Methods 2023, Cell 2022 跨平台验证最佳实践
# 新增：四象限图、一致性热图、Spearman相关、统计注释

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(Cairo)
  library(gridExtra)
})

set.seed(42)

# ==========================================
# 0. 配置
# ==========================================
base_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery"
results_dir <- file.path(base_dir, "results/L1_phenotype_anchoring")
fig_dir <- file.path(base_dir, "figures/L1")
l1_dir <- file.path(base_dir, "L1_phenotype_anchoring")

# 铜死亡基因颜色映射（按功能模块）
gene_colors <- list(
  "铜死亡核心" = "#E74C3C",
  "铜离子转运" = "#3498DB", 
  "铜伴侣蛋白" = "#2ECC71",
  "铜储存缓冲" = "#F39C12",
  "铜代谢调控" = "#9B59B6"
)

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
# 1. 加载数据（使用log2FC推导方向，避免编码问题）
# ==========================================
logfc_to_dir <- function(logfc) ifelse(!is.na(logfc), ifelse(logfc > 0, "上调", "下调"), NA)

# scRNA-seq - 从All_DEGs中提取铜死亡基因（更完整）
scrna_excel <- file.path(results_dir, "L1_scRNA_GSE174574_Summary.xlsx")
scrna_all <- read.xlsx(scrna_excel, sheet = "All_DEGs", startRow = 3, colNames = TRUE)

# 铜死亡基因列表
cuproptosis_genes <- c("FDX1", "LIAS", "DLD", "DLAT", "DLST", "PDHA1", "PDHB", "GLS", 
  "GCSH", "LIPT1", "LIPT2", "CDKN2A", "NFE2L2", "NLRP3",
  "SLC31A1", "SLC31A2", "SLC11A2", "STEAP3", "ATP7A", "ATP7B",
  "ATOX1", "CCS", "COX17", "COX11", "SCO1", "SCO2",
  "MT1A", "MT2A", "ALB", "CP", "SOD1", "SOD3",
  "COMMD1", "MTF1")

# 从All_DEGs中筛选铜死亡基因
scrna_matched <- scrna_all[toupper(scrna_all[[1]]) %in% cuproptosis_genes, ]
scrna <- data.frame(
  gene = toupper(as.character(scrna_matched[[1]])),
  log2FC = as.numeric(scrna_matched[["log2FC"]]),
  direction = logfc_to_dir(as.numeric(scrna_matched[["log2FC"]])),
  category = gene_categories[toupper(as.character(scrna_matched[[1]]))],
  stringsAsFactors = FALSE
)
scrna <- scrna[!is.na(scrna$gene) & scrna$gene != "", ]
scrna$category[is.na(scrna$category)] <- "其他"

cat("  从All_DEGs中提取铜死亡基因:", nrow(scrna), "/", length(cuproptosis_genes), "\n")

# GSE97537
gse97537 <- read.csv(file.path(l1_dir, "GSE97537_cuproptosis_DEGs.csv"), stringsAsFactors = FALSE)
gse97537$Human_Gene <- toupper(gse97537$Human_Gene)
gse97537$Direction <- logfc_to_dir(gse97537$log2FC)

# GSE61616
gse61616_excel <- file.path(results_dir, "L1_Bulk_GSE61616_Summary.xlsx")
gse61616_raw <- read.xlsx(gse61616_excel, sheet = "Cuproptosis_Genes", startRow = 3, colNames = TRUE)
gse61616 <- data.frame(
  gene = toupper(as.character(gse61616_raw[[2]])),
  log2FC = as.numeric(gse61616_raw[[4]]),
  direction = logfc_to_dir(as.numeric(gse61616_raw[[4]])),
  category = gene_categories[toupper(as.character(gse61616_raw[[2]]))],
  stringsAsFactors = FALSE
)
gse61616 <- gse61616[!is.na(gse61616$gene) & gse61616$gene != "", ]
gse61616$category[is.na(gse61616$category)] <- "其他"

cat("=== 数据加载 ===\n")
cat("  scRNA-seq:", nrow(scrna), "基因\n")
cat("  GSE97537 (24h Rat):", nrow(gse97537), "基因\n")
cat("  GSE61616 (7d Mouse):", nrow(gse61616), "基因\n")

# 构建映射
scrna_map <- setNames(scrna$log2FC, scrna$gene)
gse61616_map <- setNames(gse61616$log2FC, gse61616$gene)

# ==========================================
# 2. 横向对比：四象限图 + Spearman相关
# ==========================================
cat("\n=== 横向对比：24h 急性期 ===\n")

# 合并数据
common_genes_h <- intersect(gse97537$Human_Gene[gse97537$Status == "检出"], scrna$gene)
horizontal_df <- data.frame(
  gene = common_genes_h,
  gse97537_log2FC = gse97537$log2FC[match(common_genes_h, gse97537$Human_Gene)],
  scrna_log2FC = scrna_map[common_genes_h],
  category = gene_categories[common_genes_h],
  stringsAsFactors = FALSE
)
horizontal_df$category[is.na(horizontal_df$category)] <- "其他"

# 计算象限
horizontal_df$quadrant <- with(horizontal_df, ifelse(
  gse97537_log2FC > 0 & scrna_log2FC > 0, "Q1: 双上调",
  ifelse(gse97537_log2FC < 0 & scrna_log2FC < 0, "Q3: 双下调",
  ifelse(gse97537_log2FC > 0 & scrna_log2FC < 0, "Q2: Bulk上/scRNA下",
  "Q4: Bulk下/scRNA上"))))

# Spearman相关
cor_test <- cor.test(horizontal_df$gse97537_log2FC, horizontal_df$scrna_log2FC, method = "spearman")
cat(sprintf("  Spearman rho = %.3f, p = %.4f\n", cor_test$estimate, cor_test$p.value))

# 一致性率（排除NA）
valid_pairs <- !is.na(horizontal_df$gse97537_log2FC) & !is.na(horizontal_df$scrna_log2FC)
concordant <- sum((horizontal_df$gse97537_log2FC[valid_pairs] * horizontal_df$scrna_log2FC[valid_pairs]) > 0)
cat(sprintf("  方向一致性: %d/%d (%.1f%%)\n", concordant, sum(valid_pairs), 100*concordant/sum(valid_pairs)))

# 四象限图
p_quad <- ggplot(horizontal_df, aes(x = gse97537_log2FC, y = scrna_log2FC, color = category)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = gene), size = 2.5, vjust = -1, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dotted", size = 0.5) +
  scale_color_manual(values = gene_colors) +
  labs(
    title = "横向验证：GSE97537 (24h Rat) vs scRNA-seq (24h Mouse)",
    subtitle = sprintf("Spearman rho = %.3f, p = %.4f | 一致性 = %d/%d (%.1f%%)",
      cor_test$estimate, cor_test$p.value, concordant, nrow(horizontal_df), 100*concordant/nrow(horizontal_df)),
    x = "GSE97537 log2FC",
    y = "scRNA-seq log2FC",
    color = "功能模块"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

# 保存
CairoPDF(file.path(fig_dir, "L1_horizontal_quadrant.pdf"), width = 10, height = 8)
print(p_quad)
dev.off()

CairoPNG(file.path(fig_dir, "L1_horizontal_quadrant.png"), width = 1000, height = 800, dpi = 150)
print(p_quad)
dev.off()

cat("  四象限图已保存\n")

# ==========================================
# 3. 纵向对比：时间动态热图
# ==========================================
cat("\n=== 纵向对比：时间动态 ===\n")

# 合并三个数据集
common_genes_v <- Reduce(intersect, list(
  gse97537$Human_Gene[gse97537$Status == "检出"],
  gse61616$gene,
  scrna$gene
))

cat(sprintf("  三平台共同检出: %d 基因\n", length(common_genes_v)))

# 构建热图矩阵
heatmap_mat <- data.frame(
  gene = common_genes_v,
  GSE97537_24h = gse97537$log2FC[match(common_genes_v, gse97537$Human_Gene)],
  GSE61616_7d = gse61616$log2FC[match(common_genes_v, gse61616$gene)],
  scRNA_24h = scrna_map[common_genes_v],
  category = gene_categories[common_genes_v],
  stringsAsFactors = FALSE
)
heatmap_mat$category[is.na(heatmap_mat$category)] <- "其他"

# 按类别排序
heatmap_mat <- heatmap_mat[order(heatmap_mat$category, heatmap_mat$gene), ]

# 绘制热图
mat_for_heatmap <- as.matrix(heatmap_mat[, c("GSE97537_24h", "GSE61616_7d", "scRNA_24h")])
rownames(mat_for_heatmap) <- heatmap_mat$gene

# 注释
annotation_row <- data.frame(
  Category = heatmap_mat$category,
  row.names = heatmap_mat$gene
)

# 确保Category是因子，且水平与颜色匹配
heatmap_mat$category <- factor(heatmap_mat$category, levels = names(gene_colors))
annotation_row$Category <- factor(annotation_row$Category, levels = names(gene_colors))
annotation_colors <- list(Category = unlist(gene_colors))

# 保存热图
CairoPDF(file.path(fig_dir, "L1_timecourse_heatmap.pdf"), width = 8, height = 12)
pheatmap(mat_for_heatmap,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-2, 2, length.out = 101),
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  main = "铜死亡基因时间动态热图\n(24h Rat vs 7d Mouse vs 24h scRNA-seq)",
  fontsize_row = 8,
  fontsize_col = 10,
  cellwidth = 60,
  cellheight = 12,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  fontsize_number = 6
)
dev.off()

CairoPNG(file.path(fig_dir, "L1_timecourse_heatmap.png"), width = 800, height = 1200, dpi = 150)
pheatmap(mat_for_heatmap,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-2, 2, length.out = 101),
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  main = "铜死亡基因时间动态热图\n(24h Rat vs 7d Mouse vs 24h scRNA-seq)",
  fontsize_row = 8,
  fontsize_col = 10,
  cellwidth = 60,
  cellheight = 12,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  fontsize_number = 6
)
dev.off()

cat("  热图已保存\n")

# ==========================================
# 4. 纵向：时间动态分类条形图
# ==========================================
# 计算每个基因的时间动态
time_dynamic <- data.frame(
  gene = common_genes_v,
  dir_24h = sign(heatmap_mat$GSE97537_24h),
  dir_7d = sign(heatmap_mat$GSE61616_7d),
  stringsAsFactors = FALSE
)

time_dynamic$dynamic <- with(time_dynamic, ifelse(
  dir_24h == dir_7d & dir_24h != 0, "持续响应",
  ifelse(dir_24h != dir_7d & dir_24h != 0 & dir_7d != 0, "方向反转",
  ifelse(dir_24h != 0 & dir_7d == 0, "仅24h响应",
  ifelse(dir_24h == 0 & dir_7d != 0, "仅7d响应", "无响应")))))

# 统计
dynamic_summary <- as.data.frame(table(time_dynamic$dynamic))
colnames(dynamic_summary) <- c("Dynamic", "Count")

p_bar <- ggplot(dynamic_summary, aes(x = Dynamic, y = Count, fill = Dynamic)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c(
    "持续响应" = "#2ECC71",
    "方向反转" = "#E74C3C",
    "仅24h响应" = "#3498DB",
    "仅7d响应" = "#F39C12",
    "无响应" = "#95A5A6"
  )) +
  labs(
    title = "铜死亡基因时间动态分类",
    subtitle = sprintf("GSE97537 (24h) vs GSE61616 (7d) | n = %d", nrow(time_dynamic)),
    x = "",
    y = "基因数"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

CairoPDF(file.path(fig_dir, "L1_timecourse_dynamic_bar.pdf"), width = 8, height = 6)
print(p_bar)
dev.off()

CairoPNG(file.path(fig_dir, "L1_timecourse_dynamic_bar.png"), width = 800, height = 600, dpi = 150)
print(p_bar)
dev.off()

cat("  动态分类图已保存\n")

# ==========================================
# 5. 输出 Excel（优化版）
# ==========================================
out_file <- file.path(results_dir, "L1_Cross_Validation_v2.xlsx")
wb <- createWorkbook()

# Sheet 1: 横向对比详细
addWorksheet(wb, "Horizontal_Detail")
writeData(wb, "Horizontal_Detail", horizontal_df)

# Sheet 2: 纵向对比详细
addWorksheet(wb, "Vertical_Detail")
writeData(wb, "Vertical_Detail", heatmap_mat)

# Sheet 3: 时间动态分类
addWorksheet(wb, "Time_Dynamic")
writeData(wb, "Time_Dynamic", time_dynamic)

# Sheet 4: 统计摘要
addWorksheet(wb, "Summary")
summary_stats <- data.frame(
  Metric = c(
    "横向: 共同检出基因数",
    "横向: Spearman rho",
    "横向: Spearman p值",
    "横向: 方向一致性",
    "纵向: 三平台共同检出",
    "纵向: 持续响应",
    "纵向: 方向反转",
    "纵向: 仅24h响应",
    "纵向: 仅7d响应"
  ),
  Value = c(
    nrow(horizontal_df),
    sprintf("%.3f", cor_test$estimate),
    sprintf("%.4f", cor_test$p.value),
    sprintf("%d/%d (%.1f%%)", concordant, nrow(horizontal_df), 100*concordant/nrow(horizontal_df)),
    length(common_genes_v),
    sum(time_dynamic$dynamic == "持续响应"),
    sum(time_dynamic$dynamic == "方向反转"),
    sum(time_dynamic$dynamic == "仅24h响应"),
    sum(time_dynamic$dynamic == "仅7d响应")
  ),
  stringsAsFactors = FALSE
)
writeData(wb, "Summary", summary_stats)

saveWorkbook(wb, out_file, overwrite = TRUE)

cat("\n✅ 交叉验证 v2.0 完成！\n")
cat("  Excel:", out_file, "\n")
cat("  图表:\n")
cat("    - L1_horizontal_quadrant.png (四象限图)\n")
cat("    - L1_timecourse_heatmap.png (时间动态热图)\n")
cat("    - L1_timecourse_dynamic_bar.png (动态分类条形图)\n")
