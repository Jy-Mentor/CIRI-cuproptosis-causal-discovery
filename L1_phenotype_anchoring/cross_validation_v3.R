#!/usr/bin/env Rscript

# L1 交叉验证 v3.0：SCISSOR-inspired 优化版
# 参考：
#   1. SCISSOR (Wang et al., 2021, Cell) - 跨物种/跨平台表型评分
#   2. Nature Methods 2023 - 单细胞与Bulk一致性评估
#   3. Bland-Altman 图 - 临床测量一致性金标准

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
gene_colors <- c(
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

# 铜死亡基因列表
cuproptosis_genes <- names(gene_categories)

# ==========================================
# 1. 加载数据
# ==========================================
logfc_to_dir <- function(logfc) ifelse(!is.na(logfc), ifelse(logfc > 0, "上调", "下调"), NA)

# scRNA-seq - 从All_DEGs中提取
scrna_excel <- file.path(results_dir, "L1_scRNA_GSE174574_Summary.xlsx")
scrna_all <- read.xlsx(scrna_excel, sheet = "All_DEGs", startRow = 3, colNames = TRUE)
scrna_matched <- scrna_all[toupper(scrna_all[[1]]) %in% cuproptosis_genes, ]
scrna <- data.frame(
  gene = toupper(as.character(scrna_matched[[1]])),
  log2FC = as.numeric(scrna_matched[["log2FC"]]),
  p_value = as.numeric(scrna_matched[["P.Value"]]),
  adj_p = as.numeric(scrna_matched[["adj.P.Val"]]),
  direction = logfc_to_dir(as.numeric(scrna_matched[["log2FC"]])),
  category = gene_categories[toupper(as.character(scrna_matched[[1]]))],
  stringsAsFactors = FALSE
)
scrna <- scrna[!is.na(scrna$gene) & scrna$gene != "", ]
scrna$category[is.na(scrna$category)] <- "其他"

# GSE97537 (24h Rat)
gse97537 <- read.csv(file.path(l1_dir, "GSE97537_cuproptosis_DEGs.csv"), stringsAsFactors = FALSE)
gse97537$Human_Gene <- toupper(gse97537$Human_Gene)
gse97537$Direction <- logfc_to_dir(gse97537$log2FC)

# GSE61616 (7d Mouse)
gse61616_excel <- file.path(results_dir, "L1_Bulk_GSE61616_Summary.xlsx")
gse61616_raw <- read.xlsx(gse61616_excel, sheet = "Cuproptosis_Genes", startRow = 3, colNames = TRUE)
gse61616 <- data.frame(
  gene = toupper(as.character(gse61616_raw[[2]])),
  log2FC = as.numeric(gse61616_raw[[4]]),
  adj_p = as.numeric(gse61616_raw[[5]]),
  direction = logfc_to_dir(as.numeric(gse61616_raw[[4]])),
  category = gene_categories[toupper(as.character(gse61616_raw[[2]]))],
  stringsAsFactors = FALSE
)
gse61616 <- gse61616[!is.na(gse61616$gene) & gse61616$gene != "", ]
gse61616$category[is.na(gse61616$category)] <- "其他"

cat("=== 数据加载 ===\n")
cat("  scRNA-seq (All_DEGs筛选):", nrow(scrna), "/", length(cuproptosis_genes), "基因\n")
cat("  GSE97537 (24h Rat):", sum(gse97537$Status == "检出"), "基因\n")
cat("  GSE61616 (7d Mouse):", nrow(gse61616), "基因\n")

# 构建映射
scrna_map <- setNames(scrna$log2FC, scrna$gene)
scrna_p_map <- setNames(scrna$adj_p, scrna$gene)
gse61616_map <- setNames(gse61616$log2FC, gse61616$gene)

# ==========================================
# 2. 横向对比：四象限图 + 回归 + Bland-Altman
# ==========================================
cat("\n=== 横向对比：24h 急性期 ===\n")

# 合并数据
common_genes_h <- intersect(gse97537$Human_Gene[gse97537$Status == "检出"], scrna$gene)
horizontal_df <- data.frame(
  gene = common_genes_h,
  gse97537_log2FC = gse97537$log2FC[match(common_genes_h, gse97537$Human_Gene)],
  gse97537_p = gse97537$adj.P.Val[match(common_genes_h, gse97537$Human_Gene)],
  gse97537_sig = gse97537$Significant[match(common_genes_h, gse97537$Human_Gene)],
  scrna_log2FC = scrna_map[common_genes_h],
  scrna_p = scrna_p_map[common_genes_h],
  category = gene_categories[common_genes_h],
  stringsAsFactors = FALSE
)
horizontal_df$category[is.na(horizontal_df$category)] <- "其他"

# 计算显著性状态
horizontal_df$gse97537_status <- ifelse(horizontal_df$gse97537_p < 0.05, "显著", "不显著")
horizontal_df$scrna_status <- ifelse(horizontal_df$scrna_p < 0.05, "显著", "不显著")

# 双显著基因
both_sig <- horizontal_df$gse97537_status == "显著" & horizontal_df$scrna_status == "显著"
cat(sprintf("  双平台显著: %d/%d\n", sum(both_sig), nrow(horizontal_df)))

# Spearman相关（全部基因）
cor_test_all <- cor.test(horizontal_df$gse97537_log2FC, horizontal_df$scrna_log2FC, method = "spearman")
cat(sprintf("  Spearman rho (全部) = %.3f, p = %.4f\n", cor_test_all$estimate, cor_test_all$p.value))

# Spearman相关（双显著基因）
if(sum(both_sig) >= 3) {
  cor_test_sig <- cor.test(horizontal_df$gse97537_log2FC[both_sig], 
                           horizontal_df$scrna_log2FC[both_sig], method = "spearman")
  cat(sprintf("  Spearman rho (双显著) = %.3f, p = %.4f\n", 
    cor_test_sig$estimate, cor_test_sig$p.value))
}

# 线性回归
lm_fit <- lm(scrna_log2FC ~ gse97537_log2FC, data = horizontal_df)
lm_summary <- summary(lm_fit)
cat(sprintf("  线性回归: R² = %.3f, p = %.4f\n", lm_summary$r.squared, 
  coef(lm_summary)[2, 4]))

# 方向一致性
valid_pairs <- !is.na(horizontal_df$gse97537_log2FC) & !is.na(horizontal_df$scrna_log2FC)
concordant <- sum((horizontal_df$gse97537_log2FC[valid_pairs] * horizontal_df$scrna_log2FC[valid_pairs]) > 0)
concordance_rate <- 100 * concordant / sum(valid_pairs)
cat(sprintf("  方向一致性: %d/%d (%.1f%%)\n", concordant, sum(valid_pairs), concordance_rate))

# 四象限图（优化版：双显著基因高亮）
horizontal_df$significance <- with(horizontal_df, ifelse(
  gse97537_status == "显著" & scrna_status == "显著", "双显著",
  ifelse(gse97537_status == "显著" | scrna_status == "显著", "单显著", "均不显著")
))

p_quad <- ggplot(horizontal_df, aes(x = gse97537_log2FC, y = scrna_log2FC)) +
  geom_point(aes(color = category, size = significance, alpha = significance)) +
  geom_text(data = subset(horizontal_df, significance == "双显著"),
            aes(label = gene), size = 3, vjust = -1.2, fontface = "bold") +
  geom_text(data = subset(horizontal_df, significance != "双显著"),
            aes(label = gene), size = 2, vjust = -1, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dotted", linewidth = 0.5) +
  scale_color_manual(values = gene_colors) +
  scale_size_manual(values = c("双显著" = 4, "单显著" = 2.5, "均不显著" = 2)) +
  scale_alpha_manual(values = c("双显著" = 1, "单显著" = 0.7, "均不显著" = 0.4)) +
  annotate("text", x = min(horizontal_df$gse97537_log2FC) * 0.8, 
           y = max(horizontal_df$scrna_log2FC) * 0.9,
           label = sprintf("Spearman rho = %.3f\np = %.3f\nR² = %.3f\n一致性 = %.1f%%",
             cor_test_all$estimate, cor_test_all$p.value, lm_summary$r.squared, concordance_rate),
           hjust = 0, size = 3.5, fontface = "bold",
           color = "darkblue", alpha = 0.8) +
  labs(
    title = "横向验证：GSE97537 (24h Rat) vs scRNA-seq (24h Mouse)",
    subtitle = sprintf("n = %d 铜死亡基因 | 双显著 = %d", nrow(horizontal_df), sum(both_sig)),
    x = expression("GSE97537 log"[2]*"FC (Rat Bulk)"),
    y = expression("scRNA-seq log"[2]*"FC (Mouse)"),
    color = "功能模块",
    size = "显著性",
    alpha = "显著性"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    legend.position = "right"
  )

CairoPDF(file.path(fig_dir, "L1_horizontal_quadrant_v3.pdf"), width = 11, height = 8)
print(p_quad)
dev.off()

CairoPNG(file.path(fig_dir, "L1_horizontal_quadrant_v3.png"), width = 1100, height = 800, dpi = 150)
print(p_quad)
dev.off()

cat("  ✓ 四象限图已保存\n")

# Bland-Altman 图（评估系统偏差）
horizontal_df$mean_log2FC <- (horizontal_df$gse97537_log2FC + horizontal_df$scrna_log2FC) / 2
horizontal_df$diff_log2FC <- horizontal_df$gse97537_log2FC - horizontal_df$scrna_log2FC

mean_diff <- mean(horizontal_df$diff_log2FC, na.rm = TRUE)
sd_diff <- sd(horizontal_df$diff_log2FC, na.rm = TRUE)
loa_lower <- mean_diff - 1.96 * sd_diff
loa_upper <- mean_diff + 1.96 * sd_diff

p_bland <- ggplot(horizontal_df, aes(x = mean_log2FC, y = diff_log2FC)) +
  geom_point(aes(color = category), size = 3, alpha = 0.7) +
  geom_text(aes(label = gene), size = 2.5, vjust = -1, show.legend = FALSE) +
  geom_hline(yintercept = mean_diff, linetype = "solid", color = "blue", linewidth = 1) +
  geom_hline(yintercept = loa_lower, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = loa_upper, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = max(horizontal_df$mean_log2FC) * 0.7, y = mean_diff + 0.05,
           label = sprintf("Mean diff = %.3f", mean_diff), color = "blue", size = 3.5) +
  annotate("text", x = max(horizontal_df$mean_log2FC) * 0.7, y = loa_upper + 0.05,
           label = sprintf("+1.96 SD = %.3f", loa_upper), color = "red", size = 3) +
  annotate("text", x = max(horizontal_df$mean_log2FC) * 0.7, y = loa_lower - 0.08,
           label = sprintf("-1.96 SD = %.3f", loa_lower), color = "red", size = 3) +
  scale_color_manual(values = gene_colors) +
  labs(
    title = "Bland-Altman 图：平台间一致性评估",
    subtitle = sprintf("GSE97537 vs scRNA-seq | 95%% LoA: [%.3f, %.3f]", loa_lower, loa_upper),
    x = expression("Mean log"[2]*"FC (两个平台平均)"),
    y = expression("Difference log"[2]*"FC (GSE97537 - scRNA)"),
    color = "功能模块"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

CairoPDF(file.path(fig_dir, "L1_bland_altman.pdf"), width = 10, height = 7)
print(p_bland)
dev.off()

CairoPNG(file.path(fig_dir, "L1_bland_altman.png"), width = 1000, height = 700, dpi = 150)
print(p_bland)
dev.off()

cat("  ✓ Bland-Altman图已保存\n")

# ==========================================
# 3. 纵向对比：时间动态 + SCISSOR风格评分
# ==========================================
cat("\n=== 纵向对比：时间动态 ===\n")

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
heatmap_mat <- heatmap_mat[order(heatmap_mat$category, heatmap_mat$gene), ]

# 热图
mat_for_heatmap <- as.matrix(heatmap_mat[, c("GSE97537_24h", "GSE61616_7d", "scRNA_24h")])
rownames(mat_for_heatmap) <- heatmap_mat$gene

annotation_row <- data.frame(
  Category = factor(heatmap_mat$category, levels = names(gene_colors)),
  row.names = heatmap_mat$gene
)
annotation_colors <- list(Category = gene_colors)

CairoPDF(file.path(fig_dir, "L1_timecourse_heatmap_v3.pdf"), width = 8, height = 12)
pheatmap(mat_for_heatmap,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-2, 2, length.out = 101),
  annotation_row = annotation_row, annotation_colors = annotation_colors,
  main = "铜死亡基因时间动态热图\n(24h Rat vs 7d Mouse vs 24h scRNA-seq)",
  fontsize_row = 8, fontsize_col = 10,
  cellwidth = 60, cellheight = 12,
  display_numbers = TRUE, number_format = "%.2f",
  number_color = "black", fontsize_number = 6
)
dev.off()

CairoPNG(file.path(fig_dir, "L1_timecourse_heatmap_v3.png"), width = 800, height = 1200, dpi = 150)
pheatmap(mat_for_heatmap,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-2, 2, length.out = 101),
  annotation_row = annotation_row, annotation_colors = annotation_colors,
  main = "铜死亡基因时间动态热图\n(24h Rat vs 7d Mouse vs 24h scRNA-seq)",
  fontsize_row = 8, fontsize_col = 10,
  cellwidth = 60, cellheight = 12,
  display_numbers = TRUE, number_format = "%.2f",
  number_color = "black", fontsize_number = 6
)
dev.off()

cat("  ✓ 热图已保存\n")

# SCISSOR风格：时间动态评分
time_dynamic <- data.frame(
  gene = common_genes_v,
  log2FC_24h = heatmap_mat$GSE97537_24h,
  log2FC_7d = heatmap_mat$GSE61616_7d,
  log2FC_scRNA = heatmap_mat$scRNA_24h,
  stringsAsFactors = FALSE
)

# 计算动态评分（类似SCISSOR的表型评分）
time_dynamic$acute_score <- abs(time_dynamic$log2FC_24h)  # 急性期响应幅度
time_dynamic$recovery_score <- abs(time_dynamic$log2FC_7d)  # 恢复期响应幅度
time_dynamic$dynamic_index <- time_dynamic$acute_score - time_dynamic$recovery_score  # 动态指数
time_dynamic$direction_24h <- sign(time_dynamic$log2FC_24h)
time_dynamic$direction_7d <- sign(time_dynamic$log2FC_7d)

time_dynamic$pattern <- with(time_dynamic, ifelse(
  direction_24h == direction_7d & direction_24h != 0, "持续响应",
  ifelse(direction_24h != direction_7d & direction_24h != 0 & direction_7d != 0, "方向反转",
  ifelse(direction_24h != 0 & direction_7d == 0, "仅24h响应",
  ifelse(direction_24h == 0 & direction_7d != 0, "仅7d响应", "无响应")))))

# 动态分类条形图
dynamic_summary <- as.data.frame(table(time_dynamic$pattern))
colnames(dynamic_summary) <- c("Pattern", "Count")

p_bar <- ggplot(dynamic_summary, aes(x = Pattern, y = Count, fill = Pattern)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c(
    "持续响应" = "#2ECC71",
    "方向反转" = "#E74C3C",
    "仅24h响应" = "#3498DB",
    "仅7d响应" = "#F39C12",
    "无响应" = "#95A5A6"
  )) +
  labs(
    title = "铜死亡基因时间动态分类 (SCISSOR-inspired)",
    subtitle = sprintf("GSE97537 (24h) vs GSE61616 (7d) | n = %d", nrow(time_dynamic)),
    x = "",
    y = "基因数"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11)
  )

CairoPDF(file.path(fig_dir, "L1_timecourse_dynamic_bar_v3.pdf"), width = 9, height = 6)
print(p_bar)
dev.off()

CairoPNG(file.path(fig_dir, "L1_timecourse_dynamic_bar_v3.png"), width = 900, height = 600, dpi = 150)
print(p_bar)
dev.off()

cat("  ✓ 动态分类图已保存\n")

# ==========================================
# 4. 输出 Excel（优化版）
# ==========================================
out_file <- file.path(results_dir, "L1_Cross_Validation_v3.xlsx")
wb <- createWorkbook()

# Sheet 1: 横向对比详细
addWorksheet(wb, "Horizontal_Detail")
writeData(wb, "Horizontal_Detail", horizontal_df)

# Sheet 2: 纵向对比详细
addWorksheet(wb, "Vertical_Detail")
writeData(wb, "Vertical_Detail", heatmap_mat)

# Sheet 3: SCISSOR风格时间动态
addWorksheet(wb, "SCISSOR_Dynamic")
writeData(wb, "SCISSOR_Dynamic", time_dynamic)

# Sheet 4: 统计摘要
addWorksheet(wb, "Summary")
summary_stats <- data.frame(
  Metric = c(
    "横向: 共同检出基因数",
    "横向: 双平台显著基因数",
    "横向: Spearman rho (全部)",
    "横向: Spearman p值",
    "横向: 线性回归 R²",
    "横向: 线性回归 p值",
    "横向: 方向一致性",
    "横向: Bland-Altman Mean Diff",
    "横向: Bland-Altman 95% LoA Lower",
    "横向: Bland-Altman 95% LoA Upper",
    "纵向: 三平台共同检出",
    "纵向: 持续响应",
    "纵向: 方向反转",
    "纵向: 仅24h响应",
    "纵向: 仅7d响应"
  ),
  Value = c(
    nrow(horizontal_df),
    sum(both_sig),
    sprintf("%.3f", cor_test_all$estimate),
    sprintf("%.4f", cor_test_all$p.value),
    sprintf("%.3f", lm_summary$r.squared),
    sprintf("%.4f", coef(lm_summary)[2, 4]),
    sprintf("%d/%d (%.1f%%)", concordant, sum(valid_pairs), concordance_rate),
    sprintf("%.4f", mean_diff),
    sprintf("%.4f", loa_lower),
    sprintf("%.4f", loa_upper),
    length(common_genes_v),
    sum(time_dynamic$pattern == "持续响应"),
    sum(time_dynamic$pattern == "方向反转"),
    sum(time_dynamic$pattern == "仅24h响应"),
    sum(time_dynamic$pattern == "仅7d响应")
  ),
  stringsAsFactors = FALSE
)
writeData(wb, "Summary", summary_stats)

saveWorkbook(wb, out_file, overwrite = TRUE)

cat("\n✅ 交叉验证 v3.0 (SCISSOR-inspired) 完成！\n")
cat("  Excel:", out_file, "\n")
cat("  图表:\n")
cat("    - L1_horizontal_quadrant_v3.png (四象限图，双显著高亮)\n")
cat("    - L1_bland_altman.png (Bland-Altman一致性图)\n")
cat("    - L1_timecourse_heatmap_v3.png (时间动态热图)\n")
cat("    - L1_timecourse_dynamic_bar_v3.png (SCISSOR动态分类)\n")
