# ============================================================
# 铜死亡∩CIRI 论文图表绘制脚本
# 生成6张发表级图表
# ============================================================

# ============================================================
# 1. 环境准备
# ============================================================
install_if_missing <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs) > 0) {
    install.packages(new_pkgs, dependencies = TRUE, repos = "https://cran.r-project.org")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

pkgs <- c("ggplot2", "ggpubr", "ggsci", "dplyr", "tidyr", "reshape2",
          "scales", "RColorBrewer", "grid", "gridExtra")
install_if_missing(pkgs)

# 工作目录
WORK_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
DATA_DIR <- file.path(WORK_DIR, "figure_data")
OUTPUT_DIR <- file.path(WORK_DIR, "figures_output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 全局主题设置
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

# ============================================================
# 图1: 铜死亡GSVA评分箱线图 + 统计检验
# Stroke vs Control (GSE61616)
# ============================================================
cat("\n[1/6] 绘制图1: GSVA评分箱线图...\n")

df1 <- read.csv(file.path(DATA_DIR, "data_fig1_gsva_score.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(df1) > 0)

df1$Group <- factor(df1$Group, levels = c("Control", "Stroke"))

# 统计检验
wilcox_result <- wilcox.test(Cuproptosis_Score ~ Group, data = df1, exact = FALSE)
p_value <- wilcox_result$p.value
p_label <- ifelse(p_value < 0.001, "***",
           ifelse(p_value < 0.01, "**",
           ifelse(p_value < 0.05, "*", "ns")))

# 计算效应量 (Cohen's d)
n1 <- sum(df1$Group == "Control")
n2 <- sum(df1$Group == "Stroke")
m1 <- mean(df1$Cuproptosis_Score[df1$Group == "Control"])
m2 <- mean(df1$Cuproptosis_Score[df1$Group == "Stroke"])
s1 <- sd(df1$Cuproptosis_Score[df1$Group == "Control"])
s2 <- sd(df1$Cuproptosis_Score[df1$Group == "Stroke"])
pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
cohens_d <- (m2 - m1) / pooled_sd

p1 <- ggplot(df1, aes(x = Group, y = Cuproptosis_Score, fill = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85,
               color = "grey30", linewidth = 0.6) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6, color = "grey20") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.x.npc = "center",
                     size = 4) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Stroke" = "#E64B35"),
                    labels = c("Control (n=24)", paste0("Stroke (n=39)"))) +
  labs(
    title = "Cuproptosis GSVA Score in GSE61616",
    x = "",
    y = "Cuproptosis GSVA Enrichment Score",
    fill = "Group",
    subtitle = paste0("Wilcoxon test, Cohen's d = ", round(cohens_d, 3))
  ) +
  theme_pub +
  theme(legend.position = "top")

ggsave(file.path(OUTPUT_DIR, "Fig1_GSVA_Boxplot.pdf"), p1, width = 6, height = 5.5, dpi = 600)
ggsave(file.path(OUTPUT_DIR, "Fig1_GSVA_Boxplot.png"), p1, width = 6, height = 5.5, dpi = 600)
cat("  图1已保存\n")

# ============================================================
# 图2: GSEA富集分析结果条形图
# 7个铜死亡基因集的NES、P值、FDR
# ============================================================
cat("\n[2/6] 绘制图2: GSEA富集条形图...\n")

df2 <- read.csv(file.path(DATA_DIR, "data_fig2_gsea.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(df2) > 0)

# 简化基因集名称
df2$Label <- gsub("Cuproptosis_", "", df2$Gene_Set)
df2$Label <- gsub("_", " ", df2$Label)
df2$Label <- factor(df2$Label, levels = df2$Label[order(df2$NES)])

# 颜色映射: 激活(红) vs 抑制(蓝)
df2$Direction <- ifelse(df2$NES > 0, "Activated", "Suppressed")

p2 <- ggplot(df2, aes(x = Label, y = NES, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.65, alpha = 0.9, color = "grey30", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "grey30") +
  geom_text(aes(
    label = paste0("NES=", round(NES, 2), "\nFDR=", format.pval(FDR, digits = 2)),
    vjust = ifelse(NES > 0, -0.3, 1.3)
  ), size = 3, color = "grey20") +
  scale_fill_manual(values = c("Activated" = "#E64B35", "Suppressed" = "#4DBBD5")) +
  labs(
    title = "GSEA: Cuproptosis Gene Sets in Stroke",
    x = "",
    y = "Normalized Enrichment Score (NES)",
    fill = "Direction"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
    legend.position = "top"
  ) +
  ylim(min(df2$NES) - 0.3, max(df2$NES) + 0.3)

ggsave(file.path(OUTPUT_DIR, "Fig2_GSEA_Barplot.pdf"), p2, width = 8, height = 5.5, dpi = 600)
ggsave(file.path(OUTPUT_DIR, "Fig2_GSEA_Barplot.png"), p2, width = 8, height = 5.5, dpi = 600)
cat("  图2已保存\n")

# ============================================================
# 图3: 单细胞铜死亡评分小提琴图
# 各细胞类型 Sham vs MCAO
# ============================================================
cat("\n[3/6] 绘制图3: 单细胞小提琴图...\n")

df3 <- read.csv(file.path(DATA_DIR, "data_fig3_single_cell.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(df3) > 0)

# 转换为长格式用于分面小提琴图
df3_long <- data.frame()
for (i in 1:nrow(df3)) {
  ct <- df3$cell_type[i]
  df3_long <- rbind(df3_long, data.frame(
    cell_type = ct,
    Group = "Sham",
    mean_score = df3$mean_Sham[i],
    sd_score = df3$std_Sham[i],
    n = df3$n_Sham[i]
  ))
  df3_long <- rbind(df3_long, data.frame(
    cell_type = ct,
    Group = "MCAO",
    mean_score = df3$mean_MCAO[i],
    sd_score = df3$std_MCAO[i],
    n = df3$n_MCAO[i]
  ))
}

df3_long$cell_type <- factor(df3_long$cell_type,
  levels = df3$cell_type[order(df3$mean_total, decreasing = TRUE)])
df3_long$Group <- factor(df3_long$Group, levels = c("Sham", "MCAO"))

# 计算每个细胞类型的fold change用于标注
df3$fc_label <- paste0("FC=", round(df3$fold_change, 2))
df3$cell_type <- factor(df3$cell_type,
  levels = df3$cell_type[order(df3$mean_total, decreasing = TRUE)])

p3 <- ggplot(df3_long, aes(x = Group, y = mean_score, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6,
           alpha = 0.85, color = "grey30", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score),
                position = position_dodge(0.6), width = 0.2, linewidth = 0.5) +
  geom_text(data = df3, aes(x = 1.5, y = max(df3_long$mean_score + df3_long$sd_score) * 0.95,
                             label = fc_label),
            inherit.aes = FALSE, size = 2.8, color = "grey30") +
  facet_wrap(~ cell_type, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c("Sham" = "#4DBBD5", "MCAO" = "#E64B35")) +
  labs(
    title = "Single-Cell Cuproptosis Score by Cell Type",
    x = "",
    y = "Mean Cuproptosis Score",
    fill = "Group"
  ) +
  theme_pub +
  theme(
    strip.text = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(size = 8)
  )

ggsave(file.path(OUTPUT_DIR, "Fig3_SingleCell_Barplot.pdf"), p3, width = 10, height = 6, dpi = 600)
ggsave(file.path(OUTPUT_DIR, "Fig3_SingleCell_Barplot.png"), p3, width = 10, height = 6, dpi = 600)
cat("  图3已保存\n")

# ============================================================
# 图4: 铜死亡基因-炎症因子相关性图 (基于模板23+29改进版)
# 参考: bioR23.corplot.R (corrplot矩阵) + bioR29.bubble.R (气泡图)
# 策略: 构建基因×因子矩阵 → corrplot圆形热图 (最清晰)
# ============================================================
cat("\n[4/6] 绘制图4: 基因-炎症因子相关性图...\n")

df4 <- read.csv(file.path(DATA_DIR, "data_fig4_cytokine_corr.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(df4) > 0)

# ---- 方案C: 矩阵热图 (参考模板23 corrplot, 最清晰专业) ----
# 构建基因×因子宽矩阵
mat4 <- df4 %>%
  select(Cuproptosis_Gene, Cytokine, Correlation) %>%
  tidyr::pivot_wider(names_from = Cytokine, values_from = Correlation) %>%
  tibble::column_to_rownames("Cuproptosis_Gene")

# 只保留有数据的基因和因子
mat4 <- mat4[rowSums(!is.na(mat4)) > 0, colSums(!is.na(mat4)) > 0]
mat4[is.na(mat4)] <- 0

# 使用ggplot2绘制矩阵气泡热图 (避免corrplot依赖)
df4_mat <- df4 %>%
  select(Cuproptosis_Gene, Cytokine, Correlation, P_value) %>%
  mutate(
    neg_log10_P = -log10(pmax(P_value, 1e-300)),
    sig_label = ifelse(P_value < 0.001, "***",
                ifelse(P_value < 0.01, "**",
                ifelse(P_value < 0.05, "*", "")))
  )

# 排序: 基因按平均相关性, 因子按出现频率
gene_order <- df4_mat %>%
  group_by(Cuproptosis_Gene) %>%
  summarise(avg_r = mean(abs(Correlation))) %>%
  arrange(desc(avg_r)) %>%
  pull(Cuproptosis_Gene)

cyto_order <- df4_mat %>%
  group_by(Cytokine) %>%
  summarise(avg_r = mean(abs(Correlation))) %>%
  arrange(desc(avg_r)) %>%
  pull(Cytokine)

df4_mat$Cuproptosis_Gene <- factor(df4_mat$Cuproptosis_Gene, levels = gene_order)
df4_mat$Cytokine <- factor(df4_mat$Cytokine, levels = cyto_order)

# 矩阵气泡热图 (参考模板29气泡图 + 模板23热图配色)
p4c <- ggplot(df4_mat, aes(x = Cytokine, y = Cuproptosis_Gene)) +
  geom_tile(fill = "white", color = "grey85", linewidth = 0.3) +
  geom_point(aes(size = abs(Correlation), color = Correlation), alpha = 0.9) +
  geom_text(aes(label = sig_label), color = "black", size = 3.5, fontface = "bold") +
  scale_size_continuous(range = c(2, 10), name = "|Spearman R|", breaks = c(0.84, 0.87, 0.90)) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                        midpoint = 0.85, limits = c(0.82, 0.92), name = "Correlation") +
  labs(
    title = "Cuproptosis Gene - Cytokine Correlation Matrix",
    subtitle = paste0("n=", nrow(df4_mat), " gene-cytokine pairs, Spearman correlation"),
    x = "Inflammatory Cytokine",
    y = "Cuproptosis Gene"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    legend.position = "right",
    legend.box = "vertical"
  )

ggsave(file.path(OUTPUT_DIR, "Fig4_Cytokine_Matrix.pdf"), p4c, width = 10, height = 7, dpi = 600)
ggsave(file.path(OUTPUT_DIR, "Fig4_Cytokine_Matrix.png"), p4c, width = 10, height = 7, dpi = 600)
cat("  图4C (矩阵气泡热图) 已保存\n")

# ---- 方案D: 棒棒糖图 (参考模板30棒棒糖图, 一维排序最清晰) ----
# 只取每个基因Top-1, 共13个点, 完全不重叠
df4_top1 <- df4 %>%
  group_by(Cuproptosis_Gene) %>%
  slice_max(order_by = abs(Correlation), n = 1) %>%
  ungroup() %>%
  arrange(desc(Correlation)) %>%
  mutate(
    Pair = paste0(Cuproptosis_Gene, " - ", Cytokine),
    Pair = factor(Pair, levels = rev(Pair)),
    neg_log10_P = -log10(pmax(P_value, 1e-300))
  )

p4d <- ggplot(df4_top1, aes(x = Correlation, y = Pair)) +
  geom_segment(aes(xend = 0.80, yend = Pair), color = "grey70", linewidth = 0.6) +
  geom_point(aes(size = neg_log10_P, color = Cuproptosis_Gene), alpha = 0.9) +
  geom_text(aes(label = round(Correlation, 3), hjust = -0.25),
            size = 3.2, color = "grey20", fontface = "bold") +
  scale_size_continuous(range = c(4, 12), name = "-log10(P)") +
  scale_color_manual(values = c(
    "DLD" = "#E64B35", "PDHA2" = "#4DBBD5", "CDKN2A" = "#00A087",
    "ATP7B" = "#3C5488", "PDHA1" = "#F39B7F", "LIPT1" = "#8491B4",
    "MTF1" = "#91D1C2", "FDX1" = "#DC0000", "GCSH" = "#7E6148",
    "GLS" = "#B09C85", "LIAS" = "#3E5F8A", "DLAT" = "#9E5F8A", "PDHB" = "#5F8A6E"
  ), name = "Gene") +
  labs(
    title = "Top Cuproptosis Gene - Cytokine Correlation",
    subtitle = "Strongest correlation per gene (n=13 pairs)",
    x = "Spearman Correlation",
    y = ""
  ) +
  theme_pub +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 10, face = "bold")
  ) +
  xlim(0.80, 0.95)

ggsave(file.path(OUTPUT_DIR, "Fig4_Cytokine_Top1_Lollipop.pdf"), p4d, width = 8, height = 6.5, dpi = 600)
ggsave(file.path(OUTPUT_DIR, "Fig4_Cytokine_Top1_Lollipop.png"), p4d, width = 8, height = 6.5, dpi = 600)
cat("  图4D (Top1棒棒糖图) 已保存\n")

# ============================================================
# 图5: WGCNA模块关联 + 免疫浸润
# 组合图: 左=模块富集, 右=免疫细胞相关性
# ============================================================
cat("\n[5/6] 绘制图5: WGCNA + 免疫浸润...\n")

df5_wgcna <- read.csv(file.path(DATA_DIR, "data_fig5_wgcna_modules.csv"), stringsAsFactors = FALSE)
df5_immune <- read.csv(file.path(DATA_DIR, "data_fig5_immune.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(df5_wgcna) > 0)
stopifnot(nrow(df5_immune) > 0)

# 图5A: WGCNA模块铜死亡基因富集
df5_wgcna$Module <- factor(df5_wgcna$Module,
  levels = df5_wgcna$Module[order(df5_wgcna$Enrichment_Ratio, decreasing = TRUE)])
df5_wgcna$neg_log10_FDR <- -log10(pmax(df5_wgcna$FDR, 1e-300))

p5a <- ggplot(df5_wgcna, aes(x = Module, y = Enrichment_Ratio, fill = neg_log10_FDR)) +
  geom_bar(stat = "identity", width = 0.6, color = "grey30", linewidth = 0.4) +
  geom_text(aes(label = paste0(N_Cuproptosis, "/", Module_Size, "\n(", 
                                round(Coverage * 100, 1), "%)"),
                vjust = -0.3), size = 3.2, color = "grey20") +
  scale_fill_gradient(low = "#FDBF6F", high = "#E64B35", name = "-log10(FDR)") +
  labs(
    title = "WGCNA Module Cuproptosis Gene Enrichment",
    x = "Module",
    y = "Enrichment Ratio"
  ) +
  theme_pub +
  theme(legend.position = "right")

# 图5B: 免疫细胞与铜死亡评分相关性
df5_immune$Cell_Type <- factor(df5_immune$Cell_Type,
  levels = df5_immune$Cell_Type[order(df5_immune$Correlation)])
df5_immune$Direction <- ifelse(df5_immune$Correlation > 0, "Positive", "Negative")
df5_immune$sig_label <- ifelse(df5_immune$Significant,
  ifelse(df5_immune$P_value < 0.001, "***",
  ifelse(df5_immune$P_value < 0.01, "**", "*")), "")

p5b <- ggplot(df5_immune, aes(x = Cell_Type, y = Correlation, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.6, color = "grey30", linewidth = 0.4) +
  geom_text(aes(label = sig_label, vjust = ifelse(Correlation > 0, -0.5, 1.5)),
            size = 4, color = "grey20") +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "grey30") +
  scale_fill_manual(values = c("Positive" = "#E64B35", "Negative" = "#4DBBD5")) +
  labs(
    title = "Immune Cell Correlation with Cuproptosis Score",
    x = "",
    y = "Spearman Correlation",
    fill = "Direction"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
    legend.position = "top"
  )

# 组合图5
p5_combined <- ggarrange(p5a, p5b, ncol = 2, labels = c("A", "B"),
                         widths = c(1, 1.2),
                         font.label = list(size = 14, face = "bold"))

ggsave(file.path(OUTPUT_DIR, "Fig5_WGCNA_Immune.pdf"), p5_combined, width = 14, height = 5.5, dpi = 600)
ggsave(file.path(OUTPUT_DIR, "Fig5_WGCNA_Immune.png"), p5_combined, width = 14, height = 5.5, dpi = 600)
cat("  图5已保存\n")

# ============================================================
# 图6: PPI网络拓扑 + 多维度融合评分
# 铜死亡基因的PPI参数 + 融合评分
# ============================================================
cat("\n[6/6] 绘制图6: PPI拓扑 + 融合评分...\n")

df6 <- read.csv(file.path(DATA_DIR, "data_fig6_ppi_fusion.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(df6) > 0)

df6$Gene <- factor(df6$Gene, levels = rev(df6$Gene))

# 图6A: 棒棒糖图 - 融合评分
p6a <- ggplot(df6, aes(x = Gene, y = Fusion_Score)) +
  geom_segment(aes(xend = Gene, yend = 0), color = "grey60", linewidth = 0.8) +
  geom_point(aes(size = Degree, color = Sig_Ratio), alpha = 0.9) +
  scale_size_continuous(range = c(3, 10), name = "PPI Degree") +
  scale_color_gradient(low = "#4DBBD5", high = "#E64B35", name = "Sig. Ratio") +
  coord_flip() +
  labs(
    title = "Multi-Dimensional Fusion Score of Cuproptosis Genes",
    x = "",
    y = "Fusion Score"
  ) +
  theme_pub

# 图6B: 堆叠条形图 - 邻居基因方向
df6_stack <- df6[, c("Gene", "N_Up", "N_Down")]
df6_stack_long <- melt(df6_stack, id.vars = "Gene", 
                       variable.name = "Direction", value.name = "Count")
df6_stack_long$Direction <- factor(df6_stack_long$Direction,
  levels = c("N_Up", "N_Down"),
  labels = c("Up-regulated", "Down-regulated"))
df6_stack_long$Gene <- factor(df6_stack_long$Gene, levels = rev(levels(df6$Gene)))

p6b <- ggplot(df6_stack_long, aes(x = Gene, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7,
           alpha = 0.85, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = c("Up-regulated" = "#E64B35", "Down-regulated" = "#4DBBD5")) +
  coord_flip() +
  labs(
    title = "Significant Neighbor DEGs per Cuproptosis Gene",
    x = "",
    y = "Number of Significant Neighbors",
    fill = "Direction"
  ) +
  theme_pub +
  theme(legend.position = "top")

# 组合图6
p6_combined <- ggarrange(p6a, p6b, ncol = 2, labels = c("A", "B"),
                         font.label = list(size = 14, face = "bold"))

ggsave(file.path(OUTPUT_DIR, "Fig6_PPI_Fusion.pdf"), p6_combined, width = 14, height = 6, dpi = 600)
ggsave(file.path(OUTPUT_DIR, "Fig6_PPI_Fusion.png"), p6_combined, width = 14, height = 6, dpi = 600)
cat("  图6已保存\n")

# ============================================================
# 完成
# ============================================================
cat("\n========================================\n")
cat("所有图表已生成完毕!\n")
cat("输出目录:", OUTPUT_DIR, "\n")
cat("========================================\n")

# 列出所有输出文件
output_files <- list.files(OUTPUT_DIR, pattern = "\\.pdf$|\\.png$")
for (f in output_files) {
  cat("  ", f, "\n")
}