# ============================================================
# 跨物种铜死亡基因验证 - 4张出版级图表
# 数据集:
#   - scRNA-seq: GSE174574 (小鼠, 24h MCAO vs Sham)
#   - Bulk RNA-seq: GSE97537 (大鼠, 24h)
#   - Bulk RNA-seq: GSE61616 (小鼠, 7d)
# ============================================================

# ============================================================
# 1. 环境准备
# ============================================================
install_if_missing <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs) > 0) {
    install.packages(new_pkgs, dependencies = TRUE,
      repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

pkgs <- c("ggplot2", "ggrepel", "dplyr", "tidyr", "openxlsx",
          "scales", "RColorBrewer", "pheatmap", "grid", "gridExtra")
install_if_missing(pkgs)

# ============================================================
# 2. 定义路径
# ============================================================
base_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery"
scrna_file <- file.path(base_dir, "results/L1_phenotype_anchoring/scRNA_cuproptosis_all_genes.csv")
bulk_rat_file <- file.path(base_dir, "L1_phenotype_anchoring/GSE97537_cuproptosis_DEGs.csv")
bulk_mouse_file <- file.path(base_dir, "results/L1_phenotype_anchoring/L1_Bulk_GSE61616_Summary.xlsx")
output_dir <- file.path(base_dir, "results/L1_phenotype_anchoring/figures")

# 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 检查所有输入文件存在
stopifnot(file.exists(scrna_file))
stopifnot(file.exists(bulk_rat_file))
stopifnot(file.exists(bulk_mouse_file))

# ============================================================
# 3. 基因分类定义（小鼠基因名格式）
# ============================================================
gene_categories <- list(
  "铜死亡核心" = c("Fdx1", "Lias", "Dld", "Dlat", "Dlst", "Pdha1", "Pdhb",
                   "Gls", "Gcsh", "Lipt1", "Lipt2", "Cdkn2a", "Nfe2l2", "Nlrp3"),
  "铜离子转运" = c("Slc31a1", "Slc31a2", "Slc11a2", "Steap3", "Atp7a", "Atp7b"),
  "铜伴侣蛋白" = c("Atox1", "Ccs", "Cox17", "Cox11", "Sco1", "Sco2"),
  "铜储存缓冲" = c("Mt1a", "Mt2a", "Alb", "Cp", "Sod1", "Sod3"),
  "铜代谢调控" = c("Commd1", "Mtf1")
)

# 构建基因到分类的映射表（小鼠基因名 -> 分类）
gene_to_category <- data.frame(
  Gene = unlist(gene_categories),
  Category = rep(names(gene_categories), times = sapply(gene_categories, length)),
  stringsAsFactors = FALSE
)

# 同时构建人类基因名到分类的映射
mouse_to_human <- c(
  "Fdx1"="FDX1", "Lias"="LIAS", "Dld"="DLD", "Dlat"="DLAT", "Dlst"="DLST",
  "Pdha1"="PDHA1", "Pdhb"="PDHB", "Gls"="GLS", "Gcsh"="GCSH",
  "Lipt1"="LIPT1", "Lipt2"="LIPT2", "Cdkn2a"="CDKN2A", "Nfe2l2"="NFE2L2", "Nlrp3"="NLRP3",
  "Slc31a1"="SLC31A1", "Slc31a2"="SLC31A2", "Slc11a2"="SLC11A2", "Steap3"="STEAP3",
  "Atp7a"="ATP7A", "Atp7b"="ATP7B",
  "Atox1"="ATOX1", "Ccs"="CCS", "Cox17"="COX17", "Cox11"="COX11", "Sco1"="SCO1", "Sco2"="SCO2",
  "Mt1a"="MT1A", "Mt2a"="MT2A", "Alb"="ALB", "Cp"="CP", "Sod1"="SOD1", "Sod3"="SOD3",
  "Commd1"="COMMD1", "Mtf1"="MTF1"
)

# ============================================================
# 4. 读取数据
# ============================================================
cat("\n>>> 读取 scRNA-seq 数据...\n")
scrna <- read.csv(scrna_file, stringsAsFactors = FALSE, row.names = NULL)
stopifnot(nrow(scrna) > 0)
# 处理found列（Python的True/False可能是字符串）
scrna[["found"]] <- as.character(scrna[["found"]])
scrna <- scrna[scrna[["found"]] == "True", ]
# 只保留有表达数据的基因
scrna <- scrna[!is.na(scrna[["log2FC"]]), ]
stopifnot(nrow(scrna) > 0)
cat("  scRNA-seq 检出基因数:", nrow(scrna), "\n")

cat("\n>>> 读取 GSE97537 (大鼠 24h) 数据...\n")
rat <- read.csv(bulk_rat_file, stringsAsFactors = FALSE)
stopifnot(nrow(rat) > 0)
cat("  GSE97537 总基因数:", nrow(rat), "\n")

cat("\n>>> 读取 GSE61616 (小鼠 7d) 数据...\n")
mouse7d <- read.xlsx(bulk_mouse_file, sheet = "Cuproptosis_Genes", startRow = 3)
stopifnot(nrow(mouse7d) > 0)
cat("  GSE61616 总行数:", nrow(mouse7d), "\n")

# ============================================================
# 5. 数据预处理与标准化
# ============================================================

# 5.1 scRNA-seq: 筛选铜死亡相关基因 + 添加分类
cat("\n>>> 预处理 scRNA-seq 数据...\n")
scrna_cu <- scrna[scrna[["gene"]] %in% gene_to_category[["Gene"]], ]
colnames(scrna_cu)[colnames(scrna_cu) == "gene"] <- "Gene"
scrna_cu <- merge(scrna_cu, gene_to_category, by = "Gene", all.x = TRUE)
# 处理数值列（Python的True/False在R中是字符串）
scrna_cu[["p_value"]] <- as.numeric(scrna_cu[["p_value"]])
scrna_cu[["log2FC"]] <- as.numeric(scrna_cu[["log2FC"]])
# scRNA数据无p_adjust列，使用p_value代替（后续可计算BH校正）
scrna_cu[["p_adjust"]] <- p.adjust(scrna_cu[["p_value"]], method = "BH")
cat("  scRNA-seq 铜死亡基因数:", nrow(scrna_cu), "\n")

# 5.2 GSE97537: 添加小鼠基因名 + 分类
cat("\n>>> 预处理 GSE97537 数据...\n")
# 从Human_Gene映射回小鼠基因名
rat[["Mouse_Gene"]] <- names(mouse_to_human)[match(rat[["Human_Gene"]], mouse_to_human)]
rat[["log2FC_num"]] <- as.numeric(rat[["log2FC"]])
rat[["Significant_bool"]] <- rat[["Significant"]] == "是"
rat <- merge(rat, gene_to_category, by.x = "Mouse_Gene", by.y = "Gene", all.x = TRUE)
cat("  GSE97537 映射后行数:", nrow(rat), "\n")

# 5.3 GSE61616: 标准化列名 + 添加分类
cat("\n>>> 预处理 GSE61616 数据...\n")
colnames(mouse7d) <- c("Group", "Gene", "log2FC", "P.Value", "adj.P.Val",
                       "Direction", "Significant", "Status", 
                       "scRNA_cell_specific", "Direction_consistency", "Note")
mouse7d <- mouse7d[!is.na(mouse7d[["Gene"]]), ]
mouse7d[["log2FC"]] <- as.numeric(mouse7d[["log2FC"]])
# 从GSE61616自己的分组中提取更细的分类
mouse7d <- merge(mouse7d, gene_to_category, by = "Gene", all.x = TRUE)
cat("  GSE61616 处理后行数:", nrow(mouse7d), "\n")

# ============================================================
# 6. 合并数据用于交叉验证
# ============================================================
cat("\n>>> 合并三个数据集...\n")

# 6.1 共同基因（三个数据集都检出的铜死亡基因）
common_genes <- Reduce(intersect, list(
  scrna_cu[["Gene"]],
  rat[["Mouse_Gene"]][!is.na(rat[["Mouse_Gene"]])],
  mouse7d[["Gene"]][!is.na(mouse7d[["log2FC"]])]
))
cat("  共同基因数:", length(common_genes), "\n")

# 6.2 构建合并数据框（用于图1：scRNA-seq vs GSE97537）
merged_scrna_rat <- merge(
  scrna_cu[, c("Gene", "log2FC", "p_adjust", "Category")],
  rat[!is.na(rat[["Mouse_Gene"]]), c("Mouse_Gene", "log2FC_num", "Significant_bool", "Direction")],
  by.x = "Gene", by.y = "Mouse_Gene",
  all = FALSE
)
colnames(merged_scrna_rat)[colnames(merged_scrna_rat) == "log2FC"] <- "log2FC_scrna"
colnames(merged_scrna_rat)[colnames(merged_scrna_rat) == "log2FC_num"] <- "log2FC_rat"
colnames(merged_scrna_rat)[colnames(merged_scrna_rat) == "p_adjust"] <- "p_adj_scrna"
cat("  scRNA-seq vs GSE97537 配对基因数:", nrow(merged_scrna_rat), "\n")

# 6.3 构建三数据集log2FC矩阵（用于图2：热图）
heatmap_df <- data.frame(
  Gene = common_genes,
  stringsAsFactors = FALSE
)
# scRNA-seq log2FC
scrna_lfc <- scrna_cu[match(common_genes, scrna_cu[["Gene"]]), "log2FC"]
heatmap_df[["scRNA_24h"]] <- as.numeric(scrna_lfc)

# GSE97537 log2FC
rat_lfc <- rat[match(common_genes, rat[["Mouse_Gene"]]), "log2FC_num"]
heatmap_df[["GSE97537_24h"]] <- as.numeric(rat_lfc)

# GSE61616 log2FC
mouse_lfc <- mouse7d[match(common_genes, mouse7d[["Gene"]]), "log2FC"]
heatmap_df[["GSE61616_7d"]] <- as.numeric(mouse_lfc)

# 添加分类信息
heatmap_df <- merge(heatmap_df, gene_to_category, by = "Gene", all.x = TRUE)
heatmap_df <- heatmap_df[order(heatmap_df[["Category"]], heatmap_df[["Gene"]]), ]
rownames(heatmap_df) <- heatmap_df[["Gene"]]

cat("  热图基因数:", nrow(heatmap_df), "\n")

# 6.4 方向一致性分析（用于图3）
# 对每个基因，比较三个数据集的上调/下调方向
direction_df <- gene_to_category
direction_df[["scrna_dir"]] <- NA
direction_df[["rat_dir"]] <- NA
direction_df[["mouse7d_dir"]] <- NA

for (i in seq_len(nrow(direction_df))) {
  g <- direction_df[i, "Gene"]
  
  # scRNA方向
  idx_scrna <- which(scrna_cu[["Gene"]] == g)
  if (length(idx_scrna) > 0) {
    direction_df[i, "scrna_dir"] <- ifelse(scrna_cu[idx_scrna, "log2FC"] > 0, "上调", "下调")
  }
  
  # GSE97537方向
  idx_rat <- which(rat[["Mouse_Gene"]] == g & !is.na(rat[["Mouse_Gene"]]))
  if (length(idx_rat) > 0) {
    direction_df[i, "rat_dir"] <- rat[idx_rat[1], "Direction"]
  }
  
  # GSE61616方向
  idx_mouse <- which(mouse7d[["Gene"]] == g)
  if (length(idx_mouse) > 0) {
    direction_df[i, "mouse7d_dir"] <- mouse7d[idx_mouse[1], "Direction"]
  }
}

# 判断一致性：三个数据集方向一致（如果有两个或三个方向相同）
direction_df[["consistent"]] <- sapply(seq_len(nrow(direction_df)), function(i) {
  dirs <- na.omit(c(direction_df[i, "scrna_dir"], 
                    direction_df[i, "rat_dir"], 
                    direction_df[i, "mouse7d_dir"]))
  if (length(dirs) < 2) return(NA)
  up_count <- sum(dirs == "上调")
  down_count <- sum(dirs == "下调")
  if (up_count >= 2 || down_count >= 2) return("一致") else return("不一致")
})

# 按分类汇总
direction_summary <- direction_df %>%
  filter(!is.na(consistent)) %>%
  group_by(Category, consistent) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Category) %>%
  mutate(Total = sum(Count),
         Proportion = paste0(round(Count / Total * 100), "%"))

cat("  方向一致性汇总:\n")
print(direction_summary)

# ============================================================
# 7. 颜色方案定义
# ============================================================
category_colors <- c(
  "铜死亡核心" = "#E41A1C",
  "铜离子转运" = "#377EB8",
  "铜伴侣蛋白" = "#4DAF4A",
  "铜储存缓冲" = "#984EA3",
  "铜代谢调控" = "#FF7F00"
)

# ============================================================
# 图1: 四象限散点图 - GSE97537 vs scRNA-seq
# ============================================================
cat("\n>>> 绘制图1: 四象限散点图...\n")

# 计算Spearman相关性
cor_test <- cor.test(merged_scrna_rat[["log2FC_scrna"]], 
                      merged_scrna_rat[["log2FC_rat"]],
                      method = "spearman", use = "complete.obs",
                      exact = FALSE)
rho_val <- round(cor_test[["estimate"]], 3)
p_val_cor <- format(cor_test[["p.value"]], scientific = TRUE, digits = 3)

# 判断方向一致性
merged_scrna_rat[["direction_agree"]] <- sapply(seq_len(nrow(merged_scrna_rat)), function(i) {
  s <- merged_scrna_rat[i, "log2FC_scrna"]
  r <- merged_scrna_rat[i, "log2FC_rat"]
  if (is.na(s) || is.na(r)) return(NA)
  if ((s > 0 & r > 0) | (s < 0 & r < 0)) return("一致") else return("不一致")
})

# 剔除NA
merged_scrna_rat <- merged_scrna_rat[!is.na(merged_scrna_rat[["direction_agree"]]), ]

p1 <- ggplot(merged_scrna_rat, aes(x = log2FC_scrna, y = log2FC_rat)) +
  # 象限背景
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#FFDEDE", alpha = 0.3) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#DEE8FF", alpha = 0.3) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#E8E8E8", alpha = 0.2) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#E8E8E8", alpha = 0.2) +
  # 参考线
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  # 散点：显著实心，不显著空心
  geom_point(aes(fill = Category, shape = p_adj_scrna < 0.05),
             size = 4, color = "black", stroke = 0.8) +
  scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 21),
                     labels = c("TRUE" = "显著 (adj.p<0.05)", "FALSE" = "不显著"),
                     guide = "none") +
  # 对显著基因用实心，不显著用空心（通过fill和alpha控制）
  geom_point(data = subset(merged_scrna_rat, p_adj_scrna < 0.05),
             aes(fill = Category), shape = 21, size = 4, color = "black", stroke = 0.8) +
  geom_point(data = subset(merged_scrna_rat, p_adj_scrna >= 0.05),
             aes(fill = Category), shape = 21, size = 4, color = "black", stroke = 0.8, alpha = 0.5) +
  scale_fill_manual(values = category_colors, name = "基因分类") +
  # 基因标签
  geom_text_repel(aes(label = Gene), size = 3.2, max.overlaps = 30,
                  box.padding = 0.4, point.padding = 0.3,
                  segment.color = "grey50", segment.alpha = 0.6) +
  # 相关性注释
  annotate("label", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = paste0("Spearman r = ", rho_val, "\nP = ", p_val_cor),
           size = 4.5, fontface = "bold", family = "sans",
           fill = "white", alpha = 0.8) +
  # 象限标签
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = "方向一致", size = 4, color = "#D55E00", fontface = "bold") +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5,
           label = "方向一致", size = 4, color = "#0072B2", fontface = "bold") +
  labs(x = expression("scRNA-seq "*log[2]*"FC (小鼠 24h)"),
       y = expression("GSE97537 "*log[2]*"FC (大鼠 24h)"),
       title = "跨物种铜死亡基因表达变化一致性验证") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "grey80"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)))

# 保存
p1_file <- file.path(output_dir, "Fig1_scRNA_vs_GSE97537_scatter.png")
ggsave(p1_file, p1, width = 10, height = 8, dpi = 300)
cat("  已保存:", p1_file, "\n")

# ============================================================
# 图2: 热图 - 三数据集log2FC对比
# ============================================================
cat("\n>>> 绘制图2: 三数据集热图...\n")

# 准备热图矩阵
heat_matrix <- as.matrix(heatmap_df[, c("scRNA_24h", "GSE97537_24h", "GSE61616_7d")])
rownames(heat_matrix) <- heatmap_df[["Gene"]]
# 处理NA
heat_matrix[is.na(heat_matrix)] <- 0

# 颜色范围对称
max_abs <- max(abs(heat_matrix), na.rm = TRUE)

# 行注释 - 基因分类
annotation_row <- data.frame(
  Category = heatmap_df[["Category"]],
  row.names = rownames(heat_matrix)
)
annotation_colors_row <- list(
  Category = category_colors
)

# 列名
colnames(heat_matrix) <- c("scRNA-seq\n(24h小鼠)", "GSE97537\n(24h大鼠)", "GSE61616\n(7d小鼠)")

# 使用pheatmap
p2 <- pheatmap(
  heat_matrix,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  breaks = seq(-max_abs, max_abs, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors_row,
  display_numbers = round(heat_matrix, 2),
  number_format = "%.2f",
  number_color = "black",
  fontsize_number = 8,
  fontsize_row = 10,
  fontsize_col = 11,
  fontsize = 10,
  main = "铜死亡基因三数据集log2FC对比",
  angle_col = 45,
  na_col = "grey90",
  border_color = "grey80",
  cellwidth = 45,
  cellheight = 22,
  legend = TRUE,
  legend_breaks = round(seq(-max_abs, max_abs, length.out = 5), 1)
)

p2_file <- file.path(output_dir, "Fig2_three_dataset_log2FC_heatmap.png")
# pheatmap返回的是gtable对象，需用ggsave包装
png(p2_file, width = 10, height = max(6, nrow(heat_matrix) * 0.35 + 2), res = 300, units = "in")
grid::grid.newpage()
grid::grid.draw(p2[["gtable"]])
dev.off()
cat("  已保存:", p2_file, "\n")

# ============================================================
# 图3: 方向一致性条形图
# ============================================================
cat("\n>>> 绘制图3: 方向一致性条形图...\n")

# 准备数据
consistency_plot_data <- direction_summary %>%
  mutate(consistent = factor(consistent, levels = c("一致", "不一致")),
         Category = factor(Category, levels = names(category_colors)))

# 计算每个分类的基因总数
category_totals <- direction_df %>%
  group_by(Category) %>%
  summarise(Total_n = sum(!is.na(consistent)), .groups = "drop")

consistency_plot_data <- consistency_plot_data %>%
  left_join(category_totals, by = "Category")

p3 <- ggplot(consistency_plot_data, 
             aes(x = Category, y = Count, fill = consistent)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  # 比例标签
  geom_text(aes(label = Proportion),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold", color = "white") +
  # 总数标签在柱顶
  geom_text(data = subset(consistency_plot_data, consistent == "一致"),
            aes(y = Total, label = paste0("n=", Total)),
            vjust = -0.5, size = 3.5, color = "grey30") +
  scale_fill_manual(values = c("一致" = "#4DAF4A", "不一致" = "#E41A1C"),
                    name = "方向一致性") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "基因分类",
       y = "基因数",
       title = "铜死亡基因跨物种方向一致性分析") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "top",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

p3_file <- file.path(output_dir, "Fig3_direction_consistency_bar.png")
ggsave(p3_file, p3, width = 9, height = 7, dpi = 300)
cat("  已保存:", p3_file, "\n")

# ============================================================
# 图4: scRNA-seq铜死亡基因火山图
# ============================================================
cat("\n>>> 绘制图4: scRNA-seq铜死亡基因火山图...\n")

# 准备数据
volcano_data <- scrna_cu
volcano_data[["neg_log10_padj"]] <- -log10(pmax(volcano_data[["p_adjust"]], 1e-300))

# 确定阈值
log2fc_cutoff <- 0.25
padj_cutoff <- 0.05

# 标记上调/下调
volcano_data[["regulation"]] <- ifelse(
  volcano_data[["log2FC"]] > log2fc_cutoff & volcano_data[["p_adjust"]] < padj_cutoff, "上调",
  ifelse(volcano_data[["log2FC"]] < -log2fc_cutoff & volcano_data[["p_adjust"]] < padj_cutoff, "下调", "不显著")
)

# 统计
n_up <- sum(volcano_data[["regulation"]] == "上调", na.rm = TRUE)
n_down <- sum(volcano_data[["regulation"]] == "下调", na.rm = TRUE)
n_ns <- sum(volcano_data[["regulation"]] == "不显著", na.rm = TRUE)

cat("  上调:", n_up, " 下调:", n_down, " 不显著:", n_ns, "\n")

# 按log2FC绝对值排序，确保最显著的基因标注在最上层
volcano_data <- volcano_data[order(abs(volcano_data[["log2FC"]]), decreasing = TRUE), ]

p4 <- ggplot(volcano_data, aes(x = log2FC, y = neg_log10_padj)) +
  # 阈值线
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), 
             linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_hline(yintercept = -log10(padj_cutoff), 
             linetype = "dashed", color = "grey40", linewidth = 0.8) +
  # 散点：按分类着色
  geom_point(aes(fill = Category, shape = regulation != "不显著"),
             size = 3.5, color = "black", stroke = 0.6, alpha = 0.9) +
  scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 21),
                     guide = "none") +
  # 显著基因用实心填充，不显著用半透明
  geom_point(data = subset(volcano_data, regulation != "不显著"),
             aes(fill = Category), shape = 21, size = 4, color = "black", stroke = 0.8) +
  geom_point(data = subset(volcano_data, regulation == "不显著"),
             aes(fill = Category), shape = 21, size = 3, color = "grey60", stroke = 0.4, alpha = 0.4) +
  scale_fill_manual(values = category_colors, name = "基因分类") +
  # 基因标签 - 标注所有铜死亡基因
  geom_text_repel(aes(label = Gene, color = regulation),
                  size = 3.2, max.overlaps = 35,
                  box.padding = 0.4, point.padding = 0.3,
                  segment.color = "grey50", segment.alpha = 0.5,
                  fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("上调" = "#D55E00", "下调" = "#1B9E77", "不显著" = "grey40")) +
  # 统计文本
  annotate("label", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("上调: ", n_up, "\n下调: ", n_down, "\n不显著: ", n_ns),
           size = 4, fill = "white", alpha = 0.8,
           fontface = "bold") +
  labs(x = expression(log[2]*" Fold Change (MCAO vs Sham)"),
       y = expression(-log[10]*"(Adjusted P-value)"),
       title = "scRNA-seq铜死亡基因差异表达火山图 (GSE174574, 24h)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "grey80"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)))

p4_file <- file.path(output_dir, "Fig4_scRNA_cuproptosis_volcano.png")
ggsave(p4_file, p4, width = 10, height = 8, dpi = 300)
cat("  已保存:", p4_file, "\n")

# ============================================================
# 完成
# ============================================================
cat("\n============================================\n")
cat("所有图表已生成完毕！\n")
cat("输出目录:", output_dir, "\n")
cat("生成文件:\n")
cat("  1. Fig1_scRNA_vs_GSE97537_scatter.png\n")
cat("  2. Fig2_three_dataset_log2FC_heatmap.png\n")
cat("  3. Fig3_direction_consistency_bar.png\n")
cat("  4. Fig4_scRNA_cuproptosis_volcano.png\n")
cat("============================================\n")