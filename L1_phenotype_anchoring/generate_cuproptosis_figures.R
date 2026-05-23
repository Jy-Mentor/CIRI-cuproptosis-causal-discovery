# =============================================================================
# CIRI-铜死亡项目: 3张出版级图表生成脚本
# 数据集:
#   - Bulk RNA-seq: GSE97537 (大鼠 24h MCAO vs Sham)
#   - scRNA-seq: GSE174574 (小鼠 24h MCAO vs Sham, 按细胞类型)
#
# 图表列表:
#   Fig1: Bulk火山图 - GSE97537铜死亡基因差异表达
#   Fig2: 细胞类型特异性铜死亡差异表达热图
#   Fig3: 各亚群关键铜死亡基因差异小提琴图
# =============================================================================

# =============================================================================
# 1. 环境准备 - 智能包安装
# =============================================================================
install_if_missing <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs) > 0) {
    cat(">>> 安装缺失的包:", paste(new_pkgs, collapse = ", "), "\n")
    install.packages(new_pkgs, dependencies = TRUE,
      repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

pkgs <- c("ggplot2", "ggrepel", "dplyr", "tidyr", "pheatmap",
          "RColorBrewer", "scales", "grid", "gridExtra", "viridis")
install_if_missing(pkgs)

# =============================================================================
# 2. 定义路径
# =============================================================================
base_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery"

# 注意: GSE97537_cuproptosis_DEGs.csv 在 bulk_GSE97537_analysis.R 生成后位于
#        L1_phenotype_anchoring 目录（相对于工作目录）
#       而 celltype 数据在 results/L1_phenotype_anchoring/ 目录
bulk_deg_file     <- file.path(base_dir, "L1_phenotype_anchoring/GSE97537_GEO2R_DEGs.csv")
bulk_cupro_file   <- file.path(base_dir, "L1_phenotype_anchoring/GSE97537_cuproptosis_DEGs.csv")
celltype_deg_file <- file.path(base_dir, "results/L1_phenotype_anchoring/celltype_cuproptosis_DEGs.csv")
violin_data_file  <- file.path(base_dir, "results/L1_phenotype_anchoring/celltype_cuproptosis_violin_data.csv")
output_dir        <- file.path(base_dir, "results/L1_phenotype_anchoring/figures")

# 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(">>> 创建输出目录:", output_dir, "\n")
}

# =============================================================================
# 3. 校验输入文件
# =============================================================================
stopifnot("Bulk DEG file not found" = file.exists(bulk_deg_file))
stopifnot("Bulk cuproptosis file not found" = file.exists(bulk_cupro_file))
stopifnot("Celltype DEG file not found" = file.exists(celltype_deg_file))
stopifnot("Violin data file not found" = file.exists(violin_data_file))

cat(">>> 所有输入文件校验通过\n")

# =============================================================================
# 4. 基因分类定义（小鼠基因名格式, 与 cross_species_validation_plots.R 一致）
# =============================================================================
gene_categories <- list(
  "\u94dc\u6b7b\u4ea1\u6838\u5fc3" = c("Fdx1", "Lias", "Dld", "Dlat", "Dlst", "Pdha1", "Pdhb",
                                      "Gls", "Gcsh", "Lipt1", "Lipt2", "Cdkn2a", "Nfe2l2", "Nlrp3"),
  "\u94dc\u79bb\u5b50\u8f6c\u8fd0" = c("Slc31a1", "Slc31a2", "Slc11a2", "Steap3", "Atp7a", "Atp7b"),
  "\u94dc\u4f34\u4fa3\u86cb\u767d" = c("Atox1", "Ccs", "Cox17", "Cox11", "Sco1", "Sco2"),
  "\u94dc\u5b58\u50a8\u7f13\u51b2" = c("Mt1a", "Mt2a", "Alb", "Cp", "Sod1", "Sod3"),
  "\u94dc\u4ee3\u8c22\u8c03\u63a7" = c("Commd1", "Mtf1")
)
# 英文 fallback (避免中文编码问题)
gene_categories_en <- list(
  "Cuproptosis_Core"       = c("Fdx1", "Lias", "Dld", "Dlat", "Dlst", "Pdha1", "Pdhb",
                               "Gls", "Gcsh", "Lipt1", "Lipt2", "Cdkn2a", "Nfe2l2", "Nlrp3"),
  "Copper_Transport"       = c("Slc31a1", "Slc31a2", "Slc11a2", "Steap3", "Atp7a", "Atp7b"),
  "Copper_Chaperone"       = c("Atox1", "Ccs", "Cox17", "Cox11", "Sco1", "Sco2"),
  "Copper_Storage_Buffer"  = c("Mt1a", "Mt2a", "Alb", "Cp", "Sod1", "Sod3"),
  "Copper_Metabolic_Regulation" = c("Commd1", "Mtf1")
)

# 5. 颜色方案
category_colors <- c(
  "\u94dc\u6b7b\u4ea1\u6838\u5fc3" = "#E41A1C",
  "\u94dc\u79bb\u5b50\u8f6c\u8fd0" = "#377EB8",
  "\u94dc\u4f34\u4fa3\u86cb\u767d" = "#4DAF4A",
  "\u94dc\u5b58\u50a8\u7f13\u51b2" = "#984EA3",
  "\u94dc\u4ee3\u8c22\u8c03\u63a7" = "#FF7F00"
)

# 构建基因到分类的映射表
gene_to_category <- data.frame(
  Gene = unlist(gene_categories),
  Category = rep(names(gene_categories), times = sapply(gene_categories, length)),
  stringsAsFactors = FALSE
)
# 设置分类为有序因子，保持颜色映射顺序
gene_to_category[["Category"]] <- factor(gene_to_category[["Category"]],
  levels = names(category_colors))

# =============================================================================
# 6. 自定义学术主题 theme_academic()
# =============================================================================
theme_academic <- function(base_size = 10, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size + 2, face = "bold", colour = "black"),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 10)),
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.text = element_text(size = base_size, colour = "black"),
      legend.title = element_text(size = base_size + 1, face = "bold", colour = "black"),
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5,
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size + 1, hjust = 0.5, colour = "gray40"),
      panel.spacing = unit(1, "lines")
    )
}

# =============================================================================
# 7. 读取数据
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat(">>> 读取数据...\n")

# 7.1 读取Bulk DEG全表
cat("  读取 Bulk DEG 全表:", bulk_deg_file, "\n")
bulk_all <- read.csv(bulk_deg_file, stringsAsFactors = FALSE)
stopifnot(nrow(bulk_all) > 0)
cat("    -> 总基因数:", nrow(bulk_all), "\n")

# 确保列名一致
colnames(bulk_all)[colnames(bulk_all) == "logFC"] <- "log2FC"

# 7.2 读取铜死亡基因子集
cat("  读取铜死亡基因子集:", bulk_cupro_file, "\n")
bulk_cu <- read.csv(bulk_cupro_file, stringsAsFactors = FALSE)
stopifnot(nrow(bulk_cu) > 0)
cat("    -> 铜死亡基因数:", nrow(bulk_cu), "\n")

# 处理数值列: log2FC, P.Value, adj.P.Val 可能是字符串格式
bulk_cu[["log2FC_num"]] <- as.numeric(as.character(bulk_cu[["log2FC"]]))
bulk_cu[["P.Value_num"]] <- as.numeric(as.character(bulk_cu[["P.Value"]]))
bulk_cu[["adj.P.Val_num"]] <- as.numeric(as.character(bulk_cu[["adj.P.Val"]]))
bulk_cu[["AveExpr_num"]] <- as.numeric(as.character(bulk_cu[["AveExpr"]]))
bulk_cu[["t_value_num"]] <- as.numeric(as.character(bulk_cu[["t_value"]]))

# 添加基因分类
bulk_cu <- merge(bulk_cu, gene_to_category, by.x = "Gene", by.y = "Gene", all.x = TRUE)
# 未匹配的基因标记为"Other"
bulk_cu[["Category"]][is.na(bulk_cu[["Category"]])] <- "Other"

# 7.3 读取细胞类型DEG数据
cat("  读取细胞类型DEG数据:", celltype_deg_file, "\n")
ct_deg <- read.csv(celltype_deg_file, stringsAsFactors = FALSE)
stopifnot(nrow(ct_deg) > 0)
cat("    -> 总行数:", nrow(ct_deg), "\n")
cat("    -> 细胞类型:", paste(unique(ct_deg[["cell_type"]]), collapse = ", "), "\n")
cat("    -> 基因数:", length(unique(ct_deg[["gene"]])), "\n")

# 7.4 读取小提琴图数据
cat("  读取小提琴图数据:", violin_data_file, "\n")
vio_data <- read.csv(violin_data_file, stringsAsFactors = FALSE)
stopifnot(nrow(vio_data) > 0)
cat("    -> 总行数:", nrow(vio_data), "\n")

cat(">>> 数据读取完成\n")

# =============================================================================
# Fig 1: Bulk火山图 (GSE97537)
# =============================================================================
cat("\n>>> 绘制 Fig1: Bulk火山图 (GSE97537)...\n")

# 计算 -log10(adj.P.Val)
bulk_all[["negLog10P"]] <- -log10(pmax(bulk_all[["adj.P.Val"]], 1e-300))

# 标记显著性
log2fc_cutoff <- 0.585
adj_p_cutoff <- 0.05
bulk_all[["Significance"]] <- "Not Significant"
bulk_all[["Significance"]][bulk_all[["adj.P.Val"]] < adj_p_cutoff &
                            bulk_all[["log2FC"]] >= log2fc_cutoff] <- "Up-regulated"
bulk_all[["Significance"]][bulk_all[["adj.P.Val"]] < adj_p_cutoff &
                            bulk_all[["log2FC"]] <= -log2fc_cutoff] <- "Down-regulated"
bulk_all[["Significance"]] <- factor(bulk_all[["Significance"]],
  levels = c("Down-regulated", "Not Significant", "Up-regulated"))

# 标记铜死亡基因
bulk_all[["IsCuproptosis"]] <- bulk_all[["Gene"]] %in% gene_to_category[["Gene"]]

# 合并分类信息到bulk_all（用于铜死亡基因着色）
bulk_all <- merge(bulk_all, gene_to_category[, c("Gene", "Category")],
  by = "Gene", all.x = TRUE)
# 非铜死亡基因分类设为"Other"
bulk_all[["Category"]][is.na(bulk_all[["Category"]])] <- "Other"
# 将Category转为因子，确保颜色映射正确
all_cats <- c(names(category_colors), "Other")
bulk_all[["Category"]] <- factor(bulk_all[["Category"]], levels = all_cats)

# 统计
n_up <- sum(bulk_all[["Significance"]] == "Up-regulated", na.rm = TRUE)
n_down <- sum(bulk_all[["Significance"]] == "Down-regulated", na.rm = TRUE)
n_cupro_found <- sum(bulk_all[["IsCuproptosis"]], na.rm = TRUE)
cat("  统计: 上调=", n_up, ", 下调=", n_down, ", 铜死亡基因检出=", n_cupro_found, "\n")

# 铜死亡基因标签数据
cupro_labels <- bulk_all[bulk_all[["IsCuproptosis"]] == TRUE, ]
# 去除重复Gene
cupro_labels <- cupro_labels[!duplicated(cupro_labels[["Gene"]]), ]

# 扩展颜色映射 - 添加"Other"使用灰色
all_colors <- c(category_colors, "Other" = "#BEBEBE")

# 构建绘图颜色映射
plot_colors <- c(
  "Down-regulated" = "#1B9E77",
  "Not Significant" = "#BEBEBE",
  "Up-regulated" = "#D55E00"
)

# ====== 分层绘制火山图 ======
p1 <- ggplot() +
  # 底层: 非铜死亡基因（灰色背景）
  geom_point(
    data = bulk_all[!bulk_all[["IsCuproptosis"]], ],
    aes(x = log2FC, y = negLog10P, color = Significance),
    alpha = 0.5, size = 1.2
  ) +
  scale_color_manual(
    values = plot_colors,
    name = "Expression"
  ) +
  # 铜死亡基因光晕层
  geom_point(
    data = bulk_all[bulk_all[["IsCuproptosis"]], ],
    aes(x = log2FC, y = negLog10P),
    fill = "#0072B2", color = "#0072B2", size = 5, alpha = 0.2, shape = 21
  ) +
  # 铜死亡基因主层 - 按分类着色
  geom_point(
    data = bulk_all[bulk_all[["IsCuproptosis"]], ],
    aes(x = log2FC, y = negLog10P, fill = Category),
    size = 3.5, shape = 21, color = "black", stroke = 0.6
  ) +
  scale_fill_manual(
    values = all_colors,
    name = "Gene Category",
    breaks = names(category_colors)
  ) +
  # 阈值线
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff),
    linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_hline(yintercept = -log10(adj_p_cutoff),
    linetype = "dashed", color = "gray40", linewidth = 0.6) +
  # 基因标签
  geom_text_repel(
    data = cupro_labels,
    aes(x = log2FC, y = negLog10P, label = Gene, color = Category),
    size = 3.2, fontface = "bold", max.overlaps = 50,
    box.padding = 0.4, point.padding = 0.3,
    segment.color = "grey50", segment.alpha = 0.5,
    show.legend = FALSE
  ) +
  # 统计标注
  annotate("label",
    x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
    label = paste0(
      "Up: ", n_up, "\n",
      "Down: ", n_down, "\n",
      "Cuproptosis: ", n_cupro_found, "/", length(unlist(gene_categories))
    ),
    size = 4.2, fontface = "bold", fill = "white", alpha = 0.85
  ) +
  labs(
    x = expression(log[2] ~ "Fold Change (MCAO vs Sham)"),
    y = expression(-log[10] ~ "(Adjusted P-value)"),
    title = "GSE97537: Cuproptosis Gene Expression (Rat 24h MCAO)"
  ) +
  theme_academic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.margin = margin(0.5, 1.5, 0.5, 0.5, "cm")
  ) +
  guides(
    color = guide_legend(order = 1, title = "Significance",
      override.aes = list(size = 2.5, alpha = 0.8)),
    fill = guide_legend(order = 2, title = "Cuproptosis Category",
      override.aes = list(size = 3, shape = 21, color = "black"))
  )

# 保存
fig1_file <- file.path(output_dir, "Fig1_bulk_volcano_GSE97537.png")
ggsave(fig1_file, p1, width = 10, height = 8, dpi = 600)
cat("  已保存:", fig1_file, "\n")

# 同时保存PDF矢量版
fig1_pdf <- file.path(output_dir, "Fig1_bulk_volcano_GSE97537.pdf")
ggsave(fig1_pdf, p1, width = 10, height = 8, dpi = 300, device = cairo_pdf)
cat("  已保存:", fig1_pdf, "\n")

# =============================================================================
# Fig 2: 细胞类型特异性铜死亡差异表达热图
# =============================================================================
cat("\n>>> 绘制 Fig2: 细胞类型热图...\n")

# 准备热图矩阵: 行=基因, 列=细胞类型, 值=log2FC
ct_deg[["gene_cat"]] <- gene_to_category[["Category"]][
  match(ct_deg[["gene"]], gene_to_category[["Gene"]])]
# 未匹配基因标记
ct_deg[["gene_cat"]][is.na(ct_deg[["gene_cat"]])] <- "Other"

# 对每个基因在每个细胞类型中取唯一值（防止重复）
ct_unique <- ct_deg %>%
  group_by(gene, cell_type) %>%
  summarise(
    log2FC = mean(log2FC, na.rm = TRUE),
    p_adjust = min(p_adjust, na.rm = TRUE),
    gene_cat = dplyr::first(gene_cat),
    .groups = "drop"
  )

# 转换为宽格式矩阵
heat_matrix <- ct_unique %>%
  select(gene, cell_type, log2FC) %>%
  pivot_wider(names_from = cell_type, values_from = log2FC) %>%
  as.data.frame()

rownames(heat_matrix) <- heat_matrix[["gene"]]
heat_matrix[["gene"]] <- NULL

# 转换为数值矩阵
heat_matrix <- as.matrix(heat_matrix)
mode(heat_matrix) <- "numeric"

# 获取基因分类顺序（按分类排序）
gene_order <- gene_to_category[["Gene"]][gene_to_category[["Gene"]] %in% rownames(heat_matrix)]
# 添加不在分类列表中的基因
extra_genes <- rownames(heat_matrix)[!rownames(heat_matrix) %in% gene_order]
gene_order <- c(gene_order, extra_genes)
heat_matrix <- heat_matrix[gene_order, , drop = FALSE]

# 对称颜色范围
max_abs_lfc <- max(abs(heat_matrix), na.rm = TRUE)
if (is.finite(max_abs_lfc) && max_abs_lfc > 0) {
  color_breaks <- seq(-max_abs_lfc, max_abs_lfc, length.out = 101)
} else {
  color_breaks <- seq(-1, 1, length.out = 101)
}

# 构建显著性矩阵 (p_adjust < 0.05 且 |log2FC| > 0.25)
sig_matrix <- ct_unique %>%
  select(gene, cell_type, p_adjust, log2FC) %>%
  mutate(
    sig_up = log2FC > 0.25 & p_adjust < 0.05,
    sig_down = log2FC < -0.25 & p_adjust < 0.05,
    sig_symbol = case_when(
      sig_up ~ "*",
      sig_down ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(gene, cell_type, sig_symbol) %>%
  pivot_wider(names_from = cell_type, values_from = sig_symbol) %>%
  as.data.frame()

rownames(sig_matrix) <- sig_matrix[["gene"]]
sig_matrix[["gene"]] <- NULL
sig_matrix <- as.matrix(sig_matrix)
sig_matrix <- sig_matrix[gene_order, , drop = FALSE]
# 填充NA
sig_matrix[is.na(sig_matrix)] <- ""

# 行注释: 基因分类
annotation_row <- data.frame(
  Category = gene_to_category[["Category"]][
    match(rownames(heat_matrix), gene_to_category[["Gene"]])],
  row.names = rownames(heat_matrix)
)
annotation_row[["Category"]][is.na(annotation_row[["Category"]])] <- "Other"
annotation_row[["Category"]] <- factor(annotation_row[["Category"]],
  levels = c(names(category_colors), "Other"))

# 行注释颜色
annotation_colors <- list(
  Category = c(category_colors, "Other" = "#BEBEBE")
)

# 处理NA值
heat_matrix[is.na(heat_matrix)] <- 0

# 确保display_numbers矩阵与heat_matrix维度一致
# 使用显著性标记代替数值
display_mat <- sig_matrix

# 显示数值的函数: 显示log2FC + 显著性标记
display_numbers <- matrix("", nrow = nrow(heat_matrix), ncol = ncol(heat_matrix))
for (i in seq_len(nrow(heat_matrix))) {
  for (j in seq_len(ncol(heat_matrix))) {
    val <- round(heat_matrix[i, j], 2)
    sig <- sig_matrix[i, j]
    display_numbers[i, j] <- paste0(sprintf("%.2f", val), sig)
  }
}
rownames(display_numbers) <- rownames(heat_matrix)
colnames(display_numbers) <- colnames(heat_matrix)

# 创建热图
p2 <- pheatmap(
  heat_matrix,
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  breaks = color_breaks,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  display_numbers = display_numbers,
  number_format = "%.2f",
  number_color = "black",
  fontsize_number = 6,
  fontsize_row = 9,
  fontsize_col = 10,
  fontsize = 8,
  main = "Celltype-specific Cuproptosis Gene Expression (log2FC)",
  angle_col = 45,
  na_col = "grey90",
  border_color = "grey80",
  cellwidth = 35,
  cellheight = 18,
  legend = TRUE,
  legend_breaks = round(seq(-max_abs_lfc, max_abs_lfc, length.out = 5), 1),
  legend_labels = as.character(round(seq(-max_abs_lfc, max_abs_lfc, length.out = 5), 1))
)

# 保存热图
fig2_file <- file.path(output_dir, "Fig2_celltype_cuproptosis_heatmap.png")
# 计算合理的图高
fig_height <- max(6, nrow(heat_matrix) * 0.3 + 3)
png(fig2_file, width = 12, height = fig_height, res = 600, units = "in")
grid::grid.newpage()
grid::grid.draw(p2[["gtable"]])
dev.off()
cat("  已保存:", fig2_file, "\n")

# PDF版
fig2_pdf <- file.path(output_dir, "Fig2_celltype_cuproptosis_heatmap.pdf")
pdf(fig2_pdf, width = 12, height = fig_height)
grid::grid.newpage()
grid::grid.draw(p2[["gtable"]])
dev.off()
cat("  已保存:", fig2_pdf, "\n")

# =============================================================================
# Fig 3: 各亚群差异小提琴图 (关键铜死亡基因)
# =============================================================================
cat("\n>>> 绘制 Fig3: 关键基因小提琴图...\n")

# 确定最显著的8个基因: 在最多细胞类型中 p_adjust < 0.05 的基因
gene_sig_count <- ct_deg %>%
  filter(p_adjust < 0.05) %>%
  group_by(gene) %>%
  summarise(n_sig_celltypes = n_distinct(cell_type), .groups = "drop") %>%
  arrange(desc(n_sig_celltypes))

cat("  基因显著性排名（跨细胞类型）:\n")
print(gene_sig_count, n = 15)

# 选择前8个基因
top_n_genes <- 8
if (nrow(gene_sig_count) < top_n_genes) {
  top_n_genes <- nrow(gene_sig_count)
}
top_genes <- gene_sig_count[["gene"]][seq_len(top_n_genes)]
cat("  选择前", top_n_genes, "个基因:", paste(top_genes, collapse = ", "), "\n")

# 如果显著基因不足8个，补充其他基因
if (length(top_genes) < 8) {
  # 补充差异最大的基因（按|log2FC|排序）
  remaining <- ct_deg %>%
    filter(!gene %in% top_genes) %>%
    group_by(gene) %>%
    summarise(max_abs_lfc = max(abs(log2FC), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_lfc)) %>%
    head(8 - length(top_genes))
  top_genes <- c(top_genes, remaining[["gene"]])
  cat("  补充基因:", paste(remaining[["gene"]], collapse = ", "), "\n")
}

# 过滤小提琴数据
vio_sub <- vio_data[vio_data[["gene"]] %in% top_genes, ]
stopifnot(nrow(vio_sub) > 0)
cat("  小提琴数据行数:", nrow(vio_sub), "\n")

# 确保condition为因子
vio_sub[["condition"]] <- factor(vio_sub[["condition"]], levels = c("Sham", "MCAO"))

# 为每个基因在每个细胞类型中计算wilcox p值
stat_results <- vio_sub %>%
  group_by(gene, cell_type) %>%
  summarise(
    p_val = tryCatch({
      test <- wilcox.test(expression ~ condition, data = pick(1L))
      test[["p.value"]]
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      is.na(p_val) ~ "",
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    # 注释位置: 在每组最大值上方
    p_y = NA_real_
  )

# 计算每个面板的最大值用于注释位置
max_expr <- vio_sub %>%
  group_by(gene, cell_type) %>%
  summarise(max_val = max(expression, na.rm = TRUE), .groups = "drop")

stat_results <- stat_results %>%
  left_join(max_expr, by = c("gene", "cell_type")) %>%
  mutate(p_y = max_val * 1.1 + 0.1)

# 设置基因因子顺序（按显著性排序）
vio_sub[["gene"]] <- factor(vio_sub[["gene"]], levels = top_genes)

# 调整细胞类型顺序（按表达模式聚类，但保持可读性）
# 按第一个基因的表达均值排序
first_gene <- top_genes[1]
ct_order <- vio_sub %>%
  filter(gene == first_gene) %>%
  group_by(cell_type) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_expr))
vio_sub[["cell_type"]] <- factor(vio_sub[["cell_type"]],
  levels = ct_order[["cell_type"]])

# 自定义颜色
condition_colors <- c("Sham" = "#4DAF4A", "MCAO" = "#E41A1C")

# 计算分面中的最大y值，用于统一y轴范围
y_max_all <- max(vio_sub[["expression"]], na.rm = TRUE) * 1.3

p3 <- ggplot(vio_sub, aes(x = cell_type, y = expression, fill = condition)) +
  # 小提琴
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7,
    position = position_dodge(width = 0.8), color = "black", linewidth = 0.3) +
  # 箱线图内嵌
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.8),
    alpha = 0.9, outlier.shape = NA, color = "black", linewidth = 0.3) +
  # 散点（抖动）
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 0.3, alpha = 0.2, color = "grey30") +
  # 分面
  facet_wrap(~ gene, nrow = 2, scales = "free_y") +
  # 颜色
  scale_fill_manual(values = condition_colors, name = "Condition") +
  # 统计注释
  geom_text(
    data = stat_results,
    aes(x = cell_type, y = p_y, label = p_label),
    inherit.aes = FALSE,
    size = 3, vjust = 0.5, color = "black", fontface = "bold"
  ) +
  # 标签
  labs(
    x = "Cell Type",
    y = "Normalized Expression (log1p)",
    title = "Cuproptosis Gene Expression Across Cell Types (GSE174574, 24h MCAO)"
  ) +
  # 主题
  theme_academic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(size = 9, face = "bold.italic"),
    legend.position = "top",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    panel.spacing = unit(0.8, "lines")
  )

# 保存
fig3_file <- file.path(output_dir, "Fig3_celltype_violin.png")
ggsave(fig3_file, p3, width = 14, height = 10, dpi = 600)
cat("  已保存:", fig3_file, "\n")

# PDF版
fig3_pdf <- file.path(output_dir, "Fig3_celltype_violin.pdf")
ggsave(fig3_pdf, p3, width = 14, height = 10, dpi = 300, device = cairo_pdf)
cat("  已保存:", fig3_pdf, "\n")

# =============================================================================
# 完成
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("所有图表已生成完毕！\n")
cat("输出目录:", output_dir, "\n")
cat("生成文件:\n")
cat("  1. Fig1_bulk_volcano_GSE97537.png (10x8\", 600dpi)\n")
cat("  2. Fig2_celltype_cuproptosis_heatmap.png (12x", round(fig_height, 1), "\", 600dpi)\n", sep = "")
cat("  3. Fig3_celltype_violin.png (14x10\", 600dpi)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")