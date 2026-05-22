# ============================================================
# Volcano Plot: GSE61616 Bulk RNA-seq Differential Expression
# Highlight cuproptosis & copper homeostasis genes
# ============================================================

# 1. 环境准备 ----
install_if_missing <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs) > 0) {
    install.packages(new_pkgs, dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

pkgs <- c("ggplot2", "ggrepel", "dplyr", "scales")
install_if_missing(pkgs)

# 2. 数据读取与预处理 ----
input_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery/results/L1_phenotype_anchoring/GSE61616_GEO2R_DEGs.csv"
output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery/figures/L1"

stopifnot(file.exists(input_file))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

data <- read.csv(input_file, row.names = 1, stringsAsFactors = FALSE)
stopifnot(nrow(data) > 0)

# 验证必要列存在
required_cols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# 将行名转为Gene Symbol列
data[["Gene"]] <- rownames(data)
stopifnot(all(nchar(data[["Gene"]]) > 0))

# 计算-log10(adj.P.Val)
data[["negLog10AdjP"]] <- -log10(data[["adj.P.Val"]])
data[["negLog10AdjP"]][is.infinite(data[["negLog10AdjP"]])] <- NA

# 3. 定义显著性标签与目标基因 ----
# 显著性阈值
adj_p_cutoff <- 0.05
log2fc_cutoff <- 0.585  # |log2FC| >= 0.585 ≈ 1.5倍变化

# 标记显著基因
data[["Significance"]] <- "Not Significant"
data[["Significance"]][data[["adj.P.Val"]] < adj_p_cutoff & data[["logFC"]] >= log2fc_cutoff] <- "Up-regulated"
data[["Significance"]][data[["adj.P.Val"]] < adj_p_cutoff & data[["logFC"]] <= -log2fc_cutoff] <- "Down-regulated"

# 将Significance转为factor并设定顺序
data[["Significance"]] <- factor(
  data[["Significance"]],
  levels = c("Down-regulated", "Not Significant", "Up-regulated")
)

# 定义铜死亡与铜稳态相关基因 (35个) - 人类基因名
cuproptosis_genes_human <- c(
  # Core cuproptosis genes (20)
  "ATP7A", "ATP7B", "CDKN2A", "COX17", "DBT",
  "DLD", "DLAT", "DLST", "FDX1", "GCSH",
  "GLS", "LIAS", "LIPT1", "LIPT2", "MTF1",
  "NFE2L2", "NLRP3", "PDHA1", "PDHB", "SLC31A1",
  # Copper homeostasis genes (15)
  "SLC31A2", "SLC11A2", "STEAP3", "ATOX1", "CCS",
  "COX11", "SCO1", "SCO2", "MT1A", "MT2A",
  "ALB", "CP", "SOD1", "SOD3", "COMMD1"
)

# 将人类基因名转换为小鼠同源基因名格式
# 规则：全大写 -> 首字母大写+其余小写（小鼠命名规范）
convert_to_mouse <- function(human_genes) {
  mouse_genes <- sapply(human_genes, function(g) {
    if (nchar(g) <= 1) return(g)
    # 全大写基因名转为小鼠格式：首字母大写，其余小写
    paste0(toupper(substring(g, 1, 1)), tolower(substring(g, 2)))
  }, USE.NAMES = FALSE)
  return(mouse_genes)
}

cuproptosis_genes <- convert_to_mouse(cuproptosis_genes_human)

# 检查数据中存在的目标基因
genes_in_data <- intersect(cuproptosis_genes, data[["Gene"]])
genes_not_found <- setdiff(cuproptosis_genes, data[["Gene"]])

# 输出人类-小鼠对应关系
gene_mapping <- data.frame(
  Human_Gene = cuproptosis_genes_human,
  Mouse_Gene = cuproptosis_genes,
  Found_In_Data = cuproptosis_genes %in% data[["Gene"]],
  stringsAsFactors = FALSE
)

if (length(genes_not_found) > 0) {
  message("警告：以下目标基因未在数据中找到：", paste(genes_not_found, collapse = ", "))
}

message("成功匹配 ", length(genes_in_data), "/", length(cuproptosis_genes), " 个目标基因")

# 输出未映射基因列表
unmapped_file <- file.path(output_dir, "unmapped_cuproptosis_genes.txt")
writeLines(genes_not_found, unmapped_file)
message("未映射基因已保存至：", unmapped_file)

# 输出人类-小鼠映射关系
mapping_file <- file.path(output_dir, "gene_mapping_human_mouse.csv")
write.csv(gene_mapping, mapping_file, row.names = FALSE)
message("人类-小鼠基因映射已保存至：", mapping_file)

# 标记目标基因
data[["IsTarget"]] <- ifelse(data[["Gene"]] %in% genes_in_data, "Cuproptosis Gene", "Other")

# 提取用于标签显示的基因：目标基因 + 极显著基因 (adj.P.Val < 0.001 且 |log2FC| > 1)
label_data <- data %>%
  filter(
    IsTarget == "Cuproptosis Gene" |
    (adj.P.Val < 0.001 & abs(logFC) > 1)
  ) %>%
  distinct(Gene, .keep_all = TRUE)

# 4. 绘图 ----
# 自定义颜色
custom_colors <- c(
  "Down-regulated" = "#1B9E77",
  "Not Significant" = "#BEBEBE",
  "Up-regulated" = "#D55E00"
)

# 创建火山图
p <- ggplot(data, aes(x = logFC, y = negLog10AdjP)) +
  # 基础散点层
  geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
  # 目标基因高亮层
  geom_point(
    data = filter(data, IsTarget == "Cuproptosis Gene"),
    color = "#0072B2", size = 3, alpha = 0.9
  ) +
  # 阈值线
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_hline(yintercept = -log10(adj_p_cutoff), linetype = "dashed", color = "gray40", linewidth = 0.8) +
  # 基因标签 (使用ggrepel避免重叠)
  geom_text_repel(
    data = label_data,
    aes(label = Gene),
    color = "#0072B2",
    size = 3.5,
    fontface = "bold",
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "#0072B2",
    segment.alpha = 0.5
  ) +
  # 颜色映射
  scale_color_manual(values = custom_colors) +
  # 坐标轴标签
  labs(
    x = expression(log[2] ~ "Fold Change"),
    y = expression(-log[10] ~ "Adjusted P-value"),
    color = "Expression"
  ) +
  # 主题设置 (出版级质量)
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black")
  ) +
  # 添加标题
  ggtitle("GSE61616: Differential Expression Analysis") +
  # 添加统计信息注释
  annotate(
    "text",
    x = max(data[["logFC"]], na.rm = TRUE) * 0.7,
    y = max(data[["negLog10AdjP"]], na.rm = TRUE) * 0.95,
    label = paste0(
      "Up: ", sum(data[["Significance"]] == "Up-regulated", na.rm = TRUE), "\n",
      "Down: ", sum(data[["Significance"]] == "Down-regulated", na.rm = TRUE), "\n",
      "Targets: ", length(genes_in_data)
    ),
    hjust = 1,
    vjust = 1,
    size = 4,
    fontface = "bold",
    color = "black"
  )

# 5. 保存图像 ----
pdf_file <- file.path(output_dir, "GSE61616_volcano_plot.pdf")
png_file <- file.path(output_dir, "GSE61616_volcano_plot.png")

# 使用Cairo设备避免PDF字体问题（Windows兼容性）
if (requireNamespace("Cairo", quietly = TRUE)) {
  Cairo::CairoPDF(
    file = pdf_file,
    width = 10,
    height = 7
  )
  print(p)
  dev.off()
} else {
  install.packages("Cairo", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(Cairo)
  Cairo::CairoPDF(
    file = pdf_file,
    width = 10,
    height = 7
  )
  print(p)
  dev.off()
}

# PNG保存使用ggsave
ggsave(
  filename = png_file,
  plot = p,
  width = 10,
  height = 7,
  dpi = 600
)

message("PDF保存至：", pdf_file)
message("PNG保存至：", png_file)

# 6. 输出统计摘要 ----
summary_stats <- data.frame(
  Metric = c(
    "Total Genes",
    "Up-regulated (adj.P.Val < 0.05, log2FC >= 0.585)",
    "Down-regulated (adj.P.Val < 0.05, log2FC <= -0.585)",
    "Cuproptosis Genes Found",
    "Cuproptosis Genes Not Found"
  ),
  Count = c(
    nrow(data),
    sum(data[["Significance"]] == "Up-regulated", na.rm = TRUE),
    sum(data[["Significance"]] == "Down-regulated", na.rm = TRUE),
    length(genes_in_data),
    length(genes_not_found)
  )
)

print(summary_stats)
message("\n火山图生成完成！")
