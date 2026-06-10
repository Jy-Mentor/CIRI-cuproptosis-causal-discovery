# GSE97537 差异表达分析脚本
# 分析: MCAO vs Sham (脑缺血再灌注损伤模型)
# 筛选条件: |logFC| > 1, FDR < 0.05

# ==================== 1. 包安装与加载 ====================
cat("正在检查和安装必要的R包...\n")

packages <- c("limma", "ggplot2", "dplyr", "tidyr")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Installing package: ", pkg, "\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

for (pkg in packages) {
  install_if_missing(pkg)
}

cat("所有包加载完成!\n\n")

# ==================== 2. 设置工作目录和文件路径 ====================
work_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
setwd(work_dir)

# 文件路径
series_matrix_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/AI 代码编写/GSE23160与97537数据集合并分析/GSE97537_series_matrix.txt"
platform_file <- "C:/Users/Jy-Mentor-7/Downloads/GPL1355-10794.txt"

# 输出目录
output_dir <- file.path(work_dir, "GSE97537_limma_results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==================== 3. 读取系列矩阵数据 ====================
cat("读取 GSE97537 系列矩阵数据...\n")

# 读取表达矩阵 (跳过元数据行，从 !series_matrix_table_begin 后开始)
lines <- readLines(series_matrix_file)
data_start <- which(grepl("!series_matrix_table_begin", lines)) + 1
data_end <- which(grepl("!series_matrix_table_end", lines)) - 1

# 提取表达数据
expr_data <- read.delim(series_matrix_file, 
                        skip = data_start - 1, 
                        nrows = data_end - data_start + 1,
                        header = TRUE, 
                        check.names = FALSE,
                        row.names = 1)

# 转换为矩阵
expr_matrix <- as.matrix(expr_data)
cat(paste0("表达矩阵维度: ", nrow(expr_matrix), " probes x ", ncol(expr_matrix), " samples\n"))

# ==================== 4. 读取样本分组信息 ====================
cat("提取样本分组信息...\n")

# 从系列矩阵文件提取样本标题
sample_titles_line <- grep("^!Sample_title", lines, value = TRUE)
sample_titles <- unlist(strsplit(sample_titles_line, "\t"))[-1]
sample_titles <- gsub('"', '', sample_titles)

# 创建样本分组
sample_info <- data.frame(
  sample_id = colnames(expr_matrix),
  title = sample_titles,
  stringsAsFactors = FALSE
)

# 根据样本标题确定分组 (MCAO vs Sham)
sample_info$group <- ifelse(grepl("MCAO", sample_info$title, ignore.case = TRUE), 
                            "MCAO", "Sham")
sample_info$group <- factor(sample_info$group, levels = c("Sham", "MCAO"))

cat("样本分组情况:\n")
print(table(sample_info$group))

# ==================== 5. 读取平台注释信息 ====================
cat("读取 GPL1355 平台注释...\n")

# 读取平台注释文件 (跳过头部注释)
plat_lines <- readLines(platform_file)
header_line <- which(grepl("^ID\t", plat_lines))[1]
platform_annot <- read.delim(platform_file, 
                             skip = header_line - 1, 
                             header = TRUE,
                             check.names = FALSE,
                             stringsAsFactors = FALSE,
                             quote = "",
                             fill = TRUE)

# 提取关键列: Probe ID, Gene Symbol, Gene Title, Entrez ID
cat("平台注释列名:\n")
print(names(platform_annot)[1:15])

# 提取基因注释
annot_table <- data.frame(
  ID = platform_annot$ID,
  Gene_Symbol = platform_annot$`Gene Symbol`,
  Gene_Title = platform_annot$`Gene Title`,
  Entrez_ID = platform_annot$`ENTREZ_GENE_ID`,
  stringsAsFactors = FALSE
)

cat(paste0("平台注释基因数: ", nrow(annot_table), "\n\n"))

# ==================== 6. limma 差异表达分析 ====================
cat("开始进行 limma 差异表达分析...\n")

# 创建设计矩阵
design <- model.matrix(~0 + sample_info$group)
colnames(design) <- levels(sample_info$group)

cat("设计矩阵:\n")
print(design)

# 创建对比矩阵 (MCAO vs Sham)
contrast.matrix <- makeContrasts(MCAOvsSham = MCAO - Sham, levels = design)

cat("\n对比矩阵:\n")
print(contrast.matrix)

# 线性模型拟合
fit <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 提取所有结果
cat("\n提取差异表达结果...\n")
results <- topTable(fit2, number = Inf, adjust.method = "fdr", sort.by = "p")

# ==================== 7. 添加基因注释到结果 ====================
cat("添加基因注释信息...\n")

# 合并注释信息
results$ID <- rownames(results)
results_annotated <- merge(results, annot_table, by = "ID", all.x = TRUE)

# 重新排序
results_annotated <- results_annotated[order(results_annotated$adj.P.Val), ]

# ==================== 8. 筛选显著差异基因 ====================
cat("\n筛选显著差异基因 (|logFC| > 1, FDR < 0.05)...\n")

# 定义显著性
deg_results <- results_annotated %>%
  mutate(
    significance = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# 筛选满足条件的差异基因
sig_genes <- deg_results %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05)

# 统计结果
cat("\n===== 差异表达分析结果汇总 =====\n")
cat(paste0("总探针数: ", nrow(deg_results), "\n"))
cat(paste0("显著差异基因总数: ", nrow(sig_genes), "\n"))
cat(paste0("  - 上调基因 (logFC > 1, FDR < 0.05): ", sum(sig_genes$logFC > 0), "\n"))
cat(paste0("  - 下调基因 (logFC < -1, FDR < 0.05): ", sum(sig_genes$logFC < 0), "\n"))
cat("================================\n\n")

# ==================== 9. 保存结果文件 ====================
cat("保存结果文件...\n")

# 保存所有基因结果
write.table(deg_results, 
            file = file.path(output_dir, "GSE97537_all_genes_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 保存显著差异基因
write.table(sig_genes,
            file = file.path(output_dir, "GSE97537_DEGs_logFC1_FDR0.05.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 保存上调基因
up_genes <- sig_genes %>% filter(logFC > 0)
write.table(up_genes,
            file = file.path(output_dir, "GSE97537_upregulated_genes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 保存下调基因
down_genes <- sig_genes %>% filter(logFC < 0)
write.table(down_genes,
            file = file.path(output_dir, "GSE97537_downregulated_genes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 保存仅含Gene Symbol的列表 (用于后续分析)
sig_gene_symbols <- sig_genes$Gene_Symbol[!is.na(sig_genes$Gene_Symbol) & sig_genes$Gene_Symbol != ""]
write.table(sig_gene_symbols,
            file = file.path(output_dir, "GSE97537_DEG_gene_symbols.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(paste0("结果文件已保存到: ", output_dir, "\n"))

# ==================== 10. 生成火山图 ====================
cat("生成火山图...\n")

# 准备火山图数据
volcano_data <- deg_results %>%
  mutate(
    logP = -log10(P.Value),
    color = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# 设置颜色
colors <- c("Upregulated" = "#FF4040", "Downregulated" =="#4169E1", "Not significant" = "grey80")

# 创建火山图
p <- ggplot(volcano_data, aes(x = logFC, y = logP, color = color)) +
  geom_point(size = 1.2, alpha = 0.9) +
  scale_color_manual(values = colors, 
                     name = "Significance",
                     labels = c("Not significant" = "Not significant",
                                "Upregulated" = paste0("Up (n=", sum(volcano_data$color == "Upregulated"), ")"),
                                "Downregulated" = paste0("Down (n=", sum(volcano_data$color == "Downregulated"), ")"))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  labs(title = "GSE97537: MCAO vs Sham Differential Expression",
       subtitle = paste0("|logFC| > 1, FDR < 0.05 | DEGs: ", nrow(sig_genes)),
       x = expression("log"[2]~"Fold Change"),
       y = expression("-log"[10]~"P-value")) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

# 保存火山图
ggsave(file.path(output_dir, "GSE97537_volcano_plot.png"), 
       plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "GSE97537_volcano_plot.pdf"), 
       plot = p, width = 10, height = 8)

cat("火山图已保存!\n")

# ==================== 11. 显示Top差异基因 ====================
cat("\n===== Top 20 差异表达基因 =====\n")
top20 <- sig_genes %>%
  select(ID, Gene_Symbol, logFC, P.Value, adj.P.Val) %>%
  head(20)
print(top20, row.names = FALSE)

# ==================== 12. 总结报告 ====================
cat("\n===== 分析完成 =====\n")
cat(paste0("输出目录: ", output_dir, "\n"))
cat("\n生成文件:\n")
cat("1. GSE97537_all_genes_results.txt - 所有基因结果\n")
cat("2. GSE97537_DEGs_logFC1_FDR0.05.txt - 显著差异基因\n")
cat("3. GSE97537_upregulated_genes.txt - 上调基因\n")
cat("4. GSE97537_downregulated_genes.txt - 下调基因\n")
cat("5. GSE97537_DEG_gene_symbols.txt - 差异基因Symbol列表\n")
cat("6. GSE97537_volcano_plot.png/pdf - 火山图\n")
cat("\n分析参数:\n")
cat("- 对比组: MCAO vs Sham\n")
cat("- 筛选条件: |logFC| > 1, FDR < 0.05\n")
cat("- 差异基因总数: ", nrow(sig_genes), "\n")
cat("==================\n")
