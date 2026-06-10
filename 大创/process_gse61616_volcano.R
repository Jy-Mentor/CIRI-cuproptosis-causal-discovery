# 生物信息学分析脚本：GSE61616差异表达分析与火山图绘制

# 安装必要的包
# 设置CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

install_if_missing("tidyverse")
install_if_missing("limma")
install_if_missing("here")

# 加载包
library(tidyverse)
library(limma)
library(here)

# 设置工作目录
setwd("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创")

# 第一步：数据读取与limma质控
cat("Step 1: 数据读取与质控...\n")
data_path <- "GSE61616.top.table (1).tsv"
data <- readr::read_tsv(data_path)

# 处理列名（确保列名正确）
# 先检查原始列名
cat("  原始列名:", paste(colnames(data), collapse = ", "), "\n")

# 手动处理列名，确保关键列存在
colnames(data) <- gsub("\\.", "_", colnames(data))  # 将点替换为下划线
colnames(data) <- gsub("\\s+", "_", colnames(data))  # 移除空格

# 确保关键列名正确
if ("Gene_symbol" %in% colnames(data)) {
  # 列名已经是Gene_symbol
} else if ("Gene.symbol" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Gene.symbol"] <- "Gene_symbol"
} else if ("gene_symbol" %in% colnames(data)) {
  colnames(data)[colnames(data) == "gene_symbol"] <- "Gene_symbol"
}

if ("adj_P_Val" %in% colnames(data)) {
  # 列名已经是adj_P_Val
} else if ("adj.P.Val" %in% colnames(data)) {
  colnames(data)[colnames(data) == "adj.P.Val"] <- "adj_P_Val"
} else if ("adj_p_val" %in% colnames(data)) {
  colnames(data)[colnames(data) == "adj_p_val"] <- "adj_P_Val"
}

# 检查关键列是否存在
required_cols <- c("Gene_symbol", "logFC", "adj_P_Val")
if (!all(required_cols %in% colnames(data))) {
  stop("关键列不存在，请检查输入文件格式")
}

# 显示处理后的列名
cat("  处理后列名:", paste(colnames(data), collapse = ", "), "\n")

# 探针去重：保留|logFC|最大的记录
cat("  探针去重...\n")
data_clean <- data %>%
  filter(!is.na(Gene_symbol) & Gene_symbol != "") %>%
  group_by(Gene_symbol) %>%
  slice_max(abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup()

# 第二步：差异表达基因筛选
cat("Step 2: 差异表达基因筛选...\n")
logFC_threshold <- 0.2
pval_threshold <- 0.05

data_clean <- data_clean %>%
  mutate(Regulation = case_when(
    logFC >= logFC_threshold & adj_P_Val < pval_threshold ~ "Up",
    logFC <= -logFC_threshold & adj_P_Val < pval_threshold ~ "Down",
    TRUE ~ "NotSig"
  ))

# 统计差异基因数量
up_genes <- sum(data_clean$Regulation == "Up")
down_genes <- sum(data_clean$Regulation == "Down")
total_sig_genes <- up_genes + down_genes

cat(paste("  显著上调基因数:", up_genes, "\n"))
cat(paste("  显著下调基因数:", down_genes, "\n"))
cat(paste("  总差异基因数:", total_sig_genes, "\n"))

# 第三步：小鼠→人类同源基因映射
cat("Step 3: 小鼠→人类同源基因映射...\n")

# 使用本地映射文件
cat("  使用本地映射文件进行基因映射...\n")
mapping_file <- "大鼠 小鼠 人类映射库.txt"

# 读取映射文件
mapping_data <- readr::read_tsv(mapping_file, comment = "#")

# 提取小鼠基因到人类基因的映射
homologs <- mapping_data %>%
  select(MOUSE_ORTHOLOG_SYMBOL, HUMAN_ORTHOLOG_SYMBOL) %>%
  filter(!is.na(MOUSE_ORTHOLOG_SYMBOL) & MOUSE_ORTHOLOG_SYMBOL != "") %>%
  filter(!is.na(HUMAN_ORTHOLOG_SYMBOL) & HUMAN_ORTHOLOG_SYMBOL != "") %>%
  rename(external_gene_name = MOUSE_ORTHOLOG_SYMBOL,
         hsapiens_homolog_associated_gene_name = HUMAN_ORTHOLOG_SYMBOL) %>%
  group_by(external_gene_name) %>%
  slice(1) %>%  # 一对多保留第一个
  ungroup()

# 合并映射结果
data_mapped <- data_clean %>%
  left_join(homologs, by = c("Gene_symbol" = "external_gene_name")) %>%
  rename(Mouse_Symbol = Gene_symbol,
         Human_Symbol = hsapiens_homolog_associated_gene_name)

# 第四步：无标签火山图绘制
cat("Step 4: 绘制无标签火山图...\n")

# 计算-log10(adj.P.Val)
data_mapped <- data_mapped %>%
  mutate(neg_log10_pval = -log10(adj_P_Val))

# 绘制火山图
volcano_plot <- ggplot(data_mapped, aes(x = logFC, y = neg_log10_pval, color = Regulation)) +
  geom_point(aes(alpha = Regulation), size = 2, stroke = 0) +
  scale_alpha_manual(values = c(
    "Up" = 0.8,
    "Down" = 0.8,
    "NotSig" = 0.5
  )) +
  # 阈值线
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), 
             linetype = "dashed", color = "gray", linewidth = 0.5) +
  geom_hline(yintercept = -log10(pval_threshold), 
             color = "gray", linewidth = 0.5) +
  # 颜色映射
  scale_color_manual(values = c(
    "Up" = "#DC143C",
    "Down" = "#4169E1",
    "NotSig" = "#D3D3D3"
  )) +
  # 标题与标签
  labs(
    title = "GSE61616 DEG Volcano Plot ",
    x = "log2 Fold Change",
    y = "-log10 (Adjusted P-value)"
  ) +
  # 主题
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# 输出火山图
output_dir <- "."
pdf_file <- file.path(output_dir, "GSE61616_Volcano_NoLabel.pdf")
png_file <- file.path(output_dir, "GSE61616_Volcano_NoLabel.png")

# 保存PDF
pdf(pdf_file, width = 8, height = 6)
print(volcano_plot)
dev.off()

# 保存PNG
png(png_file, width = 2400, height = 1800, res = 300)
print(volcano_plot)
dev.off()

cat(paste("  火山图已保存至:", pdf_file, "和", png_file, "\n"))

# 第五步：输出文件清单
cat("Step 5: 输出结果文件...\n")

# 1. 完整映射结果
data_mapped %>%
  select(Mouse_Symbol, Human_Symbol, logFC, adj_P_Val, Regulation) %>%
  write_tsv("GSE61616_all_genes_mapped.tsv")

# 2. 显著上调基因列表
data_mapped %>%
  filter(Regulation == "Up" & !is.na(Human_Symbol)) %>%
  select(Human_Symbol, Mouse_Symbol, logFC, adj_P_Val) %>%
  arrange(desc(logFC)) %>%
  write_tsv("GSE61616_sig_upregulated.tsv")

# 3. 显著下调基因列表
data_mapped %>%
  filter(Regulation == "Down" & !is.na(Human_Symbol)) %>%
  select(Human_Symbol, Mouse_Symbol, logFC, adj_P_Val) %>%
  arrange(logFC) %>%
  write_tsv("GSE61616_sig_downregulated.tsv")

# 4. 显著基因合并列表
data_mapped %>%
  filter(Regulation %in% c("Up", "Down") & !is.na(Human_Symbol)) %>%
  select(Mouse_Symbol, Human_Symbol, logFC, adj_P_Val, Regulation) %>%
  rename(Direction = Regulation) %>%
  arrange(desc(abs(logFC))) %>%
  write_tsv("GSE61616_sig_DEGs_combined.tsv")

# 5. 简易文本格式
# 上调基因
data_mapped %>%
  filter(Regulation == "Up" & !is.na(Human_Symbol)) %>%
  pull(Human_Symbol) %>%
  write_lines("human_up_genes.txt")

# 下调基因
data_mapped %>%
  filter(Regulation == "Down" & !is.na(Human_Symbol)) %>%
  pull(Human_Symbol) %>%
  write_lines("human_down_genes.txt")

# 所有显著基因
data_mapped %>%
  filter(Regulation %in% c("Up", "Down") & !is.na(Human_Symbol)) %>%
  pull(Human_Symbol) %>%
  write_lines("human_all_sig_genes.txt")

# 6. 未映射基因
unmapped_genes <- data_mapped %>%
  filter(is.na(Human_Symbol)) %>%
  pull(Mouse_Symbol)

if (length(unmapped_genes) > 0) {
  write_lines(unmapped_genes, "GSE61616_unmapped_genes.txt")
  cat(paste("  未映射基因数:", length(unmapped_genes), "已保存至 GSE61616_unmapped_genes.txt\n"))
}

cat("\n分析完成！所有结果文件已输出至当前目录。\n")
