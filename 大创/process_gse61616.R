# 设置种子保证可重复性
set.seed(123)

# 加载必要的包
library(tidyverse)
library(ggrepel)

# 文件路径设置
input_file <- "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创\\GSE61616.top.table (1).tsv"
mapping_file <- "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创\\大鼠 小鼠 人类映射库.txt"
output_dir <- "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创"

# 第一步：数据读取与清洗
cat("\n=== 第一步：数据读取与清洗 ===\n")

# 读取数据
data <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE)

# 原始基因总数
original_count <- nrow(data)
cat(paste("原始基因总数：", original_count, "\n"))

# 检查关键列
required_cols <- c("Gene.symbol", "logFC", "P.Value")
if ("adj.P.Val" %in% colnames(data)) {
  required_cols <- c(required_cols, "adj.P.Val")
}

# 检查列是否存在
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop(paste("缺少必要列：", paste(missing_cols, collapse = ", ")))
}

# 处理NA值和空字符串
data_cleaned <- data %>%
  filter(!is.na(Gene.symbol) & Gene.symbol != "" & 
         !is.na(logFC) & 
         !is.na(P.Value))

# 计算丢弃的基因数
dropped_count <- original_count - nrow(data_cleaned)
cat(paste("成功保留基因数：", nrow(data_cleaned), "\n"))
cat(paste("丢弃基因数：", dropped_count, "\n"))

# 第二步：差异表达统计
cat("\n=== 第二步：差异表达统计 ===\n")

# 确定使用的P值列
p_col <- if ("adj.P.Val" %in% colnames(data_cleaned)) "adj.P.Val" else "P.Value"
cat(paste("使用的P值列：", p_col, "\n"))

# 统计显著差异基因
upregulated <- data_cleaned %>% 
  filter(logFC > 1 & !!sym(p_col) < 0.05)
downregulated <- data_cleaned %>% 
  filter(logFC < -1 & !!sym(p_col) < 0.05)
not_significant <- data_cleaned %>% 
  filter(abs(logFC) <= 1 | !!sym(p_col) >= 0.05)

cat(paste("显著上调基因数：", nrow(upregulated), "\n"))
cat(paste("显著下调基因数：", nrow(downregulated), "\n"))
cat(paste("无显著差异基因数：", nrow(not_significant), "\n"))

# 列出Top 5上下调基因
cat("\nTop 5上调基因（按logFC排序）：\n")
top5_up <- upregulated %>% 
  arrange(desc(logFC)) %>% 
  head(5)
print(top5_up[, c("Gene.symbol", "logFC", p_col)])

cat("\nTop 5下调基因（按logFC排序）：\n")
top5_down <- downregulated %>% 
  arrange(logFC) %>% 
  head(5)
print(top5_down[, c("Gene.symbol", "logFC", p_col)])

# 第三步：小鼠→人类同源基因映射
cat("\n=== 第三步：小鼠→人类同源基因映射 ===\n")

# 读取本地映射库（跳过注释行）
# 首先读取所有行
all_lines <- readLines(mapping_file)
# 找到实际的列名行（不以#开头的第一行）
header_line <- grep("^[^#]", all_lines)[1]
# 读取数据，从header_line开始
mapping_data <- read.delim(text = all_lines[header_line:length(all_lines)], 
                          header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# 检查映射库的列名
cat("映射库列名：", paste(colnames(mapping_data), collapse = ", "), "\n")

# 确定小鼠和人类基因列（根据RGD文件格式）
mouse_col <- "MOUSE_ORTHOLOG_SYMBOL"
human_col <- "HUMAN_ORTHOLOG_SYMBOL"

# 检查列是否存在
if (!mouse_col %in% colnames(mapping_data) || !human_col %in% colnames(mapping_data)) {
  stop("映射库中未找到小鼠或人类基因列")
}

cat(paste("使用的小鼠基因列：", mouse_col, "\n"))
cat(paste("使用的人类基因列：", human_col, "\n"))

# 执行映射
data_mapped <- data_cleaned %>%
  left_join(mapping_data, by = setNames(mouse_col, "Gene.symbol")) %>%
  rename(Mouse_Gene_Symbol = Gene.symbol, Human_Gene_Symbol = !!sym(human_col))

# 统计映射成功和失败的基因数
mapped_count <- sum(!is.na(data_mapped$Human_Gene_Symbol) & data_mapped$Human_Gene_Symbol != "")
unmapped_count <- nrow(data_mapped) - mapped_count

cat(paste("成功映射基因数：", mapped_count, "\n"))
cat(paste("映射失败基因数：", unmapped_count, "\n"))
cat(paste("映射成功率：", round(mapped_count / nrow(data_mapped) * 100, 2), "%\n"))

# 提取未映射的基因
unmapped_genes <- data_mapped %>%
  filter(is.na(Human_Gene_Symbol) | Human_Gene_Symbol == "") %>%
  select(Mouse_Gene_Symbol)

if (nrow(unmapped_genes) > 0) {
  cat("\n未映射的基因（前10个）：\n")
  print(head(unmapped_genes, 10))
}

# 第四步：火山图绘制
cat("\n=== 第四步：火山图绘制 ===\n")

# 准备火山图数据
volcano_data <- data_mapped %>%
  filter(!is.na(Human_Gene_Symbol) & Human_Gene_Symbol != "") %>%
  mutate(threshold = case_when(
    logFC > 1 & !!sym(p_col) < 0.05 ~ "Up",
    logFC < -1 & !!sym(p_col) < 0.05 ~ "Down",
    TRUE ~ "Not Sig"
  ))

# 选择Top 10上调和Top 10下调的基因进行标记
top10_up <- volcano_data %>%
  filter(threshold == "Up") %>%
  arrange(desc(logFC)) %>%
  head(10)
top10_down <- volcano_data %>%
  filter(threshold == "Down") %>%
  arrange(logFC) %>%
  head(10)

label_data <- bind_rows(top10_up, top10_down)

# 绘制火山图
p <- ggplot(volcano_data, aes(x = logFC, y = -log10(!!sym(p_col)), color = threshold)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_text_repel(data = label_data, aes(label = Human_Gene_Symbol), 
                 size = 3, box.padding = 0.3, point.padding = 0.5, 
                 segment.size = 0.5, segment.alpha = 0.5) +
  scale_color_manual(values = c("Up" = "#DC143C", "Down" = "#4169E1", "Not Sig" = "#D3D3D3")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  labs(title = "Human Homologs of DEGs in GSE61616 (Mouse Ischemia Model)",
       x = "log2FoldChange",
       y = "-log10(adjusted P-value)",
       color = "Differential Expression") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top")

# 保存火山图
pdf_path <- file.path(output_dir, "GSE61616_Volcano_HumanHomolog.pdf")
png_path <- file.path(output_dir, "GSE61616_Volcano_HumanHomolog.png")

ggsave(pdf_path, p, width = 8, height = 6, device = "pdf")
ggsave(png_path, p, width = 8, height = 6, dpi = 300, device = "png")

cat(paste("火山图PDF保存路径：", pdf_path, "\n"))
cat(paste("火山图PNG保存路径：", png_path, "\n"))

# 第五步：输出文件
cat("\n=== 第五步：输出文件 ===\n")

# 输出清洗后的数据
cleaned_output <- file.path(output_dir, "GSE61616_cleaned_mouse.tsv")
data_cleaned %>%
  mutate(threshold = case_when(
    logFC > 1 & !!sym(p_col) < 0.05 ~ "Up",
    logFC < -1 & !!sym(p_col) < 0.05 ~ "Down",
    TRUE ~ "Not Sig"
  )) %>%
  write_tsv(cleaned_output)

# 输出映射后的数据
mapped_output <- file.path(output_dir, "GSE61616_human_homologs.tsv")
data_mapped %>%
  mutate(threshold = case_when(
    logFC > 1 & !!sym(p_col) < 0.05 ~ "Up",
    logFC < -1 & !!sym(p_col) < 0.05 ~ "Down",
    TRUE ~ "Not Sig"
  )) %>%
  write_tsv(mapped_output)

# 输出未映射的基因
unmapped_output <- file.path(output_dir, "GSE61616_unmapped_genes.tsv")
write_tsv(unmapped_genes, unmapped_output)

cat(paste("清洗后数据输出路径：", cleaned_output, "\n"))
cat(paste("映射后数据输出路径：", mapped_output, "\n"))
cat(paste("未映射基因输出路径：", unmapped_output, "\n"))

# 生成统计报告
cat("\n=== 统计报告 ===\n")
cat(paste("1. 原始基因总数：", original_count, "\n"))
cat(paste("2. 清洗后基因数：", nrow(data_cleaned), "\n"))
cat(paste("3. 映射成功基因数：", mapped_count, "\n"))
cat(paste("4. 映射失败基因数：", unmapped_count, "\n"))
cat(paste("5. 映射成功率：", round(mapped_count / nrow(data_mapped) * 100, 2), "%\n"))
cat(paste("6. 显著上调基因数：", nrow(upregulated), "\n"))
cat(paste("7. 显著下调基因数：", nrow(downregulated), "\n"))
cat(paste("8. 无显著差异基因数：", nrow(not_significant), "\n"))

# 显示映射后的Top 5上下调人类基因
cat("\n映射后的Top 5上调人类基因：\n")
top5_up_human <- data_mapped %>%
  filter(logFC > 1 & !!sym(p_col) < 0.05 & !is.na(Human_Gene_Symbol) & Human_Gene_Symbol != "") %>%
  arrange(desc(logFC)) %>%
  head(5)
print(top5_up_human[, c("Mouse_Gene_Symbol", "Human_Gene_Symbol", "logFC", p_col)])

cat("\n映射后的Top 5下调人类基因：\n")
top5_down_human <- data_mapped %>%
  filter(logFC < -1 & !!sym(p_col) < 0.05 & !is.na(Human_Gene_Symbol) & Human_Gene_Symbol != "") %>%
  arrange(logFC) %>%
  head(5)
print(top5_down_human[, c("Mouse_Gene_Symbol", "Human_Gene_Symbol", "logFC", p_col)])

cat("\n分析完成！所有文件已输出到指定目录。\n")