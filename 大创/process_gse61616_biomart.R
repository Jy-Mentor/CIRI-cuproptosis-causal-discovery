# 设置种子保证可重复性
set.seed(123)

# 加载必要的包
library(tidyverse)
library(ggrepel)
library(limma)

# 尝试加载biomaRt包
biomart_available <- requireNamespace("biomaRt", quietly = TRUE)
if (!biomart_available) {
  cat("警告：biomaRt包未安装，将使用本地映射库作为备选方案\n")
}

# 文件路径设置
input_path <- "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创\\GSE61616.top.table (1).tsv"
output_dir <- dirname(input_path)

# 确保输出目录存在
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 第一步：数据读取与limma标准化质控
cat("\n=== 第一步：数据读取与limma标准化质控 ===\n")

# 读取数据
data <- read_tsv(input_path)

# 原始基因总数
original_count <- nrow(data)
cat(paste("原始基因总数：", original_count, "\n"))

# 检查关键列
required_cols <- c("Gene.symbol", "logFC", "adj.P.Val")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop(paste("缺少必要列：", paste(missing_cols, collapse = ", ")))
}

# 数据清洗：剔除NA和空值
data_cleaned <- data %>%
  filter(!is.na(Gene.symbol) & Gene.symbol != "" &
         !is.na(logFC) &
         !is.na(adj.P.Val))

# 处理重复探针（保留|logFC|最大的）
data_cleaned <- data_cleaned %>%
  group_by(Gene.symbol) %>%
  arrange(desc(abs(logFC))) %>%
  slice(1) %>%
  ungroup()

# 质控报告
unique_genes <- n_distinct(data_cleaned$Gene.symbol)
missing_rate <- (original_count - nrow(data_cleaned)) / original_count * 100

cat(paste("清洗后基因数：", nrow(data_cleaned), "\n"))
cat(paste("唯一基因symbol数：", unique_genes, "\n"))
cat(paste("缺失值比例：", round(missing_rate, 2), "%\n"))

# 第二步：差异表达基因筛选（limma标准）
cat("\n=== 第二步：差异表达基因筛选（limma标准） ===\n")

# 筛选阈值：|logFC|≥0.2 & adj.P.Val<0.05
upregulated <- data_cleaned %>%
  filter(logFC >= 0.2 & adj.P.Val < 0.05)
downregulated <- data_cleaned %>%
  filter(logFC <= -0.2 & adj.P.Val < 0.05)
not_significant <- data_cleaned %>%
  filter(abs(logFC) < 0.2 | adj.P.Val >= 0.05)

cat(paste("显著上调基因数：", nrow(upregulated), "\n"))
cat(paste("显著下调基因数：", nrow(downregulated), "\n"))
cat(paste("无显著差异基因数：", nrow(not_significant), "\n"))

# 列出Top 5上下调基因
cat("\nTop 5上调基因（按logFC排序）：\n")
top5_up <- upregulated %>%
  arrange(desc(logFC)) %>%
  head(5)
print(top5_up[, c("Gene.symbol", "logFC", "adj.P.Val")])

cat("\nTop 5下调基因（按logFC排序）：\n")
top5_down <- downregulated %>%
  arrange(logFC) %>%
  head(5)
print(top5_down[, c("Gene.symbol", "logFC", "adj.P.Val")])

# 第三步：小鼠→人类同源基因映射
cat("\n=== 第三步：小鼠→人类同源基因映射 ===\n")

if (biomart_available) {
  # 尝试连接biomaRt
  cat("正在连接Ensembl数据库...\n")
  tryCatch({
    # 加载biomaRt包
    library(biomaRt)
    
    # 尝试主服务器
    mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # 获取同源基因映射
    mouse_genes <- data_cleaned$Gene.symbol
    
    # 使用getBM获取同源基因信息
    mapping <- getBM(
      attributes = c(
        "external_gene_name",  # 小鼠基因symbol
        "hsapiens_homolog_associated_gene_name",  # 人类同源基因symbol
        "hsapiens_homolog_orthology_type"  # 同源关系类型
      ),
      filters = "external_gene_name",
      values = mouse_genes,
      mart = mouse_mart
    )
    
    # 处理映射结果
    # 对于一对多映射，优先选择one2one同源关系
    mapping_processed <- mapping %>%
      group_by(external_gene_name) %>%
      arrange(
        # 优先选择one2one同源关系
        desc(hsapiens_homolog_orthology_type == "ortholog_one2one"),
        # 然后选择非空的人类基因
        desc(!is.na(hsapiens_homolog_associated_gene_name) & hsapiens_homolog_associated_gene_name != "")
      ) %>%
      slice(1) %>%
      ungroup()
    
    # 执行映射
    data_mapped <- data_cleaned %>%
      left_join(mapping_processed, by = c("Gene.symbol" = "external_gene_name")) %>%
      rename(
        Mouse_Symbol = Gene.symbol,
        Human_Symbol = hsapiens_homolog_associated_gene_name,
        Orthology_Type = hsapiens_homolog_orthology_type
      ) %>%
      # 添加调控状态
      mutate(Regulation_Status = case_when(
        logFC >= 0.2 & adj.P.Val < 0.05 ~ "Up",
        logFC <= -0.2 & adj.P.Val < 0.05 ~ "Down",
        TRUE ~ "NotSig"
      ))
    
    # 统计映射结果
    mapped_count <- sum(!is.na(data_mapped$Human_Symbol) & data_mapped$Human_Symbol != "")
    unmapped_count <- nrow(data_mapped) - mapped_count
    
    cat(paste("成功映射基因数：", mapped_count, "\n"))
    cat(paste("映射失败基因数：", unmapped_count, "\n"))
    cat(paste("映射成功率：", round(mapped_count / nrow(data_mapped) * 100, 2), "%\n"))
    
    # 提取未映射的基因
    unmapped_genes <- data_mapped %>%
      filter(is.na(Human_Symbol) | Human_Symbol == "") %>%
      select(Mouse_Symbol)
    
    if (nrow(unmapped_genes) > 0) {
      cat("\n未映射的基因（前10个）：\n")
      print(head(unmapped_genes, 10))
    }
    
  }, error = function(e) {
    cat("\n警告：biomaRt连接失败，使用本地映射库作为备选方案\n")
    cat(paste("错误信息：", e$message, "\n"))
    
    # 使用本地映射库
    mapping_file <- "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创\\大鼠 小鼠 人类映射库.txt"
    
    # 读取本地映射库（跳过注释行）
    all_lines <- readLines(mapping_file)
    header_line <- grep("^[^#]", all_lines)[1]
    mapping_data <- read.delim(text = all_lines[header_line:length(all_lines)], 
                              header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    # 执行映射
    data_mapped <- data_cleaned %>%
      left_join(mapping_data, by = c("Gene.symbol" = "MOUSE_ORTHOLOG_SYMBOL")) %>%
      rename(
        Mouse_Symbol = Gene.symbol,
        Human_Symbol = HUMAN_ORTHOLOG_SYMBOL,
        Orthology_Type = "HUMAN_ORTHOLOG_SOURCE"
      ) %>%
      mutate(Regulation_Status = case_when(
        logFC >= 0.2 & adj.P.Val < 0.05 ~ "Up",
        logFC <= -0.2 & adj.P.Val < 0.05 ~ "Down",
        TRUE ~ "NotSig"
      ))
    
    # 统计映射结果
    mapped_count <- sum(!is.na(data_mapped$Human_Symbol) & data_mapped$Human_Symbol != "")
    unmapped_count <- nrow(data_mapped) - mapped_count
    
    cat(paste("成功映射基因数：", mapped_count, "\n"))
    cat(paste("映射失败基因数：", unmapped_count, "\n"))
    cat(paste("映射成功率：", round(mapped_count / nrow(data_mapped) * 100, 2), "%\n"))
    
    # 提取未映射的基因
    unmapped_genes <- data_mapped %>%
      filter(is.na(Human_Symbol) | Human_Symbol == "") %>%
      select(Mouse_Symbol)
    
    if (nrow(unmapped_genes) > 0) {
      cat("\n未映射的基因（前10个）：\n")
      print(head(unmapped_genes, 10))
    }
  })
} else {
  # 使用本地映射库
  cat("使用本地映射库进行基因映射...\n")
  mapping_file <- "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创\\大鼠 小鼠 人类映射库.txt"
  
  # 读取本地映射库（跳过注释行）
  all_lines <- readLines(mapping_file)
  header_line <- grep("^[^#]", all_lines)[1]
  mapping_data <- read.delim(text = all_lines[header_line:length(all_lines)], 
                            header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # 执行映射
  data_mapped <- data_cleaned %>%
    left_join(mapping_data, by = c("Gene.symbol" = "MOUSE_ORTHOLOG_SYMBOL")) %>%
    rename(
      Mouse_Symbol = Gene.symbol,
      Human_Symbol = HUMAN_ORTHOLOG_SYMBOL,
      Orthology_Type = "HUMAN_ORTHOLOG_SOURCE"
    ) %>%
    mutate(Regulation_Status = case_when(
      logFC >= 0.2 & adj.P.Val < 0.05 ~ "Up",
      logFC <= -0.2 & adj.P.Val < 0.05 ~ "Down",
      TRUE ~ "NotSig"
    ))
  
  # 统计映射结果
  mapped_count <- sum(!is.na(data_mapped$Human_Symbol) & data_mapped$Human_Symbol != "")
  unmapped_count <- nrow(data_mapped) - mapped_count
  
  cat(paste("成功映射基因数：", mapped_count, "\n"))
  cat(paste("映射失败基因数：", unmapped_count, "\n"))
  cat(paste("映射成功率：", round(mapped_count / nrow(data_mapped) * 100, 2), "%\n"))
  
  # 提取未映射的基因
  unmapped_genes <- data_mapped %>%
    filter(is.na(Human_Symbol) | Human_Symbol == "") %>%
    select(Mouse_Symbol)
  
  if (nrow(unmapped_genes) > 0) {
    cat("\n未映射的基因（前10个）：\n")
    print(head(unmapped_genes, 10))
  }
}

# 第四步：火山图绘制
cat("\n=== 第四步：火山图绘制 ===\n")

# 准备火山图数据
volcano_data <- data_mapped %>%
  filter(!is.na(Human_Symbol) & Human_Symbol != "") %>%
  mutate(
    neg_log10_pval = -log10(adj.P.Val),
    threshold = case_when(
      logFC >= 0.2 & adj.P.Val < 0.05 ~ "Up",
      logFC <= -0.2 & adj.P.Val < 0.05 ~ "Down",
      TRUE ~ "NotSig"
    )
  )

# 选择需要标记的基因（|logFC| > 0.5且adj.P.Val < 0.01的Top 15）
top_genes <- volcano_data %>%
  filter(abs(logFC) > 0.5 & adj.P.Val < 0.01) %>%
  arrange(desc(abs(logFC)))

# 分别取上调和下调的Top 7-8个，总共不超过15个
if (nrow(top_genes) > 15) {
  top_up <- top_genes %>% filter(logFC > 0) %>% head(8)
  top_down <- top_genes %>% filter(logFC < 0) %>% head(7)
  label_data <- bind_rows(top_up, top_down)
} else {
  label_data <- top_genes
}

# 绘制火山图
p <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_pval, color = threshold)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_text_repel(data = label_data, aes(label = Human_Symbol),
                 size = 3, box.padding = 0.3, point.padding = 0.5,
                 segment.size = 0.5, segment.alpha = 0.5) +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NotSig" = "#E0E0E0")) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Cross-species DEG Analysis: Human Homologs of GSE61616 (limma |logFC|≥0.2, BH-FDR<0.05)",
       subtitle = paste("Total DEGs (Human):", nrow(volcano_data %>% filter(threshold != "NotSig")), "(Up:", nrow(volcano_data %>% filter(threshold == "Up")), ", Down:", nrow(volcano_data %>% filter(threshold == "Down")), ")"),
       x = "log2 Fold Change (Mouse→Human Mapped)",
       y = "-log10 (BH-corrected P-value)",
       color = "Regulation Status") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "top"
  )

# 保存火山图
pdf_path <- file.path(output_dir, "GSE61616_Volcano_HumanMapped.pdf")
png_path <- file.path(output_dir, "GSE61616_Volcano_HumanMapped.png")

ggsave(pdf_path, p, width = 8, height = 6, device = "pdf")
ggsave(png_path, p, width = 8, height = 6, dpi = 300, device = "png")

cat(paste("火山图PDF保存路径：", pdf_path, "\n"))
cat(paste("火山图PNG保存路径：", png_path, "\n"))

# 第五步：结果输出
cat("\n=== 第五步：结果输出 ===\n")

# 输出清洗后小鼠数据
cleaned_output <- file.path(output_dir, "GSE61616_mouse_limma_cleaned.tsv")
data_cleaned %>%
  mutate(Regulation_Status = case_when(
    logFC >= 0.2 & adj.P.Val < 0.05 ~ "Up",
    logFC <= -0.2 & adj.P.Val < 0.05 ~ "Down",
    TRUE ~ "NotSig"
  )) %>%
  write_tsv(cleaned_output)

# 输出人类同源基因列表（仅差异表达基因）
human_deg_output <- file.path(output_dir, "GSE61616_human_homologs_DEGs.tsv")
data_mapped %>%
  filter(
    (!is.na(Human_Symbol) & Human_Symbol != "") &
    (abs(logFC) >= 0.2 & adj.P.Val < 0.05)
  ) %>%
  select(Mouse_Symbol, Human_Symbol, logFC, adj.P.Val, Regulation_Status, Orthology_Type) %>%
  write_tsv(human_deg_output)

# 输出未映射的基因
unmapped_output <- file.path(output_dir, "GSE61616_unmapped_genes_biomart.tsv")
if (exists("unmapped_genes")) {
  write_tsv(unmapped_genes, unmapped_output)
}

cat(paste("清洗后小鼠数据输出路径：", cleaned_output, "\n"))
cat(paste("人类同源基因DEG输出路径：", human_deg_output, "\n"))
if (exists("unmapped_genes")) {
  cat(paste("未映射基因输出路径：", unmapped_output, "\n"))
}

# 生成完整统计摘要
cat("\n=== 完整统计摘要 ===\n")
cat(paste("1. 原始基因总数：", original_count, "\n"))
cat(paste("2. 清洗后基因数：", nrow(data_cleaned), "\n"))
cat(paste("3. 显著上调基因数（小鼠）：", nrow(upregulated), "\n"))
cat(paste("4. 显著下调基因数（小鼠）：", nrow(downregulated), "\n"))
cat(paste("5. 成功映射基因数：", mapped_count, "\n"))
cat(paste("6. 映射失败基因数：", unmapped_count, "\n"))
cat(paste("7. 映射成功率：", round(mapped_count / nrow(data_mapped) * 100, 2), "%\n"))

# 统计人类映射后的差异基因数
human_up <- data_mapped %>%
  filter(!is.na(Human_Symbol) & Human_Symbol != "" & logFC >= 0.2 & adj.P.Val < 0.05)
human_down <- data_mapped %>%
  filter(!is.na(Human_Symbol) & Human_Symbol != "" & logFC <= -0.2 & adj.P.Val < 0.05)

cat(paste("8. 显著上调基因数（人类）：", nrow(human_up), "\n"))
cat(paste("9. 显著下调基因数（人类）：", nrow(human_down), "\n"))
cat(paste("10. 总差异基因数（人类）：", nrow(human_up) + nrow(human_down), "\n"))

# 显示映射后的Top 5上下调人类基因
cat("\n映射后的Top 5上调人类基因：\n")
top5_up_human <- data_mapped %>%
  filter(!is.na(Human_Symbol) & Human_Symbol != "" & logFC >= 0.2 & adj.P.Val < 0.05) %>%
  arrange(desc(logFC)) %>%
  head(5)
print(top5_up_human[, c("Mouse_Symbol", "Human_Symbol", "logFC", "adj.P.Val")])

cat("\n映射后的Top 5下调人类基因：\n")
top5_down_human <- data_mapped %>%
  filter(!is.na(Human_Symbol) & Human_Symbol != "" & logFC <= -0.2 & adj.P.Val < 0.05) %>%
  arrange(logFC) %>%
  head(5)
print(top5_down_human[, c("Mouse_Symbol", "Human_Symbol", "logFC", "adj.P.Val")])

# 输出会话信息
cat("\n=== 会话信息 ===\n")
print(sessionInfo())

cat("\n分析完成！所有文件已输出到指定目录。\n")