# ================================================================================
# GO 富集分析高级柱状图
# 学术发表标准：BP/CC/MF 三维度展示，重点 Term 标记
# ================================================================================

# --------------------------------------------------------------------------------
# 1. 环境配置
# --------------------------------------------------------------------------------

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(bitmapType = "cairo")

# 智能加载包
load_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("⚠  安装包:", pkg))
    if (pkg %in% c("clusterProfiler", "enrichplot", "ggplot2", "dplyr", 
                   "tidyr", "ggrepel", "scales")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
}

message("\n===== 加载必要的包 =====")
packages <- c("ggplot2", "dplyr", "tidyr", "ggrepel", "scales", "stringr")
for (pkg in packages) {
  tryCatch({
    load_package(pkg)
    message(paste("✓ 已加载:", pkg))
  }, error = function(e) {
    message(paste("❌ 加载失败:", pkg))
  })
}

# 设置工作目录
setwd("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/AI 代码编写")
message(paste("\n✓ 工作目录:", getwd()))

# 定义配色方案
colors <- list(
  bp = "#2E86AB",      # 青蓝 - BP
  cc = "#A23B72",      # 紫红 - CC
  mf = "#F18F01",      # 橙黄 - MF
  highlight = "#E63946" # 红色 - 重点标记
)

# 创建输出目录
output_dir <- "GO_Figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------------
# 2. 数据读取
# --------------------------------------------------------------------------------

message("\n===== 读取 GO 富集结果 =====")

# 读取 GO 数据
go_bp <- read.csv("GO_enrichment_results/GO_BP_enrichment_results.csv", 
                  stringsAsFactors = FALSE)
go_cc <- read.csv("GO_enrichment_results/GO_CC_enrichment_results.csv", 
                  stringsAsFactors = FALSE)
go_mf <- read.csv("GO_enrichment_results/GO_MF_enrichment_results.csv", 
                  stringsAsFactors = FALSE)

# 添加类别标签
go_bp$ontology <- "BP"
go_cc$ontology <- "CC"
go_mf$ontology <- "MF"

# 合并数据
go_all <- rbind(go_bp, go_cc, go_mf)
message(paste("✓ 总 GO 条目数:", nrow(go_all)))
message(paste("  - BP:", nrow(go_bp), "条"))
message(paste("  - CC:", nrow(go_cc), "条"))
message(paste("  - MF:", nrow(go_mf), "条"))

# --------------------------------------------------------------------------------
# 3. 数据预处理
# --------------------------------------------------------------------------------

message("\n===== 数据预处理 =====")

# 转换 p.adjust 为数值并计算 -log10
go_all$pvalue_adj <- as.numeric(go_all$p.adjust)
go_all$neg_log10_padj <- -log10(go_all$pvalue_adj)

# 转换 Count 为数值
go_all$count <- as.numeric(go_all$Count)

# 按 ontology 和 p.adjust 排序，选择每类 Top 10
go_top <- go_all %>%
  group_by(ontology) %>%
  arrange(ontology, pvalue_adj) %>%
  slice_head(n = 10) %>%
  ungroup()

# 为每个 term 添加唯一标识
go_top$term_id <- paste0(go_top$ontology, "_", 1:nrow(go_top))

# 标记与重点功能相关的 term
# 关键词：copper ion homeostasis, inflammatory response, apoptosis, etc.
keywords <- c(
  "copper ion", "inflammatory", "apoptosis", "response to metal ion",
  "cell death", "oxidative stress", "cytokine", "immune response"
)

go_top$is_highlight <- FALSE
go_top$highlight_reason <- ""

for (i in 1:nrow(go_top)) {
  desc <- tolower(go_top$Description[i])  # 注意这里是大写的 D
  
  # 检查是否包含关键词
  if (str_detect(desc, "copper ion|copper homeostasis")) {
    go_top$is_highlight[i] <- TRUE
    go_top$highlight_reason[i] <- "Copper Ion Homeostasis"
  } else if (str_detect(desc, "inflammatory")) {
    go_top$is_highlight[i] <- TRUE
    go_top$highlight_reason[i] <- "Inflammatory Response"
  } else if (str_detect(desc, "apoptosis|cell death")) {
    go_top$is_highlight[i] <- TRUE
    go_top$highlight_reason[i] <- "Apoptosis/Cell Death"
  } else if (str_detect(desc, "cytokine|immune")) {
    go_top$is_highlight[i] <- TRUE
    go_top$highlight_reason[i] <- "Immune Response"
  }
}

message(paste("✓ 高亮 Term 数:", sum(go_top$is_highlight)))

# --------------------------------------------------------------------------------
# 4. 创建 GO 柱状图
# --------------------------------------------------------------------------------

message("\n===== 创建 GO 柱状图 =====")

# 设置因子水平（按 ontology 和 p 值排序）
go_top$Description <- factor(go_top$Description, 
                              levels = rev(go_top$Description))

# 创建主图
p_go <- ggplot(go_top, aes(x = neg_log10_padj, y = Description, fill = ontology)) +
  # 绘制柱状图
  geom_bar(stat = "identity", width = 0.7, alpha = 0.9) +
  
  # 添加黑色边框
  geom_bar(stat = "identity", width = 0.7, fill = NA, 
           color = "black", linewidth = 0.5) +
  
  # 在右侧添加基因数量标签
  geom_text(aes(label = count), 
            x = max(go_top$neg_log10_padj) * 1.02,
            hjust = 0, size = 3.5, fontface = "bold", color = "gray30") +
  
  # 为高亮 term 添加红色星号标记
  geom_text(data = subset(go_top, is_highlight),
            aes(label = "*"), 
            x = max(go_top$neg_log10_padj) * 0.98,
            hjust = 1, size = 6, fontface = "bold", color = colors$highlight) +
  
  # 设置颜色
  scale_fill_manual(values = c(BP = colors$bp, CC = colors$cc, MF = colors$mf),
                    name = "Ontology",
                    labels = c("BP (Biological Process)", 
                               "CC (Cellular Component)", 
                               "MF (Molecular Function)")) +
  
  # 设置坐标轴
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.15)),
                     breaks = seq(0, ceiling(max(go_top$neg_log10_padj)), by = 5)) +
  
  # 标签和标题
  labs(
    title = "GO Enrichment Analysis",
    subtitle = "Top 10 enriched terms per category (sorted by adjusted p-value)",
    caption = "* Highlighted terms: Copper ion homeostasis | Inflammatory response",
    x = expression(-log[10]~"(adjusted p-value)"),
    y = "GO Term"
  ) +
  
  # 主题设置
  theme_bw(base_size = 11, base_family = "sans") +
  theme(
    # 背景设置
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    
    # 轴线
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    
    # 字体
    axis.text = element_text(size = 10, color = "black"),
    axis.text.y = element_text(face = "italic", size = 9.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    
    # 图例
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 9.5),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.spacing.y = unit(0.3, "cm"),
    
    # 标题
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                              margin = margin(b = 8)),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5, 
                                 color = "gray40", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = colors$highlight, 
                                face = "italic", hjust = 0.5),
    
    # 边距
    plot.margin = margin(1.5, 1.5, 1, 1, "cm")
  )

message("✓ GO 柱状图创建完成")

# --------------------------------------------------------------------------------
# 5. 输出图形
# --------------------------------------------------------------------------------

message("\n===== 保存图形 =====")

# 保存为 PDF（矢量格式，用于发表）
pdf_file <- file.path(output_dir, "GO_Barplot_Publication.pdf")
ggsave(pdf_file, p_go, width = 10, height = 8, device = cairo_pdf)
message(paste("✓ PDF 已保存:", pdf_file))

# 保存为 PNG（位图格式，用于预览）
png_file <- file.path(output_dir, "GO_Barplot_Preview.png")
ggsave(png_file, p_go, width = 10, height = 8, dpi = 300)
message(paste("✓ PNG 已保存:", png_file))

# 保存高亮 term 的详细信息
highlight_info <- go_top %>%
  filter(is_highlight) %>%
  select(ontology, Description, Count, pvalue_adj, neg_log10_padj, highlight_reason)

if (nrow(highlight_info) > 0) {
  highlight_file <- file.path(output_dir, "Highlighted_Terms.csv")
  write.csv(highlight_info, highlight_file, row.names = FALSE)
  message(paste("✓ 高亮 Term 信息已保存:", highlight_file))
}

# --------------------------------------------------------------------------------
# 6. 生成摘要报告
# --------------------------------------------------------------------------------

message("\n===== 生成摘要报告 =====")

summary_text <- c(
  "================================================================================",
  "GO 富集分析柱状图 - 摘要报告",
  "================================================================================",
  "",
  paste("生成时间:", Sys.time()),
  "",
  "【数据概览】",
  paste("  - 总 GO 条目数:", nrow(go_all)),
  paste("  - BP 条目数:", nrow(go_bp)),
  paste("  - CC 条目数:", nrow(go_cc)),
  paste("  - MF 条目数:", nrow(go_mf)),
  paste("  - 展示条目数:", nrow(go_top), "(每类 Top 10)"),
  "",
  "【高亮 Term 统计】",
  paste("  - 高亮 Term 总数:", sum(go_top$is_highlight)),
  ""
)

# 添加高亮 Term 详情
if (nrow(highlight_info) > 0) {
  summary_text <- c(summary_text, "【高亮 Term 列表】")
  for (i in 1:nrow(highlight_info)) {
    summary_text <- c(summary_text,
      paste0("  ", i, ". ", highlight_info$description[i], 
             " [", highlight_info$ontology[i], "] ",
             "(Count=", highlight_info$count[i], 
             ", p.adj=", format(highlight_info$pvalue_adj[i], digits = 3), ")"))
  }
}

summary_text <- c(summary_text,
  "",
  "【输出文件】",
  paste("  - PDF (矢量):", pdf_file),
  paste("  - PNG (300dpi):", png_file),
  paste("  - 高亮信息:", highlight_file),
  "",
  "【配色方案】",
  paste("  - BP:", colors$bp, "青蓝色"),
  paste("  - CC:", colors$cc, "紫红色"),
  paste("  - MF:", colors$mf, "橙黄色"),
  paste("  - 高亮标记:", colors$highlight, "红色星号"),
  "",
  "================================================================================"
)

# 保存摘要
summary_file <- file.path(output_dir, "GO_Barplot_Summary.txt")
writeLines(summary_text, summary_file)
message(paste("✓ 摘要报告已保存:", summary_file))

# 显示摘要
cat(paste(summary_text, collapse = "\n"))

message("\n\n✅ GO 富集柱状图完成！")
