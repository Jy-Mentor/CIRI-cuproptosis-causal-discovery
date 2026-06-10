# 差异表达基因分析和火山图绘制

# 加载必要的包
library(limma)
library(ggplot2)
library(stringr)

# 读取数据
data <- read.table("GSE61616.top.table (1).tsv", header = TRUE, sep = "\t")

# 筛选差异表达基因（绝对对数FC≥0.2，校正后P值<0.05）
degs <- data[data$adj.P.Val < 0.05 & abs(data$logFC) >= 0.2, ]

# 保存DEGs
write.table(degs, "DEGs.tsv", sep = "\t", row.names = FALSE)

# 绘制火山图
volcano_plot <- ggplot(data, aes(x = logFC, y = -log10(adj.P.Val)))
volcano_plot <- volcano_plot +
  geom_point(aes(color = ifelse(adj.P.Val < 0.05 & abs(logFC) >= 0.2, "显著差异", "无显著差异")), 
             alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("显著差异" = "red", "无显著差异" = "gray")) +
  labs(title = "差异表达基因火山图",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "表达差异") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "black")

# 保存火山图
ggsave("volcano_plot.png", volcano_plot, width = 10, height = 8, dpi = 300)

# 统计结果
cat("差异表达基因分析完成！\n")
cat("总基因数:", nrow(data), "\n")
cat("显著差异表达基因数:", nrow(degs), "\n")
cat("其中上调基因数:", sum(degs$logFC > 0), "\n")
cat("其中下调基因数:", sum(degs$logFC < 0), "\n")
