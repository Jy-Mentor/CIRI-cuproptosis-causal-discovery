# 火山图绘制脚本 - 官方代码版本 (lfc=0)
# 输入验证
stopifnot(file.exists("C:/Users/Jy-Mentor-7/Downloads/GSE97537.top.table (1).tsv"))

# 加载包 - 智能安装
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org/")
}
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cloud.r-project.org/")
}
library(ggplot2)
library(dplyr)

# 读取数据
data <- read.table("C:/Users/Jy-Mentor-7/Downloads/GSE97537.top.table (1).tsv", 
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 验证必要列存在
stopifnot("logFC" %in% colnames(data))
stopifnot("P.Value" %in% colnames(data))
stopifnot("adj.P.Val" %in% colnames(data))
stopifnot(nrow(data) > 0)

# 计算 -log10(Pvalue)
data$negLog10Pval <- -log10(data$P.Value)

# 官方代码: lfc=0, 只看padj < 0.05，不看logFC阈值
padj_threshold <- 0.05

# 分类基因 - 官方decideTests逻辑: lfc=0
# 只要padj < 0.05就着色，不管logFC大小
data <- data %>%
  mutate(
    group = case_when(
      adj.P.Val >= padj_threshold ~ "NotSig",
      logFC > 0 & adj.P.Val < padj_threshold ~ "Up",
      logFC < 0 & adj.P.Val < padj_threshold ~ "Down",
      TRUE ~ "NotSig"
    )
  )

# 设置颜色 - 匹配官方图样式
colors <- c("Down" = "#00BFFF",    # 亮青色
            "NotSig" = "#666666",   # 深灰色
            "Up" = "#FF3333")       # 鲜红色

# 创建PDF
pdf("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/volcano_plot_with_legend.pdf", width = 10, height = 10)

# 绘制火山图
p <- ggplot(data, aes(x = logFC, y = negLog10Pval, color = group)) +
  geom_point(size = 1.2, alpha = 0.9) +
  scale_color_manual(values = colors,
                     breaks = c("Down", "NotSig", "Up"),
                     labels = c("Down", "NotSig", "Up")) +
  
  # 设置标题和轴标签 - 匹配官方图格式
  labs(title = "Volcano plot\nGSE97537: Expression data from MCAO or\nSham operated rat brain [mRNA]\nSHAM vs MCAO, Padj<0.05",
       x = "log2(fold change)",
       y = "-log10(Pvalue)") +
  
  # 主题设置 - 模仿官方图样式
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.box = "horizontal",
    legend.margin = margin(b = -5),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(3, "mm")
  ) +
  
  # 图例设置
  guides(color = guide_legend(override.aes = list(size = 4),
                              nrow = 1,
                              byrow = TRUE)) +
  
  # 设置坐标轴范围
  scale_x_continuous(limits = c(-5, 1.5), breaks = seq(-5, 1, by = 1)) +
  scale_y_continuous(limits = c(0, 8.5), breaks = seq(0, 8, by = 1))

print(p)
dev.off()

cat("火山图已保存至: C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/volcano_plot_with_legend.pdf\n")

# 输出统计信息
cat("\n统计信息 (lfc=0, padj<0.05):\n")
print(table(data$group))
