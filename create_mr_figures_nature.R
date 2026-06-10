#!/usr/bin/env Rscript
# ================================================================================
# MR 分析专业图表生成 - Nature Communications 标准
# 参考权威论文和 GitHub 最佳实践
# ================================================================================

# 安装必要的包
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
  }
}

install_if_missing(c("ggplot2", "ggpubr", "cowplot", "dplyr", "tidyr", "scales", "grid", "gridExtra"))

library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(tidyr)
library(scales)
library(grid)
library(gridExtra)

# ================================================================================
# Nature 期刊风格主题
# ================================================================================

nature_theme <- function() {
  theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.key = element_rect(fill = "white", color = NA)
    )
}

# Nature 配色方案
nature_colors <- c(
  "#0072B2",  # 蓝色
  "#E69F00",  # 橙色
  "#009E73",  # 绿色
  "#CC79A7",  # 粉色
  "#D55E00",  # 红色
  "#56B4E9",  # 浅蓝色
  "#F0E442",  # 黄色
  "#999999"   # 灰色
)

# ================================================================================
# 数据加载
# ================================================================================

load_data <- function() {
  mr_file <- "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/MR_results_main_optimized.csv"
  enrich_file <- "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/functional_enrichment/Reactome_results.csv"
  
  mr_data <- NULL
  enrich_data <- NULL
  
  if (file.exists(mr_file)) {
    mr_data <- read.csv(mr_file, stringsAsFactors = FALSE)
    cat("✓ 加载 MR 结果：", nrow(mr_data), "个基因\n")
  } else {
    cat("✗ MR 结果文件不存在\n")
  }
  
  if (file.exists(enrich_file)) {
    enrich_data <- read.csv(enrich_file, stringsAsFactors = FALSE)
    cat("✓ 加载富集结果：", nrow(enrich_data), "个通路\n")
  }
  
  return(list(mr = mr_data, enrich = enrich_data))
}

# ================================================================================
# Figure 1: 森林图 (Forest Plot)
# ================================================================================

create_forest_plot <- function(mr_data, output_dir) {
  cat("\n创建 Figure 1: 森林图...\n")
  
  # 筛选显著基因和部分非显著基因
  sig_genes <- mr_data[mr_data$fdr_sig == TRUE, ]
  
  if (nrow(sig_genes) == 0) {
    cat("  ✗ 无显著基因，跳过\n")
    return()
  }
  
  # 解析置信区间
  parse_ci <- function(ci_str) {
    if (is.na(ci_str) || ci_str == "") return(c(NA, NA))
    ci_str <- gsub("[()]", "", ci_str)
    parts <- strsplit(ci_str, "-")[[1]]
    if (length(parts) == 2) {
      return(c(as.numeric(parts[1]), as.numeric(parts[2])))
    }
    return(c(NA, NA))
  }
  
  # 准备绘图数据
  plot_data <- sig_genes %>%
    mutate(
      or = as.numeric(OR_95CI),
      ci_low = sapply(OR_95CI, function(x) parse_ci(x)[1]),
      ci_high = sapply(OR_95CI, function(x) parse_ci(x)[2]),
      gene = factor(gene, levels = rev(gene))
    ) %>%
    filter(!is.na(or))
  
  if (nrow(plot_data) == 0) {
    cat("  ✗ 无有效数据，跳过\n")
    return()
  }
  
  # 创建森林图
  p <- ggplot(plot_data, aes(x = or, y = gene)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), 
                   height = 0.3, color = nature_colors[1], linewidth = 0.8) +
    geom_point(shape = 21, fill = nature_colors[1], color = "white", 
               size = 3, stroke = 1.5) +
    scale_x_log10() +
    labs(
      title = "Mendelian Randomization Results",
      subtitle = "Genetic Associations with Stroke Risk",
      x = "Odds Ratio (95% CI)",
      y = ""
    ) +
    nature_theme() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(face = "bold")
    )
  
  ggsave(file.path(output_dir, "Figure1_Forest_Plot.png"), p, 
         width = 8, height = 4, dpi = 300)
  cat("  ✓ 已保存：Figure1_Forest_Plot.png\n")
}

# ================================================================================
# Figure 2: 火山图 (Volcano Plot)
# ================================================================================

create_volcano_plot <- function(mr_data, output_dir) {
  cat("\n创建 Figure 2: 火山图...\n")
  
  # 计算 -log10(P)
  mr_data <- mr_data %>%
    mutate(
      neg_log10_pval = -log10(as.numeric(discovery_pval)),
      is_significant = ifelse(fdr_sig == TRUE, "FDR Significant", "Non-significant"),
      is_significant = factor(is_significant, levels = c("Non-significant", "FDR Significant"))
    )
  
  # 创建火山图
  p <- ggplot(mr_data, aes(x = as.numeric(discovery_b), y = neg_log10_pval)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey60", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.8) +
    geom_hline(yintercept = -log10(0.001), linetype = "dotted", color = "grey40", linewidth = 0.8) +
    geom_point(aes(fill = is_significant), shape = 21, color = "white", 
               size = 2.5, alpha = 0.7) +
    scale_fill_manual(values = c("Non-significant" = nature_colors[8], 
                                 "FDR Significant" = nature_colors[1])) +
    geom_text(data = subset(mr_data, fdr_sig == TRUE),
              aes(label = gene), vjust = -0.5, hjust = 0.5,
              size = 3, fontface = "bold") +
    labs(
      title = "Volcano Plot of MR Associations",
      x = "Beta (Effect Size)",
      y = "-log₁₀(P-value)",
      fill = ""
    ) +
    nature_theme() +
    theme(legend.position = "top")
  
  ggsave(file.path(output_dir, "Figure2_Volcano_Plot.png"), p, 
         width = 8, height = 6, dpi = 300)
  cat("  ✓ 已保存：Figure2_Volcano_Plot.png\n")
}

# ================================================================================
# Figure 3: 富集分析气泡图 (Dot Plot)
# ================================================================================

create_enrichment_plot <- function(enrich_data, output_dir) {
  cat("\n创建 Figure 3: 功能富集气泡图...\n")
  
  if (is.null(enrich_data) || nrow(enrich_data) == 0) {
    cat("  ✗ 无富集数据，跳过\n")
    return()
  }
  
  # 选择 Top 15 通路
  top_pathways <- enrich_data %>%
    arrange(pvalue) %>%
    head(15) %>%
    mutate(
      neg_log10_pval = -log10(pvalue),
      Description = factor(Description, levels = rev(Description))
    )
  
  # 创建气泡图
  p <- ggplot(top_pathways, aes(x = Count, y = Description)) +
    geom_point(aes(size = neg_log10_pval, fill = neg_log10_pval), 
               shape = 21, color = "white", stroke = 0.5) +
    scale_size(range = c(3, 10)) +
    scale_fill_gradient(low = nature_colors[8], high = nature_colors[1]) +
    labs(
      title = "Reactome Pathway Enrichment Analysis",
      x = "Gene Count",
      y = "",
      size = "-log₁₀(P)",
      fill = "-log₁₀(P)"
    ) +
    nature_theme() +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    )
  
  ggsave(file.path(output_dir, "Figure3_Enrichment_Dot_Plot.png"), p, 
         width = 10, height = 8, dpi = 300)
  cat("  ✓ 已保存：Figure3_Enrichment_Dot_Plot.png\n")
}

# ================================================================================
# Figure 4: 敏感性分析图 (Sensitivity Analysis)
# ================================================================================

create_sensitivity_plots <- function(mr_data, output_dir) {
  cat("\n创建 Figure 4: 敏感性分析图...\n")
  
  # 1. F 统计量分布
  p1 <- ggplot(mr_data, aes(x = as.numeric(F_mean))) +
    geom_histogram(fill = nature_colors[1], color = "white", bins = 20, alpha = 0.7) +
    geom_vline(xintercept = 10, linetype = "dashed", color = nature_colors[5], 
               linewidth = 1, label = "F=10") +
    labs(
      title = "Instrument Strength",
      subtitle = "F-statistic Distribution",
      x = "Mean F-statistic",
      y = "Count"
    ) +
    nature_theme()
  
  # 2. SNP 数量分布
  p2 <- ggplot(mr_data, aes(x = as.numeric(nsnp))) +
    geom_histogram(fill = nature_colors[3], color = "white", bins = 20, alpha = 0.7) +
    labs(
      title = "Number of Instruments",
      subtitle = "SNPs per Gene",
      x = "Number of SNPs",
      y = "Count"
    ) +
    nature_theme()
  
  # 3. 异质性检验
  het_data <- mr_data %>% filter(!is.na(Q_p), Q_p > 0)
  p3 <- ggplot(het_data, aes(x = -log10(as.numeric(Q_p)))) +
    geom_histogram(fill = nature_colors[2], color = "white", bins = 20, alpha = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = nature_colors[5], 
               linewidth = 1) +
    labs(
      title = "Heterogeneity Test",
      subtitle = "Cochran Q",
      x = "-log₁₀(P-value)",
      y = "Count"
    ) +
    nature_theme()
  
  # 4. 多效性检验
  pleio_data <- mr_data %>% filter(!is.na(Egger_intercept_p), Egger_intercept_p > 0)
  p4 <- ggplot(pleio_data, aes(x = -log10(as.numeric(Egger_intercept_p)))) +
    geom_histogram(fill = nature_colors[4], color = "white", bins = 20, alpha = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = nature_colors[5], 
               linewidth = 1) +
    labs(
      title = "Pleiotropy Test",
      subtitle = "MR-Egger Intercept",
      x = "-log₁₀(P-value)",
      y = "Count"
    ) +
    nature_theme()
  
  # 组合图表
  combined <- plot_grid(p1, p2, p3, p4, ncol = 2, 
                        labels = c("A", "B", "C", "D"),
                        label_size = 12, label_fontface = "bold")
  
  title <- ggdraw() + 
    draw_label("Sensitivity Analysis & Quality Control", 
               size = 14, fontface = "bold", hjust = 0.5)
  
  final_plot <- plot_grid(title, combined, ncol = 1, rel_heights = c(0.1, 1))
  
  ggsave(file.path(output_dir, "Figure4_Sensitivity_Analysis.png"), final_plot, 
         width = 12, height = 10, dpi = 300)
  cat("  ✓ 已保存：Figure4_Sensitivity_Analysis.png\n")
}

# ================================================================================
# Figure 5: 研究流程图 (Study Flowchart)
# ================================================================================

create_flowchart <- function(mr_data, output_dir) {
  cat("\n创建 Figure 5: 研究流程图...\n")
  
  # 创建简化流程图
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.95, 
             label = "Study Design: Two-Sample Mendelian Randomization",
             size = 5, fontface = "bold") +
    
    # Step 1
    annotate("rect", xmin = 0.1, xmax = 0.9, ymin = 0.80, ymax = 0.90,
             fill = nature_colors[1], color = "white", alpha = 0.8) +
    annotate("text", x = 0.5, y = 0.85, 
             label = "Step 1: Gene Selection (n = 112 genes)",
             color = "white", size = 4) +
    
    # Step 2
    annotate("rect", xmin = 0.1, xmax = 0.9, ymin = 0.60, ymax = 0.70,
             fill = nature_colors[3], color = "white", alpha = 0.8) +
    annotate("text", x = 0.5, y = 0.65, 
             label = paste("Step 2: MR Analysis (n =", nrow(mr_data), "genes)"),
             color = "white", size = 4) +
    
    # Step 3
    annotate("rect", xmin = 0.1, xmax = 0.9, ymin = 0.40, ymax = 0.50,
             fill = nature_colors[2], color = "white", alpha = 0.8) +
    annotate("text", x = 0.5, y = 0.45, 
             label = "Step 3: Sensitivity Analysis",
             color = "white", size = 4) +
    
    # Step 4
    annotate("rect", xmin = 0.1, xmax = 0.9, ymin = 0.20, ymax = 0.30,
             fill = nature_colors[4], color = "white", alpha = 0.8) +
    sig_count <- sum(mr_data$fdr_sig == TRUE, na.rm = TRUE)
    annotate("text", x = 0.5, y = 0.25, 
             label = paste("Step 4: FDR Correction (n =", sig_count, "significant genes)"),
             color = "white", size = 4) +
    
    # 箭头
    annotate("segment", x = 0.5, y = 0.78, xend = 0.5, yend = 0.72,
             arrow = arrow(type = "closed", length = unit(0.05, "inches")),
             color = "grey40", linewidth = 1) +
    annotate("segment", x = 0.5, y = 0.58, xend = 0.5, yend = 0.52,
             arrow = arrow(type = "closed", length = unit(0.05, "inches")),
             color = "grey40", linewidth = 1) +
    annotate("segment", x = 0.5, y = 0.38, xend = 0.5, yend = 0.32,
             arrow = arrow(type = "closed", length = unit(0.05, "inches")),
             color = "grey40", linewidth = 1) +
    
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
  
  ggsave(file.path(output_dir, "Figure5_Study_Flowchart.png"), p, 
         width = 10, height = 8, dpi = 300)
  cat("  ✓ 已保存：Figure5_Study_Flowchart.png\n")
}

# ================================================================================
# 主函数
# ================================================================================

main <- function() {
  cat("=" , rep("=", 69), "\n", sep = "")
  cat("MR 分析专业图表生成 - Nature Communications 标准\n")
  cat(rep("=", 70), "\n", sep = "")
  
  # 加载数据
  data <- load_data()
  mr_data <- data$mr
  enrich_data <- data$enrich
  
  if (is.null(mr_data)) {
    cat("错误：无法加载 MR 数据\n")
    return()
  }
  
  # 创建输出目录
  output_dir <- "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/figures"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 生成图表
  cat("\n生成专业图表...\n")
  create_forest_plot(mr_data, output_dir)
  create_volcano_plot(mr_data, output_dir)
  create_enrichment_plot(enrich_data, output_dir)
  create_sensitivity_plots(mr_data, output_dir)
  create_flowchart(mr_data, output_dir)
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("所有图表已生成!\n")
  cat(rep("=", 70), "\n", sep = "")
  cat("\n输出目录：", output_dir, "\n")
  cat("\n生成的图表:\n")
  cat("  Figure 1: 森林图 - MR 主要结果\n")
  cat("  Figure 2: 火山图 - 效应量 vs 显著性\n")
  cat("  Figure 3: 气泡图 - 功能富集分析\n")
  cat("  Figure 4: 敏感性分析 - 质量控制\n")
  cat("  Figure 5: 流程图 - 研究设计\n")
  cat("\n所有图表均为 Nature Communications 标准 (300 DPI)\n")
  cat(rep("=", 70), "\n", sep = "")
}

# 运行主函数
main()
