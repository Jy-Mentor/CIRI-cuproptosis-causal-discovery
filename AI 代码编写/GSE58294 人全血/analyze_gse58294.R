# 分析GSE58294人全血数据集
# 包含数据标准化和差异表达分析

# 设置CRAN镜像
options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

# 安装和加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("pheatmap")) install.packages("pheatmap")

library(tidyverse)
library(ggplot2)
library(pheatmap)

# 设置工作目录
setwd("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写\\GSE58294 人全血")

# 1. 读取数据
read_data <- function(file_path) {
  data <- read_tsv(file_path)
  # 处理可能的列名问题
  colnames(data) <- make.names(colnames(data))
  return(data)
}

# 读取所有时间点的数据
time_points <- c("3H", "5H", "24H")
data_files <- paste0("GSE58294  ", time_points, ".tsv")
data_list <- list()

for (i in 1:length(data_files)) {
  data_list[[time_points[i]]] <- read_data(data_files[i])
  cat("读取", data_files[i], "完成，行数：", nrow(data_list[[time_points[i]]]), "\n")
}

# 2. 数据标准化
# 对logFC进行标准化处理
normalize_data <- function(data) {
  # 使用z-score标准化logFC
  if ("logFC" %in% colnames(data)) {
    data$logFC_norm <- scale(data$logFC)[, 1]
  }
  return(data)
}

# 标准化所有时间点的数据
normalized_data_list <- lapply(data_list, normalize_data)

# 3. 差异表达分析
# 定义差异表达基因的阈值
padj_threshold <- 0.05
logfc_threshold <- 0.58  # 对应fold change > 1.5

# 识别差异表达基因
identify_degs <- function(data) {
  data %>%
    mutate(DEG = ifelse(adj.P.Val < padj_threshold & abs(logFC) > logfc_threshold, "Yes", "No"))
}

deg_data_list <- lapply(normalized_data_list, identify_degs)

# 4. 分析铜死亡基因和NF-κB通路基因
# 铜死亡基因列表（人类）
copper_death_genes <- c(
  "ATP7A", "ATP7B", "BACS1", "DLD", "DLST", "FDX1", "GLS", "LIAS", 
  "LIPT1", "LIPT2", "MARS1", "MTF1", "PDHA1", "PDHB", "SLC25A3", 
  "SLC25A4", "SLC25A5", "SLC31A1", "TFRC"
)

# NF-κB通路基因列表（人类）
nfkb_pathway_genes <- c(
  "NFKB1", "NFKB2", "REL", "REL1", "RELA", "RELB", "IKBKA", "IKBKB", 
  "IKBKG", "CHUK", "IKBKAP", "NFKBIA", "NFKBIB", "NFKBIE", "NFKBID", 
  "NFKBIZ", "BCL3", "TRAF2", "TRAF5", "TRAF6", "RIPK1", "IRAK1", 
  "IRAK4", "TAK1", "MAP3K7", "TGFB1", "TNF", "IL1A", "IL1B"
)

# 提取特定基因集的表达数据
extract_gene_set <- function(data, gene_list, gene_set_name) {
  if ("Gene.symbol" %in% colnames(data)) {
    data %>%
      filter(Gene.symbol %in% gene_list)
  } else {
    # 如果没有Gene.symbol列，返回空数据框
    data.frame()
  }
}

# 分析铜死亡基因
copper_death_data_list <- lapply(deg_data_list, extract_gene_set, 
                                 gene_list = copper_death_genes, 
                                 gene_set_name = "copper_death")

# 分析NF-κB通路基因
nfkb_data_list <- lapply(deg_data_list, extract_gene_set, 
                         gene_list = nfkb_pathway_genes, 
                         gene_set_name = "nfkb_pathway")

# 5. 生成分析结果
# 汇总差异表达基因统计
summarize_degs <- function(data, time_point) {
  total_genes <- nrow(data)
  upregulated <- sum(data$DEG == "Yes" & data$logFC > 0, na.rm = TRUE)
  downregulated <- sum(data$DEG == "Yes" & data$logFC < 0, na.rm = TRUE)
  
  data.frame(
    Time_point = time_point,
    Total_genes = total_genes,
    Upregulated = upregulated,
    Downregulated = downregulated,
    Total_DEGs = upregulated + downregulated
  )
}

deg_summary <- bind_rows(
  lapply(1:length(deg_data_list), function(i) {
    summarize_degs(deg_data_list[[i]], time_points[i])
  })
)

print("差异表达基因统计：")
print(deg_summary)

# 6. 可视化
# 创建输出目录
output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 6.1 差异表达基因数量条形图
png(file.path(output_dir, "deg_counts.png"), width = 800, height = 600)
ggplot(deg_summary, aes(x = Time_point, y = Total_DEGs, fill = Time_point)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Total_DEGs), vjust = -0.3) +
  labs(title = "差异表达基因数量（GSE58294）",
       x = "时间点",
       y = "差异表达基因数量") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# 6.2 铜死亡基因表达热图
# 准备热图数据
prepare_heatmap_data <- function(gene_data_list, gene_set_name) {
  # 合并所有时间点的数据
  combined_data <- bind_rows(gene_data_list, .id = "Time_point")
  
  if (nrow(combined_data) == 0) {
    return(NULL)
  }
  
  # 重塑数据为矩阵格式，处理重复值
  heatmap_data <- combined_data %>%
    select(Gene.symbol, Time_point, logFC_norm) %>%
    # 对每个基因和时间点的重复值取平均值
    group_by(Gene.symbol, Time_point) %>%
    summarise(logFC_norm = mean(logFC_norm, na.rm = TRUE)) %>%
    ungroup()
  
  return(heatmap_data)
}

# 生成铜死亡基因热图
copper_heatmap_data <- prepare_heatmap_data(copper_death_data_list, "copper_death")
if (!is.null(copper_heatmap_data)) {
  # 转换为矩阵
  copper_matrix <- copper_heatmap_data %>%
    spread(key = Time_point, value = logFC_norm)
  
  rownames(copper_matrix) <- copper_matrix$Gene.symbol
  copper_matrix <- copper_matrix[, -1]
  
  # 确保矩阵不为空
  if (nrow(copper_matrix) > 0 && ncol(copper_matrix) > 0) {
    png(file.path(output_dir, "copper_death_heatmap.png"), width = 800, height = 600)
    pheatmap(copper_matrix, 
             main = "铜死亡基因表达热图（GSE58294）",
             cluster_rows = TRUE, 
             cluster_cols = TRUE, 
             scale = "none",
             color = colorRampPalette(c("blue", "white", "red"))(100))
    dev.off()
  }
}

# 生成NF-κB通路基因热图
nfkb_heatmap_data <- prepare_heatmap_data(nfkb_data_list, "nfkb_pathway")
if (!is.null(nfkb_heatmap_data)) {
  # 转换为矩阵
  nfkb_matrix <- nfkb_heatmap_data %>%
    spread(key = Time_point, value = logFC_norm)
  
  rownames(nfkb_matrix) <- nfkb_matrix$Gene.symbol
  nfkb_matrix <- nfkb_matrix[, -1]
  
  # 确保矩阵不为空
  if (nrow(nfkb_matrix) > 0 && ncol(nfkb_matrix) > 0) {
    png(file.path(output_dir, "nfkb_pathway_heatmap.png"), width = 800, height = 600)
    pheatmap(nfkb_matrix, 
             main = "NF-κB通路基因表达热图（GSE58294）",
             cluster_rows = TRUE, 
             cluster_cols = TRUE, 
             scale = "none",
             color = colorRampPalette(c("blue", "white", "red"))(100))
    dev.off()
  }
}

# 7. 时间序列分析
# 合并所有时间点的数据进行时间序列分析
combined_all_data <- bind_rows(
  lapply(1:length(deg_data_list), function(i) {
    deg_data_list[[i]] %>% mutate(Time_point = time_points[i])
  })
)

# 分析铜死亡基因的时间序列表达
if ("Gene.symbol" %in% colnames(combined_all_data)) {
  copper_time_series <- combined_all_data %>%
    filter(Gene.symbol %in% copper_death_genes)
  
  if (nrow(copper_time_series) > 0) {
    png(file.path(output_dir, "copper_death_time_series.png"), width = 1000, height = 600)
    ggplot(copper_time_series, aes(x = Time_point, y = logFC_norm, group = Gene.symbol, color = Gene.symbol)) +
      geom_line() +
      geom_point() +
      labs(title = "铜死亡基因时间序列表达（GSE58294）",
           x = "时间点",
           y = "标准化logFC") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right") +
      scale_x_discrete(limits = time_points)
    dev.off()
  }
  
  # 分析NF-κB通路基因的时间序列表达
  nfkb_time_series <- combined_all_data %>%
    filter(Gene.symbol %in% nfkb_pathway_genes)
  
  if (nrow(nfkb_time_series) > 0) {
    png(file.path(output_dir, "nfkb_pathway_time_series.png"), width = 1000, height = 600)
    ggplot(nfkb_time_series, aes(x = Time_point, y = logFC_norm, group = Gene.symbol, color = Gene.symbol)) +
      geom_line() +
      geom_point() +
      labs(title = "NF-κB通路基因时间序列表达（GSE58294）",
           x = "时间点",
           y = "标准化logFC") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right") +
      scale_x_discrete(limits = time_points)
    dev.off()
  }
}

# 8. 相关性分析
# 分析铜死亡基因和NF-κB通路基因之间的相关性
if ("Gene.symbol" %in% colnames(combined_all_data)) {
  # 提取铜死亡基因和NF-κB通路基因的数据
  copper_nfkb_data <- combined_all_data %>%
    filter(Gene.symbol %in% c(copper_death_genes, nfkb_pathway_genes)) %>%
    select(Gene.symbol, Time_point, logFC_norm)
  
  if (nrow(copper_nfkb_data) > 0) {
    # 处理重复值
    copper_nfkb_data_dedup <- copper_nfkb_data %>%
      group_by(Gene.symbol, Time_point) %>%
      summarise(logFC_norm = mean(logFC_norm, na.rm = TRUE)) %>%
      ungroup()
    
    # 转换为宽格式
    copper_nfkb_wide <- copper_nfkb_data_dedup %>%
      spread(key = Time_point, value = logFC_norm)
    
    # 计算相关性矩阵
    if (nrow(copper_nfkb_wide) > 1) {
      correlation_matrix <- cor(copper_nfkb_wide[, -1], use = "pairwise.complete.obs")
      
      png(file.path(output_dir, "copper_nfkb_correlation.png"), width = 800, height = 600)
      pheatmap(correlation_matrix, 
               main = "铜死亡基因和NF-κB通路基因相关性矩阵（GSE58294）",
               cluster_rows = TRUE, 
               cluster_cols = TRUE, 
               scale = "none",
               color = colorRampPalette(c("blue", "white", "red"))(100))
      dev.off()
    }
  }
}

# 9. 保存分析结果
saveRDS(list(
  deg_summary = deg_summary,
  normalized_data_list = normalized_data_list,
  deg_data_list = deg_data_list,
  copper_death_data_list = copper_death_data_list,
  nfkb_data_list = nfkb_data_list
), file.path(output_dir, "analysis_results.rds"))

print("分析完成！结果保存在output目录中。")
