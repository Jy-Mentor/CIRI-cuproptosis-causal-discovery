# 分析GSE58294差异表达基因与石竹烯、脑缺血基因的三者交集

# 设置CRAN镜像
options(repos = c(CRAN = "https://mirror.lzu.edu.cn/CRAN/"))

# 安装和加载必要的包
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("pheatmap")) install.packages("pheatmap")

library(tidyverse)
library(ggplot2)
library(pheatmap)

# 设置工作目录
setwd("C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写\\GSE58294 人全血")

# 1. 读取基因列表文件
read_gene_list <- function(file_path) {
  genes <- readLines(file_path, encoding = "UTF-8")
  genes <- genes[genes != ""]  # 移除空行
  return(genes)
}

# 基因列表文件路径
bcp_file <- "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写\\石竹烯 人.txt"
brain_ischemia_file <- "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写\\脑缺血 人.txt"

# 读取基因列表
bcp_genes <- read_gene_list(bcp_file)
brain_ischemia_genes <- read_gene_list(brain_ischemia_file)

cat("读取石竹烯相关基因数：", length(bcp_genes), "\n")
cat("读取脑缺血相关基因数：", length(brain_ischemia_genes), "\n")

# 2. 读取之前的分析结果
if (file.exists("output/analysis_results.rds")) {
  analysis_results <- readRDS("output/analysis_results.rds")
  deg_data_list <- analysis_results$deg_data_list
  cat("成功读取之前的分析结果\n")
} else {
  # 如果没有之前的分析结果，重新读取原始数据
  cat("未找到之前的分析结果，重新读取原始数据\n")
  
  # 读取数据函数
  read_data <- function(file_path) {
    data <- read_tsv(file_path)
    colnames(data) <- make.names(colnames(data))
    return(data)
  }
  
  # 读取所有时间点的数据
  time_points <- c("3H", "5H", "24H")
  data_files <- paste0("GSE58294  ", time_points, ".tsv")
  data_list <- list()
  
  for (i in 1:length(data_files)) {
    data_list[[time_points[i]]] <- read_data(data_files[i])
  }
  
  # 标准化数据
  normalize_data <- function(data) {
    if ("logFC" %in% colnames(data)) {
      data$logFC_norm <- scale(data$logFC)[, 1]
    }
    return(data)
  }
  
  normalized_data_list <- lapply(data_list, normalize_data)
  
  # 识别差异表达基因
  padj_threshold <- 0.05
  logfc_threshold <- 0.58
  
  identify_degs <- function(data) {
    data %>%
      mutate(DEG = ifelse(adj.P.Val < padj_threshold & abs(logFC) > logfc_threshold, "Yes", "No"))
  }
  
  deg_data_list <- lapply(normalized_data_list, identify_degs)
}

# 3. 提取差异表达基因
extract_degs <- function(data) {
  data %>%
    filter(DEG == "Yes")
}

degs_list <- lapply(deg_data_list, extract_degs)

# 4. 计算三者交集
calculate_three_way_intersection <- function(deg_data, bcp_genes, brain_ischemia_genes) {
  if ("Gene.symbol" %in% colnames(deg_data)) {
    # 提取差异表达基因
    deg_genes <- deg_data$Gene.symbol
    
    # 计算三者交集
    three_way_genes <- intersect(intersect(deg_genes, bcp_genes), brain_ischemia_genes)
    
    # 提取交集数据
    three_way_data <- deg_data %>%
      filter(Gene.symbol %in% three_way_genes)
    
    cat("\n三者交集基因数：", length(three_way_genes), "\n")
    if (length(three_way_genes) > 0) {
      cat("交集基因：", paste(three_way_genes, collapse = ", "), "\n")
    }
    
    return(list(
      genes = three_way_genes,
      data = three_way_data
    ))
  } else {
    return(list(
      genes = character(0),
      data = data.frame()
    ))
  }
}

# 分析每个时间点的三者交集
analysis_results <- list()
time_points <- names(degs_list)

for (time_point in time_points) {
  cat("\n=== 分析时间点：", time_point, "===")
  
  # 计算三者交集
  three_way_intersect <- calculate_three_way_intersection(
    degs_list[[time_point]], 
    bcp_genes, 
    brain_ischemia_genes
  )
  
  # 保存结果
  analysis_results[[time_point]] <- list(
    three_way = three_way_intersect
  )
}

# 5. 生成交集分析的可视化
output_dir <- "output/three_way_intersection"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 5.1 交集基因数量统计
intersection_counts <- data.frame()

for (time_point in time_points) {
  counts <- data.frame(
    Time_point = time_point,
    Three_way = length(analysis_results[[time_point]]$three_way$genes)
  )
  intersection_counts <- rbind(intersection_counts, counts)
}

# 绘制交集基因数量条形图
png(file.path(output_dir, "three_way_intersection_counts.png"), width = 800, height = 600)
ggplot(intersection_counts, aes(x = Time_point, y = Three_way, fill = Time_point)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Three_way), vjust = -0.3) +
  labs(title = "GSE58294差异表达基因与石竹烯、脑缺血三者交集数量",
       x = "时间点",
       y = "交集基因数量") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# 5.2 分析交集基因的表达模式
# 合并所有时间点的数据
combined_data <- bind_rows(
  lapply(1:length(deg_data_list), function(i) {
    deg_data_list[[names(deg_data_list)[i]]] %>% 
      mutate(Time_point = names(deg_data_list)[i])
  })
)

# 分析三者交集基因的表达
if ("Gene.symbol" %in% colnames(combined_data)) {
  three_way_genes_all <- unique(unlist(lapply(analysis_results, function(x) x$three_way$genes)))
  if (length(three_way_genes_all) > 0) {
    three_way_expression <- combined_data %>%
      filter(Gene.symbol %in% three_way_genes_all)
    
    if (nrow(three_way_expression) > 0) {
      # 绘制热图
      three_way_heatmap_data <- three_way_expression %>%
        select(Gene.symbol, Time_point, logFC_norm) %>%
        group_by(Gene.symbol, Time_point) %>%
        summarise(logFC_norm = mean(logFC_norm, na.rm = TRUE)) %>%
        spread(key = Time_point, value = logFC_norm)
      
      if (nrow(three_way_heatmap_data) > 0) {
        rownames(three_way_heatmap_data) <- three_way_heatmap_data$Gene.symbol
        three_way_heatmap_matrix <- three_way_heatmap_data[, -1]
        
        if (ncol(three_way_heatmap_matrix) > 0) {
          png(file.path(output_dir, "three_way_genes_heatmap.png"), width = 800, height = 600)
          pheatmap(three_way_heatmap_matrix, 
                   main = "三者交集基因表达热图（GSE58294）",
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   scale = "none",
                   color = colorRampPalette(c("blue", "white", "red"))(100))
          dev.off()
        }
        
        # 绘制时间序列图
        png(file.path(output_dir, "three_way_genes_time_series.png"), width = 1000, height = 600)
        ggplot(three_way_expression, aes(x = Time_point, y = logFC_norm, group = Gene.symbol, color = Gene.symbol)) +
          geom_line() +
          geom_point() +
          labs(title = "三者交集基因时间序列表达（GSE58294）",
               x = "时间点",
               y = "标准化logFC") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5),
                legend.position = "right") +
          scale_x_discrete(limits = time_points)
        dev.off()
      }
    }
  }
}

# 6. 保存交集分析结果
saveRDS(list(
  analysis_results = analysis_results,
  intersection_counts = intersection_counts,
  bcp_genes = bcp_genes,
  brain_ischemia_genes = brain_ischemia_genes
), file.path(output_dir, "three_way_intersection_results.rds"))

# 7. 生成综合分析报告
report_file <- file.path(output_dir, "three_way_intersection_report.txt")
sink(report_file)

cat("GSE58294差异表达基因与石竹烯、脑缺血三者交集分析报告\n")
cat("====================================================\n\n")

# 报告基因列表信息
cat("1. 基因列表信息\n")
cat("- 石竹烯相关基因数：", length(bcp_genes), "\n")
cat("- 脑缺血相关基因数：", length(brain_ischemia_genes), "\n\n")

# 报告每个时间点的交集情况
cat("2. 各时间点差异表达基因与基因列表的三者交集\n")
for (time_point in time_points) {
  cat(paste0("\n时间点：", time_point, "\n"))
  cat("- 差异表达基因总数：", nrow(degs_list[[time_point]]), "\n")
  cat("- 三者交集基因数：", length(analysis_results[[time_point]]$three_way$genes), "\n")
  if (length(analysis_results[[time_point]]$three_way$genes) > 0) {
    cat("  交集基因：", paste(analysis_results[[time_point]]$three_way$genes, collapse = ", "), "\n")
  }
}

# 报告总体趋势
cat("\n3. 分析总结\n")
cat("- 交集基因数量统计：\n")
print(intersection_counts)

cat("\n- 可视化结果已保存到output/three_way_intersection目录\n")

# 关闭输出
cat("\n分析完成！报告已保存到output/three_way_intersection/three_way_intersection_report.txt\n")
sink()

print("分析完成！三者交集分析结果已保存到output/three_way_intersection目录")
