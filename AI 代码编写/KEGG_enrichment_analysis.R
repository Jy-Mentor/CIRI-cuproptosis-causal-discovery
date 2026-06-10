# ================================================================================
# KEGG 富集分析脚本
# 基于 clusterProfiler 包进行 ORA（Over-Representation Analysis）
# ================================================================================

# --------------------------------------------------------------------------------
# 1. 环境设置和包加载
# --------------------------------------------------------------------------------

# 设置 CRAN 镜像（使用清华镜像源，国内访问更快）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 智能安装必要的 R 包：只安装未安装的包，避免重复安装
required_packages <- c("clusterProfiler", "org.Hs.eg.db", "tidyverse", 
                       "enrichplot", "DOSE", "ggplot2")

# 检查并安装缺失的包
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("正在安装:", pkg))
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE")) {
      # Bioconductor 包
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      # CRAN 包
      install.packages(pkg)
    }
  }
}

# 加载所有需要的包
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)
library(DOSE)
library(ggplot2)

# --------------------------------------------------------------------------------
# 2. 数据读取和预处理
# --------------------------------------------------------------------------------

# 设置工作目录和输出目录
output_dir <- "KEGG_enrichment_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 读取基因列表文件
# 支持两种格式：TXT 文件（每行一个基因）或 CSV/Excel 文件（GeneSymbol 列）
gene_file <- "C:/Users/Jy-Mentor-7/Desktop/大创/80.txt"

# 智能读取文件
if (file.exists(gene_file)) {
  # 检查文件扩展名
  file_ext <- tools::file_ext(gene_file)
  
  if (file_ext == "txt") {
    # TXT 文件：每行一个基因名
    gene_data <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)
    colnames(gene_data) <- "GeneSymbol"
  } else if (file_ext %in% c("csv", "xlsx", "xls")) {
    # CSV 或 Excel 文件
    if (file_ext == "csv") {
      gene_data <- read_csv(gene_file, show_col_types = FALSE)
    } else {
      # 需要 readxl 包
      if (!requireNamespace("readxl", quietly = TRUE)) {
        install.packages("readxl")
      }
      library(readxl)
      gene_data <- read_excel(gene_file) %>% as.data.frame()
    }
    
    # 检查是否有 GeneSymbol 列
    if (!"GeneSymbol" %in% colnames(gene_data)) {
      stop("错误：CSV/Excel 文件必须包含 'GeneSymbol' 列")
    }
  } else {
    stop("错误：不支持的文件格式，请使用 TXT、CSV 或 Excel 文件")
  }
  
  message(paste("成功读取文件:", gene_file))
} else {
  stop(paste("错误：文件不存在 -", gene_file))
}

# 提取基因列表并去重
gene_list <- unique(gene_data$GeneSymbol)
gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]

message(paste("输入基因总数:", length(gene_list)))
message(paste("去重后基因数:", length(gene_list)))

# --------------------------------------------------------------------------------
# 3. 基因 ID 转换（SYMBOL -> ENTREZID）
# --------------------------------------------------------------------------------

message("\n===== 基因 ID 转换 =====")

# 使用 bitr 进行基因 ID 转换
# org.Hs.eg.db 是人类基因组注释数据库
gene_mapping <- tryCatch({
  bitr(
    gene_list,
    fromType = "SYMBOL",      # 输入 ID 类型：基因符号
    toType = "ENTREZID",      # 输出 ID 类型：Entrez ID（KEGG 需要）
    OrgDb = org.Hs.eg.db      # 物种数据库：人类
  )
}, error = function(e) {
  message(paste("基因转换警告:", e$message))
  return(NULL)
})

# 检查转换结果
if (is.null(gene_mapping) || nrow(gene_mapping) == 0) {
  stop("错误：基因 ID 转换失败，请检查基因符号是否正确")
}

# 统计转换成功率
mapped_genes <- unique(gene_mapping$SYMBOL)
unmapped_genes <- setdiff(gene_list, mapped_genes)

message(paste("成功映射基因数:", length(mapped_genes)))
message(paste("未映射基因数:", length(unmapped_genes)))

# 输出未映射基因到文件（便于后续检查）
if (length(unmapped_genes) > 0) {
  unmapped_file <- file.path(output_dir, "unmapped_genes.txt")
  write.table(unmapped_genes, unmapped_file, 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  message(paste("未映射基因已保存到:", unmapped_file))
  message("未映射基因示例:", paste(head(unmapped_genes), collapse = ", "))
}

# 去重：一个基因可能对应多个 Entrez ID，取第一个
gene_entrez <- gene_mapping %>%
  distinct(SYMBOL, .keep_all = TRUE)

# 提取 Entrez ID 向量用于富集分析
entrez_ids <- gene_entrez$ENTREZID

message(paste("用于富集分析的 Entrez ID 数:", length(entrez_ids)))

# --------------------------------------------------------------------------------
# 4. KEGG 富集分析（ORA 方法）
# --------------------------------------------------------------------------------

message("\n===== KEGG 富集分析 =====")

# 执行 KEGG ORA 分析
# enrichKEGG 使用超几何分布检验基因是否在某个通路中显著富集
kegg_result <- tryCatch({
  enrichKEGG(
    gene = entrez_ids,           # 输入基因列表（Entrez ID）
    organism = "hsa",            # 物种代码：hsa=人类，mmu=小鼠，rno=大鼠
    universe = NULL,             # 背景基因集：NULL 表示使用全基因组
    pvalueCutoff = 0.05,         # P 值阈值
    pAdjustMethod = "BH",        # 多重检验校正方法：Benjamini-Hochberg
    qvalueCutoff = 0.2,          # Q 值阈值（FDR 校正后的 P 值）
    minGSSize = 10,              # 通路最小基因数
    maxGSSize = 500,             # 通路最大基因数
    readable = FALSE             # 是否转换为基因符号（后续手动处理）
  )
}, error = function(e) {
  message(paste("KEGG 富集分析警告:", e$message))
  return(NULL)
})

# 检查分析结果
if (is.null(kegg_result)) {
  message("警告：KEGG 富集分析未获得显著结果，尝试放宽阈值...")
  
  # 放宽阈值重新分析
  kegg_result <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pvalueCutoff = 1,            # 不限制 P 值
    qvalueCutoff = 1,            # 不限制 Q 值
    minGSSize = 5,               # 降低最小基因数
    maxGSSize = 500
  )
}

# 获取富集结果数据框
if (!is.null(kegg_result)) {
  kegg_df <- as.data.frame(kegg_result)
  message(paste("显著富集的通路数:", nrow(kegg_df)))
  
  if (nrow(kegg_df) > 0) {
    message("\nTop 5 富集通路:")
    for (i in 1:min(5, nrow(kegg_df))) {
      message(paste0(i, ". ", kegg_df$Description[i], 
                    " (p.adjust = ", round(kegg_df$p.adjust[i], 4), ")"))
    }
  }
} else {
  stop("错误：KEGG 富集分析失败")
}

# --------------------------------------------------------------------------------
# 5. 结果导出
# --------------------------------------------------------------------------------

message("\n===== 导出结果 =====")

# 导出完整的富集结果表格（CSV 格式）
kegg_output_file <- file.path(output_dir, "KEGG_enrichment_results.csv")

# 整理输出列：包含所有要求的字段
kegg_export <- kegg_df %>%
  select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, 
         geneID, Count) %>%
  arrange(pvalue)  # 按 P 值排序

write.csv(kegg_export, kegg_output_file, row.names = FALSE)
message(paste("富集结果已保存:", kegg_output_file))

# 导出基因映射结果
mapping_file <- file.path(output_dir, "gene_mapping.csv")
write.csv(gene_mapping, mapping_file, row.names = FALSE)
message(paste("基因映射结果已保存:", mapping_file))

# --------------------------------------------------------------------------------
# 6. 数据可视化
# --------------------------------------------------------------------------------

message("\n===== 生成可视化图形 =====")

# 设置图形主题（tidyverse 风格）
theme_set(theme_bw(base_size = 12))

# 6.1 Dotplot（Top 20 通路）
# 气泡图：展示基因比例、富集显著性和通路大小
dotplot_file <- file.path(output_dir, "KEGG_dotplot_Top20.pdf")
pdf(dotplot_file, width = 12, height = 10)

if (nrow(kegg_df) > 0) {
  tryCatch({
    # 使用 enrichResult 对象而不是数据框
    p_dotplot <- dotplot(kegg_result, 
                         showCategory = 20,
                         title = "KEGG Pathway Enrichment (Top 20)",
                         font.size = 12) +
      theme(axis.text.y = element_text(size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
    print(p_dotplot)
    message("已生成：KEGG_dotplot_Top20.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Dotplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant KEGG Pathways", cex = 1.5)
}

dev.off()

# 6.2 Barplot（Top 20 通路）
# 柱状图：展示每个通路的基因数量和显著性
barplot_file <- file.path(output_dir, "KEGG_barplot_Top20.pdf")
pdf(barplot_file, width = 12, height = 10)

if (nrow(kegg_df) > 0) {
  tryCatch({
    p_barplot <- barplot(kegg_result, 
                         showCategory = 20,
                         title = "KEGG Pathway Enrichment (Top 20)",
                         font.size = 12,
                         colour = "p.adjust") +  # 按校正 P 值着色
      theme(axis.text.y = element_text(size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
    print(p_barplot)
    message("已生成：KEGG_barplot_Top20.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Barplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant KEGG Pathways", cex = 1.5)
}

dev.off()

# 6.3 Cnetplot（基因 - 通路网络图）
# 展示基因与通路之间的关系网络
cnetplot_file <- file.path(output_dir, "KEGG_cnetplot.pdf")
pdf(cnetplot_file, width = 14, height = 12)

if (nrow(kegg_df) > 0) {
  tryCatch({
    # 显示 Top 10 通路的基因网络
    p_cnet <- cnetplot(kegg_result,
                       showCategory = 10,      # 显示 Top 10 通路
                       foldChange = NULL,      # 无表达量数据
                       circular = FALSE,       # 线性布局
                       colorEdge = TRUE,       # 彩色边缘
                       cex_label_gene = 0.8,   # 基因标签大小
                       cex_label_category = 1) # 通路标签大小
    
    title("Gene-Pathway Network (Top 10 Pathways)", line = 2)
    message("已生成：KEGG_cnetplot.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Cnetplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant KEGG Pathways", cex = 1.5)
}

dev.off()

# 6.4 Emapplot（通路富集图谱）
# 展示通路之间的相似性和关系
emapplot_file <- file.path(output_dir, "KEGG_emapplot.pdf")
pdf(emapplot_file, width = 12, height = 10)

if (nrow(kegg_df) > 0) {
  tryCatch({
    # 基于 Kappa 统计量计算通路相似性
    p_emap <- emapplot(kegg_result,
                       showCategory = 20,      # 显示 Top 20 通路
                       layout = "kk",          # Kamada-Kawai 布局
                       node_label = "all",     # 显示所有节点标签
                       cex = 1.2,              # 节点大小
                       color = "p.adjust")     # 按校正 P 值着色
    
    title("KEGG Pathway Enrichment Map (Top 20)", line = 2)
    message("已生成：KEGG_emapplot.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Emapplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant KEGG Pathways", cex = 1.5)
}

dev.off()

# 6.5 Ridgeplot（脊线图，可选）
# 展示多个通路的基因分布
ridgeplot_file <- file.path(output_dir, "KEGG_ridgeplot.pdf")
pdf(ridgeplot_file, width = 12, height = 10)

if (nrow(kegg_df) > 0) {
  tryCatch({
    p_ridge <- ridgeplot(kegg_result,
                         showCategory = 15) +
      theme(axis.text.y = element_text(size = 10))
    
    print(p_ridge)
    message("已生成：KEGG_ridgeplot.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Ridgeplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant KEGG Pathways", cex = 1.5)
}

dev.off()

# --------------------------------------------------------------------------------
# 7. 生成分析摘要报告
# --------------------------------------------------------------------------------

message("\n===== 分析结果摘要 =====")

summary_lines <- c(
  "================================================================================",
  "KEGG 富集分析摘要报告",
  "================================================================================",
  "",
  paste("分析时间:", Sys.time()),
  paste("输入文件:", gene_file),
  "",
  "【基因统计】",
  paste("  - 输入基因总数:", length(gene_list)),
  paste("  - 成功映射基因数:", length(mapped_genes)),
  paste("  - 未映射基因数:", length(unmapped_genes)),
  paste("  - 映射成功率:", round(length(mapped_genes) / length(gene_list) * 100, 2), "%"),
  "",
  "【KEGG 富集分析参数】",
  "  - 物种：人类 (Homo sapiens)",
  "  - 分析方法：ORA (Over-Representation Analysis)",
  "  - 背景基因集：全基因组",
  "  - P 值校正方法：Benjamini-Hochberg (BH)",
  "  - P 值阈值：0.05",
  "  - Q 值阈值：0.2",
  "  - 通路大小范围：10-500 个基因",
  "",
  "【富集结果】",
  paste("  - 显著富集通路数:", nrow(kegg_df)),
  "",
  "【输出文件】",
  paste("  - 富集结果表格:", kegg_output_file),
  paste("  - 基因映射结果:", mapping_file),
  paste("  - 未映射基因:", file.path(output_dir, "unmapped_genes.txt")),
  paste("  - Dotplot:", dotplot_file),
  paste("  - Barplot:", barplot_file),
  paste("  - Cnetplot:", cnetplot_file),
  paste("  - Emapplot:", emapplot_file),
  paste("  - Ridgeplot:", ridgeplot_file),
  "",
  "================================================================================"
)

# 保存摘要报告
summary_file <- file.path(output_dir, "analysis_summary.txt")
writeLines(summary_lines, summary_file)
message(paste("分析摘要已保存:", summary_file))

# 在控制台显示摘要
cat(paste(summary_lines, collapse = "\n"))

message("\n✓ KEGG 富集分析完成！所有结果已保存到:", normalizePath(output_dir))
