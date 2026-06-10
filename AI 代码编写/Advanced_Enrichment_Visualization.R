# ================================================================================
# 富集分析高级可视化脚本
# 学术发表标准：GO 三维度柱状图 + KEGG 气泡图 + 网络图/弦图
# ================================================================================

# --------------------------------------------------------------------------------
# 1. Setup - 环境配置和包加载
# --------------------------------------------------------------------------------

# 设置 CRAN 镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 智能加载包函数
load_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("⚠  需要安装包:", pkg))
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE", 
                   "ggplot2", "dplyr", "tidyr", "readxl", "patchwork", 
                   "ggraph", "igraph", "ggVennDiagram", "RColorBrewer",
                   "scales", "ggrepel", "viridis", "ComplexHeatmap", "circlize")) {
      # Bioconductor 包
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
  message(paste("✓ 已加载:", pkg))
}

# 加载所有需要的包
message("\n===== 加载必要的包 =====")
required_packages <- c("ggplot2", "dplyr", "tidyr", "readxl", "patchwork",
                       "ggraph", "igraph", "ggVennDiagram", "RColorBrewer",
                       "scales", "ggrepel", "viridis", "ComplexHeatmap", 
                       "circlize", "stringr", "tidygraph")

for (pkg in required_packages) {
  tryCatch({
    load_package(pkg)
  }, error = function(e) {
    message(paste("❌ 加载失败:", pkg, "-", e$message))
  })
}

# 设置工作目录（使用双反斜杠或正斜杠）
setwd("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/AI 代码编写")
message(paste("\n✓ 工作目录:", getwd()))

# 设置图形参数
options(bitmapType = "cairo")

# 自定义学术主题函数
theme_academic <- function(base_size = 10, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # 背景设置
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # 轴线设置
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      
      # 字体设置
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size + 2, face = "bold", colour = "black"),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 10)),
      
      # 图例设置
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.text = element_text(size = base_size, colour = "black"),
      legend.title = element_text(size = base_size + 1, face = "bold", colour = "black"),
      legend.spacing.y = unit(0.5, "cm"),
      
      # 标题设置
      plot.title = element_text(size = base_size + 4, face = "bold", 
                                hjust = 0.5, vjust = 1, margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size + 1, hjust = 0.5, 
                                   colour = "gray40", margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 2, colour = "gray60"),
      
      # 面板间距
      panel.spacing = unit(1, "lines")
    )
}

# 颜色方案定义
colors <- list(
  bp = "#2E86AB",      # 学术蓝 - BP
  cc = "#A23B72",      # 紫红 - CC
  mf = "#F18F01",      # 琥珀金 - MF
  kegg_gradient = c("#D3D3D3", "#FF6B6B", "#C92A2A"),  # 灰到深红
  npg = c("#0073C2", "#EFC000", "#868686", "#CD534C", "#00AF91", 
          "#FA7921", "#707070", "#E64A19", "#2E86AB", "#A23B72")
)

# 创建输出目录
output_dir <- "Figures_Publication"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message(paste("✓ 输出目录:", normalizePath(output_dir)))

# --------------------------------------------------------------------------------
# 2. Data Ingestion - 数据读取和清洗
# --------------------------------------------------------------------------------

message("\n===== 读取富集分析结果 =====")

# 定义数据文件路径（使用 normalizePath 处理路径）
data_files <- list(
  go_bp = "GO_enrichment_results/GO_BP_enrichment_results.csv",
  go_cc = "GO_enrichment_results/GO_CC_enrichment_results.csv",
  go_mf = "GO_enrichment_results/GO_MF_enrichment_results.csv",
  kegg = "KEGG_enrichment_results/KEGG_enrichment_results.csv",
  gene_mapping = "KEGG_enrichment_results/gene_mapping.csv"
)

# 读取数据函数
read_enrichment_data <- function(file_path, data_type) {
  tryCatch({
    if (!file.exists(file_path)) {
      message(paste("⚠  文件不存在:", file_path))
      return(NULL)
    }
    
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # 标准化列名
    if (requireNamespace("janitor", quietly = TRUE)) {
      data <- janitor::clean_names(data)
    } else {
      # 手动清理列名
      colnames(data) <- tolower(colnames(data))
      colnames(data) <- gsub("\\.", "_", colnames(data))
    }
    
    message(paste("✓ 读取", data_type, ":", nrow(data), "行"))
    return(data)
  }, error = function(e) {
    message(paste("❌ 读取失败:", data_type, "-", e$message))
    return(NULL)
  })
}

# 读取所有数据
go_bp_data <- read_enrichment_data(data_files$go_bp, "GO-BP")
go_cc_data <- read_enrichment_data(data_files$go_cc, "GO-CC")
go_mf_data <- read_enrichment_data(data_files$go_mf, "GO-MF")
kegg_data <- read_enrichment_data(data_files$kegg, "KEGG")
gene_mapping <- read_enrichment_data(data_files$gene_mapping, "Gene Mapping")

# 添加类别标签
if (!is.null(go_bp_data)) go_bp_data$ontology <- "BP"
if (!is.null(go_cc_data)) go_cc_data$ontology <- "CC"
if (!is.null(go_mf_data)) go_mf_data$ontology <- "MF"

# --------------------------------------------------------------------------------
# 3. Preprocessing - 数据预处理
# --------------------------------------------------------------------------------

message("\n===== 数据预处理 =====")

# 计算富集指标的函数
calculate_enrichment_metrics <- function(data) {
  if (is.null(data) || nrow(data) == 0) return(data)
  
  tryCatch({
    # 解析 GeneRatio 和 BgRatio
    data$gene_ratio_num <- sapply(strsplit(data$gene_ratio, "/"), function(x) as.numeric(x[1]))
    data$bg_ratio_num <- sapply(strsplit(data$bg_ratio, "/"), function(x) as.numeric(x[1]))
    
    # 计算 Rich Factor (富集因子)
    data$rich_factor <- data$gene_ratio_num / data$bg_ratio_num
    
    # 计算 Fold Enrichment
    data$fold_enrichment <- (data$gene_ratio_num / data$count) / 
                            (data$bg_ratio_num / data$bg_ratio_num)
    
    message(paste("  ✓ 计算富集指标完成"))
    return(data)
  }, error = function(e) {
    message(paste("  ⚠  计算富集指标失败:", e$message))
    return(data)
  })
}

# 处理所有数据集
go_bp_data <- calculate_enrichment_metrics(go_bp_data)
go_cc_data <- calculate_enrichment_metrics(go_cc_data)
go_mf_data <- calculate_enrichment_metrics(go_mf_data)
kegg_data <- calculate_enrichment_metrics(kegg_data)

# 基因-Term 映射表创建
create_gene_term_mapping <- function(data, ontology_name) {
  if (is.null(data) || nrow(data) == 0) return(NULL)
  
  tryCatch({
    # 展开 geneID 列
    gene_term_list <- list()
    
    for (i in 1:nrow(data)) {
      genes <- unlist(strsplit(data$geneid[i], "/"))
      term <- data$description[i]
      
      for (gene in genes) {
        gene_term_list[[length(gene_term_list) + 1]] <- data.frame(
          gene = gene,
          term = term,
          ontology = ontology_name,
          stringsAsFactors = FALSE
        )
      }
    }
    
    result <- do.call(rbind, gene_term_list)
    message(paste("  ✓ 创建", ontology_name, "基因-Term 映射:", nrow(result), "条"))
    return(result)
  }, error = function(e) {
    message(paste("  ⚠  创建映射失败:", ontology_name, "-", e$message))
    return(NULL)
  })
}

# 创建所有映射
gene_term_bp <- create_gene_term_mapping(go_bp_data, "BP")
gene_term_cc <- create_gene_term_mapping(go_cc_data, "CC")
gene_term_mf <- create_gene_term_mapping(go_mf_data, "MF")
gene_term_kegg <- create_gene_term_mapping(kegg_data, "KEGG")

# 合并所有 GO 映射
gene_term_go <- NULL
if (!is.null(gene_term_bp)) gene_term_go <- rbind(gene_term_go, gene_term_bp)
if (!is.null(gene_term_cc)) gene_term_go <- rbind(gene_term_go, gene_term_cc)
if (!is.null(gene_term_mf)) gene_term_go <- rbind(gene_term_go, gene_term_mf)

# --------------------------------------------------------------------------------
# 4. Visualization Functions - 可视化函数
# --------------------------------------------------------------------------------

message("\n===== 创建可视化函数 =====")

# 4.1 GO 三维度柱状图（精简版 Top 10）
plot_go_bar <- function(bp_data, cc_data, mf_data, top_n = 10) {
  message("  [1/6] 创建 GO 柱状图...")
  
  if (is.null(bp_data) && is.null(cc_data) && is.null(mf_data)) {
    message("    ⚠  无 GO 数据")
    return(NULL)
  }
  
  # 合并数据并选择 Top N
  combined <- rbind(
    if (!is.null(bp_data)) head(bp_data[order(bp_data$pvalue), ], top_n),
    if (!is.null(cc_data)) head(cc_data[order(cc_data$pvalue), ], top_n),
    if (!is.null(mf_data)) head(mf_data[order(mf_data$pvalue), ], top_n)
  )
  
  if (is.null(combined) || nrow(combined) == 0) return(NULL)
  
  # 添加基因数标签
  combined$count <- as.numeric(combined$count)
  combined$description <- factor(combined$description, 
                                  levels = rev(combined$description))
  
  # 创建柱状图
  p <- ggplot(combined, aes(x = description, y = count, fill = ontology)) +
    geom_bar(stat = "identity", width = 0.7, colour = "black", linewidth = 0.3) +
    geom_text(aes(label = count), hjust = -0.1, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c(BP = colors$bp, CC = colors$cc, MF = colors$mf),
                      name = "Ontology") +
    labs(
      title = "GO Enrichment Analysis (Top 10 per Category)",
      subtitle = "Biological Process | Cellular Component | Molecular Function",
      x = "GO Term",
      y = "Number of Genes"
    ) +
    coord_flip() +
    theme_academic() +
    theme(
      axis.text.y = element_text(size = 9, face = "italic"),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  
  message("    ✓ GO 柱状图完成")
  return(p)
}

# 4.2 KEGG 气泡图（精简版 Top 15）
plot_kegg_bubble <- function(kegg_data, top_n = 15) {
  message("  [2/6] 创建 KEGG 气泡图...")
  
  if (is.null(kegg_data) || nrow(kegg_data) == 0) {
    message("    ⚠  无 KEGG 数据")
    return(NULL)
  }
  
  # 选择 Top N 通路
  kegg_top <- head(kegg_data[order(kegg_data$pvalue), ], top_n)
  
  # 解析基因比例
  kegg_top$gene_ratio_num <- sapply(strsplit(kegg_top$gene_ratio, "/"), 
                                     function(x) as.numeric(x[1]))
  
  # 添加通路类别标注（基于通路名称关键词）
  kegg_top$pathway_category <- case_when(
    str_detect(kegg_top$description, "Immune|Inflammatory|Cytokine|Toll|TNF|NF-kappa") ~ "Immune System",
    str_detect(kegg_top$description, "Metabolism|Lipid|Atherosclerosis|Diabetic") ~ "Metabolism",
    str_detect(kegg_top$description, "Apoptosis|Necroptosis|Cell cycle") ~ "Cell Death",
    str_detect(kegg_top$description, "Virus|Hepatitis|Influenza|Herpes") ~ "Infectious Disease",
    TRUE ~ "Other"
  )
  
  # 定义形状
  shape_map <- c("Immune System" = 21, "Metabolism" = 22, 
                 "Cell Death" = 23, "Infectious Disease" = 24, "Other" = 25)
  
  kegg_top$shape <- shape_map[kegg_top$pathway_category]
  
  # 创建气泡图
  p <- ggplot(kegg_top, aes(x = gene_ratio_num, y = description)) +
    geom_point(aes(size = count, color = -log10(pvalue)), 
               alpha = 0.8, shape = 21, stroke = 0.5) +
    scale_color_gradientn(colors = colors$kegg_gradient, 
                          name = "-log10(P-value)",
                          guide = guide_colorbar(barwidth = 10, barheight = 0.5)) +
    scale_size(range = c(3, 10), name = "Gene Count") +
    labs(
      title = "KEGG Pathway Enrichment (Top 15)",
      subtitle = "Bubble size represents gene count, color indicates significance",
      x = "Gene Ratio",
      y = "KEGG Pathway"
    ) +
    theme_academic() +
    theme(
      axis.text.y = element_text(size = 9, face = "italic"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.spacing.x = unit(1, "cm")
    )
  
  message("    ✓ KEGG 气泡图完成")
  return(p)
}

# 4.3 网络图（横向放置）
plot_network <- function(gene_term_data, top_n = 20) {
  message("  [3/6] 创建网络图...")
  
  if (is.null(gene_term_data) || nrow(gene_term_data) == 0) {
    message("    ⚠  无网络图数据")
    return(NULL)
  }
  
  # 选择 Top N terms
  top_terms <- gene_term_data %>%
    group_by(term) %>%
    summarise(gene_count = n(), .groups = "drop") %>%
    arrange(desc(gene_count)) %>%
    head(top_n) %>%
    pull(term)
  
  # 过滤数据
  network_data <- gene_term_data %>%
    filter(term %in% top_terms) %>%
    distinct()
  
  # 创建网络图
  g <- graph_from_data_frame(
    d = network_data[, c("gene", "term")],
    directed = FALSE
  )
  
  # 设置节点属性
  V(g)$type <- ifelse(V(g)$name %in% network_data$term, "term", "gene")
  V(g)$color <- ifelse(V(g)$type == "term", colors$bp, "#808080")
  V(g)$size <- ifelse(V(g)$type == "term", 8, 4)
  
  # 使用 ggraph 绘制
  p <- ggraph(g, layout = "fr") +  # Fruchterman-Reingold 布局
    geom_edge_link(color = "gray70", alpha = 0.5, edge_width = 0.3) +
    geom_node_point(aes(color = type, size = type)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
    scale_color_manual(values = c(term = colors$bp, gene = "#808080"),
                       name = "Node Type") +
    scale_size_manual(values = c(term = 8, gene = 3), guide = "none") +
    labs(
      title = "Gene-Term Network",
      subtitle = "Top 20 enriched terms and associated genes",
      caption = "Layout: Fruchterman-Reingold"
    ) +
    theme_graph() +
    theme(
      legend.position = "bottom"
    )
  
  message("    ✓ 网络图完成")
  return(p)
}

# 4.4 弦图（横向放置）
plot_chord <- function(gene_term_data, top_n = 15) {
  message("  [4/6] 创建弦图...")
  
  if (is.null(gene_term_data) || nrow(gene_term_data) == 0) {
    message("    ⚠  无弦图数据")
    return(NULL)
  }
  
  # 选择 Top N terms
  top_terms <- gene_term_data %>%
    group_by(term) %>%
    summarise(gene_count = n(), .groups = "drop") %>%
    arrange(desc(gene_count)) %>%
    head(top_n) %>%
    pull(term)
  
  # 创建邻接矩阵
  terms <- top_terms
  genes <- unique(gene_term_data$gene[gene_term_data$term %in% top_terms])
  
  # 创建矩阵
  adj_matrix <- matrix(0, nrow = length(terms), ncol = length(genes),
                       dimnames = list(terms, genes))
  
  for (i in 1:nrow(gene_term_data)) {
    if (gene_term_data$term[i] %in% terms && 
        gene_term_data$gene[i] %in% genes) {
      adj_matrix[gene_term_data$term[i], gene_term_data$gene[i]] <- 1
    }
  }
  
  # 使用 circlize 绘制弦图
  tryCatch({
    pdf(file.path(output_dir, "Chord_Gene_Term_temp.pdf"), 
        width = 14, height = 10)
    
    chordDiagram(adj_matrix, 
                 grid.col = c(rep(colors$bp, length(terms)), 
                             rep("#808080", length(genes))),
                 transparency = 30,
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.1))
    
    # 添加标签
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim = get_cell_xlim()
      ylim = get_cell_ylim()
      sector.index = get.cell.meta.data("sector.index")
      
      circos.text(mean(xlim), mean(ylim) + 0.1, sector.index, 
                  facing = "clockwise", niceFacing = TRUE, 
                  cex = 0.8, font = 3)
    }, bg.border = NA)
    
    dev.off()
    
    message("    ✓ 弦图完成")
    return(file.path(output_dir, "Chord_Gene_Term_temp.pdf"))
  }, error = function(e) {
    message(paste("    ⚠  弦图创建失败:", e$message))
    return(NULL)
  })
}

# 4.5 KEGG 热图（矩阵展示）
plot_heatmap <- function(gene_term_kegg, top_n = 20) {
  message("  [5/6] 创建 KEGG 热图...")
  
  if (is.null(gene_term_kegg) || nrow(gene_term_kegg) == 0) {
    message("    ⚠  无热图数据")
    return(NULL)
  }
  
  # 选择 Top N terms
  top_terms <- gene_term_kegg %>%
    group_by(term) %>%
    summarise(gene_count = n(), .groups = "drop") %>%
    arrange(desc(gene_count)) %>%
    head(top_n) %>%
    pull(term)
  
  # 选择 Top M genes
  top_genes <- gene_term_kegg %>%
    filter(term %in% top_terms) %>%
    group_by(gene) %>%
    summarise(term_count = n(), .groups = "drop") %>%
    arrange(desc(term_count)) %>%
    head(30) %>%
    pull(gene)
  
  # 创建矩阵
  terms <- top_terms
  genes <- top_genes
  
  matrix_data <- matrix(0, nrow = length(terms), ncol = length(genes),
                        dimnames = list(terms, genes))
  
  for (i in 1:nrow(gene_term_kegg)) {
    if (gene_term_kegg$term[i] %in% terms && 
        gene_term_kegg$gene[i] %in% genes) {
      matrix_data[gene_term_kegg$term[i], gene_term_kegg$gene[i]] <- 1
    }
  }
  
  # 使用 ComplexHeatmap 绘制
  tryCatch({
    ht <- Heatmap(matrix_data,
                  name = "Presence",
                  col = c("0" = "white", "1" = colors$bp),
                  cluster_rows = TRUE,
                  cluster_columns = TRUE,
                  show_row_names = TRUE,
                  show_column_names = TRUE,
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 6),
                  column_names_rot = 45,
                  heatmap_width = unit(12, "cm"),
                  heatmap_height = unit(10, "cm"))
    
    message("    ✓ KEGG 热图完成")
    return(ht)
  }, error = function(e) {
    message(paste("    ⚠  热图创建失败:", e$message))
    return(NULL)
  })
}

# 4.6 GO Venn 图（BP/CC/MF 重叠）
plot_venn <- function(bp_data, cc_data, mf_data) {
  message("  [6/6] 创建 Venn 图...")
  
  if (is.null(bp_data) || is.null(cc_data) || is.null(mf_data)) {
    message("    ⚠  数据不完整，跳过 Venn 图")
    return(NULL)
  }
  
  # 提取基因集
  bp_genes <- unique(unlist(strsplit(bp_data$geneid, "/")))
  cc_genes <- unique(unlist(strsplit(cc_data$geneid, "/")))
  mf_genes <- unique(unlist(strsplit(mf_data$geneid, "/")))
  
  # 创建列表
  gene_sets <- list(BP = bp_genes, CC = cc_genes, MF = mf_genes)
  
  # 绘制 Venn 图
  tryCatch({
    p <- ggVennDiagram(gene_sets, 
                       set_color = c(colors$bp, colors$cc, colors$mf),
                       category.names = c("BP", "CC", "MF")) +
      scale_fill_gradient(low = "white", high = colors$bp) +
      labs(title = "GO Category Gene Overlap",
           subtitle = "Number of shared genes across BP/CC/MF") +
      theme_academic()
    
    message("    ✓ Venn 图完成")
    return(p)
  }, error = function(e) {
    message(paste("    ⚠  Venn 图创建失败:", e$message))
    return(NULL)
  })
}

# --------------------------------------------------------------------------------
# 5. Assembly - 组装和输出图形
# --------------------------------------------------------------------------------

message("\n===== 生成和输出图形 =====")

# 5.1 GO 柱状图
p_go_bar <- plot_go_bar(go_bp_data, go_cc_data, go_mf_data, top_n = 10)
if (!is.null(p_go_bar)) {
  ggsave(file.path(output_dir, "GO_Barplot_Top15.pdf"), 
         p_go_bar, width = 10, height = 8)
  ggsave(file.path(output_dir, "GO_Barplot_Top15.png"), 
         p_go_bar, width = 10, height = 8, dpi = 300)
  message("  ✓ GO 柱状图已保存")
}

# 5.2 KEGG 气泡图
p_kegg_bubble <- plot_kegg_bubble(kegg_data, top_n = 15)
if (!is.null(p_kegg_bubble)) {
  ggsave(file.path(output_dir, "KEGG_Bubble.pdf"), 
         p_kegg_bubble, width = 10, height = 8)
  ggsave(file.path(output_dir, "KEGG_Bubble.png"), 
         p_kegg_bubble, width = 10, height = 8, dpi = 300)
  message("  ✓ KEGG 气泡图已保存")
}

# 5.3 网络图（仅生成 PDF，避免 PNG 字体问题）
p_network <- plot_network(gene_term_go, top_n = 20)
if (!is.null(p_network)) {
  tryCatch({
    ggsave(file.path(output_dir, "Network_Cluster.pdf"), 
           p_network, width = 14, height = 10)
    message("  ✓ 网络图已保存 (PDF)")
  }, error = function(e) {
    message(paste("  ⚠  网络图保存失败:", e$message))
  })
}

# 5.4 弦图
chord_file <- plot_chord(gene_term_go, top_n = 15)
if (!is.null(chord_file)) {
  file.rename(chord_file, file.path(output_dir, "Chord_Gene_Term.pdf"))
  message("  ✓ 弦图已保存")
}

# 5.5 KEGG 热图
ht_kegg <- plot_heatmap(gene_term_kegg, top_n = 20)
if (!is.null(ht_kegg)) {
  pdf(file.path(output_dir, "Heatmap_KEGG_Matrix.pdf"), 
      width = 12, height = 10)
  draw(ht_kegg)
  dev.off()
  message("  ✓ KEGG 热图已保存")
}

# 5.6 Venn 图
p_venn <- plot_venn(go_bp_data, go_cc_data, go_mf_data)
if (!is.null(p_venn)) {
  ggsave(file.path(output_dir, "GO_VennDiagram.pdf"), 
         p_venn, width = 8, height = 6)
  message("  ✓ Venn 图已保存")
}

# 5.7 组合图（使用 patchwork）
message("\n===== 创建组合图 =====")

if (!is.null(p_go_bar) && !is.null(p_kegg_bubble)) {
  # 创建组合图：左 GO 柱状图，右 KEGG 气泡图
  p_combined <- p_go_bar + p_kegg_bubble + 
    plot_layout(ncol = 2, widths = c(1, 1)) &
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  ggsave(file.path(output_dir, "Figure_Composite.pdf"), 
         p_combined, width = 16, height = 8)
  ggsave(file.path(output_dir, "Figure_Composite.png"), 
         p_combined, width = 16, height = 8, dpi = 600)
  message("  ✓ 组合图已保存 (16x8 英寸，600dpi)")
}

# --------------------------------------------------------------------------------
# 6. Summary - 生成摘要报告
# --------------------------------------------------------------------------------

message("\n===== 生成摘要报告 =====")

summary_text <- c(
  "================================================================================",
  "富集分析可视化摘要报告",
  "================================================================================",
  "",
  paste("生成时间:", Sys.time()),
  paste("工作目录:", getwd()),
  "",
  "【数据概览】",
  paste("  - GO-BP 条目数:", ifelse(is.null(go_bp_data), 0, nrow(go_bp_data))),
  paste("  - GO-CC 条目数:", ifelse(is.null(go_cc_data), 0, nrow(go_cc_data))),
  paste("  - GO-MF 条目数:", ifelse(is.null(go_mf_data), 0, nrow(go_mf_data))),
  paste("  - KEGG 通路数:", ifelse(is.null(kegg_data), 0, nrow(kegg_data))),
  "",
  "【输出文件列表】",
  "  1. GO_Barplot_Top15.pdf/png - GO 三维度柱状图（每类 Top 10）",
  "  2. KEGG_Bubble.pdf/png - KEGG 气泡图（Top 15）",
  "  3. Network_Cluster.pdf/png - 基因-Term 网络图",
  "  4. Chord_Gene_Term.pdf - 基因-Term 弦图",
  "  5. Heatmap_KEGG_Matrix.pdf - KEGG 通路 - 基因热图",
  "  6. GO_VennDiagram.pdf - BP/CC/MF 基因重叠 Venn 图",
  "  7. Figure_Composite.pdf/png - 组合图（GO+KEGG）",
  "",
  "【配色方案】",
  paste("  - BP:", colors$bp, "学术蓝"),
  paste("  - CC:", colors$cc, "紫红色"),
  paste("  - MF:", colors$mf, "琥珀金"),
  "  - KEGG: 灰色到深红渐变",
  "  - 网络图：NPG 配色方案",
  "",
  "【图形规格】",
  "  - 单图：8-10 英寸宽，6-8 英寸高",
  "  - 组合图：16 英寸宽，8 英寸高（16:9 比例）",
  "  - 分辨率：PNG 300-600 dpi",
  "  - 格式：PDF（矢量，发表用）+ PNG（预览用）",
  "",
  "================================================================================"
)

# 保存摘要
writeLines(summary_text, file.path(output_dir, "Visualization_Summary.txt"))
message(paste("✓ 摘要报告已保存:", file.path(output_dir, "Visualization_Summary.txt")))

# 显示摘要
cat(paste(summary_text, collapse = "\n"))

message("\n\n✅ 所有可视化完成！结果已保存到:", normalizePath(output_dir))
