# ================================================================================
# GO 富集分析脚本
# 基于 clusterProfiler 包进行基因本体论（Gene Ontology）富集分析
# ================================================================================

# --------------------------------------------------------------------------------
# 1. 环境设置和包加载
# --------------------------------------------------------------------------------

# 设置 CRAN 镜像（使用清华镜像源，国内访问更快）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 智能安装必要的 R 包：只安装未安装的包，避免重复安装
required_packages <- c("clusterProfiler", "org.Hs.eg.db", "tidyverse", 
                       "enrichplot", "DOSE", "ggplot2", "ggpubr", "cowplot")

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
library(ggpubr)
library(cowplot)

# --------------------------------------------------------------------------------
# 2. 数据读取和预处理
# --------------------------------------------------------------------------------

# 设置工作目录和输出目录
output_dir <- "GO_enrichment_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 定义基因列表（58 个基因）
gene_list <- c(
  "TP53", "IL1B", "IL6", "TNF", "STAT3", "BCL2", "NFKB1", "PTGS2", "TLR4", "SRC",
  "STAT1", "RELA", "ICAM1", "CCL2", "CCL5", "CASP8", "VCAM1", "TGFB1", "PTPRC", "IKBKB",
  "STAT5A", "CCND1", "HMOX1", "TIMP1", "NLRP3", "CDK4", "PARP1", "CCR5", "FAS", "MAPK9",
  "NFE2L2", "SREBF1", "IRF1", "IL10RA", "CXCR3", "PGR", "BID", "EGR1", "F3", "AIF1",
  "CTSS", "PTGS1", "IRAK4", "LYN", "SREBF2", "TOP2A", "GFAP", "CCNA2", "PTGES", "PTPN2",
  "ERBB4", "CTSD", "CTSB", "C3", "SQLE", "HMGCR", "LSS", "CYP51A1"
)

# 去重和清理
gene_list <- unique(gene_list)
gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]

message(paste("输入基因总数:", length(gene_list)))

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
    toType = "ENTREZID",      # 输出 ID 类型：Entrez ID（GO 分析需要）
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
} else {
  message("✓ 所有基因都成功映射！")
}

# 去重：一个基因可能对应多个 Entrez ID，取第一个
gene_entrez <- gene_mapping %>%
  distinct(SYMBOL, .keep_all = TRUE)

# 提取 Entrez ID 向量用于富集分析
entrez_ids <- gene_entrez$ENTREZID

message(paste("用于富集分析的 Entrez ID 数:", length(entrez_ids)))

# --------------------------------------------------------------------------------
# 4. GO 富集分析（三个本体：BP、CC、MF）
# --------------------------------------------------------------------------------

message("\n===== GO 富集分析 =====")

# 4.1 生物过程（Biological Process, BP）
# BP 描述基因参与的生物学过程或通路
message("\n[1/3] 进行生物过程 (BP) 富集分析...")

go_bp_result <- tryCatch({
  enrichGO(
    gene = entrez_ids,           # 输入基因列表（Entrez ID）
    OrgDb = org.Hs.eg.db,        # 物种数据库
    ont = "BP",                  # 本体类型：生物过程
    pvalueCutoff = 0.05,         # P 值阈值
    pAdjustMethod = "BH",        # 多重检验校正方法：Benjamini-Hochberg
    qvalueCutoff = 0.2,          # Q 值阈值（FDR 校正后的 P 值）
    minGSSize = 10,              # 通路最小基因数
    maxGSSize = 500,             # 通路最大基因数
    readable = TRUE              # 将 Entrez ID 转换为基因符号（便于阅读）
  )
}, error = function(e) {
  message(paste("BP 富集分析警告:", e$message))
  return(NULL)
})

if (!is.null(go_bp_result)) {
  go_bp_df <- as.data.frame(go_bp_result)
  message(paste("  ✓ BP 显著富集条目数:", nrow(go_bp_df)))
  if (nrow(go_bp_df) > 0) {
    message(paste("  Top 3 BP:", paste(head(go_bp_df$Description, 3), collapse = "; ")))
  }
} else {
  message("  ⚠ BP 富集分析未获得显著结果")
}

# 4.2 细胞组分（Cellular Component, CC）
# CC 描述基因产物在细胞中的位置
message("\n[2/3] 进行细胞组分 (CC) 富集分析...")

go_cc_result <- tryCatch({
  enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "CC",                  # 本体类型：细胞组分
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
}, error = function(e) {
  message(paste("CC 富集分析警告:", e$message))
  return(NULL)
})

if (!is.null(go_cc_result)) {
  go_cc_df <- as.data.frame(go_cc_result)
  message(paste("  ✓ CC 显著富集条目数:", nrow(go_cc_df)))
  if (nrow(go_cc_df) > 0) {
    message(paste("  Top 3 CC:", paste(head(go_cc_df$Description, 3), collapse = "; ")))
  }
} else {
  message("  ⚠ CC 富集分析未获得显著结果")
}

# 4.3 分子功能（Molecular Function, MF）
# MF 描述基因产物的分子功能或活性
message("\n[3/3] 进行分子功能 (MF) 富集分析...")

go_mf_result <- tryCatch({
  enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "MF",                  # 本体类型：分子功能
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
}, error = function(e) {
  message(paste("MF 富集分析警告:", e$message))
  return(NULL)
})

if (!is.null(go_mf_result)) {
  go_mf_df <- as.data.frame(go_mf_result)
  message(paste("  ✓ MF 显著富集条目数:", nrow(go_mf_df)))
  if (nrow(go_mf_df) > 0) {
    message(paste("  Top 3 MF:", paste(head(go_mf_df$Description, 3), collapse = "; ")))
  }
} else {
  message("  ⚠ MF 富集分析未获得显著结果")
}

# --------------------------------------------------------------------------------
# 5. 结果导出
# --------------------------------------------------------------------------------

message("\n===== 导出结果 =====")

# 5.1 导出 BP 结果
if (!is.null(go_bp_result)) {
  bp_file <- file.path(output_dir, "GO_BP_enrichment_results.csv")
  go_bp_export <- go_bp_df %>%
    select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, 
           qvalue, geneID, Count) %>%
    arrange(pvalue)
  write.csv(go_bp_export, bp_file, row.names = FALSE)
  message(paste("✓ BP 结果已保存:", bp_file))
}

# 5.2 导出 CC 结果
if (!is.null(go_cc_result)) {
  cc_file <- file.path(output_dir, "GO_CC_enrichment_results.csv")
  go_cc_export <- go_cc_df %>%
    select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, 
           qvalue, geneID, Count) %>%
    arrange(pvalue)
  write.csv(go_cc_export, cc_file, row.names = FALSE)
  message(paste("✓ CC 结果已保存:", cc_file))
}

# 5.3 导出 MF 结果
if (!is.null(go_mf_result)) {
  mf_file <- file.path(output_dir, "GO_MF_enrichment_results.csv")
  go_mf_export <- go_mf_df %>%
    select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, 
           qvalue, geneID, Count) %>%
    arrange(pvalue)
  write.csv(go_mf_export, mf_file, row.names = FALSE)
  message(paste("✓ MF 结果已保存:", mf_file))
}

# 5.4 导出基因映射结果
mapping_file <- file.path(output_dir, "gene_mapping.csv")
write.csv(gene_mapping, mapping_file, row.names = FALSE)
message(paste("✓ 基因映射结果已保存:", mapping_file))

# --------------------------------------------------------------------------------
# 6. 数据可视化
# --------------------------------------------------------------------------------

message("\n===== 生成可视化图形 =====")

# 设置图形主题（tidyverse 风格）
theme_set(theme_bw(base_size = 12))

# 6.1 GO-BP Dotplot（Top 20）
message("\n[1/8] 生成 BP Dotplot...")
bp_dotplot_file <- file.path(output_dir, "GO_BP_dotplot_Top20.pdf")
pdf(bp_dotplot_file, width = 12, height = 10)

if (!is.null(go_bp_result) && nrow(go_bp_df) > 0) {
  tryCatch({
    p_bp_dot <- dotplot(go_bp_result, 
                        showCategory = 20,
                        title = "GO Biological Process (Top 20)",
                        font.size = 12) +
      theme(axis.text.y = element_text(size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    print(p_bp_dot)
    message("  ✓ 已生成：GO_BP_dotplot_Top20.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Dotplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant BP Terms", cex = 1.5)
}

dev.off()

# 6.2 GO-BP Barplot（Top 20）
message("[2/8] 生成 BP Barplot...")
bp_barplot_file <- file.path(output_dir, "GO_BP_barplot_Top20.pdf")
pdf(bp_barplot_file, width = 12, height = 10)

if (!is.null(go_bp_result) && nrow(go_bp_df) > 0) {
  tryCatch({
    p_bp_bar <- barplot(go_bp_result, 
                        showCategory = 20,
                        title = "GO Biological Process (Top 20)",
                        font.size = 12,
                        colour = "p.adjust") +
      theme(axis.text.y = element_text(size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    print(p_bp_bar)
    message("  ✓ 已生成：GO_BP_barplot_Top20.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Barplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant BP Terms", cex = 1.5)
}

dev.off()

# 6.3 GO-CC Dotplot（Top 20）
message("[3/8] 生成 CC Dotplot...")
cc_dotplot_file <- file.path(output_dir, "GO_CC_dotplot_Top20.pdf")
pdf(cc_dotplot_file, width = 12, height = 10)

if (!is.null(go_cc_result) && nrow(go_cc_df) > 0) {
  tryCatch({
    p_cc_dot <- dotplot(go_cc_result, 
                        showCategory = 20,
                        title = "GO Cellular Component (Top 20)",
                        font.size = 12) +
      theme(axis.text.y = element_text(size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    print(p_cc_dot)
    message("  ✓ 已生成：GO_CC_dotplot_Top20.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Dotplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant CC Terms", cex = 1.5)
}

dev.off()

# 6.4 GO-CC Barplot（Top 20）
message("[4/8] 生成 CC Barplot...")
cc_barplot_file <- file.path(output_dir, "GO_CC_barplot_Top20.pdf")
pdf(cc_barplot_file, width = 12, height = 10)

if (!is.null(go_cc_result) && nrow(go_cc_df) > 0) {
  tryCatch({
    p_cc_bar <- barplot(go_cc_result, 
                        showCategory = 20,
                        title = "GO Cellular Component (Top 20)",
                        font.size = 12,
                        colour = "p.adjust") +
      theme(axis.text.y = element_text(size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    print(p_cc_bar)
    message("  ✓ 已生成：GO_CC_barplot_Top20.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Barplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant CC Terms", cex = 1.5)
}

dev.off()

# 6.5 GO-MF Dotplot（Top 20）
message("[5/8] 生成 MF Dotplot...")
mf_dotplot_file <- file.path(output_dir, "GO_MF_dotplot_Top20.pdf")
pdf(mf_dotplot_file, width = 12, height = 10)

if (!is.null(go_mf_result) && nrow(go_mf_df) > 0) {
  tryCatch({
    p_mf_dot <- dotplot(go_mf_result, 
                        showCategory = 20,
                        title = "GO Molecular Function (Top 20)",
                        font.size = 12) +
      theme(axis.text.y = element_text(size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    print(p_mf_dot)
    message("  ✓ 已生成：GO_MF_dotplot_Top20.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Dotplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant MF Terms", cex = 1.5)
}

dev.off()

# 6.6 GO-MF Barplot（Top 20）
message("[6/8] 生成 MF Barplot...")
mf_barplot_file <- file.path(output_dir, "GO_MF_barplot_Top20.pdf")
pdf(mf_barplot_file, width = 12, height = 10)

if (!is.null(go_mf_result) && nrow(go_mf_df) > 0) {
  tryCatch({
    p_mf_bar <- barplot(go_mf_result, 
                        showCategory = 20,
                        title = "GO Molecular Function (Top 20)",
                        font.size = 12,
                        colour = "p.adjust") +
      theme(axis.text.y = element_text(size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    print(p_mf_bar)
    message("  ✓ 已生成：GO_MF_barplot_Top20.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Barplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant MF Terms", cex = 1.5)
}

dev.off()

# 6.7 GO 富集网络图（Cnetplot - BP Top 10）
message("[7/8] 生成 GO 网络图...")
go_cnetplot_file <- file.path(output_dir, "GO_cnetplot.pdf")
pdf(go_cnetplot_file, width = 14, height = 12)

if (!is.null(go_bp_result) && nrow(go_bp_df) > 0) {
  tryCatch({
    p_cnet <- cnetplot(go_bp_result,
                       showCategory = 10,      # 显示 Top 10 条目
                       foldChange = NULL,      # 无表达量数据
                       circular = FALSE,       # 线性布局
                       colorEdge = TRUE,       # 彩色边缘
                       cex_label_gene = 0.8,   # 基因标签大小
                       cex_label_category = 1) # 条目标签大小
    title("Gene-GO Network (BP Top 10)", line = 2)
    message("  ✓ 已生成：GO_cnetplot.pdf")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Cnetplot 生成失败:", e$message), cex = 1.2)
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant GO Terms", cex = 1.5)
}

dev.off()

# 6.8 GO 有向无环图（DAG - BP Top 10）
message("[8/8] 生成 GO 有向无环图...")
go_dag_file <- file.path(output_dir, "GO_DAG_plot.pdf")
pdf(go_dag_file, width = 12, height = 10)

if (!is.null(go_bp_result) && nrow(go_bp_df) > 0) {
  tryCatch({
    # 使用 simplify 减少冗余条目
    go_bp_simplified <- simplify(go_bp_result, 
                                  cutoff = 0.7, 
                                  by = "p.adjust", 
                                  select_fun = min)
    
    p_dag <- dotplot(go_bp_simplified, 
                     showCategory = 15,
                     title = "GO Biological Process (Simplified)",
                     font.size = 12) +
      theme(axis.text.y = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    print(p_dag)
    message("  ✓ 已生成：GO_DAG_plot.pdf")
  }, error = function(e) {
    # 如果 simplify 失败，使用原始结果
    tryCatch({
      p_dag <- dotplot(go_bp_result, 
                       showCategory = 15,
                       title = "GO Biological Process",
                       font.size = 12) +
        theme(axis.text.y = element_text(size = 10),
              plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
      print(p_dag)
      message("  ✓ 已生成：GO_DAG_plot.pdf (未简化)")
    }, error = function(e2) {
      plot.new()
      text(0.5, 0.5, paste("DAG 生成失败:", e2$message), cex = 1.2)
    })
  })
} else {
  plot.new()
  text(0.5, 0.5, "No Significant GO Terms", cex = 1.5)
}

dev.off()

# --------------------------------------------------------------------------------
# 7. 生成分析摘要报告
# --------------------------------------------------------------------------------

message("\n===== 分析结果摘要 =====")

# 统计各本体结果
bp_count <- ifelse(is.null(go_bp_result), 0, nrow(go_bp_df))
cc_count <- ifelse(is.null(go_cc_result), 0, nrow(go_cc_df))
mf_count <- ifelse(is.null(go_mf_result), 0, nrow(go_mf_df))

summary_lines <- c(
  "================================================================================",
  "GO 富集分析摘要报告",
  "================================================================================",
  "",
  paste("分析时间:", Sys.time()),
  "",
  "【基因统计】",
  paste("  - 输入基因总数:", length(gene_list)),
  paste("  - 成功映射基因数:", length(mapped_genes)),
  paste("  - 未映射基因数:", length(unmapped_genes)),
  paste("  - 映射成功率:", round(length(mapped_genes) / length(gene_list) * 100, 2), "%"),
  "",
  "【GO 富集分析参数】",
  "  - 物种：人类 (Homo sapiens)",
  "  - 数据库：org.Hs.eg.db",
  "  - 分析方法：ORA (Over-Representation Analysis)",
  "  - P 值校正方法：Benjamini-Hochberg (BH)",
  "  - P 值阈值：0.05",
  "  - Q 值阈值：0.2",
  "  - 条目大小范围：10-500 个基因",
  "",
  "【富集结果】",
  paste("  - BP (生物过程) 显著富集条目数:", bp_count),
  paste("  - CC (细胞组分) 显著富集条目数:", cc_count),
  paste("  - MF (分子功能) 显著富集条目数:", mf_count),
  "",
  "【Top 5 BP 条目】"
)

# 添加 Top 5 BP
if (bp_count > 0) {
  bp_top5 <- head(go_bp_df, 5)
  for (i in 1:nrow(bp_top5)) {
    summary_lines <- c(summary_lines, 
      paste0("  ", i, ". ", bp_top5$Description[i], 
             " (p.adjust = ", format(bp_top5$p.adjust[i], scientific = TRUE), ")"))
  }
}

summary_lines <- c(summary_lines, "", "【Top 5 CC 条目】")

# 添加 Top 5 CC
if (cc_count > 0) {
  cc_top5 <- head(go_cc_df, 5)
  for (i in 1:nrow(cc_top5)) {
    summary_lines <- c(summary_lines, 
      paste0("  ", i, ". ", cc_top5$Description[i], 
             " (p.adjust = ", format(cc_top5$p.adjust[i], scientific = TRUE), ")"))
  }
}

summary_lines <- c(summary_lines, "", "【Top 5 MF 条目】")

# 添加 Top 5 MF
if (mf_count > 0) {
  mf_top5 <- head(go_mf_df, 5)
  for (i in 1:nrow(mf_top5)) {
    summary_lines <- c(summary_lines, 
      paste0("  ", i, ". ", mf_top5$Description[i], 
             " (p.adjust = ", format(mf_top5$p.adjust[i], scientific = TRUE), ")"))
  }
}

summary_lines <- c(summary_lines,
  "",
  "【输出文件】",
  "  - BP 结果：GO_BP_enrichment_results.csv",
  "  - CC 结果：GO_CC_enrichment_results.csv",
  "  - MF 结果：GO_MF_enrichment_results.csv",
  "  - 基因映射：gene_mapping.csv",
  "  - 未映射基因：unmapped_genes.txt",
  "  - BP Dotplot: GO_BP_dotplot_Top20.pdf",
  "  - BP Barplot: GO_BP_barplot_Top20.pdf",
  "  - CC Dotplot: GO_CC_dotplot_Top20.pdf",
  "  - CC Barplot: GO_CC_barplot_Top20.pdf",
  "  - MF Dotplot: GO_MF_dotplot_Top20.pdf",
  "  - MF Barplot: GO_MF_barplot_Top20.pdf",
  "  - 网络图：GO_cnetplot.pdf",
  "  - DAG 图：GO_DAG_plot.pdf",
  "",
  "================================================================================"
)

# 保存摘要报告
summary_file <- file.path(output_dir, "analysis_summary.txt")
writeLines(summary_lines, summary_file)
message(paste("✓ 分析摘要已保存:", summary_file))

# 在控制台显示摘要
cat(paste(summary_lines, collapse = "\n"))

message("\n✓ GO 富集分析完成！所有结果已保存到:", normalizePath(output_dir))
