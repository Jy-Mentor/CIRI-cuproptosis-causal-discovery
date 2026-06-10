#!/usr/bin/env Rscript
# ================================================================================
# 功能富集分析脚本 - MR 分析优化路线 B
# 包含：GO 富集、KEGG 通路、Reactome 通路
# ================================================================================

# 包安装与加载
install_and_load_packages <- function() {
  packages <- c(
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichplot",
    "ggplot2",
    "pathview",
    "ReactomePA",
    "dplyr",
    "readr"
  )
  
  message("正在检查并安装所需包...")
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("正在安装:", pkg))
      tryCatch({
        install.packages(pkg, repos = "https://cloud.r-project.org/")
      }, error = function(e) {
        message(paste("安装失败:", pkg, "-", e$message))
      })
    }
  }
  
  # 加载 Bioconductor 包
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  bioc_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", 
                     "pathview", "ReactomePA")
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE)
    }
  }
  
  invisible(lapply(packages, library, character.only = TRUE))
  message("所有包加载完成")
}

# 从 MR 结果中提取显著基因
extract_significant_genes <- function(mr_results_file, fdr_threshold = 0.05, 
                                       pval_threshold = 0.05) {
  message("\n=== 提取显著基因 ===")
  
  results <- read.csv(mr_results_file, stringsAsFactors = FALSE)
  
  # 筛选成功的基因
  successful_genes <- results %>%
    filter(status == "SUCCESS" | status == "HETEROGENEITY") %>%
    filter(!is.na(discovery_pval))
  
  message(paste("总成功基因数:", nrow(successful_genes)))
  
  # FDR 显著基因
  fdr_sig_genes <- successful_genes %>%
    filter(!is.na(fdr_qval)) %>%
    filter(fdr_qval < fdr_threshold) %>%
    pull(gene) %>%
    unique()
  
  message(paste("FDR 显著基因 (q < 0.05):", length(fdr_sig_genes)))
  if (length(fdr_sig_genes) > 0) {
    message("  ", paste(fdr_sig_genes, collapse = ", "))
  }
  
  # P 值显著基因 (未校正)
  pval_sig_genes <- successful_genes %>%
    filter(discovery_pval < pval_threshold) %>%
    pull(gene) %>%
    unique()
  
  message(paste("P 值显著基因 (p < 0.05):", length(pval_sig_genes)))
  
  # 所有成功基因（用于背景）
  all_genes <- successful_genes %>% pull(gene) %>% unique()
  
  message(paste("所有成功基因:", length(all_genes)))
  
  return(list(
    fdr_significant = fdr_sig_genes,
    pval_significant = pval_sig_genes,
    all_successful = all_genes,
    results_table = results
  ))
}

# GO 富集分析
run_go_enrichment <- function(gene_list, org_db = org.Hs.eg.db) {
  message("\n=== GO 富集分析 ===")
  
  # 基因 ID 转换
  gene_symbols <- gene_list
  gene_ids <- bitr(gene_symbols, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org_db)
  
  if (nrow(gene_ids) == 0) {
    message("无法转换基因 ID")
    return(NULL)
  }
  
  message(paste("成功转换基因 ID:", nrow(gene_ids), "/", length(gene_symbols)))
  
  # GO Biological Process
  message("  进行 BP 分析...")
  go_bp <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org_db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # GO Molecular Function
  message("  进行 MF 分析...")
  go_mf <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org_db,
    keyType = "SYMBOL",
    ont = "MF",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # GO Cellular Component
  message("  进行 CC 分析...")
  go_cc <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org_db,
    keyType = "SYMBOL",
    ont = "CC",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  message(paste("  BP 显著条目:", nrow(go_bp)))
  message(paste("  MF 显著条目:", nrow(go_mf)))
  message(paste("  CC 显著条目:", nrow(go_cc)))
  
  return(list(
    bp = go_bp,
    mf = go_mf,
    cc = go_cc,
    gene_mapping = gene_ids
  ))
}

# KEGG 通路富集分析
run_kegg_enrichment <- function(gene_list) {
  message("\n=== KEGG 通路富集分析 ===")
  
  # 基因 ID 转换
  gene_symbols <- gene_list
  gene_ids <- bitr(gene_symbols, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_ids) == 0) {
    message("无法转换基因 ID")
    return(NULL)
  }
  
  message(paste("成功转换基因 ID:", nrow(gene_ids)))
  
  # KEGG 分析 (使用内部数据避免网络问题)
  message("  使用内部 KEGG 数据...")
  kegg_result <- tryCatch({
    enrichKEGG(
      gene = gene_ids$ENTREZID,
      organism = "hsa",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      pvalueCutoff = 0.05,
      use_internal_data = TRUE
    )
  }, error = function(e) {
    message(paste("  KEGG 分析失败:", e$message))
    message("  尝试在线访问...")
    tryCatch({
      enrichKEGG(
        gene = gene_ids$ENTREZID,
        organism = "hsa",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        pvalueCutoff = 0.05,
        use_internal_data = FALSE
      )
    }, error = function(e2) {
      message(paste("  KEGG 在线访问也失败:", e2$message))
      return(NULL)
    })
  })
  
  if (!is.null(kegg_result)) {
    message(paste("  显著通路:", nrow(kegg_result)))
    if (nrow(kegg_result) > 0) {
      message("  Top 5 通路:")
      for (i in 1:min(5, nrow(kegg_result))) {
        message(paste("    ", i, ". ", kegg_result$Description[i], 
                      " (q=", format(kegg_result$qvalue[i], digits = 3), ")", sep = ""))
      }
    }
  } else {
    message("  无显著通路")
  }
  
  return(kegg_result)
}

# Reactome 通路富集分析
run_reactome_enrichment <- function(gene_list) {
  message("\n=== Reactome 通路富集分析 ===")
  
  # 基因 ID 转换
  gene_symbols <- gene_list
  gene_ids <- bitr(gene_symbols, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_ids) == 0) {
    message("无法转换基因 ID")
    return(NULL)
  }
  
  message(paste("成功转换基因 ID:", nrow(gene_ids)))
  
  # Reactome 分析
  reactome_result <- enrichPathway(
    gene = gene_ids$ENTREZID,
    organism = "human",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  if (!is.null(reactome_result)) {
    message(paste("  显著通路:", nrow(reactome_result)))
    if (nrow(reactome_result) > 0) {
      message("  Top 5 通路:")
      for (i in 1:min(5, nrow(reactome_result))) {
        message(paste("    ", i, ". ", reactome_result$Description[i], 
                      " (q=", format(reactome_result$qvalue[i], digits = 3), ")", sep = ""))
      }
    }
  } else {
    message("  无显著通路")
  }
  
  return(reactome_result)
}

# 可视化并保存结果
save_enrichment_results <- function(go_results, kegg_result, reactome_result, 
                                     output_dir = "./functional_enrichment") {
  message("\n=== 保存富集分析结果 ===")
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. GO BP 气泡图
  if (!is.null(go_results$bp) && nrow(go_results$bp) > 0) {
    message("  保存 GO BP 气泡图...")
    p_bp <- dotplot(go_results$bp, showCategory = 20, split = ".") +
      ggtitle("GO Biological Process Enrichment") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    ggsave(file.path(output_dir, "GO_BP_dotplot.png"), p_bp, 
           width = 10, height = 8, dpi = 300)
    
    # 保存详细结果
    write.csv(as.data.frame(go_results$bp), 
              file.path(output_dir, "GO_BP_results.csv"), row.names = FALSE)
  }
  
  # 2. GO MF 气泡图
  if (!is.null(go_results$mf) && nrow(go_results$mf) > 0) {
    message("  保存 GO MF 气泡图...")
    p_mf <- dotplot(go_results$mf, showCategory = 15) +
      ggtitle("GO Molecular Function Enrichment")
    ggsave(file.path(output_dir, "GO_MF_dotplot.png"), p_mf, 
           width = 10, height = 8, dpi = 300)
    write.csv(as.data.frame(go_results$mf), 
              file.path(output_dir, "GO_MF_results.csv"), row.names = FALSE)
  }
  
  # 3. GO CC 气泡图
  if (!is.null(go_results$cc) && nrow(go_results$cc) > 0) {
    message("  保存 GO CC 气泡图...")
    p_cc <- dotplot(go_results$cc, showCategory = 15) +
      ggtitle("GO Cellular Component Enrichment")
    ggsave(file.path(output_dir, "GO_CC_dotplot.png"), p_cc, 
           width = 10, height = 8, dpi = 300)
    write.csv(as.data.frame(go_results$cc), 
              file.path(output_dir, "GO_CC_results.csv"), row.names = FALSE)
  }
  
  # 4. KEGG 气泡图
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    message("  保存 KEGG 气泡图...")
    p_kegg <- dotplot(kegg_result, showCategory = 20) +
      ggtitle("KEGG Pathway Enrichment")
    ggsave(file.path(output_dir, "KEGG_dotplot.png"), p_kegg, 
           width = 10, height = 8, dpi = 300)
    write.csv(as.data.frame(kegg_result), 
              file.path(output_dir, "KEGG_results.csv"), row.names = FALSE)
    
    # KEGG 通路图（前 3 个）
    gene_ids <- bitr(rownames(kegg_result@result), fromType = "ENTREZID",
                     toType = "SYMBOL", OrgDb = org.Hs.eg.db)
    for (i in 1:min(3, nrow(kegg_result))) {
      pathway_id <- kegg_result$ID[i]
      message(paste("    生成 KEGG 通路图:", pathway_id))
      tryCatch({
        pathview(gene.data = gene_ids$SYMBOL, 
                 pathway.id = pathway_id,
                 species = "hsa",
                 out.suffix = "MR_significant",
                 kegg.dir = output_dir)
      }, error = function(e) {
        message(paste("    通路图生成失败:", e$message))
      })
    }
  }
  
  # 5. Reactome 气泡图
  if (!is.null(reactome_result) && nrow(reactome_result) > 0) {
    message("  保存 Reactome 气泡图...")
    p_reactome <- dotplot(reactome_result, showCategory = 20) +
      ggtitle("Reactome Pathway Enrichment")
    ggsave(file.path(output_dir, "Reactome_dotplot.png"), p_reactome, 
           width = 10, height = 8, dpi = 300)
    write.csv(as.data.frame(reactome_result), 
              file.path(output_dir, "Reactome_results.csv"), row.names = FALSE)
  }
  
  # 6. GO 有向无环图（前 10 个最显著的 BP）
  if (!is.null(go_results$bp) && nrow(go_results$bp) > 0) {
    message("  保存 GO 有向无环图...")
    go_bp_filtered <- go_results$bp %>%
      as.data.frame() %>%
      head(10)
    
    if (nrow(go_bp_filtered) > 0) {
      # 使用 simplify 去除冗余条目
      go_bp_simple <- simplify(go_results$bp, 
                               cutoff = 0.7, 
                               by = "p.adjust",
                               select_fun = min)
      
      p_dag <- heatplot(go_bp_simple, showCategory = 10) +
        ggtitle("GO BP Hierarchical Clustering")
      ggsave(file.path(output_dir, "GO_BP_heatplot.png"), p_dag, 
             width = 10, height = 8, dpi = 300)
    }
  }
  
  message(paste("所有结果已保存至:", output_dir))
}

# 生成富集分析报告
generate_enrichment_report <- function(gene_list, go_results, kegg_result, 
                                        reactome_result, output_file) {
  message("\n=== 生成富集分析报告 ===")
  
  report <- c(
    "# 功能富集分析报告",
    "",
    "## 分析概述",
    paste("分析日期:", Sys.time()),
    paste("输入基因数:", length(gene_list)),
    "",
    "## 显著基因列表",
    paste("", paste(gene_list, collapse = ", ")),
    "",
    "## GO Biological Process (Top 10)",
    ""
  )
  
  if (!is.null(go_results$bp) && nrow(go_results$bp) > 0) {
    bp_df <- as.data.frame(go_results$bp)
    for (i in 1:min(10, nrow(bp_df))) {
      report <- c(report, paste0(
        i, ". **", bp_df$Description[i], "**  \n",
        "   - 基因数：", bp_df$Count[i], "  \n",
        "   - P 值：", format(bp_df$pvalue[i], scientific = TRUE), "  \n",
        "   - q 值：", format(bp_df$p.adjust[i], digits = 3), "  \n",
        "   - 基因：", bp_df$geneID[i], "  \n"
      ))
    }
  }
  
  report <- c(report, "", "## KEGG Pathway (Top 10)", "")
  
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    kegg_df <- as.data.frame(kegg_result)
    for (i in 1:min(10, nrow(kegg_df))) {
      report <- c(report, paste0(
        i, ". **", kegg_df$Description[i], "**  \n",
        "   - 基因数：", kegg_df$Count[i], "  \n",
        "   - P 值：", format(kegg_df$pvalue[i], scientific = TRUE), "  \n",
        "   - q 值：", format(kegg_df$p.adjust[i], digits = 3), "  \n",
        "   - 基因：", kegg_df$geneID[i], "  \n"
      ))
    }
  }
  
  report <- c(report, "", "## Reactome Pathway (Top 10)", "")
  
  if (!is.null(reactome_result) && nrow(reactome_result) > 0) {
    reactome_df <- as.data.frame(reactome_result)
    for (i in 1:min(10, nrow(reactome_df))) {
      report <- c(report, paste0(
        i, ". **", reactome_df$Description[i], "**  \n",
        "   - 基因数：", reactome_df$Count[i], "  \n",
        "   - P 值：", format(reactome_df$pvalue[i], scientific = TRUE), "  \n",
        "   - q 值：", format(reactome_df$p.adjust[i], digits = 3), "  \n",
        "   - 基因：", reactome_df$geneID[i], "  \n"
      ))
    }
  }
  
  writeLines(report, output_file)
  message(paste("报告已保存:", output_file))
}

# 主函数
main <- function() {
  message(rep("=", 60))
  message("功能富集分析 - MR 分析优化路线 B")
  message(rep("=", 60))
  
  # 1. 加载包
  install_and_load_packages()
  
  # 2. 提取显著基因
  mr_results_file <- "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/MR_results_main_optimized.csv"
  
  if (!file.exists(mr_results_file)) {
    stop("找不到 MR 结果文件")
  }
  
  genes <- extract_significant_genes(mr_results_file, 
                                      fdr_threshold = 0.05,
                                      pval_threshold = 0.05)
  
  # 使用 FDR 显著基因 + P 值显著基因
  significant_gene_list <- unique(c(genes$fdr_significant, genes$pval_significant))
  
  if (length(significant_gene_list) == 0) {
    message("无显著基因，使用所有成功基因")
    significant_gene_list <- genes$all_successful
  }
  
  message(paste("\n用于富集分析的基因数:", length(significant_gene_list)))
  message("基因列表:", paste(significant_gene_list, collapse = ", "))
  
  # 3. GO 富集分析
  message("\n", rep("-", 60))
  go_results <- run_go_enrichment(significant_gene_list)
  
  if (!is.null(go_results)) {
    message("\nGO 富集分析完成")
    if (!is.null(go_results$bp) && nrow(go_results$bp) > 0) {
      message("  Top 3 BP 条目:")
      bp_df <- as.data.frame(go_results$bp)
      for (i in 1:min(3, nrow(bp_df))) {
        message(paste("    ", i, ". ", bp_df$Description[i], 
                      " (q=", format(bp_df$p.adjust[i], digits = 3), ")", sep = ""))
      }
    }
  }
  
  # 4. KEGG 富集分析
  message("\n", rep("-", 60))
  kegg_result <- run_kegg_enrichment(significant_gene_list)
  
  # 5. Reactome 富集分析
  message("\n", rep("-", 60))
  reactome_result <- run_reactome_enrichment(significant_gene_list)
  
  # 6. 保存结果
  message("\n", rep("-", 60))
  output_dir <- "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/functional_enrichment"
  save_enrichment_results(go_results, kegg_result, reactome_result, output_dir)
  
  # 7. 生成报告
  report_file <- file.path(output_dir, "enrichment_report.md")
  generate_enrichment_report(significant_gene_list, go_results, 
                              kegg_result, reactome_result, report_file)
  
  message("\n", rep("=", 60))
  message("功能富集分析完成!")
  message(paste("结果目录:", output_dir))
  message(rep("=", 60))
}

# 运行主函数
if (!interactive()) {
  main()
} else {
  message("交互模式下运行，请手动调用 main()")
}
