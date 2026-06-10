#!/usr/bin/env Rscript
# ================================================================================
# COLOC 共定位分析脚本 - MR 分析优化路线 B
# 目的：确认 eQTL 和 GWAS 共享同一因果变异
# ================================================================================

# 包安装与加载
install_and_load_packages <- function() {
  packages <- c(
    "coloc",
    "dplyr",
    "readr",
    "ggplot2",
    "data.table"
  )
  
  message("正在检查并安装所需包...")
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("正在安装:", pkg))
      tryCatch({
        if (pkg %in% c("coloc")) {
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
          }
          BiocManager::install(pkg, update = FALSE)
        } else {
          install.packages(pkg, repos = "https://cloud.r-project.org/")
        }
      }, error = function(e) {
        message(paste("安装失败:", pkg, "-", e$message))
      })
    }
  }
  
  invisible(lapply(packages, library, character.only = TRUE))
  message("所有包加载完成")
}

# 从 MR 结果中提取显著基因
extract_significant_genes <- function(mr_results_file, pval_threshold = 0.05) {
  message("\n=== 提取显著基因 ===")
  
  results <- read.csv(mr_results_file, stringsAsFactors = FALSE)
  
  # 筛选成功的基因
  successful_genes <- results %>%
    filter(status == "SUCCESS" | status == "HETEROGENEITY") %>%
    filter(!is.na(discovery_pval))
  
  # P 值显著基因
  pval_sig_genes <- successful_genes %>%
    filter(discovery_pval < pval_threshold) %>%
    pull(gene) %>%
    unique()
  
  message(paste("P 值显著基因:", length(pval_sig_genes)))
  if (length(pval_sig_genes) > 0) {
    message("  ", paste(pval_sig_genes, collapse = ", "))
  }
  
  return(list(
    pval_significant = pval_sig_genes,
    all_successful = unique(successful_genes$gene)
  ))
}

# 准备 eQTL 数据用于 COLOC
prepare_eqtl_data <- function(exposure_file) {
  message("  准备 eQTL 数据...")
  
  if (!file.exists(exposure_file)) {
    message(paste("    文件不存在:", exposure_file))
    return(NULL)
  }
  
  eqtl <- fread(exposure_file)
  
  # 提取必要列
  required_cols <- c("SNP", "BP", "BETA", "SE", "PVAL", "EAF")
  missing_cols <- setdiff(required_cols, colnames(eqtl))
  
  if (length(missing_cols) > 0) {
    message(paste("    缺少列:", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  # 转换为 coloc 格式
  dataset1 <- list(
    beta = eqtl$BETA,
    varbeta = eqtl$SE^2,
    position = eqtl$BP,
    MAF = eqtl$EAF,
    type = "quant",
    snp = eqtl$SNP
  )
  
  message(paste("    eQTL SNP 数量:", length(eqtl$SNP)))
  return(dataset1)
}

# 准备 GWAS 数据用于 COLOC
prepare_gwas_data <- function(outcome_file) {
  message("  准备 GWAS 数据...")
  
  if (!file.exists(outcome_file)) {
    message(paste("    文件不存在:", outcome_file))
    return(NULL)
  }
  
  gwas <- fread(outcome_file)
  
  # 提取必要列
  required_cols <- c("SNP", "BP", "BETA", "SE", "PVAL", "EAF")
  missing_cols <- setdiff(required_cols, colnames(gwas))
  
  if (length(missing_cols) > 0) {
    # 尝试使用替代列名
    if ("A1" %in% colnames(gwas)) {
      gwas$EAF <- gwas$A1
    }
    missing_cols <- setdiff(required_cols, colnames(gwas))
  }
  
  if (length(missing_cols) > 0) {
    message(paste("    缺少列:", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  # 转换为 coloc 格式
  dataset2 <- list(
    beta = gwas$BETA,
    varbeta = gwas$SE^2,
    position = gwas$BP,
    MAF = if("EAF" %in% colnames(gwas)) gwas$EAF else rep(0.5, nrow(gwas)),
    type = "cc",
    snp = gwas$SNP
  )
  
  message(paste("    GWAS SNP 数量:", length(gwas$SNP)))
  return(dataset2)
}

# 运行 COLOC 分析
run_coloc <- function(dataset1, dataset2, threshold = 0.8) {
  message("  运行 COLOC 分析...")
  
  tryCatch({
    # 运行 ABF 共定位
    result <- coloc.abf(dataset1 = dataset1, dataset2 = dataset2)
    
    # 提取后验概率
    pp <- result$summary
    
    # 判断是否共定位
    coloc_success <- !is.na(pp["PP.H4"]) && pp["PP.H4"] >= threshold
    
    message(paste("    PP.H0 (无关联):", format(pp["PP.H0"], digits = 3)))
    message(paste("    PP.H1 (仅 eQTL):", format(pp["PP.H1"], digits = 3)))
    message(paste("    PP.H2 (仅 GWAS):", format(pp["PP.H2"], digits = 3)))
    message(paste("    PP.H3 (不同因果):", format(pp["PP.H3"], digits = 3)))
    message(paste("    PP.H4 (共享因果):", format(pp["PP.H4"], digits = 3)))
    message(paste("    共定位成功 (PP.H4 >=", threshold, "):", coloc_success))
    
    return(list(
      success = coloc_success,
      pp = pp,
      result = result
    ))
    
  }, error = function(e) {
    message(paste("    COLOC 分析失败:", e$message))
    return(list(
      success = FALSE,
      pp = NULL,
      result = NULL,
      error = e$message
    ))
  })
}

# 可视化 COLOC 结果
plot_coloc <- function(coloc_result, gene_name, output_file) {
  message(paste("  生成共定位图:", gene_name))
  
  tryCatch({
    # 创建数据框
    result_df <- data.frame(
      Hypothesis = c("H0", "H1", "H2", "H3", "H4"),
      Description = c("无关联", "仅 eQTL", "仅 GWAS", 
                     "不同因果变异", "共享因果变异"),
      PP = coloc_result$pp
    )
    
    # 绘制柱状图
    p <- ggplot(result_df, aes(x = Hypothesis, y = PP, fill = Hypothesis)) +
      geom_bar(stat = "identity", width = 0.6) +
      geom_text(aes(label = sprintf("%.3f", PP)), vjust = -0.5, size = 5) +
      scale_fill_manual(values = c("gray", "blue", "green", "orange", "red")) +
      labs(
        title = paste("COLOC 共定位分析 -", gene_name),
        subtitle = paste("PP.H4 =", format(coloc_result$pp["PP.H4"], digits = 3)),
        x = "假设",
        y = "后验概率"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "none"
      )
    
    ggsave(output_file, p, width = 8, height = 6, dpi = 300)
    message(paste("    图已保存:", output_file))
    
  }, error = function(e) {
    message(paste("    绘图失败:", e$message))
  })
}

# 批量运行 COLOC 分析
run_batch_coloc <- function(gene_list, exposure_dir, outcome_dir, 
                             output_dir, threshold = 0.8) {
  message("\n=== 批量 COLOC 共定位分析 ===")
  
  results <- list()
  summary_table <- data.frame(
    gene = character(),
    PP_H0 = numeric(),
    PP_H1 = numeric(),
    PP_H2 = numeric(),
    PP_H3 = numeric(),
    PP_H4 = numeric(),
    coloc_success = logical(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(gene_list)) {
    gene <- gene_list[i]
    message(paste("\n[", i, "/", length(gene_list), "] 分析基因:", gene, sep = ""))
    
    # 查找暴露和结局文件
    exposure_file <- file.path(exposure_dir, paste0(gene, "_exposure.csv"))
    outcome_file <- file.path(outcome_dir, paste0(gene, "_outcome.csv"))
    
    # 检查文件是否存在
    if (!file.exists(exposure_file)) {
      message(paste("  暴露文件不存在，跳过"))
      next
    }
    
    if (!file.exists(outcome_file)) {
      message(paste("  结局文件不存在，跳过"))
      next
    }
    
    # 准备数据
    dataset1 <- prepare_eqtl_data(exposure_file)
    if (is.null(dataset1)) {
      message("  eQTL 数据准备失败，跳过")
      next
    }
    
    dataset2 <- prepare_gwas_data(outcome_file)
    if (is.null(dataset2)) {
      message("  GWAS 数据准备失败，跳过")
      next
    }
    
    # 运行 COLOC
    coloc_result <- run_coloc(dataset1, dataset2, threshold)
    
    # 保存结果
    if (!is.null(coloc_result$pp)) {
      results[[gene]] <- coloc_result
      
      # 添加到摘要表
      summary_table <- rbind(summary_table, data.frame(
        gene = gene,
        PP_H0 = coloc_result$pp["PP.H0"],
        PP_H1 = coloc_result$pp["PP.H1"],
        PP_H2 = coloc_result$pp["PP.H2"],
        PP_H3 = coloc_result$pp["PP.H3"],
        PP_H4 = coloc_result$pp["PP.H4"],
        coloc_success = coloc_result$success
      ))
      
      # 生成可视化
      plot_file <- file.path(output_dir, paste0(gene, "_coloc.png"))
      plot_coloc(coloc_result, gene, plot_file)
    }
  }
  
  return(list(
    results = results,
    summary = summary_table
  ))
}

# 保存 COLOC 结果
save_coloc_results <- function(coloc_results, output_dir) {
  message("\n=== 保存 COLOC 结果 ===")
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 保存摘要表
  summary_file <- file.path(output_dir, "coloc_summary.csv")
  write.csv(coloc_results$summary, summary_file, row.names = FALSE)
  message(paste("摘要表:", summary_file))
  
  # 统计共定位成功的基因
  success_genes <- coloc_results$summary %>%
    filter(coloc_success == TRUE) %>%
    arrange(desc(PP_H4))
  
  if (nrow(success_genes) > 0) {
    success_file <- file.path(output_dir, "coloc_success_genes.csv")
    write.csv(success_genes, success_file, row.names = FALSE)
    message(paste("共定位成功基因:", nrow(success_genes)))
    message(paste("  文件:", success_file))
    
    message("\n共定位成功的基因列表:")
    for (i in 1:nrow(success_genes)) {
      message(paste("  ", i, ". ", success_genes$gene[i], 
                    " (PP.H4 = ", format(success_genes$PP_H4[i], digits = 3), ")", sep = ""))
    }
  } else {
    message("无共定位成功的基因")
  }
  
  # 生成报告
  report_file <- file.path(output_dir, "coloc_report.md")
  generate_coloc_report(coloc_results, report_file)
  message(paste("报告:", report_file))
  
  message(paste("\n所有结果已保存至:", output_dir))
}

# 生成 COLOC 报告
generate_coloc_report <- function(coloc_results, output_file) {
  message("生成 COLOC 报告...")
  
  report <- c(
    "# COLOC 共定位分析报告",
    "",
    "## 分析概述",
    paste("分析日期:", Sys.time()),
    paste("分析基因数:", nrow(coloc_results$summary)),
    paste("共定位阈值：PP.H4 >= 0.8"),
    "",
    "## 共定位成功基因",
    ""
  )
  
  success_genes <- coloc_results$summary %>%
    filter(coloc_success == TRUE) %>%
    arrange(desc(PP_H4))
  
  if (nrow(success_genes) > 0) {
    message(paste("  共定位成功基因数:", nrow(success_genes)))
    
    for (i in 1:nrow(success_genes)) {
      gene <- success_genes$gene[i]
      report <- c(report, paste0(
        "### ", gene, "  \n",
        "- **PP.H0** (无关联): ", format(success_genes$PP_H0[i], digits = 3), "  \n",
        "- **PP.H1** (仅 eQTL): ", format(success_genes$PP_H1[i], digits = 3), "  \n",
        "- **PP.H2** (仅 GWAS): ", format(success_genes$PP_H2[i], digits = 3), "  \n",
        "- **PP.H3** (不同因果): ", format(success_genes$PP_H3[i], digits = 3), "  \n",
        "- **PP.H4** (共享因果): ", format(success_genes$PP_H4[i], digits = 3), " ✅  \n"
      ))
    }
  } else {
    report <- c(report, "无共定位成功的基因 (PP.H4 >= 0.8)")
  }
  
  report <- c(report, "", "## 所有基因结果", "")
  report <- c(report, "| Gene | PP.H0 | PP.H1 | PP.H2 | PP.H3 | PP.H4 | Success |")
  report <- c(report, "|------|-------|-------|-------|-------|-------|---------|")
  
  for (i in 1:nrow(coloc_results$summary)) {
    row <- coloc_results$summary[i, ]
    report <- c(report, paste0(
      "| ", row$gene,
      " | ", format(row$PP_H0, digits = 3),
      " | ", format(row$PP_H1, digits = 3),
      " | ", format(row$PP_H2, digits = 3),
      " | ", format(row$PP_H3, digits = 3),
      " | ", format(row$PP_H4, digits = 3),
      " | ", ifelse(row$coloc_success, "✅", "❌"),
      " |"
    ))
  }
  
  writeLines(report, output_file)
}

# 主函数
main <- function() {
  message(rep("=", 60))
  message("COLOC 共定位分析 - MR 分析优化路线 B")
  message(rep("=", 60))
  
  # 1. 加载包
  install_and_load_packages()
  
  # 2. 提取显著基因
  mr_results_file <- "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/MR_results_main_optimized.csv"
  
  if (!file.exists(mr_results_file)) {
    stop("找不到 MR 结果文件")
  }
  
  genes <- extract_significant_genes(mr_results_file, pval_threshold = 0.05)
  
  # 使用 P 值显著基因
  target_gene_list <- genes$pval_significant
  
  if (length(target_gene_list) == 0) {
    message("无 P 值显著基因，使用所有成功基因")
    target_gene_list <- head(genes$all_successful, 20)  # 限制为前 20 个
  }
  
  message(paste("\n用于 COLOC 分析的基因数:", length(target_gene_list)))
  message("基因列表:", paste(target_gene_list, collapse = ", "))
  
  # 3. 设置输入输出目录
  exposure_dir <- "D:/下载/MR_batch_results/exposure_optimized"
  outcome_dir <- "D:/下载/MR_batch_results/outcome"
  output_dir <- "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/coloc_analysis"
  
  # 4. 批量运行 COLOC
  coloc_results <- run_batch_coloc(
    gene_list = target_gene_list,
    exposure_dir = exposure_dir,
    outcome_dir = outcome_dir,
    output_dir = output_dir,
    threshold = 0.8
  )
  
  # 5. 保存结果
  save_coloc_results(coloc_results, output_dir)
  
  message("\n", rep("=", 60))
  message("COLOC 共定位分析完成!")
  message(paste("结果目录:", output_dir))
  message(rep("=", 60))
}

# 运行主函数
if (!interactive()) {
  main()
} else {
  message("交互模式下运行，请手动调用 main()")
}
