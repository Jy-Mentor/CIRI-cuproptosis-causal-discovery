#!/usr/bin/env Rscript
# ================================================================================
# 敏感性分析增强脚本 - MR 分析优化路线 B
# 包含：MR-PRESSO (异常值检测), 径向 MR 分析
# ================================================================================

# 包安装与加载
install_and_load_packages <- function() {
  packages <- c(
    "MRPRESSO",
    "RadialMR",
    "dplyr",
    "readr",
    "ggplot2",
    "data.table",
    "patchwork"
  )
  
  message("正在检查并安装所需包...")
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("正在安装:", pkg))
      tryCatch({
        if (pkg %in% c("MRPRESSO")) {
          # MRPRESSO 需要从 GitHub 安装
          if (!requireNamespace("remotes", quietly = TRUE)) {
            install.packages("remotes")
          }
          remotes::install_github("rondolab/MR-PRESSO", upgrade = "never")
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

# 加载协调后的数据
load_harmonised_data <- function(gene_name, details_dir) {
  message(paste("加载", gene_name, "的协调数据..."))
  
  detail_file <- file.path(details_dir, paste0(gene_name, "_detail.csv"))
  
  if (!file.exists(detail_file)) {
    message(paste("  详细结果文件不存在:", detail_file))
    return(NULL)
  }
  
  data <- fread(detail_file)
  
  # 检查必要列
  required_cols <- c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome", "SNP")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    message(paste("  缺少必要列:", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  message(paste("  SNP 数量:", nrow(data)))
  return(data)
}

# MR-PRESSO 异常值检测
run_mr_presso <- function(harmonised_data, gene_name, n_distribution = 1000) {
  message(paste("\n=== MR-PRESSO 分析:", gene_name, "==="))
  
  if (nrow(harmonised_data) < 3) {
    message("  SNP 数量不足 (需要>=3)，跳过")
    return(NULL)
  }
  
  tryCatch({
    # 准备数据
    presso_data <- data.frame(
      beta.exposure = harmonised_data$beta.exposure,
      beta.outcome = harmonised_data$beta.outcome,
      se.exposure = harmonised_data$se.exposure,
      se.outcome = harmonised_data$se.outcome
    )
    
    # 运行 MR-PRESSO
    message("  运行 MR-PRESSO 全局检验...")
    presso_result <- mr_presso(
      BetaOutcome = presso_data$beta.outcome,
      BetaExposure = presso_data$beta.exposure,
      SdOutcome = presso_data$se.outcome,
      SdExposure = presso_data$se.exposure,
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      NbDistribution = n_distribution,
      SignifThreshold = 0.05
    )
    
    # 提取结果
    global_test <- presso_result$`MR-PRESSO test`
    outlier_test <- presso_result$`Outlier Test`
    distortion_test <- presso_result$`Distortion test`
    
    message(paste("  全局检验 P 值:", format(global_test$`P-value`, digits = 3)))
    message(paste("  检测到异常值 SNP 数:", if(!is.null(outlier_test)) nrow(outlier_test) else 0))
    
    if (!is.null(distortion_test)) {
      message(paste("  畸变检验 P 值:", format(distortion_test$`P-value`, digits = 3)))
    }
    
    return(list(
      global_pval = global_test$`P-value`,
      outliers = outlier_test,
      distortion_pval = if(!is.null(distortion_test)) distortion_test$`P-value` else NA,
      has_outliers = !is.null(outlier_test) && nrow(outlier_test) > 0,
      result = presso_result
    ))
    
  }, error = function(e) {
    message(paste("  MR-PRESSO 分析失败:", e$message))
    return(list(
      global_pval = NA,
      outliers = NULL,
      distortion_pval = NA,
      has_outliers = FALSE,
      error = e$message
    ))
  })
}

# 径向 MR 分析
run_radial_mr <- function(harmonised_data, gene_name) {
  message(paste("\n=== 径向 MR 分析:", gene_name, "==="))
  
  if (nrow(harmonised_data) < 3) {
    message("  SNP 数量不足 (需要>=3)，跳过")
    return(NULL)
  }
  
  tryCatch({
    # 准备径向 MR 数据
    radial_data <- RadialMR::ivw_radial(
      betas = harmonised_data$beta.exposure,
      sebetas = harmonised_data$se.exposure,
      betaout = harmonised_data$beta.outcome,
      seout = harmonised_data$se.outcome
    )
    
    # 提取结果
    message(paste("  径向 IVW 估计:", format(radial_data$beta, digits = 4)))
    message(paste("  标准误:", format(radial_data$se, digits = 4)))
    message(paste("  P 值:", format(radial_data$pval, digits = 3)))
    
    # 检测异常值
    if (!is.null(radial_data$outliers)) {
      message(paste("  检测到异常值 SNP 数:", length(radial_data$outliers)))
    }
    
    return(list(
      beta = radial_data$beta,
      se = radial_data$se,
      pval = radial_data$pval,
      outliers = radial_data$outliers,
      result = radial_data
    ))
    
  }, error = function(e) {
    message(paste("  径向 MR 分析失败:", e$message))
    return(list(
      beta = NA,
      se = NA,
      pval = NA,
      outliers = NULL,
      error = e$message
    ))
  })
}

# 可视化敏感性分析结果
plot_sensitivity_analysis <- function(gene_name, harmonised_data, 
                                       presso_result, radial_result, output_dir) {
  message("  生成敏感性分析可视化...")
  
  plots <- list()
  
  # 1. 径向森林图
  tryCatch({
    if (!is.null(radial_result$result)) {
      p_radial <- ggplot2::ggplot() +
        ggplot2::geom_point(
          ggplot2::aes(x = 1/harmonised_data$se.exposure^2, 
                      y = harmonised_data$beta.outcome/harmonised_data$beta.exposure),
          size = 3, alpha = 0.6
        ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(x = 1/harmonised_data$se.exposure^2,
                      ymin = (harmonised_data$beta.outcome - 1.96*harmonised_data$se.outcome)/harmonised_data$beta.exposure,
                      ymax = (harmonised_data$beta.outcome + 1.96*harmonised_data$se.outcome)/harmonised_data$beta.exposure),
          width = 0
        ) +
        ggplot2::geom_hline(
          yintercept = radial_result$beta,
          color = "red", linetype = "dashed", linewidth = 1
        ) +
        ggplot2::labs(
          title = paste("Radial MR Forest Plot -", gene_name),
          x = "Precision (1/SE²)",
          y = "Wald Ratio"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")
        )
      
      plots$radial <- p_radial
    }
  }, error = function(e) {
    message(paste("    径向森林图失败:", e$message))
  })
  
  # 2. 漏斗图
  tryCatch({
    p_funnel <- ggplot2::ggplot(harmonised_data, 
                                ggplot2::aes(x = beta.exposure/se.exposure, 
                                            y = beta.outcome/se.outcome)) +
      ggplot2::geom_point(size = 3, alpha = 0.6) +
      ggplot2::geom_abline(
        intercept = 0, slope = radial_result$beta,
        color = "red", linetype = "dashed"
      ) +
      ggplot2::labs(
        title = paste("Funnel Plot -", gene_name),
        x = "Exposure Z-score",
        y = "Outcome Z-score"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")
      )
    
    plots$funnel <- p_funnel
  }, error = function(e) {
    message(paste("    漏斗图失败:", e$message))
  })
  
  # 保存图形
  if (length(plots) > 0) {
    for (plot_name in names(plots)) {
      output_file <- file.path(output_dir, paste0(gene_name, "_", plot_name, ".png"))
      ggsave(output_file, plots[[plot_name]], width = 8, height = 6, dpi = 300)
      message(paste("    已保存:", output_file))
    }
  }
}

# 批量运行敏感性分析
run_batch_sensitivity <- function(gene_list, details_dir, output_dir) {
  message("\n=== 批量敏感性分析增强 ===")
  
  results <- list()
  summary_table <- data.frame(
    gene = character(),
    presso_global_pval = numeric(),
    presso_outliers = integer(),
    presso_distortion_pval = numeric(),
    radial_beta = numeric(),
    radial_se = numeric(),
    radial_pval = numeric(),
    radial_outliers = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(gene_list)) {
    gene <- gene_list[i]
    message(paste("\n[", i, "/", length(gene_list), "] 分析基因:", gene, sep = ""))
    
    # 加载数据
    harmonised_data <- load_harmonised_data(gene, details_dir)
    if (is.null(harmonised_data)) {
      message("  无法加载数据，跳过")
      next
    }
    
    # MR-PRESSO 分析
    presso_result <- run_mr_presso(harmonised_data, gene, n_distribution = 1000)
    
    # 径向 MR 分析
    radial_result <- run_radial_mr(harmonised_data, gene)
    
    # 可视化
    plot_sensitivity_analysis(gene, harmonised_data, presso_result, radial_result, output_dir)
    
    # 保存结果
    results[[gene]] <- list(
      presso = presso_result,
      radial = radial_result
    )
    
    # 添加到摘要表
    summary_table <- rbind(summary_table, data.frame(
      gene = gene,
      presso_global_pval = if(!is.null(presso_result)) presso_result$global_pval else NA,
      presso_outliers = if(!is.null(presso_result) && !is.null(presso_result$outliers)) nrow(presso_result$outliers) else 0,
      presso_distortion_pval = if(!is.null(presso_result)) presso_result$distortion_pval else NA,
      radial_beta = if(!is.null(radial_result)) radial_result$beta else NA,
      radial_se = if(!is.null(radial_result)) radial_result$se else NA,
      radial_pval = if(!is.null(radial_result)) radial_result$pval else NA,
      radial_outliers = if(!is.null(radial_result) && !is.null(radial_result$outliers)) length(radial_result$outliers) else 0
    ))
  }
  
  return(list(
    results = results,
    summary = summary_table
  ))
}

# 保存敏感性分析结果
save_sensitivity_results <- function(sensitivity_results, output_dir) {
  message("\n=== 保存敏感性分析结果 ===")
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 保存摘要表
  summary_file <- file.path(output_dir, "sensitivity_summary.csv")
  write.csv(sensitivity_results$summary, summary_file, row.names = FALSE)
  message(paste("摘要表:", summary_file))
  
  # 统计显著结果
  message("\n敏感性分析结果汇总:")
  
  # MR-PRESSO 全局检验显著的基因
  presso_sig <- sensitivity_results$summary %>%
    filter(!is.na(presso_global_pval)) %>%
    filter(presso_global_pval < 0.05) %>%
    arrange(presso_global_pval)
  
  if (nrow(presso_sig) > 0) {
    message(paste("\nMR-PRESSO 全局检验显著 (P < 0.05):", nrow(presso_sig)))
    for (i in 1:nrow(presso_sig)) {
      message(paste("  ", i, ". ", presso_sig$gene[i], 
                    " (P=", format(presso_sig$presso_global_pval[i], digits = 3), 
                    ", 异常值=", presso_sig$presso_outliers[i], ")", sep = ""))
    }
  }
  
  # 径向 MR 异常值检测
  radial_outliers <- sensitivity_results$summary %>%
    filter(radial_outliers > 0) %>%
    arrange(desc(radial_outliers))
  
  if (nrow(radial_outliers) > 0) {
    message(paste("\n径向 MR 检测到异常值:", nrow(radial_outliers)))
    for (i in 1:nrow(radial_outliers)) {
      message(paste("  ", i, ". ", radial_outliers$gene[i], 
                    " (异常值数=", radial_outliers$radial_outliers[i], ")", sep = ""))
    }
  }
  
  # 生成报告
  report_file <- file.path(output_dir, "sensitivity_report.md")
  generate_sensitivity_report(sensitivity_results, report_file)
  message(paste("\n报告:", report_file))
  
  message(paste("\n所有结果已保存至:", output_dir))
}

# 生成敏感性分析报告
generate_sensitivity_report <- function(sensitivity_results, output_file) {
  message("生成敏感性分析报告...")
  
  report <- c(
    "# 敏感性分析增强报告",
    "",
    "## 分析概述",
    paste("分析日期:", Sys.time()),
    paste("分析基因数:", nrow(sensitivity_results$summary)),
    "",
    "## MR-PRESSO 全局检验",
    "",
    "检测水平多效性 (P < 0.05 表示存在多效性)",
    ""
  )
  
  presso_sig <- sensitivity_results$summary %>%
    filter(!is.na(presso_global_pval)) %>%
    filter(presso_global_pval < 0.05) %>%
    arrange(presso_global_pval)
  
  if (nrow(presso_sig) > 0) {
    for (i in 1:nrow(presso_sig)) {
      gene <- presso_sig$gene[i]
      report <- c(report, paste0(
        "### ", gene, "  \n",
        "- **全局检验 P 值**: ", format(presso_sig$presso_global_pval[i], digits = 3), " ✅  \n",
        "- **异常值 SNP 数**: ", presso_sig$presso_outliers[i], "  \n",
        "- **畸变检验 P 值**: ", format(presso_sig$presso_distortion_pval[i], digits = 3), "  \n"
      ))
    }
  } else {
    report <- c(report, "无显著多效性基因")
  }
  
  report <- c(report, "", "## 径向 MR 异常值检测", "")
  
  radial_outliers <- sensitivity_results$summary %>%
    filter(radial_outliers > 0) %>%
    arrange(desc(radial_outliers))
  
  if (nrow(radial_outliers) > 0) {
    for (i in 1:nrow(radial_outliers)) {
      gene <- radial_outliers$gene[i]
      report <- c(report, paste0(
        "### ", gene, "  \n",
        "- **异常值 SNP 数**: ", radial_outliers$radial_outliers[i], "  \n",
        "- **径向 IVW Beta**: ", format(radial_outliers$radial_beta[i], digits = 4), "  \n",
        "- **径向 IVW P 值**: ", format(radial_outliers$radial_pval[i], digits = 3), "  \n"
      ))
    }
  }
  
  report <- c(report, "", "## 所有基因结果", "")
  report <- c(report, "| Gene | PRESSO Global P | Outliers | Distortion P | Radial Beta | Radial P | Radial Outliers |")
  report <- c(report, "|------|-----------------|----------|--------------|-------------|----------|-----------------|")
  
  for (i in 1:nrow(sensitivity_results$summary)) {
    row <- sensitivity_results$summary[i, ]
    report <- c(report, paste0(
      "| ", row$gene,
      " | ", format(row$presso_global_pval, digits = 3, na.rm = TRUE),
      " | ", row$presso_outliers,
      " | ", format(row$presso_distortion_pval, digits = 3, na.rm = TRUE),
      " | ", format(row$radial_beta, digits = 4, na.rm = TRUE),
      " | ", format(row$radial_pval, digits = 3, na.rm = TRUE),
      " | ", row$radial_outliers,
      " |"
    ))
  }
  
  writeLines(report, output_file)
}

# 主函数
main <- function() {
  message(rep("=", 60))
  message("敏感性分析增强 - MR 分析优化路线 B")
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
    target_gene_list <- head(genes$all_successful, 20)
  }
  
  message(paste("\n用于敏感性分析的基因数:", length(target_gene_list)))
  message("基因列表:", paste(target_gene_list, collapse = ", "))
  
  # 3. 设置目录
  details_dir <- "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/details"
  output_dir <- "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/sensitivity_analysis"
  
  # 4. 批量运行敏感性分析
  sensitivity_results <- run_batch_sensitivity(
    gene_list = target_gene_list,
    details_dir = details_dir,
    output_dir = output_dir
  )
  
  # 5. 保存结果
  save_sensitivity_results(sensitivity_results, output_dir)
  
  message("\n", rep("=", 60))
  message("敏感性分析增强完成!")
  message(paste("结果目录:", output_dir))
  message(rep("=", 60))
}

# 运行主函数
if (!interactive()) {
  main()
} else {
  message("交互模式下运行，请手动调用 main()")
}
