#!/usr/bin/env Rscript
# ================================================================================
# 孟德尔随机化 (Mendelian Randomization) 分析脚本
# 版本：1.0.0
# R 版本要求：>= 4.1.0
# ================================================================================

# 设置随机数种子以确保可重复性
set.seed(2024)

# ================================================================================
# 第一部分：包安装与加载
# ================================================================================

install_and_load_packages <- function() {
  required_packages <- c(
    "TwoSampleMR",      # MR 分析核心包
    "ieugwasr",         # IEU GWAS 数据库接口
    "dplyr",            # 数据处理
    "ggplot2",          # 可视化
    "data.table",       # 高效数据操作
    "readxl",           # Excel 文件读取
    "writexl",          # Excel 文件写入
    "patchwork",        # 图形组合
    "forestplot",       # 森林图
    "tidyr"             # 数据整理
  )
  
  message("正在检查并安装所需包...")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("正在安装:", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
    message(paste("已加载:", pkg))
  }
  
  # 加载所有包
  invisible(lapply(required_packages, library, character.only = TRUE))
  
  # 初始化 TwoSampleMR
  tryCatch({
    TwoSampleMR::init()
    message("TwoSampleMR 初始化成功")
  }, error = function(e) {
    warning(paste("TwoSampleMR 初始化失败:", e$message))
  })
}

# ================================================================================
# 第二部分：输入验证与工具函数
# ================================================================================

validate_input <- function(data, required_cols, data_name = "数据") {
  stopifnot(
    paste(data_name, "不能为空:", deparse(substitute(data))) = !is.null(data),
    paste(data_name, "行数必须大于 0") = nrow(data) > 0
  )
  
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste(data_name, "缺少必需列:", paste(missing_cols, collapse = ", ")))
  }
  
  message(paste(data_name, "验证通过"))
  return(TRUE)
}

check_na_inf <- function(data, cols = NULL) {
  if (is.null(cols)) {
    cols <- colnames(data)
  }
  
  issues <- list()
  for (col in cols) {
    if (col %in% colnames(data)) {
      n_na <- sum(is.na(data[[col]]))
      n_inf <- sum(is.infinite(data[[col]]))
      
      if (n_na > 0 || n_inf > 0) {
        issues[[col]] <- list(na = n_na, inf = n_inf)
      }
    }
  }
  
  if (length(issues) > 0) {
    warning_msg <- paste0(
      "发现 NA/Inf 值:\n",
      paste(sapply(names(issues), function(x) {
        sprintf("  %s: NA=%d, Inf=%d", x, issues[[x]]$na, issues[[x]]$inf)
      }), collapse = "\n")
    )
    warning(warning_msg)
  }
  
  return(invisible(issues))
}

# ================================================================================
# 第三部分：数据加载与预处理
# ================================================================================

load_exposure_data <- function(file_path, format = "auto") {
  message("正在加载暴露因素数据...")
  
  if (!file.exists(file_path)) {
    stop(paste("暴露因素数据文件不存在:", file_path))
  }
  
  if (format == "auto") {
    ext <- tools::file_ext(file_path)
    format <- switch(ext,
      "xlsx" = "excel",
      "xls" = "excel",
      "csv" = "csv",
      "txt" = "csv",
      stop(paste("不支持的文件格式:", ext))
    )
  }
  
  if (format == "excel") {
    data <- readxl::read_excel(file_path)
  } else if (format == "csv") {
    data <- data.table::fread(file_path, stringsAsFactors = FALSE)
  } else {
    stop("未知的数据格式")
  }
  
  data <- as.data.frame(data)
  
  required_cols <- c("SNP", "beta.exposure", "se.exposure", 
                     "effect_allele.exposure", "other_allele.exposure",
                     "eaf.exposure", "pval.exposure")
  
  validate_input(data, required_cols, "暴露因素数据")
  check_na_inf(data, required_cols)
  
  message(paste("成功加载", nrow(data), "个 SNP"))
  return(data)
}

load_outcome_data <- function(file_path, format = "auto") {
  message("正在加载结局因素数据...")
  
  if (!file.exists(file_path)) {
    stop(paste("结局因素数据文件不存在:", file_path))
  }
  
  if (format == "auto") {
    ext <- tools::file_ext(file_path)
    format <- switch(ext,
      "xlsx" = "excel",
      "xls" = "excel",
      "csv" = "csv",
      "txt" = "csv",
      stop(paste("不支持的文件格式:", ext))
    )
  }
  
  if (format == "excel") {
    data <- readxl::read_excel(file_path)
  } else if (format == "csv") {
    data <- data.table::fread(file_path, stringsAsFactors = FALSE)
  } else {
    stop("未知的数据格式")
  }
  
  data <- as.data.frame(data)
  
  required_cols <- c("SNP", "beta.outcome", "se.outcome", 
                     "effect_allele.outcome", "other_allele.outcome",
                     "eaf.outcome", "pval.outcome")
  
  validate_input(data, required_cols, "结局因素数据")
  check_na_inf(data, required_cols)
  
  message(paste("成功加载", nrow(data), "个 SNP"))
  return(data)
}

harmonise_data <- function(exposure_dat, outcome_dat) {
  message("正在进行数据协调...")
  
  harmonised <- TwoSampleMR::harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_dat,
    action = 2,
    verbose = TRUE
  )
  
  message(paste("协调后剩余", nrow(harmonised), "个 SNP"))
  
  if (nrow(harmonised) == 0) {
    stop("数据协调后没有剩余 SNP，请检查数据格式")
  }
  
  return(harmonised)
}

# ================================================================================
# 第四部分：工具变量选择
# ================================================================================

select_instruments <- function(data, pval_threshold = 5e-8, 
                               clump_kb = 10000, 
                               clump_r2 = 0.001,
                               exposure = "exposure") {
  message("正在进行工具变量选择...")
  
  message(paste("  P 值阈值:", pval_threshold))
  message(paste("  连锁不平衡 clump 距离:", clump_kb, "kb"))
  message(paste("  连锁不平衡 r2 阈值:", clump_r2))
  
  instruments <- data %>%
    dplyr::filter(pval.exposure < pval_threshold)
  
  message(paste("  P 值筛选后剩余:", nrow(instruments), "个 SNP"))
  
  if (nrow(instruments) == 0) {
    warning("没有 SNP 满足 P 值阈值，尝试放宽阈值到 1e-5")
    instruments <- data %>%
      dplyr::filter(pval.exposure < 1e-5)
    message(paste("  放宽阈值后剩余:", nrow(instruments), "个 SNP"))
  }
  
  if (nrow(instruments) > 0) {
    tryCatch({
      instruments_clumped <- TwoSampleMR::clump_data(
        instruments,
        clump_kb = clump_kb,
        clump_r2 = clump_r2,
        clump_p1 = pval_threshold,
        clump_p2 = 1
      )
      
      message(paste("  LD clump 后剩余:", nrow(instruments_clumped), "个独立 SNP"))
      instruments <- instruments_clumped
    }, error = function(e) {
      warning(paste("LD clump 失败:", e$message, "使用原始筛选结果"))
    })
  }
  
  if (nrow(instruments) == 0) {
    stop("工具变量选择后没有剩余 SNP")
  }
  
  return(instruments)
}

# ================================================================================
# 第五部分：MR 分析核心方法
# ================================================================================

run_mr_analysis <- function(harmonised_dat) {
  message("正在运行 MR 分析...")
  
  mr_results <- TwoSampleMR::mr(harmonised_dat)
  
  message("MR 分析完成")
  message(paste("分析方法数:", nrow(mr_results)))
  
  return(mr_results)
}

run_sensitivity_analysis <- function(harmonised_dat) {
  message("正在运行敏感性分析...")
  
  sensitivity_results <- list()
  
  tryCatch({
    sensitivity_results$heterogeneity <- TwoSampleMR::mr_heterogeneity(harmonised_dat)
    message("异质性检验完成")
  }, error = function(e) {
    warning(paste("异质性检验失败:", e$message))
  })
  
  tryCatch({
    sensitivity_results$pleiotropy <- TwoSampleMR::mr_pleiotropy_test(harmonised_dat)
    message("水平多效性检验完成")
  }, error = function(e) {
    warning(paste("水平多效性检验失败:", e$message))
  })
  
  tryCatch({
    sensitivity_results$leave_one_out <- TwoSampleMR::mr_leaveoneout(harmonised_dat)
    message("留一法分析完成")
  }, error = function(e) {
    warning(paste("留一法分析失败:", e$message))
  })
  
  return(sensitivity_results)
}

# ================================================================================
# 第六部分：可视化
# ================================================================================

create_mr_plots <- function(harmonised_dat, mr_results, sensitivity_results, 
                            output_dir = "./mr_results") {
  message("正在创建可视化图表...")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  plots <- list()
  
  tryCatch({
    p1 <- TwoSampleMR::mr_scatter_plot(harmonised_dat, mr_results)
    plots$scatter <- p1
    ggsave(file.path(output_dir, "mr_scatter_plot.png"), p1[[1]], 
           width = 8, height = 6, dpi = 300)
    message("  散点图已保存")
  }, error = function(e) {
    warning(paste("散点图创建失败:", e$message))
  })
  
  tryCatch({
    p2 <- TwoSampleMR::mr_forest_plot(harmonised_dat)
    plots$forest <- p2
    ggsave(file.path(output_dir, "mr_forest_plot.png"), p2, 
           width = 10, height = 8, dpi = 300)
    message("  森林图已保存")
  }, error = function(e) {
    warning(paste("森林图创建失败:", e$message))
  })
  
  tryCatch({
    p3 <- TwoSampleMR::mr_funnel_plot(harmonised_dat)
    plots$funnel <- p3
    ggsave(file.path(output_dir, "mr_funnel_plot.png"), p3, 
           width = 8, height = 6, dpi = 300)
    message("  漏斗图已保存")
  }, error = function(e) {
    warning(paste("漏斗图创建失败:", e$message))
  })
  
  if (!is.null(sensitivity_results$leave_one_out)) {
    tryCatch({
      p4 <- TwoSampleMR::mr_leaveoneout_plot(harmonised_dat)
      plots$leave_one_out <- p4
      ggsave(file.path(output_dir, "mr_leave_one_out_plot.png"), p4[[1]], 
             width = 10, height = 8, dpi = 300)
      message("  留一法图已保存")
    }, error = function(e) {
      warning(paste("留一法图创建失败:", e$message))
    })
  }
  
  return(plots)
}

create_manhattan_plot <- function(data, output_file = "./mr_results/manhattan_plot.png") {
  message("正在创建曼哈顿图...")
  
  if (!"chr" %in% colnames(data) || !"pos" %in% colnames(data)) {
    warning("数据缺少染色体位置信息，跳过曼哈顿图")
    return(NULL)
  }
  
  p <- ggplot(data, aes(x = pos, y = -log10(pval.outcome), color = factor(chr))) +
    geom_point(alpha = 0.6) +
    facet_wrap(~chr, scales = "free_x") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.5, "lines")
    ) +
    labs(
      x = "染色体位置",
      y = "-log10(P 值)",
      title = "GWAS 关联分析曼哈顿图"
    )
  
  ggsave(output_file, p, width = 12, height = 6, dpi = 300)
  message("曼哈顿图已保存")
  
  return(p)
}

# ================================================================================
# 第七部分：结果导出
# ================================================================================

export_results <- function(mr_results, sensitivity_results, harmonised_dat, 
                           output_dir = "./mr_results") {
  message("正在导出结果...")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  mr_results_file <- file.path(output_dir, "mr_results.csv")
  write.csv(mr_results, mr_results_file, row.names = FALSE)
  message(paste("  MR 结果已保存:", mr_results_file))
  
  if (!is.null(sensitivity_results$heterogeneity)) {
    het_file <- file.path(output_dir, "heterogeneity_results.csv")
    write.csv(sensitivity_results$heterogeneity, het_file, row.names = FALSE)
    message(paste("  异质性结果已保存:", het_file))
  }
  
  if (!is.null(sensitivity_results$pleiotropy)) {
    pleio_file <- file.path(output_dir, "pleiotropy_results.csv")
    write.csv(sensitivity_results$pleiotropy, pleio_file, row.names = FALSE)
    message(paste("  多效性结果已保存:", pleio_file))
  }
  
  if (!is.null(sensitivity_results$leave_one_out)) {
    loo_file <- file.path(output_dir, "leave_one_out_results.csv")
    write.csv(sensitivity_results$leave_one_out, loo_file, row.names = FALSE)
    message(paste("  留一法结果已保存:", loo_file))
  }
  
  harmonised_file <- file.path(output_dir, "harmonised_data.csv")
  write.csv(harmonised_dat, harmonised_file, row.names = FALSE)
  message(paste("  协调后数据已保存:", harmonised_file))
  
  summary_file <- file.path(output_dir, "mr_summary.txt")
  sink(summary_file)
  cat("========================================\n")
  cat("孟德尔随机化分析总结报告\n")
  cat(paste("生成时间:", Sys.time(), "\n"))
  cat("========================================\n\n")
  
  cat("【MR 分析结果】\n")
  print(mr_results)
  cat("\n")
  
  if (!is.null(sensitivity_results$heterogeneity)) {
    cat("【异质性检验】\n")
    print(sensitivity_results$heterogeneity)
    cat("\n")
  }
  
  if (!is.null(sensitivity_results$pleiotropy)) {
    cat("【水平多效性检验】\n")
    print(sensitivity_results$pleiotropy)
    cat("\n")
  }
  
  sink()
  message(paste("  总结报告已保存:", summary_file))
  
  message("所有结果导出完成")
}

# ================================================================================
# 第八部分：主分析流程
# ================================================================================

run_mr_pipeline <- function(exposure_file = NULL,
                            outcome_file = NULL,
                            exposure_data = NULL,
                            outcome_data = NULL,
                            pval_threshold = 5e-8,
                            clump_kb = 10000,
                            clump_r2 = 0.001,
                            output_dir = "./mr_results",
                            save_plots = TRUE) {
  message("========================================")
  message("孟德尔随机化分析开始")
  message(paste("开始时间:", Sys.time()))
  message("========================================")
  
  start_time <- Sys.time()
  
  tryCatch({
    install_and_load_packages()
  }, error = function(e) {
    stop(paste("包加载失败:", e$message))
  })
  
  if (!is.null(exposure_file)) {
    exposure_dat <- load_exposure_data(exposure_file)
  } else if (!is.null(exposure_data)) {
    exposure_dat <- exposure_data
    message("使用提供的暴露因素数据框")
  } else {
    stop("必须提供暴露因素文件或数据框")
  }
  
  if (!is.null(outcome_file)) {
    outcome_dat <- load_outcome_data(outcome_file)
  } else if (!is.null(outcome_data)) {
    outcome_dat <- outcome_data
    message("使用提供的结局因素数据框")
  } else {
    stop("必须提供结局因素文件或数据框")
  }
  
  instruments <- select_instruments(
    exposure_dat,
    pval_threshold = pval_threshold,
    clump_kb = clump_kb,
    clump_r2 = clump_r2
  )
  
  outcome_snps <- outcome_dat %>%
    dplyr::filter(SNP %in% instruments$SNP)
  
  if (nrow(outcome_snps) == 0) {
    stop("暴露因素和结局因素没有共同的 SNP")
  }
  
  message(paste("共同 SNP 数量:", nrow(outcome_snps)))
  
  harmonised_dat <- harmonise_data(instruments, outcome_snps)
  
  mr_results <- run_mr_analysis(harmonised_dat)
  
  sensitivity_results <- run_sensitivity_analysis(harmonised_dat)
  
  if (save_plots) {
    create_mr_plots(harmonised_dat, mr_results, sensitivity_results, output_dir)
  }
  
  export_results(mr_results, sensitivity_results, harmonised_dat, output_dir)
  
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  message("========================================")
  message("孟德尔随机化分析完成")
  message(paste("结束时间:", end_time))
  message(paste("总耗时:", round(duration, 2), "分钟"))
  message(paste("结果保存目录:", output_dir))
  message("========================================")
  
  return(list(
    mr_results = mr_results,
    sensitivity_results = sensitivity_results,
    harmonised_data = harmonised_dat,
    instruments = instruments
  ))
}

# ================================================================================
# 第九部分：示例用法
# ================================================================================

example_usage <- function() {
  results <- run_mr_pipeline(
    exposure_file = "exposure_data.xlsx",
    outcome_file = "outcome_data.xlsx",
    pval_threshold = 5e-8,
    clump_kb = 10000,
    clump_r2 = 0.001,
    output_dir = "./mr_results",
    save_plots = TRUE
  )
  
  print(results$mr_results)
}

fetch_gwas_data <- function(exposure_id, outcome_id) {
  message("正在从 IEU GWAS 数据库获取数据...")
  
  tryCatch({
    exposure_dat <- TwoSampleMR::extract_instruments(outcome_ids = exposure_id)
    message(paste("暴露因素数据获取成功，SNP 数量:", nrow(exposure_dat)))
  }, error = function(e) {
    stop(paste("暴露因素数据获取失败:", e$message))
  })
  
  tryCatch({
    outcome_dat <- TwoSampleMR::extract_outcome_data(
      snps = exposure_dat$SNP,
      outcomes = outcome_id
    )
    message(paste("结局因素数据获取成功，SNP 数量:", nrow(outcome_dat)))
  }, error = function(e) {
    stop(paste("结局因素数据获取失败:", e$message))
  })
  
  return(list(
    exposure = exposure_dat,
    outcome = outcome_dat
  ))
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) >= 2) {
    exposure_file <- args[1]
    outcome_file <- args[2]
    output_dir <- if (length(args) >= 3) args[3] else "./mr_results"
    
    results <- run_mr_pipeline(
      exposure_file = exposure_file,
      outcome_file = outcome_file,
      output_dir = output_dir
    )
  } else {
    message("用法：Rscript mr_analysis.R <exposure_file> <outcome_file> [output_dir]")
    message("示例：Rscript mr_analysis.R exposure.csv outcome.csv ./results")
  }
}

message("MR 分析脚本加载完成")
