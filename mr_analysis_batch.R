#!/usr/bin/env Rscript
# ================================================================================
# 孟德尔随机化 (Mendelian Randomization) 批量分析脚本
# 版本：2.0.0 (适配 eQTL + FinnGen 数据格式)
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
    "tidyr"             # 数据整理
  )
  
  message("正在检查并安装所需包...")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("正在安装:", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
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
  
  message("所有包加载完成")
}

# ================================================================================
# 第二部分：输入验证与工具函数
# ================================================================================

validate_input <- function(data, required_cols, data_name = "数据") {
  if (is.null(data)) {
    stop(paste(data_name, "不能为空"))
  }
  
  if (nrow(data) == 0) {
    stop(paste(data_name, "行数必须大于 0"))
  }
  
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
# 第三部分：数据加载与格式化（适配你的数据格式）
# ================================================================================

load_exposure_data <- function(file_path) {
  message("正在加载暴露因素数据...")
  
  if (!file.exists(file_path)) {
    stop(paste("暴露因素数据文件不存在:", file_path))
  }
  
  # 读取 CSV 文件
  data <- data.table::fread(file_path, stringsAsFactors = FALSE)
  data <- as.data.frame(data)
  
  # 标准化列名（适配 eQTL 数据格式）
  # 必需列：SNP, CHR, BP, EFFECT_ALLELE, OTHER_ALLELE, BETA, SE, PVAL, EAF, GENE
  required_cols <- c("SNP", "CHR", "BP", "EFFECT_ALLELE", "OTHER_ALLELE", 
                     "BETA", "SE", "PVAL", "EAF", "GENE")
  
  # 检查必需列
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    # 尝试查找替代列名
    alt_mapping <- list(
      SNP = c("SNPID", "RSID", "RSNUMBER"),
      CHR = c("CHR.EXPOSURE", "CHROMOSOME"),
      BP = c("POS.EXPOSURE", "POSITION", "BPOSITION"),
      EFFECT_ALLELE = c("EA", "A1", "ALLELE1", "EFFECT_ALLELE.EXPOSURE"),
      OTHER_ALLELE = c("OA", "A2", "ALLELE0", "OTHER_ALLELE.EXPOSURE"),
      BETA = c("BETA.EXPOSURE", "B", "Z"),
      SE = c("SE.EXPOSURE"),
      PVAL = c("PVAL.EXPOSURE", "P", "PVALUE", "P-VALUE"),
      EAF = c("EAF.EXPOSURE", "FRQ", "FREQ", "MAF"),
      GENE = c("GENENAME", "SYMBOL")
    )
    
    for (col in missing_cols) {
      if (col %in% names(alt_mapping)) {
        for (alt_col in alt_mapping[[col]]) {
          if (alt_col %in% colnames(data)) {
            data[[col]] <- data[[alt_col]]
            message(paste("  列名映射:", alt_col, "->", col))
            break
          }
        }
      }
    }
    
    # 重新检查
    missing_cols <- setdiff(required_cols, colnames(data))
  }
  
  if (length(missing_cols) > 0) {
    stop(paste("暴露数据缺少必需列:", paste(missing_cols, collapse = ", ")))
  }
  
  # 添加 samplesize 列（Steiger 检验需要）
  # eQTLGen 全血样本量 N = 31,684
  if (!"samplesize" %in% colnames(data)) {
    data$samplesize <- 31684
  }
  
  validate_input(data, required_cols, "暴露因素数据")
  check_na_inf(data, c("BETA", "SE", "PVAL", "EAF"))
  
  message(paste("成功加载", nrow(data), "个 SNP"))
  return(data)
}

load_outcome_data <- function(file_path) {
  message("正在加载结局因素数据...")
  
  if (!file.exists(file_path)) {
    stop(paste("结局因素数据文件不存在:", file_path))
  }
  
  # 读取 CSV 文件
  data <- data.table::fread(file_path, stringsAsFactors = FALSE)
  data <- as.data.frame(data)
  
  # 标准化列名（适配 FinnGen 数据格式）
  # 必需列：SNP, BETA, SE, EFFECT_ALLELE, OTHER_ALLELE, PVAL
  required_cols <- c("SNP", "BETA", "SE", "EFFECT_ALLELE", "OTHER_ALLELE", "PVAL")
  
  # 尝试列名映射
  col_mapping <- list(
    SNP = c("rsid", "RSID", "SNPID"),
    BETA = c("beta", "BETA", "b"),
    SE = c("se", "SE", "std_err"),
    EFFECT_ALLELE = c("effect_allele", "EFFECT_ALLELE", "A1", "a1"),
    OTHER_ALLELE = c("other_allele", "OTHER_ALLELE", "A2", "a2"),
    PVAL = c("pval", "PVAL", "p_value", "P", "p")
  )
  
  # 转换列名（小写匹配）
  colnames_lower <- tolower(colnames(data))
  for (target_col in names(col_mapping)) {
    for (source_col in col_mapping[[target_col]]) {
      if (source_col %in% colnames_lower) {
        idx <- which(colnames_lower == source_col)
        colnames(data)[idx] <- target_col
        message(paste("  列名映射:", source_col, "->", target_col))
        break
      }
    }
  }
  
  validate_input(data, required_cols, "结局因素数据")
  check_na_inf(data, c("BETA", "SE", "PVAL"))
  
  message(paste("成功加载", nrow(data), "个 SNP"))
  return(data)
}

format_exposure_for_twosamplemr <- function(data) {
  message("正在格式化暴露数据为 TwoSampleMR 格式...")
  
  # 使用 TwoSampleMR 的 format_data 函数
  exp_dat <- TwoSampleMR::format_data(
    data,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "EFFECT_ALLELE",
    other_allele_col = "OTHER_ALLELE",
    eaf_col = "EAF",
    pval_col = "PVAL",
    samplesize_col = "samplesize",
    phenotype_col = "GENE"
  )
  
  message(paste("格式化完成，列名:", paste(colnames(exp_dat), collapse = ", ")))
  return(exp_dat)
}

format_outcome_for_twosamplemr <- function(data) {
  message("正在格式化结局数据为 TwoSampleMR 格式...")
  
  # 使用 TwoSampleMR 的 format_data 函数
  out_dat <- TwoSampleMR::format_data(
    data,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "EFFECT_ALLELE",
    other_allele_col = "OTHER_ALLELE",
    pval_col = "PVAL"
  )
  
  message(paste("格式化完成，列名:", paste(colnames(out_dat), collapse = ", ")))
  return(out_dat)
}

# ================================================================================
# 第四部分：数据协调与工具变量选择
# ================================================================================

harmonise_data <- function(exposure_dat, outcome_dat) {
  message("正在进行数据协调 (harmonisation)...")
  
  harmonised <- TwoSampleMR::harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_dat,
    action = 2,  # 排除 palindromic SNPs
    verbose = TRUE
  )
  
  message(paste("协调后剩余", nrow(harmonised), "个 SNP"))
  
  if (is.null(harmonised) || nrow(harmonised) == 0) {
    stop("数据协调后没有剩余 SNP，请检查数据格式")
  }
  
  # 检查是否所有 SNP 都被标记为排除
  if (all(harmonised$mr_keep == FALSE)) {
    stop("所有 SNP 因 palindromic 被排除，无法继续分析")
  }
  
  # 只保留 mr_keep 为 TRUE 的 SNP
  harmonised <- harmonised[harmonised$mr_keep == TRUE, ]
  
  if (nrow(harmonised) == 0) {
    stop("无有效 SNP 用于分析")
  }
  
  # 添加 samplesize 列（Steiger 检验需要）
  harmonised$samplesize.exposure <- 31684
  harmonised$samplesize.outcome <- 452000
  
  message(paste("有效 SNP 数量:", nrow(harmonised)))
  return(harmonised)
}

# ================================================================================
# 第五部分：MR 分析核心方法
# ================================================================================

run_mr_analysis <- function(harmonised_dat) {
  message("正在运行 MR 分析...")
  
  # 根据 SNP 数量选择分析方法
  if (nrow(harmonised_dat) == 1) {
    # 单 SNP: Wald ratio
    message("  单 SNP，使用 Wald ratio 方法")
    mr_results <- TwoSampleMR::mr(harmonised_dat, method_list = "mr_wald_ratio")
  } else {
    # 多 SNP: 多种方法
    message("  多 SNP，使用多种 MR 方法")
    mr_results <- TwoSampleMR::mr(harmonised_dat, method_list = c(
      "mr_ivw",              # 逆方差加权
      "mr_ivw_mre",          # 逆方差加权（随机效应）
      "mr_egger_regression", # MR-Egger 回归
      "mr_weighted_median",  # 加权中位数
      "mr_simple_mode",      # 简单模式
      "mr_weighted_mode"     # 加权模式
    ))
  }
  
  if (is.null(mr_results) || nrow(mr_results) == 0) {
    stop("MR 分析未返回结果")
  }
  
  message(paste("MR 分析完成，方法数:", nrow(mr_results)))
  return(mr_results)
}

run_sensitivity_analysis <- function(harmonised_dat) {
  message("正在运行敏感性分析...")
  
  sensitivity_results <- list()
  
  # 1. 异质性检验（仅多 SNP）
  if (nrow(harmonised_dat) >= 2) {
    tryCatch({
      sensitivity_results$heterogeneity <- TwoSampleMR::mr_heterogeneity(harmonised_dat)
      message("  异质性检验完成")
    }, error = function(e) {
      warning(paste("异质性检验失败:", e$message))
    })
  }
  
  # 2. 水平多效性检验（仅 MR-Egger）
  if (nrow(harmonised_dat) >= 2) {
    tryCatch({
      sensitivity_results$pleiotropy <- TwoSampleMR::mr_pleiotropy_test(harmonised_dat)
      message("  水平多效性检验完成")
    }, error = function(e) {
      warning(paste("水平多效性检验失败:", e$message))
    })
  }
  
  # 3. 留一法分析（仅 nSNP >= 3）
  if (nrow(harmonised_dat) >= 3) {
    tryCatch({
      sensitivity_results$leave_one_out <- TwoSampleMR::mr_leaveoneout(harmonised_dat)
      message("  留一法分析完成")
    }, error = function(e) {
      warning(paste("留一法分析失败:", e$message))
    })
  }
  
  # 4. Steiger 方向检验
  tryCatch({
    sensitivity_results$steiger <- TwoSampleMR::directionality_test(harmonised_dat)
    message("  Steiger 方向检验完成")
  }, error = function(e) {
    warning(paste("Steiger 方向检验失败:", e$message))
    sensitivity_results$steiger <- NULL
  })
  
  return(sensitivity_results)
}

# ================================================================================
# 第六部分：计算 F 统计量与可靠性指标
# ================================================================================

calculate_f_statistic <- function(harmonised_dat) {
  # F = beta^2 / se^2
  harmonised_dat$F <- harmonised_dat$beta.exposure^2 / harmonised_dat$se.exposure^2
  
  f_stats <- list(
    mean_f = mean(harmonised_dat$F, na.rm = TRUE),
    min_f = min(harmonised_dat$F, na.rm = TRUE),
    max_f = max(harmonised_dat$F, na.rm = TRUE),
    weak_iv = any(harmonised_dat$F < 10, na.rm = TRUE)
  )
  
  if (f_stats$weak_iv) {
    warning("存在弱工具变量 (F < 10)")
  }
  
  return(f_stats)
}

# ================================================================================
# 第七部分：可视化
# ================================================================================

create_mr_plots <- function(harmonised_dat, mr_results, sensitivity_results, 
                            gene_name, output_dir = "./mr_results") {
  message("正在创建可视化图表...")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }
  
  # 1. 森林图
  tryCatch({
    p1 <- TwoSampleMR::mr_forest_plot(mr_results)
    ggsave(file.path(plots_dir, paste0(gene_name, "_forest.png")), 
           p1[[1]], width = 8, height = 6, dpi = 300)
    message("  森林图已保存")
  }, error = function(e) {
    warning(paste("森林图创建失败:", e$message))
  })
  
  # 2. 漏斗图
  tryCatch({
    p2 <- TwoSampleMR::mr_funnel_plot(mr_results)
    ggsave(file.path(plots_dir, paste0(gene_name, "_funnel.png")), 
           p2[[1]], width = 8, height = 6, dpi = 300)
    message("  漏斗图已保存")
  }, error = function(e) {
    warning(paste("漏斗图创建失败:", e$message))
  })
  
  # 3. 散点图
  tryCatch({
    p3 <- TwoSampleMR::mr_scatter_plot(mr_results, harmonised_dat)
    ggsave(file.path(plots_dir, paste0(gene_name, "_scatter.png")), 
           p3[[1]], width = 8, height = 6, dpi = 300)
    message("  散点图已保存")
  }, error = function(e) {
    warning(paste("散点图创建失败:", e$message))
  })
  
  # 4. 留一法图（仅 nSNP >= 3）
  if (!is.null(sensitivity_results$leave_one_out)) {
    tryCatch({
      p4 <- TwoSampleMR::mr_leaveoneout_plot(sensitivity_results$leave_one_out)
      ggsave(file.path(plots_dir, paste0(gene_name, "_loo.png")), 
             p4[[1]], width = 8, height = 6, dpi = 300)
      message("  留一法图已保存")
    }, error = function(e) {
      warning(paste("留一法图创建失败:", e$message))
    })
  }
}

# ================================================================================
# 第八部分：结果整理与导出
# ================================================================================

extract_main_result <- function(mr_results, harmonised_dat, f_stats, 
                                sensitivity_results, gene_name) {
  # 提取 IVW 结果（主要方法）
  ivw_methods <- c("mr_ivw", "mr_ivw_mre", "mr_wald_ratio")
  ivw_res <- mr_results[mr_results$method %in% ivw_methods, ]
  
  if (nrow(ivw_res) == 0) {
    # 如果 IVW 结果为空，取第一行
    ivw_res <- mr_results[1, , drop = FALSE]
  }
  
  # 计算 OR 和 95%CI
  or_val <- exp(ivw_res$b)
  or_lci <- exp(ivw_res$b - 1.96 * ivw_res$se)
  or_uci <- exp(ivw_res$b + 1.96 * ivw_res$se)
  
  # Steiger 方向
  steiger_correct <- TRUE
  steiger_p <- NA
  if (!is.null(sensitivity_results$steiger)) {
    if (!is.na(sensitivity_results$steiger$correct_causal_direction)) {
      steiger_correct <- sensitivity_results$steiger$correct_causal_direction
    }
    if (!is.na(sensitivity_results$steiger$steiger_pval)) {
      steiger_p <- sensitivity_results$steiger$steiger_pval
    }
  }
  
  # 异质性检验
  q_p <- NA
  if (!is.null(sensitivity_results$heterogeneity)) {
    if (nrow(sensitivity_results$heterogeneity) > 0) {
      q_p <- sensitivity_results$heterogeneity$Q_pval[1]
    }
  }
  
  # 多效性检验
  egger_p <- NA
  if (!is.null(sensitivity_results$pleiotropy)) {
    egger_p <- sensitivity_results$pleiotropy$pval
  }
  
  # 状态标志
  status_flags <- c()
  if (f_stats$weak_iv) status_flags <- c(status_flags, "WEAK_IV")
  if (!steiger_correct) status_flags <- c(status_flags, "REVERSE")
  status <- ifelse(length(status_flags) == 0, "SUCCESS", paste(status_flags, collapse = "|"))
  
  # 构建结果行
  result_row <- data.frame(
    gene = gene_name,
    method = ivw_res$method,
    b = ivw_res$b,
    se = ivw_res$se,
    pval = ivw_res$pval,
    OR = or_val,
    OR_lci = or_lci,
    OR_uci = or_uci,
    nsnp = ivw_res$nsnp,
    F_mean = f_stats$mean_f,
    F_min = f_stats$min_f,
    Steiger_dir = ifelse(steiger_correct, "CORRECT", "REVERSE"),
    Steiger_p = steiger_p,
    Q_p = q_p,
    Egger_intercept_p = egger_p,
    status = status,
    stringsAsFactors = FALSE
  )
  
  return(result_row)
}

export_results <- function(results_list, log_list, output_dir = "./mr_results") {
  message("正在导出结果...")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. 导出主结果
  if (length(results_list) > 0) {
    final_results <- do.call(rbind, results_list)
    
    # FDR 校正（Benjamini-Hochberg）
    final_results$fdr_qval <- p.adjust(final_results$pval, method = "BH")
    final_results$fdr_sig <- final_results$fdr_qval < 0.05
    
    # 格式化 OR 显示
    final_results$OR_95CI <- sprintf("%.3f (%.3f-%.3f)", 
                                     final_results$OR, 
                                     final_results$OR_lci, 
                                     final_results$OR_uci)
    
    # 选择输出列
    output_cols <- c("gene", "method", "b", "se", "pval", "OR_95CI", 
                     "nsnp", "F_mean", "Steiger_dir", "Steiger_p", 
                     "Q_p", "Egger_intercept_p", "status", "fdr_qval", "fdr_sig")
    
    results_file <- file.path(output_dir, "MR_results_main.csv")
    write.csv(final_results[, output_cols], results_file, row.names = FALSE)
    message(paste("  主结果已保存:", results_file))
    
    # 显示显著结果
    sig_results <- final_results[final_results$pval < 0.05, ]
    if (nrow(sig_results) > 0) {
      cat(sprintf("\n显著结果 (p < 0.05): %d 个基因\n", nrow(sig_results)))
      print(sig_results[, c("gene", "pval", "OR_95CI", "status")])
    }
    
    # 显示 FDR 校正后显著结果
    fdr_sig_results <- final_results[final_results$fdr_sig, ]
    if (nrow(fdr_sig_results) > 0) {
      cat(sprintf("\nFDR 校正后显著 (q < 0.05): %d 个基因\n", nrow(fdr_sig_results)))
      print(fdr_sig_results[, c("gene", "pval", "fdr_qval", "OR_95CI")])
    } else {
      cat("\nFDR 校正后无显著基因 (q < 0.05)\n")
    }
  }
  
  # 2. 导出分析日志
  if (length(log_list) > 0) {
    log_df <- do.call(rbind, lapply(log_list, as.data.frame))
    log_file <- file.path(output_dir, "MR_analysis_log.csv")
    write.csv(log_df, log_file, row.names = FALSE)
    message(paste("  分析日志已保存:", log_file))
  }
}

# ================================================================================
# 第九部分：单基因分析流程
# ================================================================================

analyze_single_gene <- function(gene_name, exposure_file, outcome_file, 
                                output_dir = "./mr_results") {
  message(paste("\n========== 分析基因:", gene_name, "=========="))
  
  # 初始化日志
  log_entry <- list(
    gene = gene_name,
    status = "INIT",
    message = ""
  )
  
  tryCatch({
    # 1. 加载数据
    exposure_dat <- load_exposure_data(exposure_file)
    outcome_dat <- load_outcome_data(outcome_file)
    
    # 2. 格式化数据
    exp_dat <- format_exposure_for_twosamplemr(exposure_dat)
    out_dat <- format_outcome_for_twosamplemr(outcome_dat)
    
    # 3. 数据协调
    harmonised_dat <- harmonise_data(exp_dat, out_dat)
    
    # 4. 计算 F 统计量
    f_stats <- calculate_f_statistic(harmonised_dat)
    message(paste("  F 统计量：mean =", round(f_stats$mean_f, 2), 
                  ", min =", round(f_stats$min_f, 2)))
    
    # 5. MR 分析
    mr_results <- run_mr_analysis(harmonised_dat)
    
    # 6. 敏感性分析
    sensitivity_results <- run_sensitivity_analysis(harmonised_dat)
    
    # 7. 提取主要结果
    result_row <- extract_main_result(
      mr_results, harmonised_dat, f_stats, 
      sensitivity_results, gene_name
    )
    
    # 8. 可视化（仅显著基因）
    if (result_row$pval < 0.05) {
      message(paste("  显著结果 (p =", format(result_row$pval, scientific = TRUE), 
                    ")，生成可视化..."))
      create_mr_plots(harmonised_dat, mr_results, sensitivity_results, 
                      gene_name, output_dir)
    }
    
    # 9. 保存详细结果
    detail_file <- file.path(output_dir, "details", paste0(gene_name, "_detail.csv"))
    if (!dir.exists(dirname(detail_file))) {
      dir.create(dirname(detail_file), recursive = TRUE)
    }
    write.csv(mr_results, detail_file, row.names = FALSE)
    
    log_entry$status <- result_row$status
    log_entry$message <- sprintf("nSNP=%d, p=%.3e", nrow(harmonised_dat), result_row$pval)
    
    message(paste("  分析完成:", result_row$status))
    
    return(list(
      result = result_row,
      log = log_entry,
      mr_results = mr_results,
      sensitivity_results = sensitivity_results,
      harmonised_data = harmonised_dat
    ))
    
  }, error = function(e) {
    log_entry$status <- "ERROR"
    log_entry$message <- conditionMessage(e)
    message(paste("  错误:", conditionMessage(e)))
    
    return(list(
      result = NULL,
      log = log_entry,
      mr_results = NULL,
      sensitivity_results = NULL,
      harmonised_data = NULL
    ))
  })
}

# ================================================================================
# 第十部分：批量分析流程
# ================================================================================

run_batch_analysis <- function(exposure_dir, outcome_dir, gene_list_file,
                               output_dir = "./mr_results") {
  message("========================================")
  message("孟德尔随机化批量分析开始")
  message(paste("开始时间:", Sys.time()))
  message("========================================")
  
  start_time <- Sys.time()
  
  # 1. 加载包
  tryCatch({
    install_and_load_packages()
  }, error = function(e) {
    stop(paste("包加载失败:", e$message))
  })
  
  # 2. 读取基因清单
  if (!file.exists(gene_list_file)) {
    stop(paste("基因清单文件不存在:", gene_list_file))
  }
  
  gene_list <- read.table(gene_list_file, header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE)
  
  # 3. 过滤有数据的基因
  genes_to_analyze <- gene_list$GENE[gene_list$N_SNPS > 0]
  message(paste("需要分析的基因数:", length(genes_to_analyze)))
  
  # 4. 批量分析循环
  results_list <- list()
  log_list <- list()
  
  for (i in seq_along(genes_to_analyze)) {
    gene <- genes_to_analyze[i]
    message(sprintf("\n[%d/%d] 分析基因：%s", i, length(genes_to_analyze), gene))
    
    exposure_file <- file.path(exposure_dir, paste0(gene, ".exposure.csv"))
    outcome_file <- file.path(outcome_dir, paste0(gene, ".outcome.csv"))
    
    # 检查文件是否存在
    if (!file.exists(exposure_file)) {
      message("  跳过：暴露文件不存在")
      log_list[[gene]] <- list(
        gene = gene,
        status = "NO_EXPOSURE_FILE",
        message = "暴露文件不存在"
      )
      next
    }
    
    if (!file.exists(outcome_file)) {
      message("  跳过：结局文件不存在")
      log_list[[gene]] <- list(
        gene = gene,
        status = "NO_OUTCOME_FILE",
        message = "结局文件不存在"
      )
      next
    }
    
    # 分析基因
    analysis_result <- analyze_single_gene(
      gene, exposure_file, outcome_file, output_dir
    )
    
    if (!is.null(analysis_result$result)) {
      results_list[[gene]] <- analysis_result$result
    }
    
    log_list[[gene]] <- analysis_result$log
  }
  
  # 5. 导出结果
  export_results(results_list, log_list, output_dir)
  
  # 6. 分析摘要
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  message("\n========================================")
  message("孟德尔随机化批量分析完成")
  message(paste("结束时间:", end_time))
  message(paste("总耗时:", round(duration, 2), "分钟"))
  message(paste("结果保存目录:", output_dir))
  message("========================================")
  
  # 显示状态统计
  if (length(log_list) > 0) {
    log_df <- do.call(rbind, lapply(log_list, as.data.frame))
    cat("\n分析状态统计:\n")
    print(table(log_df$status))
  }
  
  return(list(
    results = if (length(results_list) > 0) do.call(rbind, results_list) else NULL,
    logs = if (length(log_list) > 0) do.call(rbind, lapply(log_list, as.data.frame)) else NULL
  ))
}

# ================================================================================
# 第十一部分：示例用法
# ================================================================================

example_usage <- function() {
  # 设置路径
  work_dir <- "D:/下载/MR_batch_results"
  exposure_dir <- file.path(work_dir, "exposure")
  outcome_dir <- file.path(work_dir, "outcome")
  gene_list_file <- file.path(work_dir, "gene_list.txt")
  output_dir <- file.path(work_dir, format(Sys.Date(), "%Y%m%d"))
  
  # 运行批量分析
  results <- run_batch_analysis(
    exposure_dir = exposure_dir,
    outcome_dir = outcome_dir,
    gene_list_file = gene_list_file,
    output_dir = output_dir
  )
  
  return(results)
}

# 命令行支持
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) >= 3) {
    exposure_dir <- args[1]
    outcome_dir <- args[2]
    gene_list_file <- args[3]
    output_dir <- if (length(args) >= 4) args[4] else "./mr_results"
    
    results <- run_batch_analysis(
      exposure_dir = exposure_dir,
      outcome_dir = outcome_dir,
      gene_list_file = gene_list_file,
      output_dir = output_dir
    )
  } else if (length(args) == 2) {
    # 单基因分析模式
    exposure_file <- args[1]
    outcome_file <- args[2]
    output_dir <- "./mr_results"
    
    gene_name <- tools::file_path_sans_ext(basename(exposure_file))
    gene_name <- sub("\\.exposure$", "", gene_name)
    
    result <- analyze_single_gene(gene_name, exposure_file, outcome_file, output_dir)
    print(result$result)
  } else {
    message("用法 1 (批量分析): Rscript mr_analysis_batch.R <exposure_dir> <outcome_dir> <gene_list.txt> [output_dir]")
    message("用法 2 (单基因): Rscript mr_analysis_batch.R <exposure.csv> <outcome.csv>")
    message("示例：Rscript mr_analysis_batch.R ./exposure ./outcome gene_list.txt ./results")
  }
}

message("MR 批量分析脚本加载完成")
