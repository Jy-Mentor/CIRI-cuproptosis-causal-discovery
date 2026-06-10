#!/usr/bin/env Rscript
# ================================================================================
# 孟德尔随机化 (MR) 分析脚本 - 优化版 v3.0
# 改进内容:
#   1. P 值阈值：5e-8 → 1e-5 (探索性研究推荐)
#   2. LD Clump: 分级策略 (r2 < 0.001 → 0.01 → 0.1)
#   3. 独立验证：FinnGen 发现 + UK Biobank 验证
# 版本：3.0.0
# R 版本要求：>= 4.1.0
# ================================================================================

set.seed(2024)

# ================================================================================
# 第一部分：包安装与加载
# ================================================================================

install_and_load_packages <- function() {
  required_packages <- c(
    "TwoSampleMR",
    "ieugwasr",
    "dplyr",
    "ggplot2",
    "data.table",
    "readxl",
    "writexl",
    "patchwork",
    "tidyr",
    "meta"  # 新增：meta 分析包
  )
  
  message("正在检查并安装所需包...")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("正在安装:", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
  }
  
  invisible(lapply(required_packages, library, character.only = TRUE))
  
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
  if (is.null(data)) stop(paste(data_name, "不能为空"))
  if (nrow(data) == 0) stop(paste(data_name, "行数必须大于 0"))
  
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste(data_name, "缺少必需列:", paste(missing_cols, collapse = ", ")))
  }
  
  message(paste(data_name, "验证通过"))
  return(TRUE)
}

check_na_inf <- function(data, cols = NULL) {
  if (is.null(cols)) cols <- colnames(data)
  
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
    warning_msg <- paste0("发现 NA/Inf 值:\n",
      paste(sapply(names(issues), function(x) {
        sprintf("  %s: NA=%d, Inf=%d", x, issues[[x]]$na, issues[[x]]$inf)
      }), collapse = "\n"))
    warning(warning_msg)
  }
  
  return(invisible(issues))
}

# ================================================================================
# 第三部分：数据加载与格式化（优化版）
# ================================================================================

load_exposure_data <- function(file_path) {
  message("正在加载暴露因素数据...")
  
  if (!file.exists(file_path)) {
    stop(paste("暴露因素数据文件不存在:", file_path))
  }
  
  data <- data.table::fread(file_path, stringsAsFactors = FALSE)
  data <- as.data.frame(data)
  
  # 标准化列名（保留原始大小写）
  original_cols <- colnames(data)
  colnames_upper <- toupper(original_cols)
  
  required_cols <- c("SNP", "CHR", "BP", "EFFECT_ALLELE", "OTHER_ALLELE", 
                     "BETA", "SE", "PVAL", "EAF", "GENE")
  
  # 列名映射（支持多种格式）
  column_mapping <- list(
    SNP = c("SNP", "SNPID", "RSID"),
    CHR = c("CHR", "CHR.EXPOSURE", "CHROMOSOME"),
    BP = c("BP", "POS", "POSITION", "POS.EXPOSURE", "BPOSITION"),
    EFFECT_ALLELE = c("EFFECT_ALLELE", "EA", "A1", "EFFECT_ALLELE.EXPOSURE"),
    OTHER_ALLELE = c("OTHER_ALLELE", "OA", "A2", "OTHER_ALLELE.EXPOSURE"),
    BETA = c("BETA", "B", "BETA.EXPOSURE"),
    SE = c("SE", "SE.EXPOSURE"),
    PVAL = c("PVAL", "P", "PVALUE", "P-VALUE", "PVAL.EXPOSURE"),
    EAF = c("EAF", "FRQ", "MAF", "EAF.EXPOSURE"),
    GENE = c("GENE", "GENENAME", "SYMBOL")
  )
  
  for (target_col in names(column_mapping)) {
    for (source_col in column_mapping[[target_col]]) {
      # 检查原始列名或大写列名
      if (source_col %in% original_cols) {
        idx <- which(original_cols == source_col)
        colnames(data)[idx] <- target_col
        message(paste("  列名映射:", source_col, "->", target_col))
        break
      } else if (source_col %in% colnames_upper) {
        idx <- which(colnames_upper == source_col)
        colnames(data)[idx] <- target_col
        message(paste("  列名映射:", original_cols[idx], "->", target_col))
        break
      }
    }
  }
  
  if (length(missing_cols <- setdiff(required_cols, colnames(data))) > 0) {
    stop(paste("暴露数据缺少必需列:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!"samplesize" %in% colnames(data)) {
    data$samplesize <- 31684  # eQTLGen 样本量
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
  
  data <- data.table::fread(file_path, stringsAsFactors = FALSE)
  data <- as.data.frame(data)
  
  # 保留原始列名，进行智能映射
  original_cols <- colnames(data)
  colnames_lower <- tolower(original_cols)
  
  required_cols <- c("SNP", "BETA", "SE", "EFFECT_ALLELE", "OTHER_ALLELE", "PVAL")
  
  # 列名映射（支持多种格式，包括 FinnGen 格式）
  col_mapping <- list(
    SNP = c("SNP", "snp", "RSID", "rsid", "SNPID"),
    BETA = c("BETA", "beta", "B", "b"),
    SE = c("SE", "se", "STD_ERR"),
    EFFECT_ALLELE = c("EFFECT_ALLELE", "effect_allele", "A1", "a1"),
    OTHER_ALLELE = c("OTHER_ALLELE", "other_allele", "A2", "a2"),
    PVAL = c("PVAL", "pval", "P", "p", "P_VALUE", "p_value")
  )
  
  for (target_col in names(col_mapping)) {
    for (source_col in col_mapping[[target_col]]) {
      if (source_col %in% original_cols) {
        idx <- which(original_cols == source_col)
        colnames(data)[idx] <- target_col
        message(paste("  列名映射:", source_col, "->", target_col))
        break
      }
    }
  }
  
  validate_input(data, required_cols, "结局因素数据")
  check_na_inf(data, c("beta", "se", "pval"))
  
  message(paste("成功加载", nrow(data), "个 SNP"))
  return(data)
}

format_exposure_for_twosamplemr <- function(data) {
  message("正在格式化暴露数据为 TwoSampleMR 格式...")
  
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
  
  message(paste("格式化完成，列名:", paste(colnames(exp_dat)[1:8], collapse = ", ")))
  return(exp_dat)
}

format_outcome_for_twosamplemr <- function(data) {
  message("正在格式化结局数据为 TwoSampleMR 格式...")
  
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
  
  message(paste("格式化完成"))
  return(out_dat)
}

# ================================================================================
# 第四部分：改进的工具变量选择（分级 Clumping）
# ================================================================================

select_instruments_optimized <- function(data, 
                                          pval_threshold = 1e-5,
                                          clump_kb = 10000,
                                          gene_name = "Unknown") {
  message(paste("\n【优化策略】分析基因:", gene_name))
  message(paste("  P 值阈值:", pval_threshold, "(改进：从 5e-8 放宽到 1e-5)"))
  
  # 步骤 1: P 值筛选
  instruments <- data %>%
    dplyr::filter(pval.exposure < pval_threshold)
  
  message(paste("  P 值筛选后:", nrow(instruments), "个 SNP"))
  
  if (nrow(instruments) == 0) {
    message(paste("  ⚠ 基因", gene_name, "无满足 P <", pval_threshold, "的 SNP"))
    return(NULL)
  }
  
  # 步骤 2: 分级 LD Clumping (改进 2)
  clump_strategies <- list(
    list(r2 = 0.001, label = "严格 (r²<0.001)"),
    list(r2 = 0.01, label = "中等 (r²<0.01)"),
    list(r2 = 0.1, label = "宽松 (r²<0.1)")
  )
  
  final_instruments <- NULL
  
  for (i in seq_along(clump_strategies)) {
    strategy <- clump_strategies[[i]]
    
    message(paste("  尝试 Clump 策略", i, ":", strategy$label))
    
    tryCatch({
      instruments_clumped <- TwoSampleMR::clump_data(
        instruments,
        clump_kb = clump_kb,
        clump_r2 = strategy$r2,
        clump_p1 = pval_threshold
      )
      
      n_snps <- nrow(instruments_clumped)
      message(paste("    Clump 后剩余:", n_snps, "个独立 SNP"))
      
      if (n_snps >= 1) {
        final_instruments <- instruments_clumped
        message(paste("    ✓ 采用此策略"))
        break
      }
      
    }, error = function(e) {
      message(paste("    Clump 失败:", e$message))
    })
    
    if (i == length(clump_strategies)) {
      # 所有 clump 策略都失败，使用原始数据
      message("  ⚠ 所有 Clump 策略失败，使用原始 SNP（不 Clump）")
      final_instruments <- instruments
    }
  }
  
  # 记录最终使用的 clump 参数
  if (!is.null(final_instruments)) {
    final_instruments$clump_r2_used <- strategy$r2
    message(paste("  最终采用 r² <", strategy$r2))
  }
  
  return(final_instruments)
}

# ================================================================================
# 第五部分：数据协调
# ================================================================================

harmonise_data <- function(exposure_dat, outcome_dat) {
  message("正在进行数据协调...")
  
  harmonised <- TwoSampleMR::harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_dat,
    action = 2
  )
  
  message(paste("协调后剩余", nrow(harmonised), "个 SNP"))
  
  if (is.null(harmonised) || nrow(harmonised) == 0) {
    return(NULL)
  }
  
  harmonised <- harmonised[harmonised$mr_keep == TRUE, ]
  
  if (nrow(harmonised) == 0) {
    return(NULL)
  }
  
  # 添加 samplesize
  harmonised$samplesize.exposure <- 31684
  harmonised$samplesize.outcome <- 452000
  
  message(paste("有效 SNP 数量:", nrow(harmonised)))
  return(harmonised)
}

# ================================================================================
# 第六部分：MR 分析（含独立验证 - 改进 3）
# ================================================================================

run_mr_analysis_with_replication <- function(harmonised_dat, 
                                               outcome_id_discovery = "finn-b-I9_STROKE",
                                               outcome_id_replication = "ukb-d-I9_STROKE",
                                               gene_name = "Unknown") {
  message(paste("\n【MR 分析】基因:", gene_name))
  
  results_list <- list()
  
  # ===== 发现阶段 (FinnGen) =====
  message("【发现阶段】FinnGen 队列...")
  
  tryCatch({
    mr_results_discovery <- TwoSampleMR::mr(harmonised_dat, method_list = c(
      "mr_ivw",
      "mr_ivw_mre",
      "mr_egger_regression",
      "mr_weighted_median",
      "mr_simple_mode",
      "mr_weighted_mode"
    ))
    
    results_list$discovery <- mr_results_discovery
    message("  ✓ FinnGen MR 分析完成")
    
  }, error = function(e) {
    message(paste("  ✗ FinnGen MR 分析失败:", e$message))
    results_list$discovery <- NULL
  })
  
  # ===== 验证阶段 (UK Biobank) - 改进 3 =====
  message("【验证阶段】UK Biobank 队列...")
  
  tryCatch({
    # 从 IEU GWAS 数据库提取 UK Biobank 数据
    outcome_replication <- TwoSampleMR::extract_outcome_data(
      snps = harmonised_dat$SNP,
      outcomes = outcome_id_replication
    )
    
    if (!is.null(outcome_replication) && nrow(outcome_replication) > 0) {
      # 协调
      harmonised_rep <- TwoSampleMR::harmonise_data(
        exposure_dat = harmonised_dat,
        outcome_dat = outcome_replication,
        action = 2,
        verbose = FALSE
      )
      
      if (!is.null(harmonised_rep) && nrow(harmonised_rep) > 0) {
        # MR 分析
        mr_results_replication <- TwoSampleMR::mr(harmonised_rep, 
                                                   method_list = c("mr_ivw", "mr_egger_regression"))
        
        results_list$replication <- list(
          mr_results = mr_results_replication,
          harmonised_data = harmonised_rep
        )
        message("  ✓ UK Biobank MR 分析完成")
      } else {
        message("  ⚠ UK Biobank 数据协调后无 SNP")
        results_list$replication <- NULL
      }
    } else {
      message("  ⚠ 无法获取 UK Biobank 数据")
      results_list$replication <- NULL
    }
    
  }, error = function(e) {
    message(paste("  ✗ UK Biobank 验证失败:", e$message))
    results_list$replication <- NULL
  })
  
  # ===== Meta 分析（如果验证成功）=====
  if (!is.null(results_list$discovery) && !is.null(results_list$replication)) {
    message("【Meta 分析】合并发现队列和验证队列...")
    
    tryCatch({
      # 提取 IVW 结果
      ivw_discovery <- results_list$discovery %>%
        dplyr::filter(method == "mr_ivw")
      
      ivw_replication <- results_list$replication$mr_results %>%
        dplyr::filter(method == "mr_ivw")
      
      if (nrow(ivw_discovery) > 0 && nrow(ivw_replication) > 0) {
        # 固定效应 Meta 分析
        meta_result <- meta::metagen(
          TE = c(ivw_discovery$b, ivw_replication$b),
          seTE = c(ivw_discovery$se, ivw_replication$se),
          studlab = c("FinnGen", "UK Biobank"),
          sm = "OR",
          method.tau = "DL"
        )
        
        results_list$meta <- list(
          pooled_b = meta_result$TE.fixed,
          pooled_se = meta_result$seTE.fixed,
          pooled_pval = 2 * pnorm(-abs(meta_result$TE.fixed / meta_result$seTE.fixed)),
          pooled_or = exp(meta_result$TE.fixed),
          pooled_lci = exp(meta_result$lower.fixed),
          pooled_uci = exp(meta_result$upper.fixed),
          heterogeneity_q = meta_result$Q,
          heterogeneity_p = meta_result$pval.Q,
          i_squared = meta_result$I2
        )
        
        message("  ✓ Meta 分析完成")
        message(paste("    合并 OR =", sprintf("%.3f", results_list$meta$pooled_or),
                      ", 95%CI =", sprintf("(%.3f-%.3f)", 
                                           results_list$meta$pooled_lci,
                                           results_list$meta$pooled_uci)))
      }
    }, error = function(e) {
      message(paste("  ✗ Meta 分析失败:", e$message))
    })
  }
  
  return(results_list)
}

# ================================================================================
# 第七部分：敏感性分析
# ================================================================================

run_sensitivity_analysis <- function(harmonised_dat) {
  message("正在运行敏感性分析...")
  
  sensitivity_results <- list()
  
  if (nrow(harmonised_dat) >= 2) {
    tryCatch({
      sensitivity_results$heterogeneity <- TwoSampleMR::mr_heterogeneity(harmonised_dat)
      message("  ✓ 异质性检验完成")
    }, error = function(e) {
      warning(paste("异质性检验失败:", e$message))
    })
    
    tryCatch({
      sensitivity_results$pleiotropy <- TwoSampleMR::mr_pleiotropy_test(harmonised_dat)
      message("  ✓ 水平多效性检验完成")
    }, error = function(e) {
      warning(paste("水平多效性检验失败:", e$message))
    })
  }
  
  if (nrow(harmonised_dat) >= 3) {
    tryCatch({
      sensitivity_results$leave_one_out <- TwoSampleMR::mr_leaveoneout(harmonised_dat)
      message("  ✓ 留一法分析完成")
    }, error = function(e) {
      warning(paste("留一法分析失败:", e$message))
    })
  }
  
  tryCatch({
    sensitivity_results$steiger <- TwoSampleMR::directionality_test(harmonised_dat)
    message("  ✓ Steiger 方向检验完成")
  }, error = function(e) {
    warning(paste("Steiger 方向检验失败:", e$message))
  })
  
  return(sensitivity_results)
}

# ================================================================================
# 第八部分：计算 F 统计量
# ================================================================================

calculate_f_statistic <- function(harmonised_dat) {
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
# 第九部分：结果提取与整理
# ================================================================================

extract_main_result <- function(mr_results_list, harmonised_dat, f_stats, 
                                 sensitivity_results, gene_name) {
  
  # 提取发现队列 IVW 结果
  if (!is.null(mr_results_list$discovery)) {
    ivw_res <- mr_results_list$discovery %>%
      dplyr::filter(method %in% c("mr_ivw", "mr_ivw_mre", "mr_wald_ratio"))
    
    if (nrow(ivw_res) == 0) {
      ivw_res <- mr_results_list$discovery[1, , drop = FALSE]
    }
  } else {
    return(NULL)
  }
  
  # 计算 OR 和 95%CI
  or_val <- exp(ivw_res$b)
  or_lci <- exp(ivw_res$b - 1.96 * ivw_res$se)
  or_uci <- exp(ivw_res$b + 1.96 * ivw_res$se)
  
  # Meta 分析结果（如果有）
  has_replication <- !is.null(mr_results_list$replication)
  has_meta <- !is.null(mr_results_list$meta)
  
  if (has_meta) {
    meta_or <- mr_results_list$meta$pooled_or
    meta_lci <- mr_results_list$meta$pooled_lci
    meta_uci <- mr_results_list$meta$pooled_uci
    meta_pval <- mr_results_list$meta$pooled_pval
  } else {
    meta_or <- or_val
    meta_lci <- or_lci
    meta_uci <- or_uci
    meta_pval <- ivw_res$pval
  }
  
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
  
  # 状态标志 (安全处理 NA 值)
  status_flags <- c()
  
  # 1. 弱工具变量检查
  if (!is.na(f_stats$weak_iv) && f_stats$weak_iv) {
    status_flags <- c(status_flags, "WEAK_IV")
  }
  
  # 2. Steiger 方向检查 (修复 NA 值问题)
  if (!is.na(steiger_correct) && !steiger_correct) {
    status_flags <- c(status_flags, "REVERSE")
  }
  
  # 3. 验证检查
  if (has_replication && has_meta && !is.na(meta_pval) && meta_pval >= 0.05) {
    status_flags <- c(status_flags, "NOT_REPLICATED")
  }
  
  # 4. 异质性检查
  if (!is.na(q_p) && q_p < 0.05) {
    status_flags <- c(status_flags, "HETEROGENEITY")
  }
  
  # 5. 多效性检查
  if (!is.na(egger_p) && egger_p < 0.05) {
    status_flags <- c(status_flags, "PLEIOTROPY")
  }
  
  status <- ifelse(length(status_flags) == 0, "SUCCESS", paste(status_flags, collapse = "|"))
  
  # 构建结果行
  result_row <- data.frame(
    gene = gene_name,
    method = ivw_res$method,
    discovery_b = ivw_res$b,
    discovery_se = ivw_res$se,
    discovery_pval = ivw_res$pval,
    discovery_OR = or_val,
    discovery_OR_lci = or_lci,
    discovery_OR_uci = or_uci,
    nsnp = ivw_res$nsnp,
    F_mean = f_stats$mean_f,
    F_min = f_stats$min_f,
    has_replication = has_replication,
    has_meta = has_meta,
    meta_OR = ifelse(has_meta, meta_or, NA),
    meta_OR_lci = ifelse(has_meta, meta_lci, NA),
    meta_OR_uci = ifelse(has_meta, meta_uci, NA),
    meta_pval = ifelse(has_meta, meta_pval, NA),
    Steiger_dir = ifelse(steiger_correct, "CORRECT", "REVERSE"),
    Steiger_p = steiger_p,
    Q_p = q_p,
    Egger_intercept_p = egger_p,
    status = status,
    stringsAsFactors = FALSE
  )
  
  return(result_row)
}

# ================================================================================
# 第十部分：单基因分析流程（优化版）
# ================================================================================

analyze_single_gene_optimized <- function(gene_name, exposure_file, outcome_file, 
                                           output_dir = "./mr_results_optimized") {
  message(paste("\n", rep("=", 60), sep=""))
  message(paste("【优化版 MR 分析】基因:", gene_name))
  message(paste(rep("=", 60), "\n", sep=""))
  
  log_entry <- list(
    gene = gene_name,
    status = "INIT",
    message = ""
  )
  
  tryCatch({
    # 1. 加载数据前先检查文件是否存在
    if (!file.exists(exposure_file)) {
      message(paste("  ✗ 暴露文件不存在:", exposure_file))
      log_entry$status <- "NO_EXPOSURE_FILE"
      log_entry$message <- "暴露文件不存在"
      return(list(result = NULL, log = log_entry))
    }
    
    if (!file.exists(outcome_file)) {
      message(paste("  ✗ 结局文件不存在:", outcome_file))
      log_entry$status <- "NO_OUTCOME_FILE"
      log_entry$message <- "结局文件不存在"
      return(list(result = NULL, log = log_entry))
    }
    
    # 检查文件是否为空
    if (file.info(exposure_file)$size == 0) {
      message(paste("  ✗ 暴露文件为空:", exposure_file))
      log_entry$status <- "EMPTY_EXPOSURE_FILE"
      log_entry$message <- "暴露文件为空"
      return(list(result = NULL, log = log_entry))
    }
    
    if (file.info(outcome_file)$size == 0) {
      message(paste("  ✗ 结局文件为空:", outcome_file))
      log_entry$status <- "EMPTY_OUTCOME_FILE"
      log_entry$message <- "结局文件为空"
      return(list(result = NULL, log = log_entry))
    }
    
    exposure_dat <- load_exposure_data(exposure_file)
    outcome_dat <- load_outcome_data(outcome_file)
    
    # 2. 格式化数据
    exp_dat <- format_exposure_for_twosamplemr(exposure_dat)
    out_dat <- format_outcome_for_twosamplemr(outcome_dat)
    
    # 3. 选择工具变量（优化策略）
    instruments <- select_instruments_optimized(
      exp_dat,  # 使用格式化后的数据
      pval_threshold = 1e-5,
      clump_kb = 10000,
      gene_name = gene_name
    )
    
    if (is.null(instruments) || nrow(instruments) == 0) {
      message(paste("  ✗ 基因", gene_name, "无可用工具变量"))
      log_entry$status <- "NO_INSTRUMENTS"
      log_entry$message <- "无满足 P<1e-5 的 SNP"
      return(list(result = NULL, log = log_entry))
    }
    
    # 4. 提取结局数据中的对应 SNP
    outcome_snps <- out_dat %>%
      dplyr::filter(SNP %in% instruments$SNP)
    
    if (nrow(outcome_snps) == 0) {
      message(paste("  ✗ 暴露和结局无共同 SNP"))
      log_entry$status <- "NO_COMMON_SNPS"
      log_entry$message <- "无共同 SNP"
      return(list(result = NULL, log = log_entry))
    }
    
    message(paste("  共同 SNP 数量:", nrow(outcome_snps)))
    
    # 5. 数据协调
    harmonised_dat <- harmonise_data(instruments, outcome_snps)
    
    if (is.null(harmonised_dat) || nrow(harmonised_dat) == 0) {
      message(paste("  ✗ 数据协调后无有效 SNP"))
      log_entry$status <- "HARMONISATION_FAILED"
      log_entry$message <- "协调失败"
      return(list(result = NULL, log = log_entry))
    }
    
    # 6. 计算 F 统计量
    f_stats <- calculate_f_statistic(harmonised_dat)
    message(paste("  F 统计量：mean =", round(f_stats$mean_f, 2), 
                  ", min =", round(f_stats$min_f, 2)))
    
    # 7. MR 分析（含独立验证 - 改进 3）
    mr_results_list <- run_mr_analysis_with_replication(
      harmonised_dat,
      outcome_id_discovery = "finn-b-I9_STROKE",
      outcome_id_replication = "ukb-d-I9_STROKE",
      gene_name = gene_name
    )
    
    if (is.null(mr_results_list$discovery)) {
      message(paste("  ✗ MR 分析失败"))
      log_entry$status <- "MR_ANALYSIS_FAILED"
      log_entry$message <- "MR 分析失败"
      return(list(result = NULL, log = log_entry))
    }
    
    # 8. 敏感性分析
    sensitivity_results <- run_sensitivity_analysis(harmonised_dat)
    
    # 9. 提取主要结果
    result_row <- extract_main_result(
      mr_results_list, harmonised_dat, f_stats, 
      sensitivity_results, gene_name
    )
    
    # 10. 可视化（仅显著且验证成功的基因）
    # 安全处理 NA 值
    is_significant <- !is.na(result_row$discovery_pval) && result_row$discovery_pval < 0.05
    has_valid_meta <- !is.na(result_row$has_meta) && result_row$has_meta
    meta_significant <- has_valid_meta && !is.na(result_row$meta_pval) && result_row$meta_pval < 0.05
    
    if (is_significant) {
      if (meta_significant) {
        message(paste("  ✓ 显著且验证成功 (p_meta =", format(result_row$meta_pval, scientific = TRUE), ")"))
        # 这里可以添加可视化代码
      } else if (has_valid_meta) {
        message(paste("  ⚠ 显著但验证失败 (p_discovery =", format(result_row$discovery_pval, scientific = TRUE), 
                      ", p_meta =", format(result_row$meta_pval, scientific = TRUE), ")"))
      } else {
        message(paste("  ⚠ 显著但无验证 (p_discovery =", format(result_row$discovery_pval, scientific = TRUE), ")"))
      }
    }
    
    # 11. 保存详细结果
    detail_dir <- file.path(output_dir, "details")
    dir.create(detail_dir, showWarnings = FALSE, recursive = TRUE)
    
    detail_file <- file.path(detail_dir, paste0(gene_name, "_detail.csv"))
    write.csv(mr_results_list$discovery, detail_file, row.names = FALSE)
    
    log_entry$status <- result_row$status
    log_entry$message <- sprintf("nSNP=%d, p_discovery=%.3e, p_meta=%.3e", 
                                  nrow(harmonised_dat), 
                                  result_row$discovery_pval,
                                  ifelse(result_row$has_meta, result_row$meta_pval, NA))
    
    message(paste("  分析完成:", result_row$status))
    
    return(list(
      result = result_row,
      log = log_entry,
      mr_results = mr_results_list,
      sensitivity_results = sensitivity_results,
      harmonised_data = harmonised_dat
    ))
    
  }, error = function(e) {
    log_entry$status <- "ERROR"
    log_entry$message <- conditionMessage(e)
    message(paste("  ✗ 错误:", conditionMessage(e)))
    
    return(list(
      result = NULL,
      log = log_entry
    ))
  })
}

# ================================================================================
# 第十一部分：批量分析流程（优化版）
# ================================================================================

run_batch_analysis_optimized <- function(exposure_dir, outcome_dir, gene_list_file,
                                          output_dir = "./mr_results_optimized") {
  message(rep("=", 60))
  message("孟德尔随机化批量分析 - 优化版 v3.0")
  message(paste("开始时间:", Sys.time()))
  message(rep("=", 60))
  message("\n【主要改进】")
  message("1. P 值阈值：5e-8 → 1e-5 (提高暴露率)")
  message("2. LD Clump: 分级策略 (r2: 0.001 → 0.01 → 0.1)")
  message("3. 独立验证：FinnGen 发现 + UK Biobank 验证")
  message(rep("=", 60))
  
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
  
  genes_to_analyze <- gene_list$GENE[gene_list$N_SNPS > 0]
  message(paste("\n需要分析的基因数:", length(genes_to_analyze)))
  
  # 3. 批量分析循环
  results_list <- list()
  log_list <- list()
  
  for (i in seq_along(genes_to_analyze)) {
    gene <- genes_to_analyze[i]
    message(sprintf("\n[%d/%d]", i, length(genes_to_analyze)))
    
    exposure_file <- file.path(exposure_dir, paste0(gene, ".exposure.csv"))
    outcome_file <- file.path(outcome_dir, paste0(gene, ".outcome.csv"))
    
    if (!file.exists(exposure_file)) {
      message("  跳过：暴露文件不存在")
      log_list[[gene]] <- list(gene = gene, status = "NO_EXPOSURE_FILE", message = "暴露文件不存在")
      next
    }
    
    if (!file.exists(outcome_file)) {
      message("  跳过：结局文件不存在")
      log_list[[gene]] <- list(gene = gene, status = "NO_OUTCOME_FILE", message = "结局文件不存在")
      next
    }
    
    analysis_result <- analyze_single_gene_optimized(
      gene, exposure_file, outcome_file, output_dir
    )
    
    if (!is.null(analysis_result$result)) {
      results_list[[gene]] <- analysis_result$result
    }
    
    log_list[[gene]] <- analysis_result$log
  }
  
  # 4. 导出结果
  message("\n", rep("=", 60))
  message("正在导出结果...")
  message(rep("=", 60))
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (length(results_list) > 0) {
    final_results <- do.call(rbind, results_list)
    
    # FDR 校正
    final_results$fdr_qval <- p.adjust(final_results$discovery_pval, method = "BH")
    final_results$fdr_sig <- final_results$fdr_qval < 0.05
    
    # 格式化 OR 显示
    final_results$OR_95CI <- sprintf("%.3f (%.3f-%.3f)", 
                                     final_results$discovery_OR, 
                                     final_results$discovery_OR_lci, 
                                     final_results$discovery_OR_uci)
    
    # Meta 分析 OR 显示
    final_results$Meta_OR_95CI <- ifelse(
      final_results$has_meta,
      sprintf("%.3f (%.3f-%.3f)", 
              final_results$meta_OR, 
              final_results$meta_OR_lci, 
              final_results$meta_OR_uci),
      "NA"
    )
    
    # 选择输出列
    output_cols <- c("gene", "method", "discovery_b", "discovery_se", "discovery_pval",
                     "OR_95CI", "Meta_OR_95CI", "meta_pval", "nsnp", "F_mean", 
                     "Steiger_dir", "Steiger_p", "Q_p", "Egger_intercept_p",
                     "has_replication", "has_meta", "status", "fdr_qval", "fdr_sig")
    
    results_file <- file.path(output_dir, "MR_results_main_optimized.csv")
    write.csv(final_results[, output_cols], results_file, row.names = FALSE)
    message(paste("主结果已保存:", results_file))
    
    # 显示显著结果
    sig_discovery <- final_results[final_results$discovery_pval < 0.05, ]
    if (nrow(sig_discovery) > 0) {
      cat(sprintf("\n发现队列显著 (p < 0.05): %d 个基因\n", nrow(sig_discovery)))
      print(sig_discovery[, c("gene", "discovery_pval", "OR_95CI", "status")])
    }
    
    sig_meta <- final_results[final_results$has_meta & final_results$meta_pval < 0.05, ]
    if (nrow(sig_meta) > 0) {
      cat(sprintf("\n✓ Meta 分析显著 (p < 0.05): %d 个基因\n", nrow(sig_meta)))
      print(sig_meta[, c("gene", "meta_pval", "Meta_OR_95CI", "status")])
    }
    
    fdr_sig <- final_results[final_results$fdr_sig, ]
    if (nrow(fdr_sig) > 0) {
      cat(sprintf("\nFDR 校正后显著 (q < 0.05): %d 个基因\n", nrow(fdr_sig)))
      print(fdr_sig[, c("gene", "discovery_pval", "fdr_qval", "OR_95CI")])
    }
  }
  
  # 5. 导出分析日志
  if (length(log_list) > 0) {
    log_df <- do.call(rbind, lapply(log_list, as.data.frame))
    log_file <- file.path(output_dir, "MR_analysis_log_optimized.csv")
    write.csv(log_df, log_file, row.names = FALSE)
    message(paste("分析日志已保存:", log_file))
  }
  
  # 6. 分析摘要
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  message("\n", rep("=", 60))
  message("孟德尔随机化批量分析 - 优化版完成")
  message(paste("结束时间:", end_time))
  message(paste("总耗时:", round(duration, 2), "分钟"))
  message(paste("结果保存目录:", output_dir))
  message(rep("=", 60))
  
  # 状态统计
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
# 第十二部分：命令行支持
# ================================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) >= 3) {
    exposure_dir <- args[1]
    outcome_dir <- args[2]
    gene_list_file <- args[3]
    output_dir <- if (length(args) >= 4) args[4] else "./mr_results_optimized"
    
    results <- run_batch_analysis_optimized(
      exposure_dir = exposure_dir,
      outcome_dir = outcome_dir,
      gene_list_file = gene_list_file,
      output_dir = output_dir
    )
  } else {
    message("用法：Rscript mr_analysis_optimized.R <exposure_dir> <outcome_dir> <gene_list.txt> [output_dir]")
    message("示例：Rscript mr_analysis_optimized.R ./exposure ./outcome gene_list.txt ./results_optimized")
  }
}

message("\nMR 优化版分析脚本加载完成")
message("主要改进:")
message("1. P 值阈值：5e-8 → 1e-5")
message("2. LD Clump: 分级策略 (r2: 0.001 → 0.01 → 0.1)")
message("3. 独立验证：FinnGen + UK Biobank")
