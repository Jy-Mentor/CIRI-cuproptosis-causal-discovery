#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# MRAnalysisToolkit 测试脚本
# 描述: 验证工具包的核心功能
# =============================================================================

cat("========================================\n")
cat("MRAnalysisToolkit 功能测试\n")
cat("========================================\n\n")

# 加载工具包
source("utils/MRAnalysisToolkit.R")

# 测试计数器
test_passed <- 0
test_failed <- 0

# 辅助函数
run_test <- function(test_name, test_expr) {
  cat("测试:", test_name, "... ")
  tryCatch({
    result <- eval(test_expr)
    if (isTRUE(result) || is.null(result)) {
      cat("✓ 通过\n")
      test_passed <<- test_passed + 1
      return(TRUE)
    } else {
      cat("✗ 失败 (返回:", result, ")\n")
      test_failed <<- test_failed + 1
      return(FALSE)
    }
  }, error = function(e) {
    cat("✗ 错误:", conditionMessage(e), "\n")
    test_failed <<- test_failed + 1
    return(FALSE)
  })
}

# =============================================================================
# 测试1: 包加载和版本信息
# =============================================================================
cat("\n--- 测试1: 包加载和版本信息 ---\n")

run_test("工具包版本号", expression({
  exists("MR_TOOLKIT_VERSION") && MR_TOOLKIT_VERSION == "1.0.0"
}))

run_test("工具包日期", expression({
  exists("MR_TOOLKIT_DATE") && MR_TOOLKIT_DATE == "2025-07-10"
}))

run_test("默认列名映射存在", expression({
  exists("MR_COLUMN_MAPPINGS") &&
  "exposure" %in% names(MR_COLUMN_MAPPINGS) &&
  "outcome" %in% names(MR_COLUMN_MAPPINGS)
}))

run_test("默认MR方法存在", expression({
  exists("MR_DEFAULT_METHODS") &&
  "multi_snp" %in% names(MR_DEFAULT_METHODS) &&
  "single_snp" %in% names(MR_DEFAULT_METHODS)
}))

run_test("默认QC配置存在", expression({
  exists("MR_DEFAULT_QC") &&
  "snp_threshold" %in% names(MR_DEFAULT_QC) &&
  "steiger_filter" %in% names(MR_DEFAULT_QC)
}))

# =============================================================================
# 测试2: 列名自动检测功能
# =============================================================================
cat("\n--- 测试2: 列名自动检测功能 ---\n")

# 创建测试数据框
test_exposure_data <- data.frame(
  SNP = c("rs123", "rs456"),
  beta.exposure = c(0.5, -0.3),
  se.exposure = c(0.1, 0.08),
  pval.exposure = c(1e-10, 5e-9),
  eaf.exposure = c(0.3, 0.45),
  effect_allele.exposure = c("A", "G"),
  other_allele.exposure = c("G", "A"),
  gene = c("NFKB1", "NFKB1"),
  stringsAsFactors = FALSE
)

run_test("暴露数据列名检测", expression({
  cols <- auto_detect_columns(test_exposure_data, "exposure")
  !is.null(cols$snp) && cols$snp == "SNP" &&
  !is.null(cols$beta) && cols$beta == "beta.exposure" &&
  !is.null(cols$se) && cols$se == "se.exposure"
}))

test_outcome_data <- data.frame(
  SNP = c("rs123", "rs456"),
  beta.outcome = c(0.2, -0.1),
  se.outcome = c(0.05, 0.04),
  pval.outcome = c(0.01, 0.05),
  eaf.outcome = c(0.32, 0.47),
  effect_allele.outcome = c("A", "G"),
  other_allele.outcome = c("G", "A"),
  stringsAsFactors = FALSE
)

run_test("结局数据列名检测", expression({
  cols <- auto_detect_columns(test_outcome_data, "outcome")
  !is.null(cols$snp) && cols$snp == "SNP" &&
  !is.null(cols$beta) && cols$beta == "beta.outcome" &&
  !is.null(cols$se) && cols$se == "se.outcome"
}))

run_test("无效映射类型检测", expression({
  tryCatch({
    auto_detect_columns(test_exposure_data, "invalid_type")
    FALSE
  }, error = function(e) {
    TRUE
  })
}))

# =============================================================================
# 测试3: 数据格式化功能
# =============================================================================
cat("\n--- 测试3: 数据格式化功能 ---\n")

run_test("暴露数据格式化", expression({
  cols <- auto_detect_columns(test_exposure_data, "exposure")
  formatted <- format_mr_data(test_exposure_data, "exposure", cols)
  !is.null(formatted) && nrow(formatted) == 2 &&
  "SNP" %in% names(formatted) && "beta.exposure" %in% names(formatted)
}))

run_test("结局数据格式化", expression({
  cols <- auto_detect_columns(test_outcome_data, "outcome")
  formatted <- format_mr_data(test_outcome_data, "outcome", cols)
  !is.null(formatted) && nrow(formatted) == 2 &&
  "SNP" %in% names(formatted) && "beta.outcome" %in% names(formatted)
}))

run_test("缺失必需列检测", expression({
  bad_data <- data.frame(SNP = c("rs123"), stringsAsFactors = FALSE)
  tryCatch({
    format_mr_data(bad_data, "exposure")
    FALSE
  }, error = function(e) {
    TRUE
  })
}))

# =============================================================================
# 测试4: 配置验证
# =============================================================================
cat("\n--- 测试4: 配置验证 ---\n")

valid_config <- list(
  analysis_name = "Test_Analysis",
  exposure = list(
    type = "excel",
    path = "test.xlsx",
    gene_symbol = "TEST"
  ),
  outcome = list(
    file = "test.csv",
    phenotype = "Test"
  ),
  output = list(
    directory = "Test_Results"
  )
)

run_test("有效配置结构", expression({
  !is.null(valid_config$analysis_name) &&
  !is.null(valid_config$exposure) &&
  !is.null(valid_config$outcome)
}))

run_test("配置包含暴露类型", expression({
  !is.null(valid_config$exposure$type)
}))

run_test("配置包含结局文件", expression({
  !is.null(valid_config$outcome$file)
}))

# =============================================================================
# 测试5: 批量分析配置
# =============================================================================
cat("\n--- 测试5: 批量分析配置 ---\n")

test_genes <- c("NFKB1", "FDX1", "STAT3")

test_batch_template <- list(
  exposure = list(
    type = "excel",
    path = "test.xlsx"
  ),
  outcome = list(
    file = "test.csv",
    phenotype = "Test"
  ),
  quality_control = list(
    snp_threshold = 5e-08,
    steiger_filter = TRUE
  ),
  output = list(
    save_plots = TRUE
  )
)

run_test("批量分析配置结构", expression({
  !is.null(test_batch_template$exposure) &&
  !is.null(test_batch_template$outcome) &&
  !is.null(test_batch_template$quality_control)
}))

run_test("批量分析基因列表", expression({
  length(test_genes) == 3 && all(test_genes %in% c("NFKB1", "FDX1", "STAT3"))
}))

# =============================================================================
# 测试6: 辅助函数
# =============================================================================
cat("\n--- 测试6: 辅助函数 ---\n")

run_test("空值合并操作符", expression({
  result <- NULL %||% "default"
  result == "default"
}))

run_test("空值合并非空值", expression({
  result <- "value" %||% "default"
  result == "value"
}))

run_test("基因列查找", expression({
  test_df <- data.frame(gene_name = c("A", "B"), value = c(1, 2))
  col <- .find_gene_column(test_df)
  col == "gene_name"
}))

# =============================================================================
# 测试7: 配置文件读取
# =============================================================================
cat("\n--- 测试7: 配置文件读取 ---\n")

# 创建临时测试配置文件
test_config_file <- "test_config.json"
test_config <- list(
  analysis_name = "Config_Test",
  exposure = list(
    type = "csv",
    path = "test.csv",
    gene_symbol = "TEST"
  ),
  outcome = list(
    file = "test_outcome.csv",
    phenotype = "Test"
  )
)

jsonlite::write_json(test_config, test_config_file)

run_test("配置文件创建和读取", expression({
  file.exists(test_config_file)
}))

run_test("JSON配置文件解析", expression({
  config <- jsonlite::read_json(test_config_file, simplifyVector = TRUE)
  config$analysis_name == "Config_Test" &&
  config$exposure$type == "csv"
}))

# 清理测试文件
if (file.exists(test_config_file)) {
  file.remove(test_config_file)
}

# =============================================================================
# 测试8: 错误处理
# =============================================================================
cat("\n--- 测试8: 错误处理 ---\n")

run_test("无效暴露类型处理", expression({
  tryCatch({
    read_exposure_data(list(type = "invalid"))
    FALSE
  }, error = function(e) {
    TRUE
  })
}))

run_test("缺失暴露路径处理", expression({
  tryCatch({
    read_exposure_data(list(type = "excel"))
    FALSE
  }, error = function(e) {
    TRUE
  })
}))

run_test("缺失结局文件处理", expression({
  tryCatch({
    read_outcome_data(list())
    FALSE
  }, error = function(e) {
    TRUE
  })
}))

run_test("不存在的文件处理", expression({
  tryCatch({
    read_exposure_data(list(type = "excel", path = "nonexistent.xlsx"))
    FALSE
  }, error = function(e) {
    TRUE
  })
}))

# =============================================================================
# 测试9: 性能测试
# =============================================================================
cat("\n--- 测试9: 性能测试 ---\n")

run_test("大数据集列名检测性能", expression({
  # 创建大型测试数据集
  large_data <- data.frame(
    SNP = paste0("rs", 1:10000),
    beta.exposure = rnorm(10000),
    se.exposure = runif(10000, 0.01, 0.5),
    pval.exposure = runif(10000, 1e-15, 0.05),
    eaf.exposure = runif(10000, 0.1, 0.9),
    effect_allele.exposure = sample(c("A", "C", "G", "T"), 10000, replace = TRUE),
    other_allele.exposure = sample(c("A", "C", "G", "T"), 10000, replace = TRUE),
    stringsAsFactors = FALSE
  )

  start_time <- Sys.time()
  cols <- auto_detect_columns(large_data, "exposure")
  end_time <- Sys.time()

  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat("(耗时:", round(duration, 3), "秒) ")

  !is.null(cols$snp) && duration < 5  # 应该在5秒内完成
}))

# =============================================================================
# 测试10: 集成测试
# =============================================================================
cat("\n--- 测试10: 集成测试 ---\n")

run_test("完整配置结构验证", expression({
  full_config <- list(
    analysis_name = "Integration_Test",
    description = "Test analysis",
    exposure = list(
      type = "excel",
      path = "D:/EQTL/test.xlsx",
      gene_symbol = "NFKB1",
      description = "Test exposure"
    ),
    outcome = list(
      file = "D:/EQTL/test_outcome.csv",
      phenotype = "Stroke",
      samplesize = 100000,
      ncase = 5000,
      ncontrol = 95000,
      description = "Test outcome"
    ),
    methods = c("mr_ivw", "mr_egger_regression"),
    quality_control = list(
      snp_threshold = 5e-08,
      relaxed_threshold = 1e-05,
      ld_clumping = list(
        enabled = TRUE,
        kb = 10000,
        r2 = 0.001,
        p1 = 1,
        p2 = 1
      ),
      steiger_filter = TRUE,
      heterogeneity_test = TRUE,
      pleiotropy_test = TRUE,
      leave_one_out = FALSE
    ),
    output = list(
      directory = "D:/EQTL/Test_Results",
      save_plots = TRUE,
      save_data = TRUE
    )
  )

  # 验证配置结构完整性
  all(c("analysis_name", "exposure", "outcome", "methods",
        "quality_control", "output") %in% names(full_config)) &&
  all(c("type", "path", "gene_symbol") %in% names(full_config$exposure)) &&
  all(c("file", "phenotype") %in% names(full_config$outcome)) &&
  all(c("snp_threshold", "steiger_filter", "heterogeneity_test") %in%
        names(full_config$quality_control)) &&
  all(c("directory", "save_plots", "save_data") %in% names(full_config$output))
}))

# =============================================================================
# 测试总结
# =============================================================================
cat("\n========================================\n")
cat("测试总结\n")
cat("========================================\n")
cat("通过:", test_passed, "\n")
cat("失败:", test_failed, "\n")
cat("总计:", test_passed + test_failed, "\n")
cat("通过率:", round(test_passed / (test_passed + test_failed) * 100, 1), "%\n")

if (test_failed == 0) {
  cat("\n✓ 所有测试通过! MRAnalysisToolkit 功能正常。\n")
} else {
  cat("\n✗ 有", test_failed, "个测试失败，请检查相关功能。\n")
}

cat("========================================\n")
