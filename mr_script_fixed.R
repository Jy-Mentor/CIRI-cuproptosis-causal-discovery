#!/usr/bin/env Rscript

# 铜死亡与卒中Mendelian Randomization分析
# 关键修复版本：解决变量名冲突和错误追踪问题

library(data.table)
library(TwoSampleMR)
library(coloc)
library(ggplot2)
# library(fst)  # 未使用，注释掉
# library(pingr)  # 未使用，注释掉

# 1. 配置和日志
config <- list(
  eqtl_file = "D:/EQTL/共定位顺式/eqtlgen_ieu_100kbcis.rds",  # 使用100kb窗口的eQTL数据
  gwas_file = "D:/EQTL/finngen_R12_I9_STR",  # Finngen卒中GWAS数据
  output_dir = "D:/EQTL/MR_100kbcis_Results",  # 新输出目录避免覆盖
  eqtl_sample_size = 31684,  # eQTLGen实际样本量
  n_cores = 2,
  cache_dir = "D:/EQTL/MR_Cache",  # 缓存目录
  cache_enabled = TRUE,  # 是否启用缓存
  cache_max_size = 10 * 1024 * 1024 * 1024  # 缓存最大大小（10GB）
)

# 2. 缓存管理模块
# 使用全局环境存储可变状态（解决R的复制-on-修改语义问题）
cache_env <- new.env(hash = TRUE)

# 使用全局变量存储缓存统计
.cache_stats <- new.env(parent = emptyenv())
.cache_stats$size <- 0
.cache_stats$max_size <- config$cache_max_size
.cache_stats$entries <- 0
.cache_stats$max_entries <- 1000
.cache_stats$hits <- 0
.cache_stats$misses <- 0

# 元数据访问器函数
get_cache_stats <- function() {
  c(
    size = .cache_stats$size,
    max_size = .cache_stats$max_size,
    entries = .cache_stats$entries,
    hits = .cache_stats$hits,
    misses = .cache_stats$misses
  )
}

# 创建缓存目录
dir.create(config$cache_dir, recursive = TRUE, showWarnings = FALSE)

# 缓存键生成函数
generate_cache_key <- function(type, ...) {
  params <- list(...)
  params_str <- paste(sapply(params, function(x) {
    if (is.character(x)) return(x)
    if (is.numeric(x)) return(as.character(x))
    if (is.logical(x)) return(as.character(x))
    if (is.data.frame(x)) return(paste(nrow(x), ncol(x), sep = "x"))
    return(toString(x))
  }), collapse = "_")
  key <- paste(type, params_str, sep = "_")
  key <- gsub("[^a-zA-Z0-9_]", "_", key)
  return(key)
}

# 缓存读写函数
cache_get <- function(key) {
  if (!config$cache_enabled) return(NULL)
  if (exists(key, envir = cache_env)) {
    # 更新访问时间
    attr(cache_env[[key]], "last_accessed") <- Sys.time()
    # 增加命中计数
    .cache_stats$hits <- .cache_stats$hits + 1
    return(cache_env[[key]])
  }
  # 尝试从磁盘加载
  disk_path <- file.path(config$cache_dir, paste0(key, ".rds"))
  if (file.exists(disk_path)) {
    tryCatch({
      data <- readRDS(disk_path)
      cache_env[[key]] <- data
      attr(cache_env[[key]], "last_accessed") <- Sys.time()
      # 增加命中计数
      .cache_stats$hits <- .cache_stats$hits + 1
      return(data)
    }, error = function(e) {
      # 增加未命中计数
      .cache_stats$misses <- .cache_stats$misses + 1
      return(NULL)
    })
  }
  # 增加未命中计数
  .cache_stats$misses <- .cache_stats$misses + 1
  return(NULL)
}

cache_set <- function(key, value) {
  if (!config$cache_enabled) return(FALSE)
  
  # 计算数据大小
  data_size <- object.size(value)
  
  # 检查缓存大小
  if (.cache_stats$size + data_size > .cache_stats$max_size) {
    # 清理缓存
    cleanup_cache(data_size)
  }
  
  # 存储到内存
  cache_env[[key]] <- value
  attr(cache_env[[key]], "last_accessed") <- Sys.time()
  attr(cache_env[[key]], "size") <- data_size
  
  # 更新缓存元数据
  .cache_stats$size <- .cache_stats$size + data_size
  .cache_stats$entries <- .cache_stats$entries + 1
  
  # 存储到磁盘
  disk_path <- file.path(config$cache_dir, paste0(key, ".rds"))
  tryCatch({
    saveRDS(value, disk_path)
  }, error = function(e) {
    # 磁盘存储失败，仅保存在内存
  })
  
  return(TRUE)
}

# 缓存清理函数
cleanup_cache <- function(required_size) {
  # 获取所有缓存项
  keys <- ls(envir = cache_env)
  if (length(keys) == 0) return()
  
  # 按访问时间排序
  access_times <- sapply(keys, function(key) {
    attr(cache_env[[key]], "last_accessed", exact = TRUE)
  })
  keys_sorted <- keys[order(access_times)]
  
  # 清理最久未使用的项
  for (key in keys_sorted) {
    if (.cache_stats$size <= .cache_stats$max_size - required_size) break
    
    # 移除内存中的缓存
    item_size <- attr(cache_env[[key]], "size", exact = TRUE)
    .cache_stats$size <- .cache_stats$size - item_size
    .cache_stats$entries <- .cache_stats$entries - 1
    rm(list = key, envir = cache_env)
    
    # 移除磁盘中的缓存
    disk_path <- file.path(config$cache_dir, paste0(key, ".rds"))
    if (file.exists(disk_path)) {
      tryCatch({
        file.remove(disk_path)
      }, error = function(e) {})
    }
  }
}

# 缓存失效函数
cache_invalidate <- function(pattern = NULL) {
  if (!config$cache_enabled) return()
  
  if (is.null(pattern)) {
    # 清除所有缓存
    keys <- ls(envir = cache_env)
    for (key in keys) {
      rm(list = key, envir = cache_env)
    }
    # 清除磁盘缓存
    cache_files <- list.files(config$cache_dir, pattern = "\\.rds$", full.names = TRUE)
    for (file in cache_files) {
      tryCatch({
        file.remove(file)
      }, error = function(e) {})
    }
    # 重置元数据
    .cache_stats$size <- 0
    .cache_stats$entries <- 0
  } else {
    # 清除匹配的缓存
    keys <- ls(envir = cache_env, pattern = pattern)
    for (key in keys) {
      item_size <- attr(cache_env[[key]], "size", exact = TRUE)
      .cache_stats$size <- .cache_stats$size - item_size
      .cache_stats$entries <- .cache_stats$entries - 1
      rm(list = key, envir = cache_env)
      
      # 移除磁盘中的缓存
      disk_path <- file.path(config$cache_dir, paste0(key, ".rds"))
      if (file.exists(disk_path)) {
        tryCatch({
          file.remove(disk_path)
        }, error = function(e) {})
      }
    }
  }
}

# 缓存统计函数
cache_stats <- function() {
  cat("=== 缓存统计 ===\n")
  cat(sprintf("缓存条目数: %d\n", .cache_stats$entries))
  cat(sprintf("缓存大小: %.2f MB\n", .cache_stats$size / 1024 / 1024))
  cat(sprintf("缓存最大大小: %.2f MB\n", .cache_stats$max_size / 1024 / 1024))
  cat(sprintf("缓存命中率: %.2f%%\n", .cache_stats$hits / (.cache_stats$hits + .cache_stats$misses) * 100))
  cat("================\n")
}

# 核心基因列表（整合PC因果推断结果）
core_genes <- c(
  # 上游调控（PC根节点）
  "TLR4", "STAT1", "TGFB1", "NFE2L2",
  # 中间传递（PC中间节点）
  "JAK1", "STAT3", "CCL2", "ICAM1", "HMOX1",
  # 铜死亡执行（PC未覆盖，机制终点）
  "FDX1", "LIAS", "SLC31A1", "ATOX1", "DLAT", "DLD", "PDHA1", "PDHB"
)

# 阈值设置
thresholds <- list(
  pval_exposure = 5.0e-08,
  pval_exposure_relaxed = 5.0e-05,
  f_stat = 10.0,  # 传统的F>10标准，确保工具变量强度
  eaf_diff = 0.2,
  steiger_pval = 0.05,
  ivw_pval = 0.05,
  coloc_pph4_strong = 0.8
)

# 创建输出目录
dir.create(config$output_dir, recursive = TRUE)
log_file <- file.path(config$output_dir, "mr_analysis.log")
sink(log_file, append = TRUE, split = TRUE)

# 清除旧缓存（因为修改了F统计量阈值）
if (config$cache_enabled) {
  cat("清除旧缓存...\n")
  cache_invalidate()
}

cat("\n=== 铜死亡与卒中MR分析（关键修复版）===", "\n")
cat("开始时间:", Sys.time(), "\n")
cat("核心基因数:", length(core_genes), "\n")

# 2. 读取eQTL数据
cat("\n=== 阶段一: 读取eQTL数据 ===", "\n")

# 加载clump后的eQTL数据（支持Excel和RDS文件）
load_eqtl_data_clump <- function(eqtl_file, core_genes) {
  # 生成缓存键
  cache_key <- generate_cache_key("eqtl_data_clump", eqtl_file, paste(sort(core_genes), collapse = ","))
  
  # 尝试从缓存加载
  cached_data <- cache_get(cache_key)
  if (!is.null(cached_data)) {
    cat("从缓存加载eQTL数据...\n")
    return(cached_data)
  }
  
  cat("加载clump后的eQTL数据...\n")
  
  # 根据文件扩展名选择读取方法
  file_ext <- tolower(tools::file_ext(eqtl_file))
  
  if (file_ext == "rds") {
    # 读取RDS文件
    eqtl_data <- as.data.table(readRDS(eqtl_file))
  } else if (file_ext %in% c("xlsx", "xls")) {
    # 加载readxl包
    if (!requireNamespace("readxl", quietly = TRUE)) {
      install.packages("readxl", repos = "https://cran.r-project.org")
    }
    library(readxl)
    
    # 读取Excel文件
    eqtl_data <- as.data.table(read_excel(eqtl_file))
  } else {
    stop(sprintf("不支持的文件格式: %s", file_ext))
  }
  
  cat(sprintf("eQTL文件读取完成，共 %d 条记录\n", nrow(eqtl_data)))
  
  # 动态识别列名
  col_names <- colnames(eqtl_data)
  cat(sprintf("eQTL数据列名: %s\n", paste(col_names, collapse = ", ")))
  
  # 确定关键列名
  pval_col <- if ("pval.exposure" %in% col_names) "pval.exposure" else if ("P" %in% col_names) "P" else stop("无法找到p值列")
  beta_col <- if ("beta.exposure" %in% col_names) "beta.exposure" else if ("beta" %in% col_names) "beta" else stop("无法找到beta列")
  se_col <- if ("se.exposure" %in% col_names) "se.exposure" else if ("se" %in% col_names) "se" else stop("无法找到se列")
  gene_col <- if ("gene" %in% col_names) "gene" else if ("Gene" %in% col_names) "Gene" else stop("无法找到基因列")
  snp_col <- if ("SNP" %in% col_names) "SNP" else if ("snp" %in% col_names) "snp" else stop("无法找到SNP列")
  
  # 识别等位基因相关列
  effect_allele_col <- if ("effect_allele.exposure" %in% col_names) "effect_allele.exposure" else if ("A1" %in% col_names) "A1" else if ("EA" %in% col_names) "EA" else if ("effect_allele" %in% col_names) "effect_allele" else ""
  other_allele_col <- if ("other_allele.exposure" %in% col_names) "other_allele.exposure" else if ("A2" %in% col_names) "A2" else if ("NEA" %in% col_names) "NEA" else if ("other_allele" %in% col_names) "other_allele" else ""
  eaf_col <- if ("eaf.exposure" %in% col_names) "eaf.exposure" else if ("EAF" %in% col_names) "EAF" else if ("ea_freq" %in% col_names) "ea_freq" else if ("effect_allele_frequency" %in% col_names) "effect_allele_frequency" else ""
  
  cat(sprintf("使用列名 - P值: %s, Beta: %s, SE: %s, 基因: %s, SNP: %s\n", 
                  pval_col, beta_col, se_col, gene_col, snp_col))
  
  # 统一列名
  old_names <- c(pval_col, beta_col, se_col, gene_col, snp_col)
  new_names <- c("pval.exposure", "beta.exposure", "se.exposure", "gene", "SNP")
  
  # 添加等位基因相关列的映射
  if (effect_allele_col != "") {
    old_names <- c(old_names, effect_allele_col)
    new_names <- c(new_names, "effect_allele.exposure")
  }
  if (other_allele_col != "") {
    old_names <- c(old_names, other_allele_col)
    new_names <- c(new_names, "other_allele.exposure")
  }
  if (eaf_col != "") {
    old_names <- c(old_names, eaf_col)
    new_names <- c(new_names, "eaf.exposure")
  }
  
  setnames(eqtl_data, old = old_names, new = new_names)
  
  # 检查前10行数据，看看基因列的值
  cat(sprintf("eQTL数据前10行基因: %s\n", paste(head(eqtl_data$gene), collapse = ", ")))
  
  # 筛选核心基因
  eqtl_data_core <- eqtl_data[toupper(gene) %in% toupper(core_genes), ]
  cat(sprintf("核心基因相关记录: %d 条\n", nrow(eqtl_data_core)))
  
  # 强制标准化列名，解决F统计量NA问题
  if (!"beta.exposure" %in% colnames(eqtl_data_core)) {
    if ("beta" %in% colnames(eqtl_data_core)) {
      setnames(eqtl_data_core, "beta", "beta.exposure")
    }
  }
  if (!"se.exposure" %in% colnames(eqtl_data_core)) {
    if ("se" %in% colnames(eqtl_data_core)) {
      setnames(eqtl_data_core, "se", "se.exposure")
    }
  }
  
  # 计算F统计量，看看数据质量
  if (nrow(eqtl_data_core) > 0 && "beta.exposure" %in% colnames(eqtl_data_core) && "se.exposure" %in% colnames(eqtl_data_core)) {
    N <- config$eqtl_sample_size  # 使用配置文件中的实际样本量
    eqtl_data_core[, R2 := beta.exposure^2 / (beta.exposure^2 + se.exposure^2 * N)]
    eqtl_data_core[, F_stat := (N - 1 - 1) / 1 * R2 / (1 - R2)]
    # 处理NA值
    f_stat_values <- eqtl_data_core$F_stat
    f_stat_values <- f_stat_values[!is.na(f_stat_values) & is.finite(f_stat_values)]
    if (length(f_stat_values) > 0) {
      cat(sprintf("核心基因F统计量 - 最小值: %.2f, 中位数: %.2f, 最大值: %.2f\n", 
                    min(f_stat_values), median(f_stat_values), max(f_stat_values)))
    } else {
      cat("核心基因F统计量计算失败，可能存在NA值或无效数据\n")
    }
  } else {
    cat("核心基因数据缺少必要列，无法计算F统计量\n")
  }
  
  # 如果没有核心基因相关记录，打印核心基因列表
  if (nrow(eqtl_data_core) == 0) {
    cat(sprintf("没有找到核心基因相关记录。核心基因列表: %s\n", paste(core_genes, collapse = ", ")))
    # 为了测试，暂时返回所有数据
    # 缓存数据
    cache_set(cache_key, eqtl_data)
    return(eqtl_data)
  }
  
  # 缓存数据
  cache_set(cache_key, eqtl_data_core)
  return(eqtl_data_core)
}

eqtl_data <- load_eqtl_data_clump(config$eqtl_file, core_genes)

# 检查基因列
if (!"gene" %in% colnames(eqtl_data)) {
  stop("eQTL数据缺少'gene'列")
}

# 3. 读取Finngen GWAS数据
cat("\n=== 阶段二: 读取Finngen GWAS数据 ===", "\n")
parse_finngen_gwas <- function(gwas_file) {
  # 生成缓存键
  cache_key <- generate_cache_key("gwas_data", gwas_file)
  
  # 尝试从缓存加载
  cached_data <- cache_get(cache_key)
  if (!is.null(cached_data)) {
    cat("从缓存加载Finngen GWAS数据...\n")
    return(cached_data)
  }
  
  cat("读取Finngen GWAS数据...\n")
  
  # 读取数据
  gwas_data <- fread(gwas_file, sep = "\t", header = TRUE)
  
  cat(sprintf("Finngen数据维度: %d x %d\n", nrow(gwas_data), ncol(gwas_data)))
  cat(sprintf("Finngen数据列名: %s\n", paste(colnames(gwas_data), collapse = ", ")))
  
  # 显示前10行数据
  cat("Finngen数据前10行:\n")
  print(head(gwas_data))
  
  # 重命名列名以符合TwoSampleMR要求
  setnames(gwas_data, 
           old = c("rsids", "beta", "sebeta", "pval", "alt", "ref", "af_alt"), 
           new = c("SNP", "beta.outcome", "se.outcome", "pval.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome"))
  
  # 添加必要的列
  gwas_data[, id.outcome := "finngen_R12_I9_STR"]
  gwas_data[, outcome := "Ischemic Stroke"]
  
  # 清理数据
  gwas_data <- gwas_data[!is.na(SNP) & !is.na(beta.outcome) & !is.na(se.outcome) & !is.na(pval.outcome),
                       .(SNP, effect_allele.outcome, other_allele.outcome, beta.outcome, se.outcome, pval.outcome, eaf.outcome, id.outcome, outcome)]
  
  cat(sprintf("Finngen数据清理完成，共 %d 行有效数据\n", nrow(gwas_data)))
  
  # 缓存数据
  cache_set(cache_key, gwas_data)
  return(gwas_data)
}

outcome_data <- parse_finngen_gwas(config$gwas_file)

# 4. 计算F统计量
calculate_f_stat <- function(dat, N) {
  dat <- as.data.table(dat)
  k <- nrow(dat)
  if (k == 0) return(dat)
  
  # 生成缓存键
  cache_key <- generate_cache_key("f_stat", nrow(dat), N, paste(head(dat$SNP, 5), collapse = ","))
  
  # 尝试从缓存加载
  cached_result <- cache_get(cache_key)
  if (!is.null(cached_result)) {
    cat("从缓存加载F统计量计算结果...\n")
    return(cached_result)
  }
  
  # 确保列存在
  if (!all(c("beta.exposure", "se.exposure") %in% colnames(dat))) {
    cat("缺少 beta.exposure 或 se.exposure 列，无法计算F统计量\n")
    dat[, F_stat := NA_real_]
    return(dat)
  }
  
  # 检查数据是否有NA值
  if (anyNA(dat$beta.exposure) || anyNA(dat$se.exposure)) {
    cat("beta.exposure 或 se.exposure 列存在NA值，移除这些行\n")
    dat <- dat[!is.na(beta.exposure) & !is.na(se.exposure), ]
    if (nrow(dat) == 0) return(dat)
  }
  
  # 计算F统计量（使用更稳健的方法）
  # 简化但稳健的近似：F = (beta/se)^2
  dat[, F_stat := (beta.exposure / se.exposure)^2]
  
  # 保留原公式作为参考（注释掉）
  # dat[, R2 := beta.exposure^2 / (beta.exposure^2 + se.exposure^2 * N)]
  # dat[R2 >= 1, R2 := 0.9999]
  # if (k >= N - 1) {
  #   cat(sprintf("SNP数量 %d 接近或超过样本量 %d，使用近似F统计量\n", k, N))
  #   dat[, F_stat := R2 / (1 - R2) * 1000]
  # } else {
  #   dat[, F_stat := (N - k - 1) / k * R2 / (1 - R2)]
  # }
  
  # 移除 F 统计量为 NA 或无穷大的行
  dat <- dat[is.finite(F_stat) & !is.na(F_stat), ]
  cat(sprintf("计算F统计量后剩余 %d 行数据\n", nrow(dat)))
  
  # 缓存结果
  cache_set(cache_key, dat)
  return(dat)
}

# 5. MR分析函数
run_gene_mr <- function(gene_name, eqtl_data, outcome_data, thresholds, config) {
  # 生成缓存键
  cache_key <- generate_cache_key("mr_analysis", gene_name, paste(thresholds, collapse = ","))
  
  # 尝试从缓存加载
  cached_result <- cache_get(cache_key)
  if (!is.null(cached_result)) {
    cat(sprintf("\n=== 分析基因: %s (从缓存加载) ===\n", gene_name))
    return(cached_result)
  }
  
  cat(sprintf("\n=== 分析基因: %s ===\n", gene_name))
  
  # 初始化结果
  results <- list(
    main = data.table(),
    steiger = data.table(),
    coloc = data.table(),
    sensitivity = data.table(),
    snp_level = data.table(),
    diagnostics = data.table(
      gene = gene_name,
      status = "initialized",
      raw_snps = 0,
      fstat_snps = 0,
      common_snps = 0,
      harmonised_snps = 0
    )
  )
  
  tryCatch({
    # 1. 精确提取目标基因数据（修复变量名冲突）
    cat(sprintf("提取基因 %s 的数据...\n", gene_name))
    rows_match <- which(toupper(eqtl_data$gene) == toupper(gene_name))
    gene_exposure <- eqtl_data[rows_match, ]
    
    if (nrow(gene_exposure) == 0) {
      cat(sprintf("%s 没有可用的eQTL数据，跳过\n", gene_name))
      results$diagnostics$status <- "no_eqtl_data"
      return(results)
    }
    
    results$diagnostics$raw_snps <- nrow(gene_exposure)
    cat(sprintf("提取到 %d 行eQTL数据\n", nrow(gene_exposure)))
    
    # 2. 计算F统计量并过滤弱工具变量
    gene_exposure <- calculate_f_stat(gene_exposure, config$eqtl_sample_size)
    gene_exposure <- gene_exposure[F_stat > thresholds$f_stat, ]
    results$diagnostics$fstat_snps <- nrow(gene_exposure)
    
    if (nrow(gene_exposure) == 0) {
      cat(sprintf("%s 计算F统计量后没有数据，跳过\n", gene_name))
      results$diagnostics$status <- "no_fstat_data"
      return(results)
    }
    
    # 3. 添加TwoSampleMR必需列
    gene_exposure[, `:=`(
      id.exposure = gene_name,
      exposure = gene_name,
      samplesize.exposure = config$eqtl_sample_size  # Steiger必需
    )]
    
    # 4. 检查SNP列
    if (!"SNP" %in% colnames(gene_exposure)) {
      cat(sprintf("%s 数据中没有SNP列，跳过\n", gene_name))
      results$diagnostics$status <- "no_snp_column"
      return(results)
    }
    
    # 5. 匹配SNP
    cat("计算共同SNP...\n")
    common_snps <- intersect(gene_exposure$SNP, outcome_data$SNP)
    results$diagnostics$common_snps <- length(common_snps)
    
    if (length(common_snps) < 3) {
      cat(sprintf("%s 共同SNP数量不足（需要至少3个），跳过\n", gene_name))
      results$diagnostics$status <- "insufficient_snps"
      return(results)
    }
    
    cat(sprintf("找到 %d 个共同SNP\n", length(common_snps)))
    
    # 6. 过滤数据并添加结局样本量
    gene_exposure <- gene_exposure[SNP %in% common_snps, ]
    outcome_data_subset <- outcome_data[SNP %in% common_snps, ]
    outcome_data_subset[, samplesize.outcome := 208200]  # FinnGen实际样本量
    
    # 7. Harmonise数据
    cat("Harmonising数据...\n")
    dat <- harmonise_data(
      exposure_dat = gene_exposure,
      outcome_dat = outcome_data_subset,
      action = 2
    )
    
    if (nrow(dat) < 3) {
      cat(sprintf("%s Harmonise后没有足够数据，跳过\n", gene_name))
      results$diagnostics$status <- "harmonise_empty"
      return(results)
    }
    
    results$diagnostics$harmonised_snps <- nrow(dat)
    cat(sprintf("Harmonise后剩余 %d 个SNP\n", nrow(dat)))
    
    # 检查harmonise后的列名并标准化
    cat("Harmonise后列名:", paste(colnames(dat), collapse = ", "), "\n")
    
    # 标准化列名（解决列名不一致问题）
    if (!"effect_allele.exposure" %in% colnames(dat) && "effect_allele" %in% colnames(dat)) {
      setnames(dat, "effect_allele", "effect_allele.exposure")
    }
    if (!"other_allele.exposure" %in% colnames(dat) && "other_allele" %in% colnames(dat)) {
      setnames(dat, "other_allele", "other_allele.exposure")
    }
    if (!"eaf.exposure" %in% colnames(dat) && "eaf" %in% colnames(dat)) {
      setnames(dat, "eaf", "eaf.exposure")
    }
    if (!"effect_allele.outcome" %in% colnames(dat) && "effect_allele_outcome" %in% colnames(dat)) {
      setnames(dat, "effect_allele_outcome", "effect_allele.outcome")
    }
    if (!"eaf.outcome" %in% colnames(dat) && "eaf_outcome" %in% colnames(dat)) {
      setnames(dat, "eaf_outcome", "eaf.outcome")
    }
    
    # 防御性回文SNP过滤（仅当列存在时）
    if (all(c("effect_allele.exposure", "other_allele.exposure") %in% colnames(dat))) {
      cat("检查回文SNP...\n")
      # 定义回文等位基因对
      palindromic <- (dat$effect_allele.exposure == "A" & dat$other_allele.exposure == "T") | 
                     (dat$effect_allele.exposure == "T" & dat$other_allele.exposure == "A") | 
                     (dat$effect_allele.exposure == "C" & dat$other_allele.exposure == "G") | 
                     (dat$effect_allele.exposure == "G" & dat$other_allele.exposure == "C")
      
      n_palindromic <- sum(palindromic, na.rm = TRUE)
      if (n_palindromic > 0) {
        cat(sprintf("移除 %d 个回文SNP\n", n_palindromic))
        dat <- dat[!palindromic, ]
      }
      
      if (nrow(dat) < 3) {
        cat(sprintf("%s 移除回文SNP后样本不足，跳过\n", gene_name))
        results$diagnostics$status <- "palindromic_removed"
        return(results)
      }
    } else {
      cat("警告：缺少等位基因列，跳过回文SNP检查\n")
    }
    
    # 8. 进行MR分析
    cat("进行MR分析...\n")
    res <- mr(dat)
    results$main <- as.data.table(res)
    
    # 9. Steiger过滤
    cat("进行Steiger过滤...\n")
    steiger_res <- steiger_filtering(dat)
    # 输出Steiger结果（高优先级修复）
    steiger_dir <- NA
    steiger_p <- NA
    if (!is.null(steiger_res) && length(steiger_res) > 0) {
      if ("steiger_dir" %in% names(steiger_res)) {
        steiger_dir <- steiger_res$steiger_dir[1]  # TRUE表示方向正确
      }
      if ("steiger_pval" %in% names(steiger_res)) {
        steiger_p <- steiger_res$steiger_pval[1]
      }
    }
    cat(sprintf("Steiger方向检验: 方向正确=%s, p值=%.4f\n", steiger_dir, steiger_p))
    
    # Steiger过滤（强化版：只有当方向正确且p值显著时才保留）
    if (!is.na(steiger_dir) && !is.na(steiger_p)) {
      if (!steiger_dir && steiger_p < thresholds$steiger_pval) {
        # 方向错误且p值显著，排除
        cat(sprintf("排除: %s 存在显著反向因果 (Steiger P=%.4f)\n", gene_name, steiger_p))
        results$diagnostics$status <- "reverse_causation"
        return(results)
      } else if (steiger_dir && steiger_p < thresholds$steiger_pval) {
        # 方向正确且p值显著，通过
        cat(sprintf("通过: %s 方向正确 (Steiger P=%.4f)\n", gene_name, steiger_p))
      } else {
        # 方向不确定或p值不显著，发出警告
        cat(sprintf("警告: %s 方向不确定或p值不显著 (Steiger P=%.4f)\n", gene_name, steiger_p))
      }
    }
    
    # 保存Steiger结果
    steiger_dt <- as.data.table(steiger_res)
    if (nrow(steiger_dt) == 0) {
      steiger_dt <- data.table(
        gene = gene_name,
        steiger_dir = steiger_dir,
        steiger_pval = steiger_p
      )
    } else {
      steiger_dt[, gene := gene_name]
    }
    results$steiger <- steiger_dt
    
    # 10. 敏感性分析（强化版）
    cat("进行敏感性分析...\n")
    tryCatch({
      pleio <- mr_pleiotropy_test(dat)
      hetero <- mr_heterogeneity(dat)
      
      # 敏感性分析防御（高优先级修复）
      egger_p <- NA_real_
      q_ivw_p <- NA_real_
      
      if (!is.null(pleio) && nrow(pleio) > 0 && "pval" %in% names(pleio)) {
        egger_p <- pleio$pval[1]
      } else {
        cat("多效性检验失败: 结果为空或格式错误\n")
      }
      
      if (!is.null(hetero) && nrow(hetero) > 0) {
        # 查找IVW行
        ivw_idx <- grep("Inverse variance weighted", hetero$method, ignore.case = TRUE)
        if (length(ivw_idx) > 0) {
          q_ivw_p <- hetero$pval[ivw_idx[1]]
        } else {
          cat("异质性检验失败: 未找到IVW结果\n")
        }
      } else {
        cat("异质性检验失败: 结果为空或格式错误\n")
      }
      
      # 添加MR-PRESSO（中优先级）
      presso_p <- NA_real_
      tryCatch({
        if (nrow(dat) >= 3) {
          # 根据SNP数量动态调整分布数
          n_dist <- min(5000, max(1000, nrow(dat) * 100))
          presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                             SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                             data = dat, NbDistribution = n_dist)
          if (!is.null(presso) && !is.null(presso$`MR-PRESSO results`)) {
            presso_p <- presso$`MR-PRESSO results`$`Global test`$`P-value`
          }
        }
      }, error = function(e) {
        cat(sprintf("MR-PRESSO失败: %s，跳过\n", conditionMessage(e)))
      })
      
      sensitivity_res <- data.table(
        gene = gene_name,
        egger_p = egger_p,
        Q_ivw_p = q_ivw_p,
        presso_p = presso_p,
        n_snps = nrow(dat)  # 添加SNP数量用于诊断
      )
      results$sensitivity <- sensitivity_res
    }, error = function(e) {
      cat(sprintf("敏感性分析失败: %s，跳过\n", conditionMessage(e)))
      results$sensitivity <- data.table(
        gene = gene_name,
        egger_p = NA_real_,
        Q_ivw_p = NA_real_,
        presso_p = NA_real_,
        n_snps = nrow(dat)
      )
    })
    
    # 11. SNP水平结果
    results$snp_level <- dat
    
    # 12. 共定位分析（如果有足够SNP）
    if (nrow(dat) >= 15 && !anyNA(dat$eaf.exposure) && !anyNA(dat$eaf.outcome)) {
      cat("进行共定位分析...\n")
      tryCatch({
        # 准备coloc数据
        dat <- as.data.table(dat)
        
        # 共定位数据清洗（高优先级修复）
        dat_clean <- dat[complete.cases(pval.exposure, pval.outcome, eaf.exposure, eaf.outcome, SNP), ]
        dat_clean <- dat_clean[is.finite(pval.exposure) & pval.exposure > .Machine$double.eps & pval.exposure < 1, ]
        dat_clean <- dat_clean[is.finite(pval.outcome) & pval.outcome > .Machine$double.eps & pval.outcome < 1, ]
        
        # 确保长度一致
        if (nrow(dat_clean) < 15) {
          cat(sprintf("SNP数量不足(%d<15)，跳过共定位分析（统计效力不足）\n", nrow(dat_clean)))
          return(results)
        }
        
        # 确保MAF在0-0.5范围内
        dat_clean[, maf1 := ifelse(eaf.exposure > 0.5, 1-eaf.exposure, eaf.exposure)]
        dat_clean[, maf2 := ifelse(eaf.outcome > 0.5, 1-eaf.outcome, eaf.outcome)]
        
        # 确保MAF有效
        dat_clean <- dat_clean[is.finite(maf1) & maf1 > 0 & maf1 < 0.5, ]
        dat_clean <- dat_clean[is.finite(maf2) & maf2 > 0 & maf2 < 0.5, ]
        
        if (nrow(dat_clean) < 15) {
          cat(sprintf("SNP数量不足(%d<15)，跳过共定位分析（统计效力不足）\n", nrow(dat_clean)))
          return(results)
        }
        
        # 准备数据集
        dataset1 <- list(
          pvalues = dat_clean$pval.exposure,
          N = config$eqtl_sample_size,
          type = "quant",
          MAF = dat_clean$maf1,
          snp = dat_clean$SNP
        )
        
        dataset2 <- list(
          pvalues = dat_clean$pval.outcome,
          N = 208200,
          type = "cc",
          MAF = dat_clean$maf2,
          snp = dat_clean$SNP,
          s = 0.02  # 缺血性卒中患病率约2%
        )
        
        # 运行共定位分析
        coloc_res <- coloc.abf(
          dataset1 = dataset1,
          dataset2 = dataset2
        )
        
        # 提取结果
        if (!is.null(coloc_res) && !is.null(coloc_res$summary)) {
          results$coloc <- data.table(
            gene = gene_name,
            pph0 = coloc_res$summary[1,2],
            pph1 = coloc_res$summary[2,2],
            pph2 = coloc_res$summary[3,2],
            pph3 = coloc_res$summary[4,2],
            pph4 = coloc_res$summary[5,2]
          )
          cat(sprintf("共定位分析完成，PP.H4值: %.3f%%\n", coloc_res$summary[5,2] * 100))
        } else {
          cat("共定位分析失败: 结果格式错误，跳过\n")
        }
      }, error = function(e) {
        cat(sprintf("共定位分析失败: %s，跳过\n", conditionMessage(e)))
      })
    }
    
    results$diagnostics$status <- "completed"
    cat(sprintf("%s 分析完成\n", gene_name))
    
  }, error = function(e) {
    cat(sprintf("%s 分析失败: %s\n", gene_name, conditionMessage(e)))
    results$diagnostics$status <- paste0("error: ", conditionMessage(e))
  })
  
  # 缓存结果
  cache_set(cache_key, results)
  return(results)
}

# 6. 主分析
cat("\n=== 阶段三: 基因水平MR分析 ===", "\n")
all_results <- list()
all_diagnostics <- data.table()

for (gene in core_genes) {
  res <- run_gene_mr(gene, eqtl_data, outcome_data, thresholds, config)
  all_results[[gene]] <- res
  all_diagnostics <- rbind(all_diagnostics, res$diagnostics)
}

# 7. 结果汇总和保存
cat("\n=== 阶段四: 结果汇总 ===", "\n")

# 保存诊断信息
diagnostics_file <- file.path(config$output_dir, "diagnostics.csv")
fwrite(all_diagnostics, diagnostics_file)
cat("诊断信息已保存到:", diagnostics_file, "\n")

# 汇总主结果
main_results <- data.table()
for (gene in core_genes) {
  if (nrow(all_results[[gene]]$main) > 0) {
    tmp <- all_results[[gene]]$main
    tmp$gene <- gene
    main_results <- rbind(main_results, tmp)
  }
}

# 保存主结果
if (nrow(main_results) > 0) {
  main_file <- file.path(config$output_dir, "mr_main_results.csv")
  fwrite(main_results, main_file)
  cat("主结果已保存到:", main_file, "\n")
  
  # 筛选显著结果
  significant_results <- main_results[pval < thresholds$ivw_pval, ]
  if (nrow(significant_results) > 0) {
    sig_file <- file.path(config$output_dir, "mr_significant_results.csv")
    fwrite(significant_results, sig_file)
    cat("显著结果已保存到:", sig_file, "\n")
  }
  
  # 汇总敏感性分析结果（中优先级修复）
  sens_results <- data.table()
  for (gene in core_genes) {
    if (nrow(all_results[[gene]]$sensitivity) > 0) {
      sens_results <- rbind(sens_results, all_results[[gene]]$sensitivity)
    }
  }
  if (nrow(sens_results) > 0) {
    sens_file <- file.path(config$output_dir, "mr_sensitivity.csv")
    fwrite(sens_results, sens_file)
    cat("敏感性分析结果已保存到:", sens_file, "\n")
  }
  
  # 汇总Steiger结果
  steiger_results <- data.table()
  for (gene in core_genes) {
    if (nrow(all_results[[gene]]$steiger) > 0) {
      steiger_results <- rbind(steiger_results, all_results[[gene]]$steiger)
    }
  }
  if (nrow(steiger_results) > 0) {
    steiger_file <- file.path(config$output_dir, "mr_steiger.csv")
    fwrite(steiger_results, steiger_file)
    cat("Steiger方向检验结果已保存到:", steiger_file, "\n")
  }
  
  # 汇总共定位结果
  coloc_results <- data.table()
  for (gene in core_genes) {
    if (nrow(all_results[[gene]]$coloc) > 0) {
      coloc_results <- rbind(coloc_results, all_results[[gene]]$coloc)
    }
  }
  if (nrow(coloc_results) > 0) {
    coloc_file <- file.path(config$output_dir, "mr_coloc.csv")
    fwrite(coloc_results, coloc_file)
    cat("共定位分析结果已保存到:", coloc_file, "\n")
  }
}

# 8. 分析摘要
cat("\n=== 分析摘要 ===\n")
cat("总基因数:", length(core_genes), "\n")
cat("成功分析数:", sum(all_diagnostics$status == "completed"), "\n")

if (exists("main_results")) {
  significant_count <- ifelse(exists("significant_results"), nrow(significant_results), 0)
  cat("显著基因数 (IVW p<0.05):", significant_count, "\n")
}

# 输出缓存统计
if (config$cache_enabled) {
  cat("\n=== 缓存统计 ===\n")
  cat(sprintf("缓存条目数: %d\n", .cache_stats$entries))
  cat(sprintf("缓存大小: %.2f MB\n", .cache_stats$size / 1024 / 1024))
  cat(sprintf("缓存最大大小: %.2f MB\n", .cache_stats$max_size / 1024 / 1024))
  if (.cache_stats$hits + .cache_stats$misses > 0) {
    cat(sprintf("缓存命中率: %.2f%%\n", .cache_stats$hits / (.cache_stats$hits + .cache_stats$misses) * 100))
  } else {
    cat("缓存命中率: 0.00%\n")
  }
  cat("================\n")
}

cat("\n结果已保存到:", config$output_dir, "\n")
cat("完成时间:", Sys.time(), "\n")

# 关闭日志
sink()
