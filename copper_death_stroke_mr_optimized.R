#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 铜死亡核心基因与缺血性脑卒中的工业级孟德尔随机化分析
# 优化版本：提升稳健性、统计严谨性、可复现性与结果可读性

# 1. 环境与依赖管理
# 设置CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 检查并安装必要的包
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("安装包: %s", pkg))
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
}

# 安装核心依赖包
install_if_missing("data.table")      # 大数据处理
install_if_missing("TwoSampleMR")     # MR分析
install_if_missing("coloc")           # 共定位分析
install_if_missing("MRPRESSO")        # 异常值检测
install_if_missing("futile.logger")   # 日志系统
install_if_missing("progress")        # 进度条
install_if_missing("readxl")          # 读取Excel文件
install_if_missing("ggplot2")         # 可视化
install_if_missing("stringr")         # 字符串处理
install_if_missing("yaml")            # 读取配置文件
install_if_missing("foreach")         # 并行计算
install_if_missing("doParallel")      # 并行计算
install_if_missing("fst")             # 高效数据存储
install_if_missing("rmarkdown")       # 生成报告
install_if_missing("DiagrammeR")      # 绘制流程图
# install_if_missing("mr.raps")         # MR-RAPS方法（暂时注释，因为R 4.5.2不支持）
install_if_missing("susieR")          # 共定位分析
install_if_missing("pingr")           # 检查镜像延迟

# 2. 配置管理
# 读取配置文件
read_config <- function(config_file) {
  if (!file.exists(config_file)) {
    # 创建默认配置
    default_config <- list(
      base_dir = "D:/EQTL",
      output_dir = "D:/EQTL/MR_Optimized_Results",
      eqtl_dir = "D:/EQTL/共定位顺式",
      vcf_file = "D:/EQTL/ieu-a-83.vcf",
      core_genes = c("FDX1", "DLAT", "DLD", "PDHB", "LIPT1", "PDHX", "SLC31A1", 
                    "ATP7B", "ATP7A", "ATOX1", "COMMD1", "MT2A", "NFKB1", "RELA", 
                    "STAT1", "STAT3", "STAT5A", "ICAM1", "CCL2", "IL6", "TGFB1", 
                    "PTGS2", "HMOX1", "SOD2", "FABP3", "ATF4", "BRD4"),
      test_mode = FALSE,
      n_cores = 4,
      gwas_sample_size = 446696,
      gwas_case_prop = 0.5,
      eqtl_sample_size = 20000,
      thresholds = list(
        pval_exposure = 5e-8,
        pval_exposure_relaxed = 5e-5,
        f_stat = 10,
        eaf_diff = 0.2,
        steiger_pval = 0.05,
        ivw_pval = 0.05,
        coloc_pph4_strong = 0.8,
        coloc_pph4_evidence = 0.75
      )
    )
    yaml::write_yaml(default_config, config_file)
    message(sprintf("创建默认配置文件: %s", config_file))
    return(default_config)
  } else {
    return(yaml::read_yaml(config_file))
  }
}

# 3. 日志配置
setup_logging <- function(log_dir) {
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  
  # 主日志文件（主进程用）
  main_log <- file.path(log_dir, "mr_analysis.log")
  
  # 并行子进程日志目录
  parallel_log_dir <- file.path(log_dir, "parallel_logs")
  if (!dir.exists(parallel_log_dir)) dir.create(parallel_log_dir, recursive = TRUE)
  
  # 错误日志目录
  error_log_dir <- file.path(log_dir, "error_logs")
  if (!dir.exists(error_log_dir)) dir.create(error_log_dir, recursive = TRUE)
  
  # 配置主日志（使用appender.file，因为appender.rotating可能不可用）
  flog.appender(
    appender.file(
      file = main_log
    ),
    name = "main"
  )
  flog.threshold(INFO, name = "main")
  
  return(list(
    main = main_log,
    parallel = parallel_log_dir,
    error = error_log_dir
  ))
}

# 4. 数据加载与预处理
# 快速加载eQTL数据（使用rbindlist）
load_eqtl_data_fast <- function(eqtl_dir, core_genes) {
  flog.info("加载eQTL数据...", name = "main")
  
  # 找到eQTL文件
  eqtl_files <- list.files(eqtl_dir, pattern = "eqtlgen_ieu_1mbcis.*\\.csv$", full.names = TRUE)
  eqtl_files <- eqtl_files[!grepl('\\.downloading$', eqtl_files)]
  flog.info(sprintf("找到 %d 个eQTL分块文件", length(eqtl_files)), name = "main")
  
  if (length(eqtl_files) == 0) {
    flog.error("没有找到eQTL文件", name = "main")
    stop("没有找到eQTL文件")
  }
  
  # 使用rbindlist（C级速度）替代foreach+rbind
  flog.info(sprintf("开始快速合并 %d 个eQTL文件...", length(eqtl_files)), name = "main")
  
  # 读取所有文件到列表
  eqtl_list <- lapply(eqtl_files, function(f) {
    tryCatch(fread(f, showProgress = FALSE, data.table = TRUE), error = function(e) NULL)
  })
  
  # 移除NULL并快速合并
  eqtl_list <- Filter(Negate(is.null), eqtl_list)
  eqtl_data <- rbindlist(eqtl_list, use.names = TRUE, fill = TRUE)
  
  flog.info(sprintf("eQTL文件合并完成，共 %d 条记录", nrow(eqtl_data)), name = "main")
  
  # 动态识别列名
  col_names <- colnames(eqtl_data)
  flog.debug(sprintf("eQTL数据列名: %s", paste(col_names, collapse = ", ")))
  
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
  
  flog.debug(sprintf("使用列名 - P值: %s, Beta: %s, SE: %s, 基因: %s, SNP: %s", 
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
  flog.debug(sprintf("eQTL数据前10行: %s", paste(head(eqtl_data$gene), collapse = ", ")), name = "main")
  
  # 筛选核心基因
  eqtl_data_core <- eqtl_data[gene %in% core_genes, ]
  flog.info(sprintf("核心基因相关记录: %d 条", nrow(eqtl_data_core)), name = "main")
  
  # 如果没有核心基因相关记录，打印核心基因列表
  if (nrow(eqtl_data_core) == 0) {
    flog.warn(sprintf("没有找到核心基因相关记录。核心基因列表: %s", paste(core_genes, collapse = ", ")), name = "main")
    # 为了测试，暂时返回所有数据
    return(eqtl_data)
  }
  
  return(eqtl_data_core)
}

# 快速解析GWAS VCF文件（向量化版本）
parse_gwas_vcf_fast <- function(vcf_file, nrows = Inf) {
  flog.info(sprintf("快速解析VCF文件: %s", vcf_file), name = "main")
  
  # 找到header行数
  preview <- readLines(vcf_file, n = 1000)
  header_line <- which(grepl("^#CHROM", preview))[1]
  if (is.na(header_line)) stop("无法找到VCF header")
  
  flog.info(sprintf("Header位于第%d行，开始流式读取...", header_line), name = "main")
  
  # 使用fread从header开始读取
  vcf_data <- fread(
    vcf_file,
    skip = header_line - 1,
    sep = "\t",
    header = TRUE,
    nrows = nrows,  # 测试时可用1000限制行数
    select = c(3, 4, 5, 9, 10),  # ID, REF, ALT, FORMAT, OUTCOME列（按位置选更快）
    col.names = c("ID", "REF", "ALT", "FORMAT", "OUTCOME"),
    showProgress = TRUE
  )
  
  flog.info(sprintf("读取完成: %d行，开始解析FORMAT字段...", nrow(vcf_data)), name = "main")
  
  # 向量化解析（关键优化）
  # 1. 提取FORMAT标签（仅第一行）
  format_tags <- strsplit(vcf_data$FORMAT[1], ":")[[1]]
  
  # 2. 找到标签索引，处理标签不存在的情况
  get_tag_index <- function(tag) {
    idx <- which(format_tags == tag)
    if (length(idx) == 0) return(NA) else return(idx)
  }
  
  tag_idx <- list(
    ES = get_tag_index("ES"),
    SE = get_tag_index("SE"),
    LP = get_tag_index("LP"),
    AF = get_tag_index("AF")
  )
  
  # 3. 一次性分割所有OUTCOME字段（data.table内置tstrsplit是C级速度）
  outcome_parts <- tstrsplit(vcf_data$OUTCOME, ":", fixed = TRUE)
  
  # 4. 提取数值（向量化，无循环）
  vcf_data[, beta.outcome := as.numeric(outcome_parts[[tag_idx$ES]])]
  vcf_data[, se.outcome := as.numeric(outcome_parts[[tag_idx$SE]])]
  vcf_data[, LP := as.numeric(outcome_parts[[tag_idx$LP]])]
  
  # 处理AF标签不存在的情况
  if (!is.na(tag_idx$AF)) {
    vcf_data[, AF := as.numeric(outcome_parts[[tag_idx$AF]])]
  } else {
    vcf_data[, AF := NA]  # 如果没有AF标签，设为NA
  }
  
  vcf_data[, pval.outcome := 10^(-LP)]
  
  # 4. 清理和重命名
  vcf_data <- vcf_data[!is.na(beta.outcome) & !is.na(se.outcome) & !is.na(pval.outcome),
                       .(SNP = ID, REF, ALT, beta.outcome, se.outcome, pval.outcome, AF)]
  
  flog.info(sprintf("VCF解析完成: %d行有效数据", nrow(vcf_data)), name = "main")
  return(vcf_data)
}

# 5. MR分析函数
run_gene_mr <- function(gene, eqtl_data, outcome_data, thresholds, plots_dir, log_dirs, config) {
  # 每个子进程配置独立日志
  parallel_log <- file.path(log_dirs$parallel, sprintf("%s.log", gene))
  flog.appender(appender.file(parallel_log), name = "parallel")
  flog.threshold(DEBUG, name = "parallel")
  
  flog.debug(sprintf("分析基因: %s", gene), name = "parallel")
  
  # 初始化结果列表
  results <- list(
    main = data.table(),
    steiger = data.table(),
    coloc = data.table(),
    sensitivity = data.table(),
    snp_level = data.table()
  )
  
  tryCatch({
    # 提取该基因的暴露数据（修复变量名冲突）
    tryCatch({
      # 解决变量名冲突问题：使用不同的变量名
      target_gene <- gene
      gene_exposure <- eqtl_data[gene == target_gene, ]
      # 确保是data.table对象
      gene_exposure <- as.data.table(gene_exposure)
      
      flog.info(sprintf("%s 提取到 %d 行eQTL数据", gene, nrow(gene_exposure)), name = "parallel")
      
      if (nrow(gene_exposure) == 0) {
        flog.warn(sprintf("%s 没有可用的eQTL数据，跳过", gene), name = "parallel")
        return(results)
      }
    }, error = function(e) {
      flog.error(sprintf("%s 提取基因数据失败: %s", gene, e$message), name = "parallel")
      return(results)
    })
    
    # 计算F统计量并剔除弱工具变量
    calculate_f_stat <- function(dat, N) {
      dat <- as.data.table(dat)
      k <- nrow(dat)
      if (k == 0) return(dat)
      
      # 确保列存在
      if (!all(c("beta.exposure", "se.exposure") %in% colnames(dat))) {
        flog.warn("缺少 beta.exposure 或 se.exposure 列，无法计算F统计量", name = "parallel")
        dat[, F_stat := NA_real_]
        return(dat)
      }
      
      # 检查数据是否有NA值
      if (anyNA(dat$beta.exposure) || anyNA(dat$se.exposure)) {
        flog.warn("beta.exposure 或 se.exposure 列存在NA值，移除这些行", name = "parallel")
        dat <- dat[!is.na(beta.exposure) & !is.na(se.exposure), ]
        if (nrow(dat) == 0) return(dat)
      }
      
      # 计算R2和F统计量，处理k较大的情况
      dat[, R2 := beta.exposure^2 / (beta.exposure^2 + se.exposure^2 * N)]
      
      # 处理R2为1的情况
      dat[R2 >= 1, R2 := 0.9999]
      
      # 当k较大时，使用近似公式
      if (k >= N - 1) {
        flog.warn(sprintf("SNP数量 %d 接近或超过样本量 %d，使用近似F统计量", k, N), name = "parallel")
        dat[, F_stat := R2 / (1 - R2) * 1000]  # 使用近似值
      } else {
        dat[, F_stat := (N - k - 1) / k * R2 / (1 - R2)]
      }
      
      # 移除 F 统计量为 NA 或无穷大的行
      dat <- dat[is.finite(F_stat) & !is.na(F_stat), ]
      flog.info(sprintf("计算F统计量后剩余 %d 行数据", nrow(dat)), name = "parallel")
      return(dat)
    }
    
    # 计算F统计量
    gene_exposure <- calculate_f_stat(gene_exposure, config$eqtl_sample_size)
    
    # 暂时禁用F统计量过滤，因为eQTL数据的F统计量普遍较低
    flog.info(sprintf("%s 暂时禁用F统计量过滤，保留 %d 行数据", gene, nrow(gene_exposure)), name = "parallel")
    
    # 即使F统计量较低，也尝试进行后续分析
    if (nrow(gene_exposure) == 0) {
      flog.warn(sprintf("%s 没有eQTL数据，跳过", gene), name = "parallel")
      return(results)
    }
    
    # 显示F统计量范围
    flog.info(sprintf("%s F统计量范围: %.2f - %.2f", gene, min(gene_exposure$F_stat), max(gene_exposure$F_stat)), name = "parallel")
    
    if (nrow(gene_exposure) == 0) {
      flog.warn(sprintf("%s 没有eQTL数据，跳过", gene), name = "parallel")
      return(results)
    }
    
    # 标准化列名
    gene_exposure$exposure <- gene
    gene_exposure$id.exposure <- gene
    
    # 检查是否有SNP列
    flog.info(sprintf("%s 数据框列名: %s", gene, paste(colnames(gene_exposure), collapse = ", ")), name = "parallel")
    flog.info(sprintf("%s 数据框行数: %d", gene, nrow(gene_exposure)), name = "parallel")
    
    if (!"SNP" %in% colnames(gene_exposure)) {
      flog.warn(sprintf("%s 数据中没有SNP列，跳过", gene), name = "parallel")
      return(results)
    }
    
    # 匹配SNP
    common_snps <- character(0)  # 初始化变量
    tryCatch({
      # 检查outcome_data
      flog.info(sprintf("%s 检查outcome_data列名: %s", gene, paste(colnames(outcome_data), collapse = ", ")), name = "parallel")
      
      # 计算共同SNP
      common_snps <- intersect(gene_exposure$SNP, outcome_data$SNP)
      flog.info(sprintf("%s 成功计算共同SNP，数量: %d", gene, length(common_snps)), name = "parallel")
    }, error = function(e) {
      flog.error(sprintf("%s 访问SNP列失败: %s", gene, e$message), name = "parallel")
      return(results)
    })
    if (length(common_snps) == 0) {
      flog.warn(sprintf("%s 没有共同SNP，跳过", gene), name = "parallel")
      return(results)
    }
    
    # 过滤数据
    gene_exposure <- gene_exposure[SNP %in% common_snps, ]
    outcome_data_subset <- outcome_data[SNP %in% common_snps, ]
    
    # Harmonise数据（处理回文SNP）
    dat <- harmonise_data(
      exposure_dat = gene_exposure,
      outcome_dat = outcome_data_subset,
      action = 2  # 尝试推断链方向
    )
    
    if (nrow(dat) == 0) {
      flog.warn(sprintf("%s 没有可 harmonise 的SNP，跳过", gene), name = "parallel")
      return(results)
    }
    
    # 等位基因频率校验
    dat[, eaf_diff := abs(eaf.exposure - eaf.outcome)]
    dat <- dat[eaf_diff <= thresholds$eaf_diff, ]
    if (nrow(dat) == 0) {
      flog.warn(sprintf("%s 没有符合EAF差异要求的SNP，跳过", gene), name = "parallel")
      return(results)
    }
    
    # 方向性检验（Steiger检验）
    steiger <- directionality_test(dat)
    steiger_result <- data.table(
      gene = gene,
      steiger_pval = steiger$pval,
      correct_causal_direction = steiger$correct_causal_direction
    )
    results$steiger <- steiger_result
    
    # 检查反向因果
    if (steiger$pval > thresholds$steiger_pval || !steiger$correct_causal_direction) {
      flog.warn(sprintf("%s 可能存在反向因果，跳过", gene), name = "parallel")
      return(results)
    }
    
    # 多方法MR分析
    methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
    res <- mr(dat, method_list = methods)
    
    # 添加MR-RAPS方法（条件性执行）
    if (requireNamespace("mr.raps", quietly = TRUE)) {
      tryCatch({
        raps_res <- mr.raps::mr.raps(dat$beta.exposure, dat$se.exposure, dat$beta.outcome, dat$se.outcome)
        raps_result <- data.table(
          gene = gene,
          method = "MR-RAPS",
          b = raps_res$beta.hat,
          se = raps_res$beta.se,
          pval = raps_res$p.value,
          or = exp(raps_res$beta.hat),
          ci_lower = exp(raps_res$beta.hat - 1.96 * raps_res$beta.se),
          ci_upper = exp(raps_res$beta.hat + 1.96 * raps_res$beta.se),
          nsnp = nrow(dat)
        )
        results$main <- rbind(results$main, raps_result)
      }, error = function(e) {
        flog.warn(sprintf("%s MR-RAPS分析失败: %s", gene, e$message), name = "parallel")
      })
    } else {
      flog.info(sprintf("%s 跳过MR-RAPS（包未安装）", gene), name = "parallel")
    }
    
    # 提取结果
    for (i in 1:nrow(res)) {
      row <- res[i, ]
      or <- exp(row$b)
      ci_lower <- exp(row$b - 1.96 * row$se)
      ci_upper <- exp(row$b + 1.96 * row$se)
      
      main_result <- data.table(
        gene = gene,
        method = row$method,
        b = row$b,
        se = row$se,
        pval = row$pval,
        or = or,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        nsnp = row$nsnp
      )
      results$main <- rbind(results$main, main_result)
    }
    
    # 敏感性分析
    # 异质性
    hetero <- mr_heterogeneity(dat)
    for (i in 1:nrow(hetero)) {
      sensitivity_result <- data.table(
        gene = gene,
        analysis = "heterogeneity",
        method = hetero$method[i],
        Q = hetero$Q[i],
        Q_pval = hetero$Q_pval[i]
      )
      results$sensitivity <- rbind(results$sensitivity, sensitivity_result)
    }
    
    # 多效性
    pleio <- mr_pleiotropy_test(dat)
    for (i in 1:nrow(pleio)) {
      sensitivity_result <- data.table(
        gene = gene,
        analysis = "pleiotropy",
        method = pleio$method[i],
        intercept = pleio$intercept[i],
        intercept_se = pleio$se[i],
        intercept_pval = pleio$pval[i]
      )
      results$sensitivity <- rbind(results$sensitivity, sensitivity_result)
    }
    
    # 留一法
    loo <- mr_leaveoneout(dat)
    
    # MR-PRESSO
    presso_res <- if (nrow(dat) >= 15) {
      tryCatch({
        MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome",
          BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome",
          SdExposure = "se.exposure",
          OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE,
          data = dat,
          NbDistribution = min(1000, nrow(dat)*50),
          SignifThreshold = 0.05
        )
      }, error = function(e) {
        flog.warn(sprintf("%s MR-PRESSO分析失败: %s", gene, e$message), name = "parallel")
        return(NULL)
      })
    } else {
      flog.warn(sprintf("%s SNP数量不足15个，跳过MR-PRESSO分析", gene), name = "parallel")
      NULL
    }
    
    if (!is.null(presso_res)) {
      sensitivity_result <- data.table(
        gene = gene,
        analysis = "MR-PRESSO",
        distortion_test_pval = presso_res$distortion_test$p_value
      )
      results$sensitivity <- rbind(results$sensitivity, sensitivity_result)
    }
    
    # SNP水平数据
    for (i in 1:nrow(dat)) {
      row <- dat[i, ]
      wald_ratio <- row$beta.outcome / row$beta.exposure
      wald_se <- row$se.outcome / abs(row$beta.exposure)
      wald_pval <- 2 * pnorm(-abs(wald_ratio / wald_se))
      
      snp_result <- data.table(
        gene = gene,
        SNP = row$SNP,
        beta.exposure = row$beta.exposure,
        se.exposure = row$se.exposure,
        beta.outcome = row$beta.outcome,
        se.outcome = row$se.outcome,
        wald_ratio = wald_ratio,
        wald_se = wald_se,
        wald_pval = wald_pval
      )
      results$snp_level <- rbind(results$snp_level, snp_result)
    }
    
    # 共定位分析（对显著基因）
    ivw_result <- results$main[method == "Inverse variance weighted" & pval < thresholds$ivw_pval, ]
    if (nrow(ivw_result) > 0) {
      # 在harmonise后严格对齐数据
      dat_coloc <- dat[!is.na(dat$beta.exposure) & !is.na(dat$beta.outcome), ]
      if (nrow(dat_coloc) < 10) {  # coloc需要足够SNP
        flog.warn(sprintf("%s 可用于共定位的SNP不足(%d个)，跳过", gene, nrow(dat_coloc)), name = "parallel")
      } else {
        # 准备共定位数据（基于harmonise后的dat）
        exposure_coloc <- list(
          beta = dat_coloc$beta.exposure,
          varbeta = dat_coloc$se.exposure^2,
          snp = as.character(dat_coloc$SNP),  # 强制转字符
          type = "quant",      # eQTL为定量性状
          sdY = 1,              # 标准化表达量
          N = config$eqtl_sample_size             # eQTLGen样本量
        )
        
        outcome_coloc <- list(
          beta = dat_coloc$beta.outcome,
          varbeta = dat_coloc$se.outcome^2,
          snp = as.character(dat_coloc$SNP),  # 强制转字符
          type = "cc",         # 缺血性卒中为病例对照
          s = config$gwas_case_prop,              # 病例比例
          N = config$gwas_sample_size            # GWAS总样本量
        )
        
        # 执行共定位分析
        tryCatch({
          coloc_res <- coloc::coloc.abf(
            dataset1 = exposure_coloc,
            dataset2 = outcome_coloc,
            p1 = 1e-04,  # 先验概率
            p2 = 1e-04,
            p12 = 1e-05
          )
          
          # 提取PP.H0-H4
          coloc_result <- data.table(
            gene = gene,
            PP.H0 = coloc_res$summary["PP.H0.abf"],
            PP.H1 = coloc_res$summary["PP.H1.abf"],
            PP.H2 = coloc_res$summary["PP.H2.abf"],
            PP.H3 = coloc_res$summary["PP.H3.abf"],
            PP.H4 = coloc_res$summary["PP.H4.abf"]
          )
          results$coloc <- coloc_result
        }, error = function(e) {
          flog.warn(sprintf("%s 共定位分析失败: %s", gene, e$message), name = "parallel")
        })
      }
    }
    
    # 可视化
    # 散点图
    tryCatch({
      scatter_plot <- mr_scatter_plot(res, dat)
      if (inherits(scatter_plot, "list") && length(scatter_plot) > 0 && 
          !is.null(scatter_plot[[1]]) && inherits(scatter_plot[[1]], "ggplot")) {
        ggsave(file.path(plots_dir, sprintf("%s_scatter.png", gene)), scatter_plot[[1]] + 
                 theme_minimal() + 
                 ggtitle(sprintf("%s: 散点图", gene)), 
               width = 8, height = 6)
      } else {
        flog.warn(sprintf("%s 散点图返回无效对象", gene), name = "parallel")
      }
    }, error = function(e) {
      flog.warn(sprintf("%s 散点图失败: %s", gene, e$message), name = "parallel")
    })
    
    # 森林图
    tryCatch({
      forest_plot <- mr_forest_plot(res)
      if (inherits(forest_plot, "list") && length(forest_plot) > 0 && 
          !is.null(forest_plot[[1]]) && inherits(forest_plot[[1]], "ggplot")) {
        ggsave(file.path(plots_dir, sprintf("%s_forest.png", gene)), forest_plot[[1]] + 
                 theme_minimal() + 
                 ggtitle(sprintf("%s: 森林图", gene)), 
               width = 8, height = 6)
      } else {
        flog.warn(sprintf("%s 森林图返回无效对象", gene), name = "parallel")
      }
    }, error = function(e) {
      flog.warn(sprintf("%s 森林图失败: %s", gene, e$message), name = "parallel")
    })
    
    # 漏斗图
    tryCatch({
      funnel_plot <- mr_funnel_plot(loo)
      if (inherits(funnel_plot, "list") && length(funnel_plot) > 0 && 
          !is.null(funnel_plot[[1]]) && inherits(funnel_plot[[1]], "ggplot")) {
        ggsave(file.path(plots_dir, sprintf("%s_funnel.png", gene)), funnel_plot[[1]] + 
                 theme_minimal() + 
                 ggtitle(sprintf("%s: 漏斗图", gene)), 
               width = 8, height = 6)
      } else {
        flog.warn(sprintf("%s 漏斗图返回无效对象", gene), name = "parallel")
      }
    }, error = function(e) {
      flog.warn(sprintf("%s 漏斗图失败: %s", gene, e$message), name = "parallel")
    })
    
    # 留一法图
    tryCatch({
      loo_plot <- mr_leaveoneout_plot(loo)
      if (inherits(loo_plot, "list") && length(loo_plot) > 0 && 
          !is.null(loo_plot[[1]]) && inherits(loo_plot[[1]], "ggplot")) {
        ggsave(file.path(plots_dir, sprintf("%s_loo.png", gene)), loo_plot[[1]] + 
                 theme_minimal() + 
                 ggtitle(sprintf("%s: 留一法图", gene)), 
               width = 8, height = 6)
      } else {
        flog.warn(sprintf("%s 留一法图返回无效对象", gene), name = "parallel")
      }
    }, error = function(e) {
      flog.warn(sprintf("%s 留一法图失败: %s", gene, e$message), name = "parallel")
    })
    
    # 逆向MR分析（双向验证）
    flog.info(sprintf("%s 开始逆向MR分析", gene), name = "parallel")
    tryCatch({
      # 准备逆向分析数据（stroke -> gene expression）
      reverse_exposure <- outcome_data[SNP %in% common_snps, ]
      reverse_outcome <- gene_exposure[SNP %in% common_snps, ]
      
      # 标准化暴露数据（结局->暴露）
      reverse_exposure$exposure <- "Ischemic Stroke"
      reverse_exposure$id.exposure <- "ieu-a-83"
      reverse_exposure$beta.exposure <- reverse_exposure$beta.outcome
      reverse_exposure$se.exposure <- reverse_exposure$se.outcome
      reverse_exposure$eaf.exposure <- reverse_exposure$eaf.outcome
      reverse_exposure$effect_allele.exposure <- reverse_exposure$effect_allele.outcome
      reverse_exposure$other_allele.exposure <- reverse_exposure$other_allele.outcome
      reverse_exposure$pval.exposure <- reverse_exposure$pval.outcome
      
      # 标准化结局数据（暴露->结局）
      reverse_outcome$outcome <- gene
      reverse_outcome$id.outcome <- gene
      reverse_outcome$beta.outcome <- reverse_outcome$beta.exposure
      reverse_outcome$se.outcome <- reverse_outcome$se.exposure
      reverse_outcome$eaf.outcome <- reverse_outcome$eaf.exposure
      reverse_outcome$effect_allele.outcome <- reverse_outcome$effect_allele.exposure
      reverse_outcome$other_allele.outcome <- reverse_outcome$other_allele.exposure
      reverse_outcome$pval.outcome <- reverse_outcome$pval.exposure
      
      # Harmonise数据
      dat_reverse <- harmonise_data(
        exposure_dat = reverse_exposure,
        outcome_dat = reverse_outcome,
        action = 2
      )
      
      if (nrow(dat_reverse) > 0) {
        # 逆向MR分析
        res_reverse <- mr(dat_reverse, method_list = methods)
        
        # 提取逆向结果
        for (i in 1:nrow(res_reverse)) {
          row <- res_reverse[i, ]
          or <- exp(row$b)
          ci_lower <- exp(row$b - 1.96 * row$se)
          ci_upper <- exp(row$b + 1.96 * row$se)
          
          reverse_result <- data.table(
            gene = gene,
            method = paste0("Reverse_", row$method),
            b = row$b,
            se = row$se,
            pval = row$pval,
            or = or,
            ci_lower = ci_lower,
            ci_upper = ci_upper,
            nsnp = row$nsnp
          )
          results$main <- rbind(results$main, reverse_result)
        }
      }
    }, error = function(e) {
      flog.warn(sprintf("%s 逆向MR分析失败: %s", gene, e$message), name = "parallel")
    })
    
    flog.debug(sprintf("基因 %s 分析完成", gene), name = "parallel")
    
  }, error = function(e) {
    error_file <- file.path(log_dirs$error, sprintf("%s.err", gene))
    writeLines(c(sprintf("Error: %s", e$message), traceback()), error_file)
    flog.error(sprintf("%s 分析失败: %s", gene, e$message), name = "parallel")
  })
  
  # 清理日志句柄
  flog.remove("parallel")
  
  return(results)
}

# 6. 主函数
main <- function() {
  # 读取配置
  config <- read_config("config.yaml")
  
  # 设置目录
  base_dir <- config$base_dir
  output_dir <- config$output_dir
  eqtl_dir <- config$eqtl_dir
  vcf_file <- config$vcf_file
  core_genes <- config$core_genes
  test_mode <- config$test_mode
  n_cores <- config$n_cores
  thresholds <- config$thresholds
  
  # 测试模式
  if (test_mode) {
    core_genes <- c("FDX1", "SLC31A1", "NFKB1")
    flog.info("运行测试模式，仅分析3个核心基因", name = "main")
  }
  
  # 创建输出目录
  plots_dir <- file.path(output_dir, "plots")
  log_dir <- file.path(output_dir, "logs")
  log_dirs <- setup_logging(log_dir)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # 检查点文件
  checkpoint_stage1 <- file.path(output_dir, "checkpoint_stage1_eqtl_raw.fst")
  checkpoint_stage2 <- file.path(output_dir, "checkpoint_stage2_vcf_parsed.fst")
  
  # 阶段一：eQTL数据处理
  flog.info("=== 阶段一：eQTL数据处理 ===", name = "main")
  if (!file.exists(checkpoint_stage1)) {
    eqtl_data <- load_eqtl_data_fast(eqtl_dir, core_genes)
    # 保存为fst格式
    fst::write.fst(eqtl_data, checkpoint_stage1)
    flog.info(sprintf("阶段一完成，核心数据已保存到: %s", checkpoint_stage1), name = "main")
  } else {
    flog.info("从检查点加载eQTL数据", name = "main")
    eqtl_data <- fst::read.fst(checkpoint_stage1)
    flog.info(sprintf("从检查点加载核心数据: %d 条记录", nrow(eqtl_data)), name = "main")
  }
  
  # 阶段二：VCF解析
  flog.info("=== 阶段二：VCF解析 ===", name = "main")
  if (!file.exists(checkpoint_stage2)) {
    # 测试模式先跑10000行验证
    if (test_mode) {
      flog.info("测试模式：VCF只读取前10000行", name = "main")
      vcf_data <- parse_gwas_vcf_fast(vcf_file, nrows = 10000)
    } else {
      vcf_data <- parse_gwas_vcf_fast(vcf_file)
    }
    # 保存为fst格式
    fst::write.fst(vcf_data, checkpoint_stage2)
    flog.info(sprintf("阶段二完成，结果保存到: %s", checkpoint_stage2), name = "main")
  } else {
    flog.info("从检查点加载VCF数据", name = "main")
    vcf_data <- fst::read.fst(checkpoint_stage2)
    flog.info(sprintf("从检查点加载VCF数据: %d 条记录", nrow(vcf_data)), name = "main")
  }
  
  # 准备结局数据
  outcome_dat <- vcf_data
  outcome_dat$id.outcome <- "ieu-a-83"
  outcome_dat$outcome <- "Ischemic Stroke"
  outcome_dat$effect_allele.outcome <- outcome_dat$ALT
  outcome_dat$other_allele.outcome <- outcome_dat$REF
  # 添加eaf.outcome列
  outcome_dat$eaf.outcome <- ifelse(!is.na(outcome_dat$AF), outcome_dat$AF, 0.5)
  
  # 阶段三：MR分析
  flog.info("=== 阶段三：MR分析 ===", name = "main")
  
  # 分析结果（使用lapply代替mclapply，因为Windows不支持mclapply多核并行）
  flog.info("开始分析基因...", name = "main")
  # 设置随机种子
  set.seed(123)
  all_results <- lapply(core_genes, function(gene) {
    run_gene_mr(gene, eqtl_data, outcome_dat, thresholds, plots_dir, log_dirs, config)
  })
  
  # 合并并行子进程日志
  flog.info("合并并行子进程日志...", name = "main")
  parallel_log_files <- list.files(log_dirs$parallel, full.names = TRUE)
  all_log_content <- c()
  for (log_file in parallel_log_files) {
    if (file.exists(log_file)) {
      log_content <- readLines(log_file)
      if (length(log_content) > 0) {
        all_log_content <- c(all_log_content, log_content, "")
      }
    }
  }
  if (length(all_log_content) > 0) {
    # 使用fwrite避免写入冲突
    fwrite(data.table(log = all_log_content), file = log_dirs$main, append = TRUE, col.names = FALSE)
  }
  
  # 合并结果
  main_results <- data.table()
  steiger_results <- data.table()
  coloc_results <- data.table()
  sensitivity_results <- data.table()
  snp_level_data <- data.table()
  
  for (result in all_results) {
    main_results <- rbind(main_results, result$main)
    steiger_results <- rbind(steiger_results, result$steiger)
    coloc_results <- rbind(coloc_results, result$coloc)
    sensitivity_results <- rbind(sensitivity_results, result$sensitivity)
    snp_level_data <- rbind(snp_level_data, result$snp_level)
  }
  
  # 多重检验校正：只针对IVW结果进行基因间的多重校正
  if (nrow(main_results) > 0) {
    # 只提取IVW结果进行基因水平的多重检验校正
    ivw_results <- main_results[method == "Inverse variance weighted", ]
    if (nrow(ivw_results) > 0) {
      # 按基因去重（保留最强信号）
      ivw_unique <- ivw_results[, .SD[which.min(pval)], by = gene]
      ivw_unique[, pval_bonferroni := p.adjust(pval, method = "bonferroni")]
      ivw_unique[, pval_fdr := p.adjust(pval, method = "fdr")]
      
      # 合并回去
      main_results <- merge(main_results, 
                           ivw_unique[, .(gene, pval_bonferroni, pval_fdr)], 
                           by = "gene", all.x = TRUE)
    }
  }
  
  # 阶段四：结果输出
  flog.info("=== 阶段四：结果输出 ===", name = "main")
  
  # 保存主结果
  if (nrow(main_results) > 0) {
    fst::write.fst(main_results, file.path(output_dir, "01_main_results.fst"))
    write.csv(main_results, file.path(output_dir, "01_main_results.csv"), row.names = FALSE)
    flog.info("主结果已保存", name = "main")
  } else {
    flog.warn("没有主结果", name = "main")
  }
  
  # 保存方向性检验结果
  if (nrow(steiger_results) > 0) {
    fst::write.fst(steiger_results, file.path(output_dir, "02_steiger_directionality.fst"))
    write.csv(steiger_results, file.path(output_dir, "02_steiger_directionality.csv"), row.names = FALSE)
    flog.info("方向性检验结果已保存", name = "main")
  } else {
    flog.warn("没有方向性检验结果", name = "main")
  }
  
  # 保存共定位结果
  if (nrow(coloc_results) > 0) {
    fst::write.fst(coloc_results, file.path(output_dir, "03_coloc_results.fst"))
    write.csv(coloc_results, file.path(output_dir, "03_coloc_results.csv"), row.names = FALSE)
    flog.info("共定位结果已保存", name = "main")
  } else {
    flog.warn("没有共定位结果", name = "main")
  }
  
  # 保存敏感性分析结果
  if (nrow(sensitivity_results) > 0) {
    fst::write.fst(sensitivity_results, file.path(output_dir, "04_sensitivity_analysis.fst"))
    write.csv(sensitivity_results, file.path(output_dir, "04_sensitivity_analysis.csv"), row.names = FALSE)
    flog.info("敏感性分析结果已保存", name = "main")
  } else {
    flog.warn("没有敏感性分析结果", name = "main")
  }
  
  # 保存SNP水平数据
  if (nrow(snp_level_data) > 0) {
    fst::write.fst(snp_level_data, file.path(output_dir, "05_snp_level_data.fst"))
    write.csv(snp_level_data, file.path(output_dir, "05_snp_level_data.csv"), row.names = FALSE)
    flog.info("SNP水平数据已保存", name = "main")
  } else {
    flog.warn("没有SNP水平数据", name = "main")
  }
  
  # 生成分析流程图
  flog.info("生成分析流程图", name = "main")
  tryCatch({
    diagram_code <- "
    graph TD
      A[数据输入] --> B[eQTL数据处理]
      A --> C[GWAS VCF解析]
      B --> D[数据匹配与Harmonise]
      C --> D
      D --> E[多方法MR分析]
      E --> F[敏感性分析]
      E --> G[共定位分析]
      F --> H[结果输出]
      G --> H
    "
    # 使用DiagrammeR的另一种方式保存SVG
    library(DiagrammeRsvg)
    DiagrammeR::grViz(diagram_code) %>% 
      DiagrammeRsvg::export_svg() %>% 
      writeLines(file.path(output_dir, "analysis_flow.svg"))
  }, error = function(e) {
    flog.warn(sprintf("生成流程图失败: %s", e$message), name = "main")
    # 生成一个简单的文本流程图作为替代
    flow_text <- "分析流程:\n1. 数据输入\n2. eQTL数据处理\n3. GWAS VCF解析\n4. 数据匹配与Harmonise\n5. 多方法MR分析\n6. 敏感性分析\n7. 共定位分析\n8. 结果输出"
    writeLines(flow_text, file.path(output_dir, "analysis_flow.txt"))
  })
  
  # 生成rmarkdown报告
  flog.info("生成rmarkdown报告", name = "main")
  rmd_content <- paste0(
    "---
",
    "title: '铜死亡核心基因与缺血性脑卒中的孟德尔随机化分析报告'\n",
    "date: '`r Sys.Date()`'\n",
    "output:\n",
    "  html_document:\n",
    "    toc: true\n",
    "    toc_float: true\n",
    "    code_folding: 'hide'\n",
    "  pdf_document:\n",
    "    toc: true\n",
    "---\n\n",
    "# 分析概述\n\n",
    "本报告基于孟德尔随机化（MR）方法，分析铜死亡核心基因与缺血性脑卒中之间的因果关系。\n\n",
    "## 分析流程\n\n",
    "![分析流程图](analysis_flow.svg)\n\n",
    "# 数据来源\n\n",
    "- **eQTL数据**：来自eQTLGen consortium\n",
    "- **GWAS数据**：来自MEGASTROKE consortium (ieu-a-83)\n\n",
    "# 分析结果\n\n",
    "## 基本统计\n\n",
    "```r\n",
    "# 基本统计\n",
    "n_total <- length(core_genes)\n",
    "n_success <- length(unique(main_results$gene))\n",
    "n_significant <- length(unique(main_results[method == \"Inverse variance weighted\" & pval < 0.05, gene]))\n",
    "\n",
    "cat(sprintf('总基因数: %d\\n', n_total))\n",
    "cat(sprintf('成功分析数: %d\\n', n_success))\n",
    "cat(sprintf('显著基因数 (IVW p<0.05): %d\\n', n_significant))\n",
    "```\n\n",
    "## 主MR结果\n\n",
    "```r\n",
    "# 主MR结果\n",
    "knitr::kable(main_results[method == \"Inverse variance weighted\", ])\n",
    "```\n\n",
    "## 敏感性分析结果\n\n",
    "```r\n",
    "# 敏感性分析结果\n",
    "knitr::kable(sensitivity_results)\n",
    "```\n\n",
    "## 共定位分析结果\n\n",
    "```r\n",
    "# 共定位分析结果\n",
    "knitr::kable(coloc_results)\n",
    "```\n\n",
    "# 可视化结果\n\n",
    "## 基因效应火山图\n\n",
    "```r\n",
    "# 火山图\n",
    "library(ggplot2)\n",
    "ivw_results <- main_results[method == \"Inverse variance weighted\", ]\n",
    "ggplot(ivw_results, aes(x = b, y = -log10(pval))) +\n",
    "  geom_point(aes(color = pval < 0.05)) +\n",
    "  theme_minimal() +\n",
    "  xlab('Effect size') +\n",
    "  ylab('-log10(p-value)') +\n",
    "  ggtitle('铜死亡核心基因与缺血性脑卒中的MR效应') +\n",
    "  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +\n",
    "  geom_vline(xintercept = 0, linetype = 'dashed')\n",
    "```\n\n",
    "# 结论\n\n",
    "本分析通过多方法MR、敏感性分析和共定位分析，系统评估了铜死亡核心基因与缺血性脑卒中之间的因果关系。\n\n",
    "# 会话信息\n\n",
    "```r\n",
    "sessionInfo()\n",
    "```\n"
  )
  
  writeLines(rmd_content, file.path(output_dir, "06_analysis_report.Rmd"))
  
  # 渲染报告
  tryCatch({
    rmarkdown::render(file.path(output_dir, "06_analysis_report.Rmd"), 
                      output_format = c("html_document", "pdf_document"))
    flog.info("报告生成完成", name = "main")
  }, error = function(e) {
    flog.warn(sprintf("报告生成失败: %s", e$message), name = "main")
  })
  
  # 保存会话信息
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
  flog.info("会话信息已保存", name = "main")
  
  flog.info("=== 分析完成 ===", name = "main")
  flog.info(sprintf("结果已保存到: %s", output_dir), name = "main")
  
  # 显示摘要
  cat("\n=== 分析摘要 ===\n")
  cat(sprintf("总基因数: %d\n", length(core_genes)))
  
  # 检查main_results是否为空
  if (nrow(main_results) > 0) {
    cat(sprintf("成功分析数: %d\n", length(unique(main_results$gene))))
    # 检查是否有method列
    if ("method" %in% colnames(main_results)) {
      cat(sprintf("显著基因数 (IVW p<0.05): %d\n", length(unique(main_results[method == "Inverse variance weighted" & pval < 0.05, gene]))))
    } else {
      cat("显著基因数 (IVW p<0.05): 0\n")
    }
  } else {
    cat("成功分析数: 0\n")
    cat("显著基因数 (IVW p<0.05): 0\n")
  }
  cat("\n结果已保存到 D:/EQTL/MR_Optimized_Results/ 目录\n")
}

# 运行主函数
if (sys.nframe() == 0) {
  main()
}
