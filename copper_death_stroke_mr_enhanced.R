#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 铜死亡核心基因与缺血性脑卒中的工业级孟德尔随机化分析
# 实现多方法MR分析、方向性验证、共定位分析、敏感性分析等功能

# 设置CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 包加载和安装
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("安装包: %s", pkg))
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
}

# 安装必要的包
install_if_missing("data.table")      # 大数据处理
install_if_missing("TwoSampleMR")     # MR分析
install_if_missing("coloc")           # 共定位分析
install_if_missing("MRPRESSO")        # 异常值检测
install_if_missing("futile.logger")   # 日志系统
install_if_missing("progress")        # 进度条
install_if_missing("readxl")          # 读取Excel文件
install_if_missing("ggplot2")         # 可视化
install_if_missing("stringr")         # 字符串处理

# 设置随机种子
set.seed(123)

# 目录设置
base_dir <- "D:/EQTL"
output_dir <- "D:/EQTL/MR_Enhanced_Results"
plots_dir <- file.path(output_dir, "plots")
log_dir <- file.path(output_dir, "logs")

# 创建输出目录
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

# 核心基因列表
core_genes <- c("FDX1", "DLAT", "DLD", "PDHB", "LIPT1", "PDHX", "SLC31A1", 
                "ATP7B", "ATP7A", "ATOX1", "COMMD1", "MT2A", "NFKB1", "RELA", 
                "STAT1", "STAT3", "STAT5A", "ICAM1", "CCL2", "IL6", "TGFB1", 
                "PTGS2", "HMOX1", "SOD2", "FABP3", "ATF4", "BRD4")

# 测试模式
test_mode <- FALSE
if (test_mode) {
  core_genes <- c("FDX1", "SLC31A1", "NFKB1")
  flog.info("运行测试模式，仅分析3个核心基因")
}

# 检查点文件
checkpoint_stage1 <- file.path(output_dir, "checkpoint_stage1_eqtl_raw.RData")
checkpoint_stage2 <- file.path(output_dir, "checkpoint_stage2_vcf_parsed.RData")
checkpoint_stage3 <- file.path(output_dir, "checkpoint_stage3_harmonised.RData")

# 日志配置
flog.appender(appender.file(file.path(log_dir, "mr_analysis.log")))
flog.threshold(INFO)
flog.info("开始孟德尔随机化分析")
flog.info(sprintf("分析基因数量: %d", length(core_genes)))

# 阶段一：eQTL数据处理
flog.info("=== 阶段一：eQTL数据处理 ===")
eqtl_dir <- file.path(base_dir, "共定位顺式")
eqtl_files <- list.files(eqtl_dir, pattern = "eqtlgen_ieu_1mbcis.*\\.csv$", full.names = TRUE)
eqtl_files <- eqtl_files[!grepl('\\.downloading$', eqtl_files)]
flog.info(sprintf("找到 %d 个eQTL分块文件", length(eqtl_files)))

# 合并eQTL数据 - 使用lapply + rbindlist提高效率
flog.info("合并eQTL数据...")
eqtl_data <- rbindlist(lapply(eqtl_files, function(file) {
  flog.info(sprintf("读取文件: %s", basename(file)))
  fread(file, showProgress = FALSE, data.table = TRUE)
}))

flog.info(sprintf("共读取 %d 条eQTL记录", nrow(eqtl_data)))

# 动态识别列名
col_names <- colnames(eqtl_data)
flog.info(sprintf("eQTL数据列名: %s", paste(col_names, collapse = ", ")))

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

flog.info(sprintf("使用列名 - P值: %s, Beta: %s, SE: %s, 基因: %s, SNP: %s", 
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

# 筛选核心基因（保留所有p值记录，供敏感性分析复用）
eqtl_data_core <- eqtl_data[gene %in% core_genes, ]
flog.info(sprintf("核心基因相关记录: %d 条", nrow(eqtl_data_core)))

# 保存阶段一结果
save(eqtl_data_core, file = checkpoint_stage1)
flog.info(sprintf("阶段一完成，核心数据已保存到: %s", checkpoint_stage1))

# 阶段二：VCF解析
flog.info("=== 阶段二：VCF解析 ===")
vcf_file <- "D:/EQTL/ieu-a-83.vcf"
flog.info(sprintf("读取VCF文件: %s", vcf_file))

# 使用readLines读取VCF文件
vcf_lines <- readLines(vcf_file, encoding = "UTF-8")

# 找到header行
header_idx <- which(grepl("^#CHROM", vcf_lines))
if (length(header_idx) == 0) stop("无法找到VCF header行")

# 提取数据行（去除注释行）
data_lines <- vcf_lines[(header_idx + 1):length(vcf_lines)]

# 提取header
header_line <- vcf_lines[header_idx]
header <- strsplit(header_line, "\t")[[1]]

# 找到FORMAT列和ieu-a-83列
format_idx <- which(header == "FORMAT")
outcome_idx <- which(header == "ieu-a-83")

# 提取FORMAT标签
format_line <- data_lines[1]
format_parts <- strsplit(format_line, "\t")[[1]]
format_tags <- strsplit(format_parts[format_idx], ":")[[1]]
flog.info(sprintf("FORMAT标签: %s", paste(format_tags, collapse = ", ")))

# 创建标签到索引的映射
tag_map <- setNames(seq_along(format_tags), format_tags)

# 解析VCF数据 - 修正：提取AF作为eaf.outcome
parse_vcf_field <- function(field) {
  parts <- strsplit(field, ":")[[1]]
  result <- list()
  for (tag in c("ES", "SE", "LP", "AF")) {
    if (tag %in% names(tag_map)) {
      idx <- tag_map[[tag]]
      if (length(parts) >= idx) {
        result[[tag]] <- as.numeric(parts[idx])
      } else {
        result[[tag]] <- NA
      }
    } else {
      result[[tag]] <- NA
    }
  }
  # 计算P值：pval = 10^(-LP)
  if (!is.na(result$LP)) {
    result$P <- 10^(-result$LP)
  } else {
    result$P <- NA
  }
  return(result)
}

# 提取SNP和解析ieu-a-83列 - 修正：使用lapply + rbindlist提高效率
flog.info("解析VCF数据...")
vcf_data <- rbindlist(lapply(data_lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= outcome_idx) {
    snp <- parts[3]  # ID列
    ref <- parts[4]  # REF列
    alt <- parts[5]  # ALT列
    outcome_field <- parts[outcome_idx]
    parsed <- parse_vcf_field(outcome_field)
    
    data.table(
      SNP = snp,
      REF = ref,
      ALT = alt,
      BETA = parsed$ES,
      SE = parsed$SE,
      P = parsed$P,
      AF = parsed$AF  # 提取AF作为eaf.outcome
    )
  }
}))

# 剔除NA值
vcf_data <- vcf_data[!is.na(vcf_data$BETA) & !is.na(vcf_data$SE) & !is.na(vcf_data$P), ]
flog.info(sprintf("VCF文件读取完成，共 %d 行数据", nrow(vcf_data)))

# 统一列名
setnames(vcf_data, old = c("BETA", "SE", "P"), new = c("beta.outcome", "se.outcome", "pval.outcome"))

# 保存阶段二结果
save(vcf_data, file = checkpoint_stage2)
flog.info(sprintf("阶段二完成，结果保存到: %s", checkpoint_stage2))

# 阶段三：数据匹配和Harmonise
flog.info("=== 阶段三：数据匹配和Harmonise ===")

# 检查暴露数据中的SNP
snp_exposure <- unique(eqtl_data_core$SNP)
snp_exposure <- snp_exposure[!is.na(snp_exposure) & snp_exposure != ""]
flog.info(sprintf("暴露数据中有 %d 个唯一SNP", length(snp_exposure)))

# 检查结局数据中的SNP
snp_outcome <- unique(vcf_data$SNP)
snp_outcome <- snp_outcome[!is.na(snp_outcome) & snp_outcome != ""]
flog.info(sprintf("结局数据中有 %d 个唯一SNP", length(snp_outcome)))

# 匹配SNP
common_snps <- intersect(snp_exposure, snp_outcome)
flog.info(sprintf("共同SNP数量: %d", length(common_snps)))

if (length(common_snps) == 0) {
  flog.error("没有找到共同SNP，分析无法继续")
  stop("没有找到共同SNP，分析无法继续")
}

# 准备暴露数据
exposure_dat <- eqtl_data_core[eqtl_data_core$SNP %in% common_snps, ]

# 计算F统计量并剔除弱工具变量
exposure_dat[, F_stat := (beta.exposure / se.exposure)^2]
exposure_dat <- exposure_dat[F_stat >= 10, ]
flog.info(sprintf("剔除F<10后剩余 %d 条记录", nrow(exposure_dat)))

# 准备结局数据
outcome_dat <- vcf_data[vcf_data$SNP %in% common_snps, ]

# 添加必要的列
outcome_dat$id.outcome <- "ieu-a-83"
outcome_dat$outcome <- "Ischemic Stroke"
outcome_dat$effect_allele.outcome <- outcome_dat$ALT
outcome_dat$other_allele.outcome <- outcome_dat$REF
# 添加eaf.outcome列 - 修正：使用从VCF中提取的AF值
outcome_dat$eaf.outcome <- ifelse(!is.na(outcome_dat$AF), outcome_dat$AF, 0.5)  # AF缺失时使用0.5作为默认值

# 保存阶段三结果
save(exposure_dat, outcome_dat, common_snps, file = checkpoint_stage3)
flog.info(sprintf("阶段三完成，结果保存到: %s", checkpoint_stage3))

# 阶段四：MR分析
flog.info("=== 阶段四：MR分析 ===")

# 结果存储
main_results <- data.table()
steiger_results <- data.table()
coloc_results <- data.table()
sensitivity_results <- data.table()
snp_level_data <- data.table()

# 进度条
pb <- progress_bar$new(
  format = "[:bar] :percent (:current/:total) :elapsed",
  total = length(core_genes),
  clear = FALSE
)

for (gene in core_genes) {
  pb$tick()
  flog.info(sprintf("分析基因: %s", gene))
  
  tryCatch({
    # 提取该基因的暴露数据
    gene_exposure <- exposure_dat[exposure_dat$gene == gene, ]
    if (nrow(gene_exposure) == 0) {
      flog.warn(sprintf("%s 没有可用的eQTL数据，跳过", gene))
      next
    }
    
    # 标准化列名
    gene_exposure$exposure <- gene
    gene_exposure$id.exposure <- gene
    
    # Harmonise数据（处理回文SNP）
    dat <- harmonise_data(
      exposure_dat = gene_exposure,
      outcome_dat = outcome_dat,
      action = 2  # 尝试推断链方向
    )
    
    if (nrow(dat) == 0) {
      flog.warn(sprintf("%s 没有可 harmonise 的SNP，跳过", gene))
      next
    }
    
    # 等位基因频率校验
    dat[, eaf_diff := abs(eaf.exposure - eaf.outcome)]
    dat <- dat[eaf_diff <= 0.2, ]
    if (nrow(dat) == 0) {
      flog.warn(sprintf("%s 没有符合EAF差异要求的SNP，跳过", gene))
      next
    }
    
    # 方向性检验（Steiger检验）
    steiger <- directionality_test(dat)
    steiger_result <- data.table(
      gene = gene,
      steiger_pval = steiger$pval,
      correct_causal_direction = steiger$correct_causal_direction
    )
    steiger_results <- rbind(steiger_results, steiger_result)
    
    # 检查反向因果
    if (steiger$pval > 0.05 || !steiger$correct_causal_direction) {
      flog.warn(sprintf("%s 可能存在反向因果，跳过", gene))
      next
    }
    
    # 多方法MR分析
    methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
    res <- mr(dat, method_list = methods)
    
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
      main_results <- rbind(main_results, main_result)
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
      sensitivity_results <- rbind(sensitivity_results, sensitivity_result)
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
      sensitivity_results <- rbind(sensitivity_results, sensitivity_result)
    }
    
    # 留一法
    loo <- mr_leaveoneout(dat)
    
    # MR-PRESSO - 修正：使用命名参数并设置正确的参数
    presso_res <- tryCatch({
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = dat,
        NbDistribution = 1000,
        SignifThreshold = 0.05
      )
    }, error = function(e) {
      flog.warn(sprintf("%s MR-PRESSO分析失败: %s", gene, e$message))
      return(NULL)
    })
    
    if (!is.null(presso_res)) {
      sensitivity_result <- data.table(
        gene = gene,
        analysis = "MR-PRESSO",
        distortion_test_pval = presso_res$distortion_test$p_value
      )
      sensitivity_results <- rbind(sensitivity_results, sensitivity_result)
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
      snp_level_data <- rbind(snp_level_data, snp_result)
    }
    
    # 共定位分析（对显著基因） - 修正：使用harmonise后的数据并指定trait类型
if (nrow(res) > 0 && any(res$method == "Inverse variance weighted" & res$pval < 0.05)) {
  # 准备共定位数据（基于harmonise后的dat） - 文献依据：Giambartolomei et al. (2014, PLoS Genet)
  exposure_coloc <- list(
    beta = dat$beta.exposure,
    varbeta = dat$se.exposure^2,
    snp = dat$SNP,
    type = "quant",      # eQTL为定量性状
    sdY = 1              # 标准化表达量
  )
  
  outcome_coloc <- list(
    beta = dat$beta.outcome,
    varbeta = dat$se.outcome^2,
    snp = dat$SNP,
    type = "cc",         # 缺血性卒中为病例对照
    s = 0.5              # ieu-a-83病例比例 - 文献依据：Wang et al. (2020, Nat Commun)
  )
  
  # 执行共定位分析
  coloc_res <- coloc::coloc.abf(
    dataset1 = exposure_coloc,
    dataset2 = outcome_coloc,
    p1 = 1e-04,  # 先验概率 - 标准设置
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
  coloc_results <- rbind(coloc_results, coloc_result)
}
    
    # 可视化
    # 散点图
    tryCatch({
      scatter_plot <- mr_scatter_plot(res, dat)
      if (!is.null(scatter_plot) && length(scatter_plot) > 0 && !is.null(scatter_plot[[1]])) {
        ggsave(file.path(plots_dir, sprintf("%s_scatter.png", gene)), scatter_plot[[1]], width = 8, height = 6)
      }
    }, error = function(e) {
      flog.warn(sprintf("%s 散点图失败: %s", gene, e$message))
    })
    
    # 森林图
    tryCatch({
      forest_plot <- mr_forest_plot(res)
      if (!is.null(forest_plot) && length(forest_plot) > 0 && !is.null(forest_plot[[1]])) {
        ggsave(file.path(plots_dir, sprintf("%s_forest.png", gene)), forest_plot[[1]], width = 8, height = 6)
      }
    }, error = function(e) {
      flog.warn(sprintf("%s 森林图失败: %s", gene, e$message))
    })
    
    # 漏斗图
    tryCatch({
      funnel_plot <- mr_funnel_plot(loo)
      if (!is.null(funnel_plot) && length(funnel_plot) > 0 && !is.null(funnel_plot[[1]])) {
        ggsave(file.path(plots_dir, sprintf("%s_funnel.png", gene)), funnel_plot[[1]], width = 8, height = 6)
      }
    }, error = function(e) {
      flog.warn(sprintf("%s 漏斗图失败: %s", gene, e$message))
    })
    
    # 留一法图
    tryCatch({
      loo_plot <- mr_leaveoneout_plot(loo)
      if (!is.null(loo_plot) && length(loo_plot) > 0 && !is.null(loo_plot[[1]])) {
        ggsave(file.path(plots_dir, sprintf("%s_loo.png", gene)), loo_plot[[1]], width = 8, height = 6)
      }
    }, error = function(e) {
      flog.warn(sprintf("%s 留一法图失败: %s", gene, e$message))
    })
    
    flog.info(sprintf("%s 分析完成", gene))
    
  }, error = function(e) {
    flog.error(sprintf("%s 分析失败: %s", gene, e$message))
  })
}

# 阶段五：结果输出
flog.info("=== 阶段五：结果输出 ===")

# 保存主结果
if (nrow(main_results) > 0) {
  write.csv(main_results, file.path(output_dir, "01_main_results.csv"), row.names = FALSE)
  flog.info("主结果已保存")
} else {
  flog.warn("没有主结果")
}

# 保存方向性检验结果
if (nrow(steiger_results) > 0) {
  write.csv(steiger_results, file.path(output_dir, "02_steiger_directionality.csv"), row.names = FALSE)
  flog.info("方向性检验结果已保存")
} else {
  flog.warn("没有方向性检验结果")
}

# 保存共定位结果
if (nrow(coloc_results) > 0) {
  write.csv(coloc_results, file.path(output_dir, "03_coloc_results.csv"), row.names = FALSE)
  flog.info("共定位结果已保存")
} else {
  flog.warn("没有共定位结果")
}

# 保存敏感性分析结果
if (nrow(sensitivity_results) > 0) {
  write.csv(sensitivity_results, file.path(output_dir, "04_sensitivity_analysis.csv"), row.names = FALSE)
  flog.info("敏感性分析结果已保存")
} else {
  flog.warn("没有敏感性分析结果")
}

# 保存SNP水平数据
if (nrow(snp_level_data) > 0) {
  write.csv(snp_level_data, file.path(output_dir, "05_snp_level_data.csv"), row.names = FALSE)
  flog.info("SNP水平数据已保存")
} else {
  flog.warn("没有SNP水平数据")
}

# 生成统计报告
flog.info("=== 生成统计报告 ===")

# 统计摘要
n_total <- length(core_genes)
n_success <- length(unique(main_results$gene))
n_significant <- length(unique(main_results[method == "Inverse variance weighted" & pval < 0.05, gene]))

# 异质性显著基因
heterogeneity_genes <- unique(sensitivity_results[analysis == "heterogeneity" & Q_pval < 0.05, gene])

# 多效性显著基因
pleiotropy_genes <- unique(sensitivity_results[analysis == "pleiotropy" & intercept_pval < 0.05, gene])

# 共定位强证据基因
coloc_strong_genes <- unique(coloc_results[PP.H4 > 0.75, gene])

# 建议剔除的基因
rejected_genes <- unique(steiger_results[steiger_pval > 0.05 | !correct_causal_direction, gene])

# 生成Markdown报告
report_content <- paste0(
  "# 孟德尔随机化分析统计报告\n\n",
  "## 样本重叠警告\n",
  "- ieu-a-83（欧洲人群）与eQTLGen（欧洲人群）存在样本重叠\n",
  "- 已通过Steiger检验过滤可能的反向因果\n",
  "- 对Steiger检验通过的基因，steiger_pval可作为敏感性指标\n\n",
  "## 基本统计\n",
  sprintf("- 总基因数: %d\n", n_total),
  sprintf("- 成功分析数: %d\n", n_success),
  sprintf("- 显著基因数 (IVW p<0.05): %d\n\n", n_significant),
  "## 异质性显著基因 (Q_pval<0.05)\n",
  if (length(heterogeneity_genes) > 0) {
    paste0("- ", paste(heterogeneity_genes, collapse = "\n- "), "\n\n")
  } else {
    "- 无\n\n"
  },
  "## 多效性显著基因 (Egger_pval<0.05)\n",
  if (length(pleiotropy_genes) > 0) {
    paste0("- ", paste(pleiotropy_genes, collapse = "\n- "), "\n\n")
  } else {
    "- 无\n\n"
  },
  "## 共定位强证据基因\n",
  "- PP.H4 > 0.8 (强证据):\n",
  if (length(coloc_strong_genes) > 0) {
    paste0("  - ", paste(coloc_strong_genes[coloc_results[gene %in% coloc_strong_genes, PP.H4 > 0.8], ], collapse = "\n  - "), "\n")
  } else {
    "  - 无\n"
  },
  "- PP.H4 > 0.75 (证据):\n",
  if (length(coloc_strong_genes) > 0) {
    paste0("  - ", paste(coloc_strong_genes, collapse = "\n  - "), "\n\n")
  } else {
    "  - 无\n\n"
  },
  "## 建议剔除的基因 (Steiger检验失败或反向因果)\n",
  if (length(rejected_genes) > 0) {
    paste0("- ", paste(rejected_genes, collapse = "\n- "), "\n\n")
  } else {
    "- 无\n\n"
  },
  "## 分析时间\n",
  sprintf("- 分析完成时间: %s\n", Sys.time())
)

# 保存报告
writeLines(report_content, file.path(output_dir, "06_statistical_report.md"))
flog.info("统计报告已生成")

flog.info("=== 分析完成 ===")
flog.info(sprintf("结果已保存到: %s", output_dir))

# 显示摘要
cat("\n=== 分析摘要 ===\n")
cat(sprintf("总基因数: %d\n", n_total))
cat(sprintf("成功分析数: %d\n", n_success))
cat(sprintf("显著基因数 (IVW p<0.05): %d\n", n_significant))
cat("\n结果已保存到 D:/EQTL/MR_Enhanced_Results/ 目录\n")
