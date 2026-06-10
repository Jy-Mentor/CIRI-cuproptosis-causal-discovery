#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 铜死亡核心基因与缺血性脑卒中的孟德尔随机化分析
# 解决Windows下fread读取VCF失败、控制台滚动数字、列名识别错误等问题

# 包安装和加载
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("安装包: %s", pkg))
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
  if (!library(pkg, character.only = TRUE, quietly = TRUE, logical.return = TRUE)) {
    stop(sprintf("无法加载包: %s", pkg))
  }
}

# 安装必要的包
install_if_missing("data.table")
install_if_missing("dplyr")
install_if_missing("TwoSampleMR")
install_if_missing("MRPRESSO")
install_if_missing("ggplot2")

# 命令行参数解析
args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args
resume_mode <- "--resume" %in% args

# 目录设置
base_dir <- "D:/EQTL"
output_dir <- "D:/EQTL/MR_Results"
gene_results_dir <- file.path(output_dir, "gene_results")
ld_dir <- "D:/EQTL/clump"

# 创建输出目录
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(gene_results_dir)) dir.create(gene_results_dir, recursive = TRUE)

# 核心基因列表
core_genes <- c("FDX1", "DLAT", "DLD", "PDHB", "LIPT1", "PDHX", "SLC31A1", 
                "ATP7B", "ATP7A", "ATOX1", "COMMD1", "MT2A", "NFKB1", "RELA", 
                "STAT1", "STAT3", "STAT5A", "ICAM1", "CCL2", "IL6", "TGFB1", 
                "PTGS2", "HMOX1", "SOD2", "FABP3", "ATF4", "BRD4")

# 测试模式下仅使用3个核心基因
if (test_mode) {
  core_genes <- c("FDX1", "SLC31A1", "NFKB1")
  message("运行测试模式，仅分析3个核心基因")
}

# 检查点文件
checkpoint_stage1 <- file.path(output_dir, "checkpoint_stage1.RData")
checkpoint_stage2 <- file.path(output_dir, "checkpoint_stage2.RData")

# 检查点控制
skip_stage1 <- resume_mode && file.exists(checkpoint_stage1)
skip_stage2 <- resume_mode && file.exists(checkpoint_stage2)

# 阶段一：eQTL数据处理（修复版）
if (!skip_stage1) {
  message("\n=== 阶段一：eQTL数据处理 ===")
  
  # 读取eQTL分块文件
  eqtl_dir <- file.path(base_dir, "共定位顺式")
  eqtl_files <- list.files(eqtl_dir, pattern = "eqtlgen_ieu_1mbcis.*\\.csv$", full.names = TRUE)
  # 过滤掉正在下载的文件
  eqtl_files <- eqtl_files[!grepl('\\.downloading$', eqtl_files)]
  message(sprintf("找到 %d 个eQTL分块文件", length(eqtl_files)))
  
  if (length(eqtl_files) == 0) {
    # 尝试其他可能的文件模式
    eqtl_files <- list.files(eqtl_dir, pattern = "\\.csv$", full.names = TRUE)
    eqtl_files <- eqtl_files[!grepl('\\.downloading$', eqtl_files)]
    message(sprintf("尝试其他文件模式，找到 %d 个CSV文件", length(eqtl_files)))
    
    if (length(eqtl_files) == 0) {
      stop("没有找到有效的eQTL分块文件")
    }
  }
  
  # 合并eQTL数据
  eqtl_data <- data.table()
  for (file in eqtl_files) {
    message(sprintf("读取文件: %s", basename(file)))
    temp_data <- fread(file, showProgress = FALSE, data.table = TRUE)
    eqtl_data <- rbind(eqtl_data, temp_data)
  }
  
  message(sprintf("共读取 %d 条eQTL记录", nrow(eqtl_data)))
  
  # 动态识别列名
  col_names <- colnames(eqtl_data)
  message(sprintf("eQTL数据列名: %s", paste(col_names, collapse = ", ")))
  
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
  
  message(sprintf("使用列名 - P值: %s, Beta: %s, SE: %s, 基因: %s, SNP: %s", 
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
  message(sprintf("核心基因相关记录: %d 条", nrow(eqtl_data_core)))
  
  # 保存阶段一完整数据（关键修改：保存core而非filtered）
  save(eqtl_data_core, file = checkpoint_stage1)
  message(sprintf("阶段一完成，核心数据已保存到: %s", checkpoint_stage1))
  
  # 主分析筛选：p<5e-05（宽松参数）
  eqtl_filtered <- eqtl_data_core[pval.exposure < 5e-05, ]
  message(sprintf("主分析p<5e-05筛选后剩余 %d 条记录", nrow(eqtl_filtered)))
  
  # 计算F统计量并剔除弱工具变量
  eqtl_filtered[, F_stat := (beta.exposure / se.exposure)^2]
  eqtl_filtered <- eqtl_filtered[F_stat >= 10, ]
  message(sprintf("主分析剔除F<10后剩余 %d 条记录", nrow(eqtl_filtered)))
  
  # 使用预clumping结果（r²=0.01）
  clump_files <- list.files("D:/EQTL/clump", 
                             pattern = "eQTLgen_allgene_p_5e-[08]+_kb_10000_r2_0\\.001\\.xlsx$", 
                             full.names = TRUE)
  
  if(length(clump_files) == 0) {
    # 尝试其他模式
    clump_files <- list.files("D:/EQTL/clump", 
                               pattern = "eQTLgen.*5e.*\\.xlsx$", 
                               full.names = TRUE)
  }
  
  if(length(clump_files) > 0) {
    clump_file <- clump_files[1]  # 使用第一个匹配的文件
    message(sprintf("使用预clumping文件: %s", basename(clump_file)))
    
    if(!require("readxl")) install.packages("readxl")
    clumped_snps <- readxl::read_excel(clump_file)$SNP
    message(sprintf("预clumping文件包含 %d 个SNP", length(clumped_snps)))
    
    # 应用预clumping结果
    eqtl_filtered <- eqtl_filtered[SNP %in% clumped_snps, ]
    message(sprintf("预clumping后剩余 %d 条记录", nrow(eqtl_filtered)))
  } else {
    message("预clumping文件不存在，跳过clumping筛选")
  }
  
} else {
  message("\n=== 阶段一：跳过，从检查点恢复 ===")
  load(checkpoint_stage1)  # 加载eqtl_data_core
  message(sprintf("从检查点加载核心数据: %d 条记录", nrow(eqtl_data_core)))
  
  # 重新派生主分析数据
  eqtl_filtered <- eqtl_data_core[pval.exposure < 5e-05, ]
  eqtl_filtered[, F_stat := (beta.exposure / se.exposure)^2]
  eqtl_filtered <- eqtl_filtered[F_stat >= 10, ]
  message(sprintf("主分析数据派生完成: %d 条记录", nrow(eqtl_filtered)))
  
  # 使用预clumping结果（r²=0.01）
  clump_files <- list.files("D:/EQTL/clump", 
                             pattern = "eQTLgen_allgene_p_5e-[08]+_kb_10000_r2_0\\.001\\.xlsx$", 
                             full.names = TRUE)
  
  if(length(clump_files) > 0) {
    clump_file <- clump_files[1]  # 使用第一个匹配的文件
    message(sprintf("使用预clumping文件: %s", basename(clump_file)))
    
    if(!require("readxl")) install.packages("readxl")
    clumped_snps <- readxl::read_excel(clump_file)$SNP
    
    # 应用预clumping结果
    eqtl_filtered <- eqtl_filtered[SNP %in% clumped_snps, ]
    message(sprintf("预clumping后剩余 %d 条记录", nrow(eqtl_filtered)))
  } else {
    message("预clumping文件不存在，跳过clumping筛选")
  }
}

# 阶段二：VCF解析
if (!skip_stage2) {
  message("\n=== 阶段二：VCF解析 ===")
  
  vcf_file <- "D:/EQTL/ieu-a-83.vcf"
  message(sprintf("读取VCF文件: %s", vcf_file))
  
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
  
  # 创建临时文件
  temp_vcf <- tempfile(fileext = ".vcf")
  writeLines(c(header_line, data_lines), temp_vcf)
  
  # 使用read.table读取VCF数据
  outcome_data <- read.table(
    temp_vcf,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",      # 关键：禁用注释符号，防止#CHROM被截断
    check.names = FALSE,    # 关键：保留ieu-a-83连字符
    stringsAsFactors = FALSE,
    na.strings = c("NA", ".", ""),
    encoding = "UTF-8"
  )
  
  # 清理临时文件
  unlink(temp_vcf)
  
  message(sprintf("VCF文件读取完成，共 %d 行数据", nrow(outcome_data)))
  
  # 识别ieu-a-83列
  outcome_col <- if ("ieu-a-83" %in% colnames(outcome_data)) "ieu-a-83" else stop("无法找到ieu-a-83列")
  
  # 向量化解析ieu-a-83列
  parse_vcf_field <- function(field) {
    parts <- strsplit(field, ":")[[1]]
    if (length(parts) >= 4) {
      es <- as.numeric(parts[1])
      se <- as.numeric(parts[2])
      lp <- as.numeric(parts[3])
      ss <- as.numeric(parts[4])
      pval <- exp(-lp * log(10))
      return(c(BETA = es, SE = se, PVAL = pval, N = ss))
    } else {
      return(c(BETA = NA, SE = NA, PVAL = NA, N = NA))
    }
  }
  
  # 应用解析函数
  parsed_data <- t(vapply(outcome_data[[outcome_col]], parse_vcf_field, numeric(4)))
  outcome_data <- cbind(outcome_data, parsed_data)
  
  # 剔除NA值
  outcome_data <- outcome_data[!is.na(outcome_data$BETA) & !is.na(outcome_data$SE) & !is.na(outcome_data$PVAL), ]
  message(sprintf("剔除NA后剩余 %d 行数据", nrow(outcome_data)))
  
  # 查看VCF文件的列名
  message(sprintf("VCF文件列名: %s", paste(colnames(outcome_data), collapse = ", ")))
  
  # 统一列名
  if ("ID" %in% colnames(outcome_data)) {
    setnames(outcome_data, old = "ID", new = "SNP")
  } else if ("rsid" %in% colnames(outcome_data)) {
    setnames(outcome_data, old = "rsid", new = "SNP")
  } else if ("snp" %in% colnames(outcome_data)) {
    setnames(outcome_data, old = "snp", new = "SNP")
  }
  
  # 确保BETA、SE、PVAL列存在
  if ("BETA" %in% colnames(outcome_data)) {
    setnames(outcome_data, old = "BETA", new = "beta.outcome")
  }
  if ("SE" %in% colnames(outcome_data)) {
    setnames(outcome_data, old = "SE", new = "se.outcome")
  }
  if ("PVAL" %in% colnames(outcome_data)) {
    setnames(outcome_data, old = "PVAL", new = "pval.outcome")
  }
  if ("N" %in% colnames(outcome_data)) {
    setnames(outcome_data, old = "N", new = "samples.outcome")
  }
  
  # 再次查看处理后的列名
  message(sprintf("处理后VCF文件列名: %s", paste(colnames(outcome_data), collapse = ", ")))
  
  # 保存阶段二结果
  save(outcome_data, file = checkpoint_stage2)
  message(sprintf("阶段二完成，结果保存到: %s", checkpoint_stage2))
} else {
  message("\n=== 阶段二：跳过，从检查点恢复 ===")
  load(checkpoint_stage2)
  message(sprintf("从检查点加载 %d 行VCF数据", nrow(outcome_data)))
}

# 匹配SNP
message("\n=== 匹配暴露与结局SNP ===")

# 检查暴露数据中的SNP
message("检查暴露数据中的SNP...")
snp_exposure <- unique(eqtl_filtered$SNP)
snp_exposure <- snp_exposure[!is.na(snp_exposure) & snp_exposure != ""]
message(sprintf("暴露数据中有 %d 个唯一SNP", length(snp_exposure)))
if (length(snp_exposure) > 0) {
  message(sprintf("前5个暴露SNP: %s", paste(head(snp_exposure), collapse = ", ")))
}

# 检查结局数据中的SNP
message("检查结局数据中的SNP...")
snp_outcome <- unique(outcome_data$SNP)
snp_outcome <- snp_outcome[!is.na(snp_outcome) & snp_outcome != ""]
message(sprintf("结局数据中有 %d 个唯一SNP", length(snp_outcome)))
if (length(snp_outcome) > 0) {
  message(sprintf("前5个结局SNP: %s", paste(head(snp_outcome), collapse = ", ")))
} else {
  # 查看结局数据的前几行，了解SNP列的内容
  message("查看结局数据的前几行:")
  print(head(outcome_data[, c("SNP", "beta.outcome", "se.outcome", "pval.outcome")]))
}

# 匹配SNP
common_snps <- intersect(snp_exposure, snp_outcome)
message(sprintf("共同SNP数量: %d", length(common_snps)))
if (length(common_snps) > 0) {
  message(sprintf("前5个共同SNP: %s", paste(head(common_snps), collapse = ", ")))
}

if (length(common_snps) == 0) {
  message("警告: 没有找到共同SNP，分析无法继续")
  message("可能的原因:")
  message("1. 暴露数据和结局数据中的SNP ID格式不同")
  message("2. 没有重叠的SNP")
  message("3. 结局数据中的SNP列包含空值或无效值")
  # 不停止脚本，而是继续执行，这样至少可以看到其他错误信息
  # stop("没有找到共同SNP，分析无法继续")
}

# 阶段三：MR分析
message("\n=== 阶段三：MR分析 ===")

# 准备暴露数据
exposure_dat <- eqtl_filtered[eqtl_filtered$SNP %in% common_snps, ]
# 确保包含所有必要的列
exposure_dat <- exposure_dat[, c("SNP", "gene", "beta.exposure", "se.exposure", "pval.exposure", 
                               "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")]
# 标准化列名
exposure_dat$exposure <- exposure_dat$gene
exposure_dat$id.exposure <- exposure_dat$gene
# 移除不需要的列
exposure_dat <- exposure_dat[, c("SNP", "exposure", "id.exposure", "beta.exposure", "se.exposure", "pval.exposure", 
                               "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")]

# 准备结局数据
outcome_dat <- outcome_data[outcome_data$SNP %in% common_snps, ]
# 确保包含所有必要的列
outcome_dat <- outcome_dat[, c("SNP", "beta.outcome", "se.outcome", "pval.outcome", "REF", "ALT")]
# 添加必要的列
outcome_dat$id.outcome <- "ieu-a-83"
outcome_dat$outcome <- "Ischemic Stroke"
outcome_dat$effect_allele.outcome <- outcome_dat$ALT
outcome_dat$other_allele.outcome <- outcome_dat$REF
# 添加eaf.outcome列（假设等位基因频率为0.5作为默认值）
outcome_dat$eaf.outcome <- 0.5
# 移除不需要的列
outcome_dat <- outcome_dat[, c("SNP", "outcome", "id.outcome", "beta.outcome", "se.outcome", "pval.outcome", 
                             "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")]

# 按基因分组分析
mr_results <- list()
sensitivity_results <- list()

for (gene in core_genes) {
  message(sprintf("\n分析基因: %s", gene))
  
  # 检查是否已分析过
  result_file <- file.path(gene_results_dir, sprintf("%s.RData", gene))
  if(file.exists(result_file)){
    load(result_file)
    if(exists("res") && !is.null(res) && nrow(res)>0){
      message(gene, " 已完成，跳过")
      mr_results[[gene]] <- res
      sensitivity_results[[gene]] <- list(heterogeneity = hetero, pleiotropy = pleio, presso = presso_res, leaveoneout = loo)
      next
    } else {
      message(gene, " 之前失败，重新分析")
      file.remove(result_file)
    }
  }
  
  # 提取该基因的暴露数据
  gene_exposure <- exposure_dat[exposure_dat$exposure == gene, ]
  if (nrow(gene_exposure) == 0) {
    message(sprintf("%s 没有可用的eQTL数据，跳过", gene))
    next
  }
  
  # Harmonise数据
  dat <- NULL
  tryCatch({
    dat <- harmonise_data(
      exposure_dat = gene_exposure,
      outcome_dat = outcome_dat
    )
  }, error = function(e) {
    message(sprintf("%s Harmonise失败: %s", gene, e$message))
    next
  })
  
  if (is.null(dat) || nrow(dat) == 0) {
    message(sprintf("%s 没有可 harmonise 的SNP，跳过", gene))
    next
  }
  
  # MR分析：仅使用逆方差加权（IVW）法
  res <- NULL
  tryCatch({
    res <- mr(dat, method_list = c("mr_ivw"))
  }, error = function(e) {
    message(sprintf("%s MR分析失败: %s", gene, e$message))
    next
  })
  
  if (is.null(res) || nrow(res) == 0) {
    message(sprintf("%s MR分析无结果，跳过", gene))
    next
  }
  
  # 敏感性分析
  hetero <- NULL
  tryCatch({
    hetero <- mr_heterogeneity(dat)
  }, error = function(e) {
    message(sprintf("%s 异质性分析失败: %s", gene, e$message))
  })
  
  pleio <- NULL
  tryCatch({
    pleio <- mr_pleiotropy_test(dat)
  }, error = function(e) {
    message(sprintf("%s 多效性分析失败: %s", gene, e$message))
  })
  
  # MR-PRESSO
  presso_res <- NULL
  tryCatch({
    # 使用MendelianRandomization::mr_presso
    if (!require("MendelianRandomization")) install.packages("MendelianRandomization")
    presso_res <- MendelianRandomization::mr_presso(
      BetaOutcome="beta.outcome", BetaExposure="beta.exposure",
      SdOutcome="se.outcome", SdExposure="se.exposure",
      OUTLIERtest=TRUE, DISTORTIONtest=TRUE, data=dat, nboot=1000
    )
  }, error=function(e) message("PRESSO跳过: ", conditionMessage(e)))
  
  # 留一法分析
  loo <- NULL
  tryCatch({
    loo <- mr_leaveoneout(dat)
  }, error = function(e) {
    message(sprintf("%s 留一法分析失败: %s", gene, e$message))
    next
  })
  
  # 保存结果
  save(res, hetero, pleio, presso_res, loo, file = result_file)
  mr_results[[gene]] <- res
  sensitivity_results[[gene]] <- list(heterogeneity = hetero, pleiotropy = pleio, presso = presso_res, leaveoneout = loo)
  
  # 生成图形
  if (gene %in% c("FDX1", "SLC31A1", "NFKB1")) {
    # 散点图
    tryCatch({
      scatter_plot <- mr_scatter_plot(res, dat)
      if(!is.null(scatter_plot) && length(scatter_plot)>0 && !is.null(scatter_plot[[1]])){
        ggsave(file.path(output_dir, sprintf("%s_scatter.png", gene)), scatter_plot[[1]], width = 8, height = 6)
      }
    }, error=function(e) message(gene, " 散点图失败: ", conditionMessage(e)))
  }
  
  # 森林图
  tryCatch({
    forest_plot <- mr_forest_plot(res)
    if(!is.null(forest_plot) && length(forest_plot)>0 && !is.null(forest_plot[[1]])){
      ggsave(file.path(output_dir, sprintf("%s_forest.png", gene)), forest_plot[[1]], width = 8, height = 6)
    }
  }, error=function(e) message(gene, " 森林图失败: ", conditionMessage(e)))
  
  # 漏斗图
  tryCatch({
    funnel_plot <- mr_funnel_plot(loo)
    if(!is.null(funnel_plot) && length(funnel_plot)>0 && !is.null(funnel_plot[[1]])){
      ggsave(file.path(output_dir, sprintf("%s_funnel.png", gene)), funnel_plot[[1]], width = 8, height = 6)
    }
  }, error=function(e) message(gene, " 漏斗图失败: ", conditionMessage(e)))
  
  # 留一法图
  tryCatch({
    loo_plot <- mr_leaveoneout_plot(loo)
    if(!is.null(loo_plot) && length(loo_plot)>0 && !is.null(loo_plot[[1]])){
      ggsave(file.path(output_dir, sprintf("%s_leaveoneout.png", gene)), loo_plot[[1]], width = 8, height = 6)
    }
  }, error=function(e) message(gene, " 留一法图失败: ", conditionMessage(e)))
  
  message(sprintf("%s 分析完成", gene))
}

# 阶段四：结果输出
message("\n=== 阶段四：结果输出 ===")

# 合并主结果
results_list <- list()
for(gene in core_genes){
  f <- file.path(gene_results_dir, paste0(gene,".RData"))
  if(file.exists(f)){
    load(f)
    if(exists("res") && !is.null(res) && is.data.frame(res) && nrow(res)>0){
      res$gene <- gene
      results_list[[gene]] <- res
    }
  }
}
if(length(results_list)>0){
  # 使用rbindlist处理列数不匹配的情况
  final <- data.table::rbindlist(results_list, fill = TRUE)
  write.csv(final, file.path(output_dir,"MR_main_results.csv"), row.names = FALSE)
  message("主结果已保存")
} else {
  message("警告：无有效MR结果")
}

# 敏感性分析：使用严格参数（p<5×10⁻⁸，r²=0.001，窗口10000kb）
message("\n=== 敏感性分析：使用严格参数 ===")
sensitivity_output_dir <- file.path(output_dir, "sensitivity_strict")
if (!dir.exists(sensitivity_output_dir)) dir.create(sensitivity_output_dir, recursive = TRUE)

# 直接从内存或checkpoint获取核心数据（不再重新fread所有CSV）
if (!exists("eqtl_data_core")) {
  if (file.exists(checkpoint_stage1)) {
    load(checkpoint_stage1)
    message(sprintf("从阶段一加载核心数据: %d 条记录", nrow(eqtl_data_core)))
  } else {
    stop("阶段一数据不存在，请先运行阶段一")
  }
}

# 使用严格参数筛选：p<5e-08
eqtl_filtered_sensitivity <- eqtl_data_core[pval.exposure < 5e-08, ]
message(sprintf("严格参数p<5e-08筛选后剩余 %d 条记录", nrow(eqtl_filtered_sensitivity)))

# 计算F统计量
eqtl_filtered_sensitivity[, F_stat := (beta.exposure / se.exposure)^2]

# 剔除F<10的记录
eqtl_filtered_sensitivity <- eqtl_filtered_sensitivity[F_stat >= 10, ]
message(sprintf("严格参数剔除F<10后剩余 %d 条记录", nrow(eqtl_filtered_sensitivity)))

# 使用严格预clumping结果（r²=0.001）
clump_file_strict <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-8_kb_10000_r2_0.001.xlsx"
if(file.exists(clump_file_strict)){
  if(!require("readxl")) install.packages("readxl")
  clumped_snps_strict <- readxl::read_excel(clump_file_strict)$SNP
  eqtl_filtered_sensitivity <- eqtl_filtered_sensitivity[SNP %in% clumped_snps_strict, ]
  message(sprintf("严格参数预clumping后剩余 %d 条记录", nrow(eqtl_filtered_sensitivity)))
} else {
  message("警告：严格参数预clumping文件不存在，跳过clumping筛选")
}

# 准备暴露数据（复用阶段二的outcome_dat和common_snps）
exposure_dat_sensitivity <- eqtl_filtered_sensitivity[SNP %in% common_snps, ]
if(nrow(exposure_dat_sensitivity) == 0){
  message("警告：严格参数下没有共同SNP，跳过敏感性分析")
} else {
  exposure_dat_sensitivity$exposure <- exposure_dat_sensitivity$gene
  exposure_dat_sensitivity$id.exposure <- exposure_dat_sensitivity$gene
  
  # 确保包含所有必要的列
  exposure_dat_sensitivity <- exposure_dat_sensitivity[, c("SNP", "exposure", "id.exposure", "beta.exposure", "se.exposure", "pval.exposure", 
                                                         "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")]
}

# 按基因分组分析敏感性结果
sensitivity_results_list <- list()
if(nrow(exposure_dat_sensitivity) > 0) {
  for (gene in core_genes) {
    message(sprintf("\n敏感性分析基因: %s", gene))
    
    # 提取该基因的暴露数据
    gene_exposure <- exposure_dat_sensitivity[exposure_dat_sensitivity$exposure == gene, ]
    if (nrow(gene_exposure) == 0) {
      message(sprintf("%s 没有可用的eQTL数据，跳过", gene))
      next
    }
    
    # Harmonise数据
    dat_sensitivity <- NULL
    tryCatch({
      dat_sensitivity <- harmonise_data(
        exposure_dat = gene_exposure,
        outcome_dat = outcome_dat
      )
    }, error = function(e) {
      message(sprintf("%s Harmonise失败: %s", gene, e$message))
      next
    })
    
    if (is.null(dat_sensitivity) || nrow(dat_sensitivity) == 0) {
      message(sprintf("%s 没有可 harmonise 的SNP，跳过", gene))
      next
    }
    
    # MR分析：仅使用逆方差加权（IVW）法
    res_sensitivity <- NULL
    tryCatch({
      res_sensitivity <- mr(dat_sensitivity, method_list = c("mr_ivw"))
    }, error = function(e) {
      message(sprintf("%s MR分析失败: %s", gene, e$message))
      next
    })
    
    if (is.null(res_sensitivity) || nrow(res_sensitivity) == 0) {
      message(sprintf("%s MR分析无结果，跳过", gene))
      next
    }
    
    res_sensitivity$gene <- gene
    sensitivity_results_list[[gene]] <- res_sensitivity
    message(sprintf("%s 敏感性分析完成", gene))
  }
  
  # 保存敏感性分析结果
  if(length(sensitivity_results_list)>0){
    final_sensitivity <- data.table::rbindlist(sensitivity_results_list, fill = TRUE)
    write.csv(final_sensitivity, file.path(sensitivity_output_dir,"MR_sensitivity_strict.csv"), row.names = FALSE)
    message("敏感性分析结果已保存")
  } else {
    message("没有敏感性分析结果")
  }
} else {
  message("跳过敏感性分析结果保存")
}

# 合并敏感性结果
sensitivity_list <- list()
for(gene in core_genes){
  f <- file.path(gene_results_dir, paste0(gene,".RData"))
  if(file.exists(f)){
    load(f)
    if(exists("hetero") && !is.null(hetero) && is.data.frame(hetero) && nrow(hetero)>0){
      hetero$gene <- gene
      sensitivity_list[[paste0(gene, "_hetero")]] <- hetero
    }
    if(exists("pleio") && !is.null(pleio) && is.data.frame(pleio) && nrow(pleio)>0){
      pleio$gene <- gene
      sensitivity_list[[paste0(gene, "_pleio")]] <- pleio
    }
  }
}
if(length(sensitivity_list)>0){
  # 使用rbindlist处理列数不匹配的情况
  final_sensitivity <- data.table::rbindlist(sensitivity_list, fill = TRUE)
  write.csv(final_sensitivity, file.path(output_dir,"MR_sensitivity.csv"), row.names = FALSE)
  message("敏感性分析结果已保存")
} else {
  message("没有敏感性分析结果")
}

# 按基因输出结果
if(length(results_list)>0){
  # 使用rbindlist处理列数不匹配的情况
  final_by_gene <- data.table::rbindlist(results_list, fill = TRUE)
  write.csv(final_by_gene, file.path(output_dir,"MR_by_gene.csv"), row.names = FALSE)
  message("按基因分组结果已保存")
} else {
  message("没有按基因分组结果")
}

# Bonferroni校正
bonferroni_threshold <- 0.05 / length(core_genes)
message(sprintf("Bonferroni校正阈值: %.6f", bonferroni_threshold))

# 输出重点基因的OR(95%CI)
message("\n=== 重点基因分析结果 ===")
key_genes <- c("FDX1", "SLC31A1", "NFKB1")
for (gene in key_genes) {
  if (gene %in% names(results_list)) {
    res <- results_list[[gene]]
    ivw_res <- res[res$method == "Inverse variance weighted", ]
    if (nrow(ivw_res) > 0) {
      or <- exp(ivw_res$b)
      ci_lower <- exp(ivw_res$b - 1.96 * ivw_res$se)
      ci_upper <- exp(ivw_res$b + 1.96 * ivw_res$se)
      pval <- ivw_res$pval
      message(sprintf("%s: OR = %.3f (95%%CI: %.3f-%.3f), P = %.6f", 
                      gene, or, ci_lower, ci_upper, pval))
    }
  }
}

message("\n=== 分析完成 ===")
message(sprintf("结果已保存到: %s", output_dir))