#!/usr/bin/env Rscript
# ================================================================================
# 整合 eQTLGen 全血与 GTEx 脑组织双源 eQTL 数据的 MR 分析
# 使用本地 GTEx v11 parquet 文件
# 目标：提高工具变量数量，发现组织特异性效应
# ================================================================================

# 包管理
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("正在安装:", pkg))
      tryCatch({
        install.packages(pkg, repos = "https://cloud.r-project.org/")
      }, error = function(e) {
        message(paste("安装失败:", pkg, "-", e$message))
      })
    }
  }
}

# 安装 TwoSampleMR
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("MRCIEU/TwoSampleMR", upgrade = "never")
}

install_if_missing(c("dplyr", "readr", "data.table", "stringr", "tidyr", "arrow"))

library(TwoSampleMR)
library(dplyr)
library(readr)
library(data.table)
library(stringr)
library(tidyr)
library(arrow)

# ================================================================================
# 配置参数
# ================================================================================

# 本地 GTEx v11 数据路径
GTEx_BRAIN_FILE <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Brain_Cortex.v11.eQTLs.signif_pairs.parquet"
GTEx_BLOOD_FILE <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eQTLs.signif_pairs.parquet"

# 结果数据 - MEGASTROKE
OUTCOME_FILE <- "D:/EQTL/finngen_R12_C_STROKE"

# 基因列表
GENE_LIST_FILE <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/gene_list_optimized.txt"

# 输出目录
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/dual_source_mr_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# MR 参数
PVAL_THRESHOLD <- 1e-5  # P 值阈值
CLUMP_R2 <- 0.01        # LD clump r²阈值
CLUMP_KB <- 10000       # LD clump 窗口 (kb)
F_THRESHOLD <- 10       # F 统计量阈值

cat("======================================================================\n")
cat("双源 eQTL 整合 MR 分析\n")
cat("======================================================================\n")
cat(sprintf("\n输出目录：%s\n", OUTPUT_DIR))

# ================================================================================
# 加载基因列表
# ================================================================================

gene_list <- character(0)
if (file.exists(GENE_LIST_FILE)) {
  gene_list <- readLines(GENE_LIST_FILE, warn = FALSE)
  gene_list <- gene_list[gene_list != ""]
  cat(sprintf("✓ 加载基因列表：%d 个基因\n", length(gene_list)))
} else {
  cat(sprintf("✗ 基因列表文件不存在：%s\n", GENE_LIST_FILE))
}

# ================================================================================
# 加载 GTEx v11 eQTL 数据
# ================================================================================

cat("\n======================================================================\n")
cat("加载 GTEx v11 eQTL 数据\n")
cat("======================================================================\n")

# 加载脑组织数据
brain_eqtl <- NULL
if (file.exists(GTEx_BRAIN_FILE)) {
  cat("\n加载脑皮层 eQTL 数据...\n")
  tryCatch({
    brain_eqtl <- read_parquet(GTEx_BRAIN_FILE)
    cat(sprintf("  ✓ 脑皮层数据：%d 个 eQTL 对\n", nrow(brain_eqtl)))
    brain_genes <- unique(brain_eqtl$gene_id)
    cat(sprintf("    包含：%d 个基因\n", length(brain_genes)))
  }, error = function(e) {
    cat(sprintf("  ✗ 加载失败：%s\n", e$message))
    brain_eqtl <- NULL
  })
} else {
  cat(sprintf("  ✗ 文件不存在：%s\n", GTEx_BRAIN_FILE))
}

# 加载全血数据
blood_eqtl <- NULL
if (file.exists(GTEx_BLOOD_FILE)) {
  cat("\n加载全血 eQTL 数据...\n")
  tryCatch({
    blood_eqtl <- read_parquet(GTEx_BLOOD_FILE)
    cat(sprintf("  ✓ 全血数据：%d 个 eQTL 对\n", nrow(blood_eqtl)))
    blood_genes <- unique(blood_eqtl$gene_id)
    cat(sprintf("    包含：%d 个基因\n", length(blood_genes)))
  }, error = function(e) {
    cat(sprintf("  ✗ 加载失败：%s\n", e$message))
    blood_eqtl <- NULL
  })
} else {
  cat(sprintf("  ✗ 文件不存在：%s\n", GTEx_BLOOD_FILE))
}

# ================================================================================
# 格式化 eQTL 数据为 MR 输入
# ================================================================================

format_eqtl_for_mr <- function(eqtl_data, gene_name, tissue_type) {
  """格式化 eQTL 数据为 MR 分析输入格式"""
  
  if (is.null(eqtl_data) || nrow(eqtl_data) == 0) {
    return(NULL)
  }
  
  # 筛选该基因的 eQTL
  gene_data <- eqtl_data %>% filter(gene_id == gene_name)
  
  if (nrow(gene_data) == 0) {
    return(NULL)
  }
  
  # 选择最强的 eQTL (最低 P 值) - 限制 Top 50
  gene_data <- gene_data %>% arrange(pval_nominal) %>% head(50)
  
  # 解析 variant_id 获取染色体和位置
  # variant_id 格式：chr_pos_ref_alt_b38
  parse_variant_id <- function(variant_id) {
    tryCatch({
      parts <- strsplit(as.character(variant_id), "_")[[1]]
      if (length(parts) >= 4) {
        chr_ <- gsub("chr", "", parts[1])
        pos <- as.integer(parts[2])
        ref <- parts[3]
        alt <- parts[4]
        return(list(chr=chr_, pos=pos, ref=ref, alt=alt))
      }
    }, error = function(e) {})
    return(list(chr="NA", pos=0, ref="NA", alt="NA"))
  }
  
  # 解析变异信息
  variant_info <- lapply(gene_data$variant_id, parse_variant_id)
  gene_data$CHR <- sapply(variant_info, function(x) x$chr)
  gene_data$BP <- sapply(variant_info, function(x) x$pos)
  gene_data$REF <- sapply(variant_info, function(x) x$ref)
  gene_data$ALT <- sapply(variant_info, function(x) x$alt)
  
  # 格式化输出
  mr_format <- gene_data %>%
    transmute(
      SNP = variant_id,
      CHR = CHR,
      BP = BP,
      EFFECT_ALLELE = ALT,
      OTHER_ALLELE = REF,
      BETA = slope,
      SE = slope_se,
      PVAL = pval_nominal,
      EAF = af,
      GENE = gene_name,
      TISSUE = tissue_type,
      TSS_DISTANCE = tss_distance
    )
  
  return(mr_format)
}

# ================================================================================
# 为每个基因创建暴露数据
# ================================================================================

cat("\n======================================================================\n")
cat("准备每个基因的暴露数据\n")
cat("======================================================================\n")

# 获取所有可用基因
brain_genes <- character(0)
if (!is.null(brain_eqtl)) {
  brain_genes <- unique(brain_eqtl$gene_id)
}

blood_genes <- character(0)
if (!is.null(blood_eqtl)) {
  blood_genes <- unique(blood_eqtl$gene_id)
}

# 使用基因列表
if (length(gene_list) > 0) {
  target_genes <- gene_list
  cat(sprintf("\n目标基因：%d 个（来自基因列表）\n", length(target_genes)))
} else {
  # 使用所有在任一组织中表达的基因
  target_genes <- unique(c(brain_genes, blood_genes))
  cat(sprintf("\n目标基因：%d 个（来自 eQTL 数据）\n", length(target_genes)))
}

# 统计
stats <- list(
  brain_only = 0,
  blood_only = 0,
  both = 0,
  neither = 0
)

# 创建暴露列表
exposure_list <- list()

cat(sprintf("\n开始处理 %d 个基因...\n", length(target_genes)))

for (i in seq_along(target_genes)) {
  gene <- target_genes[i]
  
  # 检查基因是否在数据中
  in_brain <- gene %in% brain_genes
  in_blood <- gene %in% blood_genes
  
  all_data <- list()
  
  # 添加脑组织数据
  if (in_brain && !is.null(brain_eqtl)) {
    brain_data <- format_eqtl_for_mr(brain_eqtl, gene, "Brain_Cortex")
    if (!is.null(brain_data) && nrow(brain_data) > 0) {
      all_data[[1]] <- brain_data
    }
  }
  
  # 添加全血数据
  if (in_blood && !is.null(blood_eqtl)) {
    blood_data <- format_eqtl_for_mr(blood_eqtl, gene, "Whole_Blood")
    if (!is.null(blood_data) && nrow(blood_data) > 0) {
      all_data[[length(all_data) + 1]] <- blood_data
    }
  }
  
  # 保存
  if (length(all_data) > 0) {
    # 合并所有组织的数据
    combined <- bind_rows(all_data)
    exposure_list[[gene]] <- combined
    
    # 更新统计
    if (in_brain && in_blood) {
      stats$both