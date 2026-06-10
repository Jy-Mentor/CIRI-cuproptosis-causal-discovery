#!/usr/bin/env Rscript
# MR分析脚本 - GTEx v11 Whole Blood数据集
# 使用分块读取处理50GB结局文件

library(readxl)
library(data.table)
library(readr)
library(TwoSampleMR)
library(dplyr)
library(purrr)

# ============================================
# 参数设置
# ============================================
exposure_file <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eGenes.txt"
signif_pairs_file <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eQTLs.signif_pairs.parquet"
outcome_file  <- "D:/EQTL/eqtlgen_ieu_outcome.csv"
genes         <- c("NFKB1", "FDX1", "STAT3", "HIF1A", "HMOX1",
                   "GPX4", "HSPA5", "AGER", "DLAT")
output_dir    <- "D:/EQTL/mr_results_gtex_v11"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== GTEx v11 Whole Blood MR 分析 ==========\n")
cat("暴露文件:", exposure_file, "\n")
cat("结局文件:", outcome_file, "\n")
cat("分析基因:", paste(genes, collapse = ", "), "\n\n")

# ============================================
# 1. 读取GTEx eGenes数据并筛选目标基因
# ============================================
cat("步骤1: 读取GTEx eGenes数据...\n")

# 读取eGenes文件
egenes_data <- fread(exposure_file, header = TRUE)
cat("eGenes数据行数:", nrow(egenes_data), "\n")
cat("eGenes列名:", paste(names(egenes_data), collapse = ", "), "\n\n")

# 筛选目标基因（gene_name列）
target_egenes <- egenes_data[gene_name %in% genes, ]
cat("目标基因eGenes数:", nrow(target_egenes), "\n")

if (nrow(target_egenes) == 0) {
  cat("错误: 未找到目标基因，请检查基因名是否正确\n")
  quit(status = 1)
}

# 获取基因ID列表
target_gene_ids <- unique(target_egenes$gene_id)
cat("目标基因ID:", paste(target_gene_ids, collapse = ", "), "\n\n")

# ============================================
# 2. 读取signif_pairs.parquet文件获取SNP信息
# ============================================
cat("步骤2: 读取signif_pairs数据...\n")

# 检查arrow包是否安装
if (!requireNamespace("arrow", quietly = TRUE)) {
  cat("安装arrow包...\n")
  install.packages("arrow", repos = "https://cran.r-project.org")
}
library(arrow)

# 读取parquet文件
signif_pairs <- read_parquet(signif_pairs_file)
cat("signif_pairs总行数:", nrow(signif_pairs), "\n")
cat("signif_pairs列名:", paste(names(signif_pairs), collapse = ", "), "\n")

# 检查gene_id列是否存在
gene_id_col <- NULL
if ("gene_id" %in% names(signif_pairs)) {
  gene_id_col <- "gene_id"
} else if ("phenotype_id" %in% names(signif_pairs)) {
  gene_id_col <- "phenotype_id"
}

if (is.null(gene_id_col)) {
  cat("错误: 找不到gene_id列\n")
  quit(status = 1)
}

# 筛选目标基因的SNP
exposure_snps <- signif_pairs[signif_pairs[[gene_id_col]] %in% target_gene_ids, ]
cat("目标基因SNP数:", nrow(exposure_snps), "\n\n")

if (nrow(exposure_snps) == 0) {
  cat("错误: 未找到目标基因的SNP\n")
  quit(status = 1)
}

# 添加基因名
gene_id_to_name <- setNames(target_egenes$gene_name, target_egenes$gene_id)
exposure_snps$gene_name <- gene_id_to_name[exposure_snps[[gene_id_col]]]

# 格式化暴露数据
# GTEx v11 signif_pairs列名：phenotype_id, variant_id, start_distance, af, ma_samples, ma_count, 
# pval_nominal, slope, slope_se, pval_nominal_threshold, min_pval_nominal, pval_beta
exposure_dat <- data.frame(
  SNP = exposure_snps$variant_id,
  beta.exposure = exposure_snps$slope,
  se.exposure = exposure_snps$slope_se,
  pval.exposure = exposure_snps$pval_nominal,
  eaf.exposure = exposure_snps$af,
  gene.exposure = exposure_snps$gene_name,
  stringsAsFactors = FALSE
)

cat("格式化后暴露数据行数:", nrow(exposure_dat), "\n")
cat("涉及基因:", paste(unique(exposure_dat$gene.exposure), collapse = ", "), "\n\n")

# 提取所有暴露SNP（用于过滤50GB结局文件）
all_snps <- unique(exposure_dat$SNP)
cat("暴露SNP总数:", length(all_snps), "\n")

# 保存SNP列表（备用）
writeLines(all_snps, file.path(output_dir, "exposure_snplist.txt"))

# ============================================
# 3. 读取结局数据（优先使用Python预过滤的文件）
# ============================================
fast_filtered_file <- "D:/EQTL/mr_results_p5e-06/outcome_filtered_fast.csv"

if (file.exists(fast_filtered_file)) {
  cat("\n步骤3: 检测到Python预过滤的结局文件，直接读取...\n")
  cat("文件:", fast_filtered_file, "\n")
  outcome_small <- fread(fast_filtered_file)
  cat("过滤后结局数据行数:", nrow(outcome_small), "\n\n")
} else {
  cat("\n步骤3: 未找到预过滤文件，进行分块过滤50GB结局文件（这可能需要10-30分钟）...\n")
  cat("Chunk size: 50万行/块\n")
  
  # 先读100行看列名
  peek <- read_csv(outcome_file, n_max = 5, show_col_types = FALSE)
  cat("结局文件列名:", paste(names(peek), collapse = ", "), "\n")
  
  # 分块过滤函数
  filter_callback <- function(chunk, pos) {
    snp_col <- NULL
    if ("SNP" %in% names(chunk)) snp_col <- "SNP"
    else if ("rsids" %in% names(chunk)) snp_col <- "rsids"
    
    if (!is.null(snp_col)) {
      sub <- chunk[chunk[[snp_col]] %in% all_snps, ]
      if (nrow(sub) > 0) return(sub)
    }
    return(NULL)
  }
  
  # 执行分块读取
  start_time <- Sys.time()
  
  outcome_filtered <- read_csv_chunked(
    file        = outcome_file,
    callback    = DataFrameCallback$new(filter_callback),
    chunk_size  = 500000,
    col_types   = cols(.default = "c"),
    progress    = TRUE
  )
  
  end_time <- Sys.time()
  cat("\n过滤完成！耗时:", round(difftime(end_time, start_time, units = "min"), 1), "分钟\n")
  cat("过滤后结局数据行数:", nrow(outcome_filtered), "\n")
  
  if (nrow(outcome_filtered) == 0) {
    cat("错误: 没有匹配到任何SNP\n")
    quit(status = 1)
  }
  
  # 保存过滤后的结局
  fwrite(as.data.table(outcome_filtered),
         file.path(output_dir, "outcome_filtered.csv"))
  
  rm(outcome_filtered)
  gc()
  
  cat("步骤4: 读入过滤后的结局数据...\n")
  outcome_small <- fread(file.path(output_dir, "outcome_filtered.csv"))
}
cat("过滤后结局数据行数:", nrow(outcome_small), "\n\n")

# 确定结局数据列名
outcome_cols <- names(outcome_small)
snp_col <- ifelse("SNP" %in% outcome_cols, "SNP",
                  ifelse("rsids" %in% outcome_cols, "rsids", outcome_cols[1]))
beta_col <- ifelse("BETA" %in% outcome_cols, "BETA",
                   ifelse("beta" %in% outcome_cols, "beta", outcome_cols[2]))
se_col <- ifelse("SE" %in% outcome_cols, "SE",
                 ifelse("se" %in% outcome_cols, "se", outcome_cols[3]))
ea_col <- ifelse("A1" %in% outcome_cols, "A1",
                 ifelse("effect_allele" %in% outcome_cols, "effect_allele", outcome_cols[4]))
oa_col <- ifelse("A2" %in% outcome_cols, "A2",
                 ifelse("other_allele" %in% outcome_cols, "other_allele", outcome_cols[5]))
pval_col <- ifelse("P" %in% outcome_cols, "P",
                   ifelse("pval" %in% outcome_cols, "pval", outcome_cols[6]))

# 将data.table转换为data.frame
outcome_small_df <- as.data.frame(outcome_small)

# 格式化结局数据
outcome_dat <- format_data(
  dat          = outcome_small_df,
  type         = "outcome",
  snp_col      = snp_col,
  beta_col     = beta_col,
  se_col       = se_col,
  effect_allele_col = ea_col,
  other_allele_col  = oa_col,
  pval_col     = pval_col
)

# ============================================
# 5. 循环每个基因跑MR
# ============================================
cat("步骤5: 逐个基因跑MR...\n\n")

results_list <- list()

for (gene in genes) {
  cat("==========", gene, "==========\n")
  
  # 提取该基因暴露
  exp_gene <- exposure_dat[exposure_dat$gene.exposure == gene, ]
  cat("  暴露SNP:", nrow(exp_gene), "个\n")
  
  if (nrow(exp_gene) == 0) {
    cat("  跳过: 暴露数据中无该基因\n\n")
    next
  }
  
  if (nrow(exp_gene) < 2) {
    cat("  跳过: SNP不足2个，无法做MR\n\n")
    next
  }
  
  # 使用TwoSampleMR的format_data格式化
  exp_gene_fmt <- format_data(
    dat = exp_gene,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    pval_col = "pval.exposure"
  )
  
  # Harmonise
  dat <- tryCatch({
    harmonise_data(exposure_dat = exp_gene_fmt, outcome_dat = outcome_dat)
  }, error = function(e) {
    cat("  harmonise失败:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(dat) || nrow(dat) == 0) {
    cat("  跳过: harmonise后无匹配SNP\n\n")
    next
  }
  
  cat("  Harmonise后SNP:", nrow(dat), "个\n")
  
  # 计算F统计量并过滤弱工具变量
  dat$F_stat <- (dat$beta.exposure / dat$se.exposure)^2
  dat <- dat[dat$F_stat > 10, ]
  
  if (nrow(dat) == 0) {
    cat("  跳过: F统计量均<10（弱工具变量）\n\n")
    next
  }
  
  cat("  F>10的SNP:", nrow(dat), "个\n")
  
  # 跑MR
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                 "mr_weighted_median", "mr_weighted_mode"))
  
  # 敏感性分析
  hetero <- mr_heterogeneity(dat)
  pleio  <- mr_pleiotropy_test(dat)
  
  # 保存结果
  results_list[[gene]] <- list(
    gene    = gene,
    dat     = dat,
    res     = res,
    hetero  = hetero,
    pleio   = pleio
  )
  
  # 打印关键结果
  ivw <- res[res$method == "Inverse variance weighted", ]
  if (nrow(ivw) > 0 && !is.na(ivw$b[1])) {
    cat("  IVW Beta:", round(ivw$b[1], 4),
        "| OR:", round(exp(ivw$b[1]), 3),
        "| P:", format(ivw$pval[1], digits = 3, scientific = TRUE), "\n\n")
  } else {
    cat("  IVW结果无效或缺失\n\n")
  }
  
  # 保存该基因结果
  write.csv(res, file.path(output_dir, paste0(gene, "_MR_results.csv")), row.names = FALSE)
  write.csv(dat, file.path(output_dir, paste0(gene, "_harmonised_data.csv")), row.names = FALSE)
}

# ============================================
# 6. 汇总结果表
# ============================================
cat("步骤6: 生成汇总结果...\n")

summary_df <- data.frame(
  Gene = character(),
  Method = character(),
  Beta = numeric(),
  SE = numeric(),
  OR = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (gene in names(results_list)) {
  res <- results_list[[gene]]$res
  for (i in 1:nrow(res)) {
    # 检查值是否有效
    beta_val <- res$b[i]
    if (!is.na(beta_val) && is.numeric(beta_val)) {
      summary_df <- rbind(summary_df, data.frame(
        Gene = gene,
        Method = res$method[i],
        Beta = beta_val,
        SE = res$se[i],
        OR = exp(beta_val),
        P_value = res$pval[i]
      ))
    }
  }
}

write.csv(summary_df, file.path(output_dir, "MR_summary_all_genes.csv"), row.names = FALSE)

# 只保留IVW结果
ivw_summary <- summary_df[summary_df$Method == "Inverse variance weighted", ]
write.csv(ivw_summary, file.path(output_dir, "MR_summary_IVW_only.csv"), row.names = FALSE)

cat("\n========== 分析完成 ==========\n")
cat("输出目录:", output_dir, "\n")
cat("分析基因数:", length(results_list), "\n")
cat("汇总结果已保存\n")
