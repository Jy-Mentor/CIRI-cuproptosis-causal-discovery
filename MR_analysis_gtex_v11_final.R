#!/usr/bin/env Rscript
# MR分析脚本 - GTEx v11 Whole Blood数据集（最终版）
# 结合signif_pairs（所有SNP）和eGenes（ref/alt信息）

library(readxl)
library(data.table)
library(readr)
library(TwoSampleMR)
library(dplyr)
library(purrr)

# ============================================
# 参数设置
# ============================================
egenes_file <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eGenes.txt"
signif_pairs_file <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eQTLs.signif_pairs.parquet"
outcome_file  <- "D:/EQTL/eqtlgen_ieu_outcome.csv"
genes         <- c("NFKB1", "FDX1", "STAT3", "HIF1A", "HMOX1",
                   "GPX4", "HSPA5", "AGER", "DLAT")
output_dir    <- "D:/EQTL/mr_results_gtex_v11"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== GTEx v11 Whole Blood MR 分析（最终版） ==========\n")
cat("eGenes文件:", egenes_file, "\n")
cat("signif_pairs文件:", signif_pairs_file, "\n")
cat("结局文件:", outcome_file, "\n")
cat("分析基因:", paste(genes, collapse = ", "), "\n\n")

# ============================================
# 1. 读取GTEx eGenes数据并创建ref/alt映射
# ============================================
cat("步骤1: 读取GTEx eGenes数据...\n")

egenes_data <- fread(egenes_file, header = TRUE)
cat("eGenes数据行数:", nrow(egenes_data), "\n")

# 筛选目标基因
target_egenes <- egenes_data[gene_name %in% genes, ]
cat("目标基因eGenes数:", nrow(target_egenes), "\n")

if (nrow(target_egenes) == 0) {
  cat("错误: 未找到目标基因\n")
  quit(status = 1)
}

# 创建variant_id到ref/alt的映射
variant_to_alleles <- setNames(
  paste(target_egenes$ref, target_egenes$alt, sep = "/"),
  target_egenes$variant_id
)

# 创建gene_id到gene_name的映射
gene_id_to_name <- setNames(target_egenes$gene_name, target_egenes$gene_id)

cat("完成: 创建", length(variant_to_alleles), "个SNP的等位基因映射\n\n")

# ============================================
# 2. 读取signif_pairs.parquet文件
# ============================================
cat("步骤2: 读取signif_pairs数据...\n")

if (!requireNamespace("arrow", quietly = TRUE)) {
  cat("安装arrow包...\n")
  install.packages("arrow", repos = "https://cran.r-project.org")
}
library(arrow)

signif_pairs <- read_parquet(signif_pairs_file)
cat("signif_pairs总行数:", nrow(signif_pairs), "\n")

# 获取目标基因ID
target_gene_ids <- unique(target_egenes$gene_id)

# 筛选目标基因的SNP（使用phenotype_id列）
exposure_snps <- signif_pairs[signif_pairs$phenotype_id %in% target_gene_ids, ]
cat("目标基因SNP数:", nrow(exposure_snps), "\n\n")

if (nrow(exposure_snps) == 0) {
  cat("错误: 未找到目标基因的SNP\n")
  quit(status = 1)
}

# ============================================
# 3. 格式化暴露数据（添加ref/alt信息）
# ============================================
cat("步骤3: 格式化暴露数据...\n")

# 添加基因名
exposure_snps$gene_name <- gene_id_to_name[exposure_snps$phenotype_id]

# 添加ref/alt信息（从eGenes映射）
alleles <- variant_to_alleles[exposure_snps$variant_id]
alleles_split <- strsplit(as.character(alleles), "/")
exposure_snps$ref <- sapply(alleles_split, function(x) ifelse(length(x) >= 2, x[1], NA))
exposure_snps$alt <- sapply(alleles_split, function(x) ifelse(length(x) >= 2, x[2], NA))

# 移除无法映射的SNP
valid_snps <- !is.na(exposure_snps$ref) & !is.na(exposure_snps$alt)
exposure_snps <- exposure_snps[valid_snps, ]
cat("有效SNP数（有ref/alt）:", nrow(exposure_snps), "\n")

# 格式化暴露数据
exposure_dat <- data.frame(
  SNP = exposure_snps$variant_id,
  beta.exposure = exposure_snps$slope,
  se.exposure = exposure_snps$slope_se,
  pval.exposure = exposure_snps$pval_nominal,
  eaf.exposure = exposure_snps$af,
  effect_allele.exposure = exposure_snps$alt,
  other_allele.exposure = exposure_snps$ref,
  gene.exposure = exposure_snps$gene_name,
  stringsAsFactors = FALSE
)

cat("格式化后暴露数据行数:", nrow(exposure_dat), "\n")
cat("涉及基因:", paste(unique(exposure_dat$gene.exposure), collapse = ", "), "\n\n")

# 提取所有暴露SNP
all_snps <- unique(exposure_dat$SNP)
cat("暴露SNP总数:", length(all_snps), "\n")

# 保存SNP列表
writeLines(all_snps, file.path(output_dir, "exposure_snplist.txt"))

# ============================================
# 4. 读取结局数据
# ============================================
fast_filtered_file <- "D:/EQTL/mr_results_p5e-06/outcome_filtered_fast.csv"

if (file.exists(fast_filtered_file)) {
  cat("\n步骤4: 读取预过滤的结局文件...\n")
  outcome_small <- fread(fast_filtered_file)
} else {
  cat("\n步骤4: 错误 - 未找到预过滤文件\n")
  quit(status = 1)
}

cat("过滤后结局数据行数:", nrow(outcome_small), "\n")
cat("结局数据列名:", paste(names(outcome_small), collapse = ", "), "\n\n")

# 正确识别结局数据列（根据实际列名）
outcome_cols <- names(outcome_small)

# 查找正确的列名
find_col <- function(patterns, cols) {
  for (p in patterns) {
    if (p %in% cols) return(p)
  }
  return(NULL)
}

snp_col <- find_col(c("SNP", "rsids", "rsid"), outcome_cols)
beta_col <- find_col(c("beta.outcome", "BETA", "beta"), outcome_cols)
se_col <- find_col(c("se.outcome", "SE", "se"), outcome_cols)
ea_col <- find_col(c("effect_allele.outcome", "A1", "effect_allele", "alt"), outcome_cols)
oa_col <- find_col(c("other_allele.outcome", "A2", "other_allele", "ref"), outcome_cols)
pval_col <- find_col(c("pval.outcome", "P", "pval", "p"), outcome_cols)

if (any(sapply(list(snp_col, beta_col, se_col, ea_col, oa_col), is.null))) {
  cat("错误: 无法识别结局数据列\n")
  cat("请检查列名:\n")
  cat("  SNP:", snp_col, "\n")
  cat("  Beta:", beta_col, "\n")
  cat("  SE:", se_col, "\n")
  cat("  Effect allele:", ea_col, "\n")
  cat("  Other allele:", oa_col, "\n")
  quit(status = 1)
}

cat("结局数据列映射:\n")
cat("  SNP列:", snp_col, "\n")
cat("  Beta列:", beta_col, "\n")
cat("  SE列:", se_col, "\n")
cat("  Effect allele列:", ea_col, "\n")
cat("  Other allele列:", oa_col, "\n")
cat("  P-value列:", pval_col, "\n\n")

# 将data.table转换为data.frame并格式化
outcome_small_df <- as.data.frame(outcome_small)

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
cat("\n步骤5: 逐个基因跑MR...\n\n")

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
  
  # 格式化暴露数据
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
  dat_filtered <- dat[dat$F_stat > 10, ]
  
  if (nrow(dat_filtered) == 0) {
    cat("  跳过: F统计量均<10（弱工具变量）\n\n")
    next
  }
  
  cat("  F>10的SNP:", nrow(dat_filtered), "个\n")
  
  # 跑MR
  res <- tryCatch({
    mr(dat_filtered, method_list = c("mr_ivw", "mr_egger_regression",
                                     "mr_weighted_median", "mr_weighted_mode"))
  }, error = function(e) {
    cat("  MR分析失败:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(res) || nrow(res) == 0) {
    cat("  MR分析无结果\n\n")
    next
  }
  
  # 敏感性分析
  hetero <- tryCatch(mr_heterogeneity(dat_filtered), error = function(e) NULL)
  pleio  <- tryCatch(mr_pleiotropy_test(dat_filtered), error = function(e) NULL)
  
  # 保存结果
  results_list[[gene]] <- list(
    gene    = gene,
    dat     = dat_filtered,
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
  write.csv(dat_filtered, file.path(output_dir, paste0(gene, "_harmonised_data.csv")), row.names = FALSE)
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
if (nrow(ivw_summary) > 0) {
  write.csv(ivw_summary, file.path(output_dir, "MR_summary_IVW_only.csv"), row.names = FALSE)
}

cat("\n========== 分析完成 ==========\n")
cat("输出目录:", output_dir, "\n")
cat("成功分析基因数:", length(results_list), "\n")
if (length(results_list) > 0) {
  cat("成功分析的基因:", paste(names(results_list), collapse = ", "), "\n")
}
cat("汇总结果已保存\n")
