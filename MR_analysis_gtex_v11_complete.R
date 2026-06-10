#!/usr/bin/env Rscript
# MR分析脚本 - GTEx v11 Whole Blood数据集（完整版）
# 使用eGenes文件获取完整SNP信息（含ref/alt）

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
outcome_file  <- "D:/EQTL/eqtlgen_ieu_outcome.csv"
genes         <- c("NFKB1", "FDX1", "STAT3", "HIF1A", "HMOX1",
                   "GPX4", "HSPA5", "AGER", "DLAT")
output_dir    <- "D:/EQTL/mr_results_gtex_v11"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== GTEx v11 Whole Blood MR 分析（完整版） ==========\n")
cat("eGenes文件:", egenes_file, "\n")
cat("结局文件:", outcome_file, "\n")
cat("分析基因:", paste(genes, collapse = ", "), "\n\n")

# ============================================
# 1. 读取GTEx eGenes数据并格式化
# ============================================
cat("步骤1: 读取GTEx eGenes数据...\n")

# 读取eGenes文件
egenes_data <- fread(egenes_file, header = TRUE)
cat("eGenes数据行数:", nrow(egenes_data), "\n")
cat("eGenes列名:", paste(names(egenes_data), collapse = ", "), "\n\n")

# 筛选目标基因（使用gene_name列）
target_egenes <- egenes_data[gene_name %in% genes, ]
cat("目标基因eGenes数:", nrow(target_egenes), "\n")

if (nrow(target_egenes) == 0) {
  cat("错误: 未找到目标基因\n")
  quit(status = 1)
}

# 格式化暴露数据（使用eGenes中的完整信息）
# eGenes列名包含: variant_id, ref, alt, af, pval_nominal, slope, slope_se
exposure_dat <- data.frame(
  SNP = target_egenes$variant_id,
  beta.exposure = target_egenes$slope,
  se.exposure = target_egenes$slope_se,
  pval.exposure = target_egenes$pval_nominal,
  eaf.exposure = target_egenes$af,
  effect_allele.exposure = target_egenes$alt,
  other_allele.exposure = target_egenes$ref,
  gene.exposure = target_egenes$gene_name,
  stringsAsFactors = FALSE
)

cat("格式化后暴露数据行数:", nrow(exposure_dat), "\n")
cat("涉及基因:", paste(unique(exposure_dat$gene.exposure), collapse = ", "), "\n\n")

# 检查是否有缺失值
missing_info <- exposure_dat[is.na(exposure_dat$effect_allele.exposure) | 
                              is.na(exposure_dat$other_allele.exposure), ]
if (nrow(missing_info) > 0) {
  cat("警告: 有", nrow(missing_info), "个SNP缺少等位基因信息\n")
  # 移除这些SNP
  exposure_dat <- exposure_dat[!(is.na(exposure_dat$effect_allele.exposure) | 
                                  is.na(exposure_dat$other_allele.exposure)), ]
  cat("移除后暴露数据行数:", nrow(exposure_dat), "\n\n")
}

# 提取所有暴露SNP
all_snps <- unique(exposure_dat$SNP)
cat("暴露SNP总数:", length(all_snps), "\n")

# 保存SNP列表
writeLines(all_snps, file.path(output_dir, "exposure_snplist.txt"))

# ============================================
# 2. 读取结局数据（优先使用Python预过滤的文件）
# ============================================
fast_filtered_file <- "D:/EQTL/mr_results_p5e-06/outcome_filtered_fast.csv"

if (file.exists(fast_filtered_file)) {
  cat("\n步骤2: 检测到Python预过滤的结局文件，直接读取...\n")
  cat("文件:", fast_filtered_file, "\n")
  outcome_small <- fread(fast_filtered_file)
  cat("过滤后结局数据行数:", nrow(outcome_small), "\n\n")
} else {
  cat("\n步骤2: 未找到预过滤文件...\n")
  quit(status = 1)
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

cat("结局数据列映射:\n")
cat("  SNP列:", snp_col, "\n")
cat("  Beta列:", beta_col, "\n")
cat("  SE列:", se_col, "\n")
cat("  Effect allele列:", ea_col, "\n")
cat("  Other allele列:", oa_col, "\n\n")

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
# 3. 循环每个基因跑MR
# ============================================
cat("步骤3: 逐个基因跑MR...\n\n")

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
# 4. 汇总结果表
# ============================================
cat("步骤4: 生成汇总结果...\n")

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
cat("成功分析基因数:", length(results_list), "\n")
cat("汇总结果已保存\n")
