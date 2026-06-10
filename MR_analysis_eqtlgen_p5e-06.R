#!/usr/bin/env Rscript
# MR分析脚本 - eQTLgen p=5e-06数据集
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
exposure_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-06_kb_1000_r2_0.01.xlsx"
outcome_file  <- "D:/EQTL/eqtlgen_ieu_outcome.csv"
genes         <- c("NFKB1", "FDX1", "STAT3", "HIF1A", "HMOX1",
                   "GPX4", "HSPA5", "AGER", "DLAT")
output_dir    <- "D:/EQTL/mr_results_p5e-06"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== eQTLGen MR 分析 (p=5e-06) ==========\n")
cat("暴露文件:", exposure_file, "\n")
cat("结局文件:", outcome_file, "\n")
cat("分析基因:", paste(genes, collapse = ", "), "\n\n")

# ============================================
# 1. 读取暴露数据（Excel，较小，直接读入）
# ============================================
cat("步骤1: 读取暴露数据...\n")

# 先查看Excel的sheet名和列名
sheets <- excel_sheets(exposure_file)
cat("Excel sheets:", paste(sheets, collapse = ", "), "\n")

# 读取第一个sheet
exp_raw <- read_excel(exposure_file, sheet = 1)

cat("暴露数据原始列名:", paste(names(exp_raw), collapse = ", "), "\n")
cat("暴露数据行数:", nrow(exp_raw), "\n\n")

# 检查必要的列是否存在
required_cols <- c("SNP", "beta.exposure", "se.exposure", "pval.exposure",
                   "effect_allele.exposure", "other_allele.exposure", "gene")
missing_cols <- setdiff(required_cols, names(exp_raw))
if (length(missing_cols) > 0) {
  cat("警告: 缺失列:", paste(missing_cols, collapse = ", "), "\n")
  cat("使用备选列名...\n")
  
  # 尝试使用备选列名
  col_mapping <- list(
    SNP = c("SNP", "rsid", "RSID"),
    beta.exposure = c("beta.exposure", "BETA", "beta", "Beta"),
    se.exposure = c("se.exposure", "SE", "se"),
    pval.exposure = c("pval.exposure", "P", "pval", "p.value"),
    effect_allele.exposure = c("effect_allele.exposure", "A1", "EA", "effect_allele"),
    other_allele.exposure = c("other_allele.exposure", "A2", "OA", "other_allele"),
    eaf.exposure = c("eaf.exposure", "FRQ", "EAF", "eaf", "maf"),
    gene = c("gene", "GENE", "Gene")
  )
  
  # 找到实际使用的列名
  actual_cols <- list()
  for (col in names(col_mapping)) {
    for (candidate in col_mapping[[col]]) {
      if (candidate %in% names(exp_raw)) {
        actual_cols[[col]] <- candidate
        break
      }
    }
  }
  
  cat("实际列名映射:\n")
  for (col in names(actual_cols)) {
    cat("  ", col, "->", actual_cols[[col]], "\n")
  }
  cat("\n")
}

# 格式化暴露数据
exposure_dat <- format_data(
  dat          = exp_raw,
  type         = "exposure",
  snp_col      = ifelse("SNP" %in% names(exp_raw), "SNP", actual_cols[["SNP"]]),
  beta_col     = ifelse("beta.exposure" %in% names(exp_raw), "beta.exposure", actual_cols[["beta.exposure"]]),
  se_col       = ifelse("se.exposure" %in% names(exp_raw), "se.exposure", actual_cols[["se.exposure"]]),
  effect_allele_col = ifelse("effect_allele.exposure" %in% names(exp_raw), "effect_allele.exposure", actual_cols[["effect_allele.exposure"]]),
  other_allele_col  = ifelse("other_allele.exposure" %in% names(exp_raw), "other_allele.exposure", actual_cols[["other_allele.exposure"]]),
  eaf_col      = ifelse("eaf.exposure" %in% names(exp_raw), "eaf.exposure", actual_cols[["eaf.exposure"]]),
  pval_col     = ifelse("pval.exposure" %in% names(exp_raw), "pval.exposure", actual_cols[["pval.exposure"]]),
  gene_col     = ifelse("gene" %in% names(exp_raw), "gene", actual_cols[["gene"]])
)

# 只保留目标基因
exposure_dat <- exposure_dat[exposure_dat$gene.exposure %in% genes, ]
cat("过滤后暴露数据行数:", nrow(exposure_dat), "\n")
cat("涉及基因:", paste(unique(exposure_dat$gene.exposure), collapse = ", "), "\n\n")

# 提取所有暴露SNP（用于过滤50GB结局文件）
all_snps <- unique(exposure_dat$SNP)
cat("暴露SNP总数:", length(all_snps), "\n")

# 保存SNP列表（备用）
writeLines(all_snps, file.path(output_dir, "exposure_snplist.txt"))

# ============================================
# 2. 分块读取50GB结局文件，只保留暴露SNP
# ============================================
cat("\n步骤2: 分块过滤50GB结局文件（这可能需要10-30分钟）...\n")
cat("Chunk size: 50万行/块\n")

# 先读100行看列名
peek <- read_csv(outcome_file, n_max = 5, show_col_types = FALSE)
cat("结局文件列名:", paste(names(peek), collapse = ", "), "\n")

# 分块过滤函数
filter_callback <- function(chunk, pos) {
  # chunk 是 data.frame，带正确列名
  snp_col <- NULL
  if ("SNP" %in% names(chunk)) snp_col <- "SNP"
  else if ("rsids" %in% names(chunk)) snp_col <- "rsids"
  else if ("RSID" %in% names(chunk)) snp_col <- "RSID"
  
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
  cat("错误: 没有匹配到任何SNP，请检查SNP ID格式是否一致\n")
  cat("暴露数据SNP示例:", paste(head(all_snps, 5), collapse = ", "), "\n")
  if ("SNP" %in% names(outcome_filtered)) {
    cat("结局数据SNP示例:", paste(head(outcome_filtered$SNP, 5), collapse = ", "), "\n")
  }
  quit(status = 1)
}

cat("匹配到的SNP数:", length(unique(outcome_filtered[[1]])), "\n\n")

# 保存过滤后的结局
fwrite(as.data.table(outcome_filtered),
       file.path(output_dir, "outcome_filtered.csv"))

# 释放大对象
rm(outcome_filtered)
gc()

# ============================================
# 3. 重新读入过滤后的小结局文件（快）
# ============================================
cat("步骤3: 读入过滤后的结局数据...\n")
outcome_small <- fread(file.path(output_dir, "outcome_filtered.csv"))
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
# 4. 循环每个基因跑MR
# ============================================
cat("步骤4: 逐个基因跑MR...\n\n")

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

  # Harmonise
  dat <- tryCatch({
    harmonise_data(exposure_dat = exp_gene, outcome_dat = outcome_dat)
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
# 5. 汇总结果表
# ============================================
cat("步骤5: 生成汇总结果...\n")

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
