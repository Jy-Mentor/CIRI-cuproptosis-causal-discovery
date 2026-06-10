#!/usr/bin/env Rscript
# MR分析脚本 - eQTLGen p=5e-08数据集 (修复版)
# 使用chr/pos/allele匹配，MEGASTROKE结局数据

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(TwoSampleMR)
  library(dplyr)
})

# ============================================
# 参数设置
# ============================================
exposure_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-08_kb_1000_r2_0.01.xlsx"
outcome_file  <- "D:/EQTL/mr_results_megastroke/megastroke_outcome.csv"
genes         <- c("NFKB1", "STAT3", "HIF1A", "HSPA5", "HMOX1",
                   "RELA", "NFE2L2", "CP", "LIAS", "IKBKB",
                   "JAK1", "PARP1", "CASP8", "MTOR", "PTPRC")
output_dir    <- "D:/EQTL/mr_results_p5e-08_megastroke"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== eQTLGen p=5e-08 MR 分析（MEGASTROKE结局） ==========\n")
cat("暴露文件:", exposure_file, "\n")
cat("结局文件:", outcome_file, "\n\n")

# ============================================
# 1. 读取并处理暴露数据
# ============================================
cat("步骤1: 读取暴露数据...\n")

exposure_data <- read_excel(exposure_file)
exposure_data <- as.data.frame(exposure_data)

cat("暴露数据总行数:", nrow(exposure_data), "\n")

# 筛选目标基因并创建标准列
exposure_snps <- exposure_data[toupper(exposure_data$gene) %in% toupper(genes), ]
exposure_snps$chr <- as.integer(exposure_snps$chr.exposure)
exposure_snps$pos <- as.numeric(exposure_snps$pos.exposure)
exposure_snps$gene_name <- toupper(as.character(exposure_snps$gene))
exposure_snps$beta <- as.numeric(exposure_snps$beta.exposure)
exposure_snps$se <- as.numeric(exposure_snps$se.exposure)
exposure_snps$pval <- as.numeric(exposure_snps$pval.exposure)
exposure_snps$eaf <- as.numeric(exposure_snps$eaf.exposure)
exposure_snps$effect_allele <- exposure_snps$effect_allele.exposure
exposure_snps$other_allele <- exposure_snps$other_allele.exposure
exposure_snps$SNP <- exposure_snps$SNP

# 移除NA
valid <- !is.na(exposure_snps$chr) & !is.na(exposure_snps$pos) & 
         !is.na(exposure_snps$beta) & !is.na(exposure_snps$se) &
         !is.na(exposure_snps$effect_allele) & !is.na(exposure_snps$other_allele)
exposure_snps <- exposure_snps[valid, ]

# 创建chr_pos
exposure_snps$chr_pos <- paste(exposure_snps$chr, exposure_snps$pos, sep = ":")

cat("有效SNP数:", nrow(exposure_snps), "\n")
cat("唯一chr:pos数:", length(unique(exposure_snps$chr_pos)), "\n\n")

# ============================================
# 2. 读取并处理结局数据（带列名自动检测）
# ============================================
cat("步骤2: 读取结局数据...\n")

outcome_raw <- fread(outcome_file)
cat("结局数据行数:", nrow(outcome_raw), "\n")
cat("原始列名:", paste(names(outcome_raw), collapse = ", "), "\n")

# 自动检测列名（支持标准TwoSampleMR格式和自定义格式）
chr_col_out <- intersect(c("chr", "Chr", "CHR", "chromosome", "chrom", "chr.outcome"), names(outcome_raw))[1]
pos_col_out <- intersect(c("pos", "Pos", "POS", "position", "bp", "BP", "base_pair_location", "pos.outcome"), names(outcome_raw))[1]

if (is.na(chr_col_out) || is.na(pos_col_out)) {
  cat("可用列名:", paste(names(outcome_raw), collapse = ", "), "\n")
  stop("错误：无法识别结局数据的染色体/位置列名")
}

cat("检测到chr列:", chr_col_out, ", pos列:", pos_col_out, "\n")

outcome_raw$chr.outcome <- as.integer(outcome_raw[[chr_col_out]])
outcome_raw$pos.outcome <- as.numeric(outcome_raw[[pos_col_out]])
outcome_raw$chr_pos <- paste(outcome_raw$chr.outcome, outcome_raw$pos.outcome, sep = ":")

cat("唯一chr:pos数:", length(unique(outcome_raw$chr_pos)), "\n\n")

# ============================================
# 3. 循环跑MR
# ============================================
cat("步骤3: 逐个基因跑MR...\n\n")

results_list <- list()

for (gene in genes) {
  cat("==========", gene, "==========\n")
  
  # 提取该基因暴露数据
  exp_gene <- exposure_snps[exposure_snps$gene_name == toupper(gene), ]
  cat("  暴露SNP:", nrow(exp_gene), "个\n")
  
  if (nrow(exp_gene) == 0) {
    cat("  跳过: 无暴露数据\n\n")
    next
  }
  
  # chr:pos匹配
  exp_chr_pos <- unique(exp_gene$chr_pos)
  cat("  暴露chr:pos:", length(exp_chr_pos), "个\n")
  
  out_matched <- outcome_raw[outcome_raw$chr_pos %in% exp_chr_pos, ]
  cat("  匹配到结局SNP:", nrow(out_matched), "个\n")
  
  if (nrow(out_matched) == 0) {
    cat("  跳过: 无匹配结局\n\n")
    next
  }
  
  # 准备合并数据框
  exp_df <- data.frame(
    chr_pos = exp_gene$chr_pos,
    chr = exp_gene$chr,
    pos = exp_gene$pos,
    beta.exp = exp_gene$beta,
    se.exp = exp_gene$se,
    pval.exp = exp_gene$pval,
    eaf.exp = exp_gene$eaf,
    effect_allele.exp = toupper(exp_gene$effect_allele),
    other_allele.exp = toupper(exp_gene$other_allele),
    SNP.exp = exp_gene$SNP,
    gene = exp_gene$gene_name,
    stringsAsFactors = FALSE
  )
  
  out_df <- data.frame(
    chr_pos = out_matched$chr_pos,
    chr = out_matched$chr.outcome,
    pos = out_matched$pos.outcome,
    beta.out = out_matched$beta.outcome,
    se.out = out_matched$se.outcome,
    pval.out = out_matched$pval.outcome,
    eaf.out = out_matched$eaf.outcome,
    effect_allele.out = toupper(out_matched$effect_allele.outcome),
    other_allele.out = toupper(out_matched$other_allele.outcome),
    SNP.out = out_matched$SNP,
    samplesize.out = ifelse("samplesize.outcome" %in% names(out_matched), 
                            out_matched$samplesize.outcome, NA),
    outcome = out_matched$outcome,
    stringsAsFactors = FALSE
  )
  
  # 🔴 去重：防止笛卡尔积
  exp_df <- exp_df %>% group_by(chr, pos) %>% slice_min(pval.exp, n = 1, with_ties = FALSE) %>% ungroup()
  out_df <- out_df %>% group_by(chr, pos) %>% slice_min(pval.out, n = 1, with_ties = FALSE) %>% ungroup()
  cat("  去重后暴露SNP:", nrow(exp_df), "个, 结局SNP:", nrow(out_df), "个\n")
  
  # 合并
  merged <- merge(exp_df, out_df, by = "chr_pos", suffixes = c(".exp", ".out"))
  cat("  合并后:", nrow(merged), "行\n")
  
  if (nrow(merged) == 0) {
    cat("  跳过: 合并失败\n\n")
    next
  }
  
  # 等位基因匹配检查
  merged$allele_match <- with(merged,
    (effect_allele.exp == effect_allele.out & other_allele.exp == other_allele.out) |
    (effect_allele.exp == other_allele.out & other_allele.exp == effect_allele.out)
  )
  
  cat("  等位基因匹配:", sum(merged$allele_match), "/", nrow(merged), "\n")
  
  merged_clean <- merged[merged$allele_match, ]
  
  if (nrow(merged_clean) == 0) {
    cat("  跳过: 无等位基因匹配\n\n")
    next
  }
  
  # 处理等位基因翻转
  need_flip <- merged_clean$effect_allele.exp != merged_clean$effect_allele.out
  
  merged_clean$beta.outcome <- ifelse(need_flip, -merged_clean$beta.out, merged_clean$beta.out)
  merged_clean$eaf.outcome <- ifelse(need_flip, 1 - merged_clean$eaf.out, merged_clean$eaf.out)
  merged_clean$effect_allele.outcome <- merged_clean$effect_allele.exp
  merged_clean$other_allele.outcome <- merged_clean$other_allele.exp
  
  # 构建harmonised数据框
  dat <- data.frame(
    SNP = merged_clean$SNP.out,
    beta.exposure = merged_clean$beta.exp,
    beta.outcome = merged_clean$beta.outcome,
    se.exposure = merged_clean$se.exp,
    se.outcome = merged_clean$se.out,
    effect_allele.exposure = merged_clean$effect_allele.exp,
    other_allele.exposure = merged_clean$other_allele.exp,
    effect_allele.outcome = merged_clean$effect_allele.outcome,
    other_allele.outcome = merged_clean$other_allele.outcome,
    eaf.exposure = merged_clean$eaf.exp,
    eaf.outcome = merged_clean$eaf.outcome,
    pval.exposure = merged_clean$pval.exp,
    pval.outcome = merged_clean$pval.out,
    samplesize.outcome = merged_clean$samplesize.out,
    gene.exposure = merged_clean$gene,
    id.exposure = merged_clean$gene,
    exposure = merged_clean$gene,
    outcome = "Ischemic Stroke",
    id.outcome = "ebi-a-GCST006908",
    mr_keep = TRUE,
    stringsAsFactors = FALSE
  )
  
  # 清洗NA
  dat <- dat[complete.cases(dat[, c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")]), ]
  cat("  Harmonised SNP:", nrow(dat), "个\n")
  
  if (nrow(dat) < 1) {
    cat("  跳过: SNP不足\n\n")
    next
  }
  
  # F统计量过滤
  dat$F_stat <- (dat$beta.exposure / dat$se.exposure)^2
  dat <- dat[dat$F_stat > 10, ]
  
  if (nrow(dat) == 0) {
    cat("  跳过: 弱工具变量\n\n")
    next
  }
  
  cat("  F>10 SNP:", nrow(dat), "个")
  cat(" (F范围:", round(min(dat$F_stat), 1), "-", round(max(dat$F_stat), 1), ")\n")
  
  # 跑MR
  res <- tryCatch({
    mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                             "mr_weighted_median", "mr_weighted_mode"))
  }, error = function(e) {
    cat("  MR错误:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(res) || nrow(res) == 0) {
    cat("  MR无结果\n\n")
    next
  }
  
  # 敏感性分析
  hetero <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  pleio  <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  loo    <- tryCatch(mr_leaveoneout(dat), error = function(e) NULL)
  
  # 保存结果
  results_list[[gene]] <- list(
    gene = gene,
    dat = dat,
    res = res,
    hetero = hetero,
    pleio = pleio,
    loo = loo
  )
  
  # 打印IVW结果
  ivw <- res[res$method == "Inverse variance weighted", ]
  if (nrow(ivw) > 0 && !is.na(ivw$b[1])) {
    cat("  IVW Beta:", round(ivw$b[1], 4),
        "| OR:", round(exp(ivw$b[1]), 3),
        "| P:", format(ivw$pval[1], digits = 3, scientific = TRUE), "\n")
  }
  cat("\n")
  
  # 保存文件
  write.csv(res, file.path(output_dir, paste0(gene, "_MR_results.csv")), row.names = FALSE)
  write.csv(dat, file.path(output_dir, paste0(gene, "_harmonised_data.csv")), row.names = FALSE)
  if (!is.null(hetero)) write.csv(hetero, file.path(output_dir, paste0(gene, "_hetero.csv")), row.names = FALSE)
  if (!is.null(pleio)) write.csv(pleio, file.path(output_dir, paste0(gene, "_pleio.csv")), row.names = FALSE)
  if (!is.null(loo)) write.csv(loo, file.path(output_dir, paste0(gene, "_loo.csv")), row.names = FALSE)
}

# ============================================
# 4. 汇总结果
# ============================================
cat("\n========== 汇总结果 ==========\n")

if (length(results_list) > 0) {
  summary_df <- data.frame(
    Gene = character(),
    Method = character(),
    Beta = numeric(),
    SE = numeric(),
    OR = numeric(),
    P_value = numeric(),
    F_stat_mean = numeric(),
    SNP_n = integer(),
    stringsAsFactors = FALSE
  )
  
  for (gene in names(results_list)) {
    res <- results_list[[gene]]$res
    dat <- results_list[[gene]]$dat
    
    for (i in seq_len(nrow(res))) {
      summary_df <- rbind(summary_df, data.frame(
        Gene = gene,
        Method = res$method[i],
        Beta = res$b[i],
        SE = res$se[i],
        OR = exp(res$b[i]),
        P_value = res$pval[i],
        F_stat_mean = mean(dat$F_stat),
        SNP_n = nrow(dat)
      ))
    }
  }
  
  write.csv(summary_df, file.path(output_dir, "MR_summary_all_genes.csv"), row.names = FALSE)
  
  # IVW汇总
  ivw_summary <- summary_df[summary_df$Method == "Inverse variance weighted", ]
  if (nrow(ivw_summary) > 0) {
    write.csv(ivw_summary, file.path(output_dir, "MR_summary_IVW_only.csv"), row.names = FALSE)
    cat("\nIVW结果:\n")
    print(ivw_summary[, c("Gene", "Beta", "OR", "P_value", "SNP_n")])
  }
  
  cat("\n成功分析基因数:", length(results_list), "\n")
  cat("结果保存至:", output_dir, "\n")
} else {
  cat("无成功分析的基因\n")
}
