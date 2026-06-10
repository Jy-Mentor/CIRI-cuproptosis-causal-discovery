#!/usr/bin/env Rscript
# MR分析脚本 - eQTLGen p=5e-06数据集
suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(TwoSampleMR)
  library(dplyr)
})

exposure_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-06_kb_1000_r2_0.01.xlsx"
outcome_file  <- "D:/EQTL/mr_results_megastroke/megastroke_outcome.csv"
genes         <- c("NFKB1", "STAT3", "HIF1A", "HSPA5", "HMOX1",
                   "RELA", "NFE2L2", "CP", "LIAS", "IKBKB",
                   "JAK1", "PARP1", "CASP8", "MTOR", "PTPRC")
output_dir    <- "D:/EQTL/mr_results_p5e-06_megastroke"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== eQTLGen p=5e-06 MR 分析（MEGASTROKE结局） ==========\n\n")

# 读取暴露数据
cat("步骤1: 读取暴露数据...\n")
exposure_data <- read_excel(exposure_file)
exposure_data <- as.data.frame(exposure_data)

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

valid <- !is.na(exposure_snps$chr) & !is.na(exposure_snps$pos) & 
         !is.na(exposure_snps$beta) & !is.na(exposure_snps$se) &
         !is.na(exposure_snps$effect_allele) & !is.na(exposure_snps$other_allele)
exposure_snps <- exposure_snps[valid, ]
exposure_snps$chr_pos <- paste(exposure_snps$chr, exposure_snps$pos, sep = ":")

cat("有效SNP数:", nrow(exposure_snps), "\n\n")

# 读取结局数据
cat("步骤2: 读取结局数据...\n")
outcome_raw <- fread(outcome_file)
outcome_raw$chr.outcome <- as.integer(outcome_raw$chr)
outcome_raw$pos.outcome <- as.numeric(outcome_raw$pos.outcome)
outcome_raw$chr_pos <- paste(outcome_raw$chr.outcome, outcome_raw$pos.outcome, sep = ":")
cat("结局数据行数:", nrow(outcome_raw), "\n\n")

# 跑MR
cat("步骤3: 逐个基因跑MR...\n\n")
results_list <- list()

for (gene in genes) {
  cat("==========", gene, "==========\n")
  
  exp_gene <- exposure_snps[exposure_snps$gene_name == toupper(gene), ]
  cat("  暴露SNP:", nrow(exp_gene), "个\n")
  
  if (nrow(exp_gene) == 0) {
    cat("  跳过: 无暴露数据\n\n")
    next
  }
  
  exp_chr_pos <- unique(exp_gene$chr_pos)
  out_matched <- outcome_raw[outcome_raw$chr_pos %in% exp_chr_pos, ]
  cat("  匹配到结局SNP:", nrow(out_matched), "个\n")
  
  if (nrow(out_matched) == 0) {
    cat("  跳过: 无匹配结局\n\n")
    next
  }
  
  # 准备数据框
  exp_df <- data.frame(
    chr_pos = exp_gene$chr_pos,
    beta.exp = exp_gene$beta, se.exp = exp_gene$se,
    pval.exp = exp_gene$pval, eaf.exp = exp_gene$eaf,
    effect_allele.exp = toupper(exp_gene$effect_allele),
    other_allele.exp = toupper(exp_gene$other_allele),
    SNP.exp = exp_gene$SNP, gene = exp_gene$gene_name,
    stringsAsFactors = FALSE
  )
  
  out_df <- data.frame(
    chr_pos = out_matched$chr_pos,
    beta.out = out_matched$beta.outcome, se.out = out_matched$se.outcome,
    pval.out = out_matched$pval.outcome, eaf.out = out_matched$eaf.outcome,
    effect_allele.out = toupper(out_matched$effect_allele.outcome),
    other_allele.out = toupper(out_matched$other_allele.outcome),
    SNP.out = out_matched$SNP, samplesize.out = out_matched$samplesize.outcome,
    outcome = out_matched$outcome, stringsAsFactors = FALSE
  )
  
  merged <- merge(exp_df, out_df, by = "chr_pos")
  cat("  合并后:", nrow(merged), "行\n")
  
  if (nrow(merged) == 0) next
  
  merged$allele_match <- with(merged,
    (effect_allele.exp == effect_allele.out & other_allele.exp == other_allele.out) |
    (effect_allele.exp == other_allele.out & other_allele.exp == effect_allele.out)
  )
  merged_clean <- merged[merged$allele_match, ]
  cat("  等位基因匹配:", nrow(merged_clean), "个\n")
  
  if (nrow(merged_clean) == 0) next
  
  need_flip <- merged_clean$effect_allele.exp != merged_clean$effect_allele.out
  merged_clean$beta.outcome <- ifelse(need_flip, -merged_clean$beta.out, merged_clean$beta.out)
  merged_clean$eaf.outcome <- ifelse(need_flip, 1 - merged_clean$eaf.out, merged_clean$eaf.out)
  
  dat <- data.frame(
    SNP = merged_clean$SNP.out,
    beta.exposure = merged_clean$beta.exp,
    beta.outcome = merged_clean$beta.outcome,
    se.exposure = merged_clean$se.exp,
    se.outcome = merged_clean$se.out,
    effect_allele.exposure = merged_clean$effect_allele.exp,
    other_allele.exposure = merged_clean$other_allele.exp,
    effect_allele.outcome = merged_clean$effect_allele.exp,
    other_allele.outcome = merged_clean$other_allele.exp,
    eaf.exposure = merged_clean$eaf.exp,
    eaf.outcome = merged_clean$eaf.outcome,
    pval.exposure = merged_clean$pval.exp,
    pval.outcome = merged_clean$pval.out,
    samplesize.outcome = merged_clean$samplesize.out,
    id.exposure = merged_clean$gene,
    exposure = merged_clean$gene,
    outcome = "Ischemic Stroke",
    id.outcome = "ebi-a-GCST006908",
    mr_keep = TRUE, stringsAsFactors = FALSE
  )
  
  dat <- dat[complete.cases(dat[, c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")]), ]
  cat("  Harmonised SNP:", nrow(dat), "个\n")
  
  if (nrow(dat) < 1) next
  
  dat$F_stat <- (dat$beta.exposure / dat$se.exposure)^2
  dat <- dat[dat$F_stat > 10, ]
  cat("  F>10 SNP:", nrow(dat), "个\n")
  
  if (nrow(dat) == 0) next
  
  res <- tryCatch({
    mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                             "mr_weighted_median", "mr_weighted_mode"))
  }, error = function(e) NULL)
  
  if (is.null(res) || nrow(res) == 0) {
    cat("  MR无结果\n\n")
    next
  }
  
  hetero <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  pleio  <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  
  results_list[[gene]] <- list(gene = gene, dat = dat, res = res, hetero = hetero, pleio = pleio)
  
  ivw <- res[res$method == "Inverse variance weighted", ]
  if (nrow(ivw) > 0 && !is.na(ivw$b[1])) {
    cat("  IVW Beta:", round(ivw$b[1], 4), "| OR:", round(exp(ivw$b[1]), 3), 
        "| P:", format(ivw$pval[1], digits = 3, scientific = TRUE), "\n")
  }
  cat("\n")
  
  write.csv(res, file.path(output_dir, paste0(gene, "_MR_results.csv")), row.names = FALSE)
  write.csv(dat, file.path(output_dir, paste0(gene, "_harmonised_data.csv")), row.names = FALSE)
  if (!is.null(hetero)) write.csv(hetero, file.path(output_dir, paste0(gene, "_hetero.csv")), row.names = FALSE)
  if (!is.null(pleio)) write.csv(pleio, file.path(output_dir, paste0(gene, "_pleio.csv")), row.names = FALSE)
}

# 汇总
cat("\n========== 汇总结果 ==========\n")
if (length(results_list) > 0) {
  summary_df <- data.frame(Gene = character(), Method = character(), Beta = numeric(),
                           SE = numeric(), OR = numeric(), P_value = numeric(),
                           SNP_n = integer(), stringsAsFactors = FALSE)
  for (gene in names(results_list)) {
    res <- results_list[[gene]]$res
    dat <- results_list[[gene]]$dat
    for (i in seq_len(nrow(res))) {
      summary_df <- rbind(summary_df, data.frame(
        Gene = gene, Method = res$method[i], Beta = res$b[i], SE = res$se[i],
        OR = exp(res$b[i]), P_value = res$pval[i], SNP_n = nrow(dat)
      ))
    }
  }
  write.csv(summary_df, file.path(output_dir, "MR_summary_all_genes.csv"), row.names = FALSE)
  ivw_summary <- summary_df[summary_df$Method == "Inverse variance weighted", ]
  if (nrow(ivw_summary) > 0) {
    write.csv(ivw_summary, file.path(output_dir, "MR_summary_IVW_only.csv"), row.names = FALSE)
    cat("\nIVW结果:\n")
    print(ivw_summary[, c("Gene", "Beta", "OR", "P_value", "SNP_n")])
  }
  cat("\n成功分析基因数:", length(results_list), "\n")
}
cat("结果保存至:", output_dir, "\n")
