#!/usr/bin/env Rscript
# MR分析脚本 - 146个基因完整版, eQTLGen p=5e-08, MEGASTROKE结局
# 包含完整敏感性分析

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(TwoSampleMR)
  library(dplyr)
})

# 146个基因列表
genes <- c("LYN", "FABP2", "MB", "F3", "PTPRF", "ALDH9A1", "GCH1", "RENBP", "OAZ1", "PA2G4",
           "PRKCQ", "ATP7B", "B2M", "POLR2D", "C3", "HPGDS", "S100A6", "ACADVL", "DLAT", "CCNA2",
           "STK4", "HSD17B4", "NMT1", "RBM39", "HSD17B10", "ITGA1", "PTPRJ", "DDC", "STARD13", "PTGR1",
           "ZHX2", "ACAD11", "TDP1", "HBS1L", "MAPKAPK2", "CITED2", "CASK", "CUL4B", "CTSD", "LEF1",
           "MKNK2", "PDCD6", "MAN2B1", "CHFR", "SEC13", "HIBADH", "MGAT1", "BST1", "PDCD6IP", "SAT1",
           "SERPINB10", "PARP12", "IL10RA", "NUDCD2", "PCTP", "SAT2", "IGFBP2", "CNDP2", "PTPRC", "XRCC6",
           "ACADM", "RHOC", "TCN2", "ZEB1", "TSPO", "TOP2A", "TNF", "TGFB1", "TBXAS1", "STAT3",
           "STAT1", "SREBF1", "SCN9A", "RELA", "PTGS1", "PPARG", "PARP1", "PABPC1", "NR1H3", "NFKB1",
           "NFE2L2", "MTOR", "KCNA5", "JAK1", "IRF1", "IMPDH2", "IL6", "IKBKB", "HTR2C", "HTR2B",
           "HSPA5", "HMOX1", "HIF1A", "GPX1", "GFAP", "FLT4", "EPHX1", "EGFR", "CTSS", "CTSL",
           "CTSK", "CTSF", "CTSC", "CTSB", "CPT2", "CPT1A", "CP", "COL1A1", "CNR2", "CDK4",
           "CCND1", "CCL2", "CAT", "CASP8", "BRD3", "AKT1", "AIF1", "ADRB1", "ACTA2", "CXCR3",
           "PTPN2", "MAOB", "FABP4", "FABP5", "NR3C1", "CCR5", "PLA2G4A", "SPHK1", "FNTA", "TIMP1",
           "PTPN6", "ICAM1", "STAT5A", "XDH")

exposure_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-08_kb_1000_r2_0.01.xlsx"
outcome_file  <- "D:/EQTL/mr_results_megastroke/megastroke_outcome_146genes.csv"
output_dir    <- "D:/EQTL/mr_results_146genes_p5e-08"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=" , paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("       146基因 MR完整分析 (p=5e-08, MEGASTROKE)\n")
cat("=" , paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("基因数:", length(genes), "\n\n")

# 读取数据
cat("【步骤1】读取暴露数据...\n")
exposure_data <- read_excel(exposure_file)
exposure_data <- as.data.frame(exposure_data)

cat("  暴露数据总行数:", nrow(exposure_data), "\n")

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

cat("  目标基因SNP数:", nrow(exposure_snps), "\n\n")

# 读取结局数据
cat("【步骤2】读取结局数据...\n")
outcome_raw <- fread(outcome_file)
cat("  结局数据行数:", nrow(outcome_raw), "\n\n")

chr_col_out <- intersect(c("chr", "Chr", "CHR", "chromosome", "chrom", "chr.outcome"), names(outcome_raw))[1]
pos_col_out <- intersect(c("pos", "Pos", "POS", "position", "bp", "BP", "base_pair_location", "pos.outcome"), names(outcome_raw))[1]

outcome_raw$chr.outcome <- as.integer(outcome_raw[[chr_col_out]])
outcome_raw$pos.outcome <- as.numeric(outcome_raw[[pos_col_out]])
outcome_raw$chr_pos <- paste(outcome_raw$chr.outcome, outcome_raw$pos.outcome, sep = ":")

# 跑MR
cat("【步骤3】MR分析...\n\n")
results_list <- list()

for (gene in genes) {
  cat("----------", gene, "----------\n")
  
  exp_gene <- exposure_snps[exposure_snps$gene_name == toupper(gene), ]
  cat("  暴露SNP:", nrow(exp_gene), "个\n")
  
  if (nrow(exp_gene) == 0) {
    cat("  ⚠️ 跳过\n\n")
    next
  }
  
  exp_chr_pos <- unique(exp_gene$chr_pos)
  out_matched <- outcome_raw[outcome_raw$chr_pos %in% exp_chr_pos, ]
  cat("  匹配结局SNP:", nrow(out_matched), "个\n")
  
  if (nrow(out_matched) == 0) {
    cat("  ⚠️ 跳过\n\n")
    next
  }
  
  exp_df <- data.frame(
    chr_pos = exp_gene$chr_pos, chr = exp_gene$chr, pos = exp_gene$pos,
    beta.exp = exp_gene$beta, se.exp = exp_gene$se, pval.exp = exp_gene$pval, eaf.exp = exp_gene$eaf,
    effect_allele.exp = toupper(exp_gene$effect_allele), other_allele.exp = toupper(exp_gene$other_allele),
    SNP.exp = exp_gene$SNP, gene = exp_gene$gene_name, stringsAsFactors = FALSE
  )
  
  out_df <- data.frame(
    chr_pos = out_matched$chr_pos, chr = out_matched$chr.outcome, pos = out_matched$pos.outcome,
    beta.out = out_matched$beta.outcome, se.out = out_matched$se.outcome,
    pval.out = out_matched$pval.outcome, eaf.out = out_matched$eaf.outcome,
    effect_allele.out = toupper(out_matched$effect_allele.outcome),
    other_allele.out = toupper(out_matched$other_allele.outcome),
    SNP.out = out_matched$SNP,
    samplesize.out = ifelse("samplesize.outcome" %in% names(out_matched), out_matched$samplesize.outcome, NA),
    outcome = out_matched$outcome, stringsAsFactors = FALSE
  )
  
  exp_df <- exp_df %>% group_by(chr, pos) %>% slice_min(pval.exp, n = 1, with_ties = FALSE) %>% ungroup()
  out_df <- out_df %>% group_by(chr, pos) %>% slice_min(pval.out, n = 1, with_ties = FALSE) %>% ungroup()
  
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
    beta.exposure = merged_clean$beta.exp, beta.outcome = merged_clean$beta.outcome,
    se.exposure = merged_clean$se.exp, se.outcome = merged_clean$se.out,
    effect_allele.exposure = merged_clean$effect_allele.exp, other_allele.exposure = merged_clean$other_allele.exp,
    effect_allele.outcome = merged_clean$effect_allele.exp, other_allele.outcome = merged_clean$other_allele.exp,
    eaf.exposure = merged_clean$eaf.exp, eaf.outcome = merged_clean$eaf.outcome,
    pval.exposure = merged_clean$pval.exp, pval.outcome = merged_clean$pval.out,
    samplesize.outcome = merged_clean$samplesize.out,
    id.exposure = merged_clean$gene, exposure = merged_clean$gene,
    outcome = "Ischemic Stroke", id.outcome = "ebi-a-GCST006908",
    mr_keep = TRUE, stringsAsFactors = FALSE
  )
  
  dat <- dat[complete.cases(dat[, c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")]), ]
  cat("  Harmonised SNP:", nrow(dat), "个\n")
  
  if (nrow(dat) < 1) next
  
  dat$F_stat <- (dat$beta.exposure / dat$se.exposure)^2
  dat <- dat[dat$F_stat > 10, ]
  if (nrow(dat) > 0) {
    cat("  F>10 SNP:", nrow(dat), "个\n")
  } else {
    cat("  ⚠️ 跳过\n\n")
    next
  }
  
  res <- tryCatch({
    mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  }, error = function(e) NULL)
  
  if (is.null(res) || nrow(res) == 0) {
    cat("  ⚠️ MR无结果\n\n")
    next
  }
  
  # 敏感性分析
  cat("  敏感性分析:\n")
  hetero <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  loo <- tryCatch(mr_leaveoneout(dat), error = function(e) NULL)
  
  if (!is.null(hetero) && nrow(hetero) > 0) {
    ivw_q <- hetero$Q_pval[hetero$method == "Inverse variance weighted"]
    if (length(ivw_q) > 0 && !is.na(ivw_q[1])) {
      hetero_flag <- ifelse(ivw_q[1] < 0.05, "⚠️ 显著异质性", "✓ 无异质性")
      cat("    - 异质性:", hetero_flag, "\n")
    }
  }
  
  if (!is.null(pleio) && nrow(pleio) > 0) {
    pleio_p <- pleio$pval[1]
    pleio_flag <- ifelse(pleio_p < 0.05, "⚠️ 多效性", "✓ 无多效性")
    cat("    - 多效性:", pleio_flag, "\n")
  }
  
  results_list[[gene]] <- list(gene = gene, dat = dat, res = res, hetero = hetero, pleio = pleio, loo = loo)
  
  ivw <- res[res$method == "Inverse variance weighted", ]
  if (nrow(ivw) > 0 && !is.na(ivw$b[1])) {
    sig_marker <- ifelse(ivw$pval[1] < 0.05, " ***", "")
    cat("  IVW: OR =", round(exp(ivw$b[1]), 3), ", P =", format(ivw$pval[1], digits = 3, scientific = TRUE), sig_marker, "\n")
  }
  cat("\n")
  
  write.csv(res, file.path(output_dir, paste0(gene, "_MR_results.csv")), row.names = FALSE)
  write.csv(dat, file.path(output_dir, paste0(gene, "_harmonised_data.csv")), row.names = FALSE)
  if (!is.null(hetero)) write.csv(hetero, file.path(output_dir, paste0(gene, "_hetero.csv")), row.names = FALSE)
  if (!is.null(pleio)) write.csv(pleio, file.path(output_dir, paste0(gene, "_pleio.csv")), row.names = FALSE)
}

# 汇总
cat("=" , paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("                    汇总结果\n")
cat("=" , paste(rep("=", 80), collapse = ""), "\n", sep = "")

if (length(results_list) > 0) {
  summary_df <- data.frame(
    Gene = character(), Method = character(), Beta = numeric(), SE = numeric(),
    OR = numeric(), P_value = numeric(), SNP_n = integer(), stringsAsFactors = FALSE
  )
  
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
    
    cat("\nIVW结果 (按P值排序, 前20个):\n")
    ivw_sorted <- ivw_summary[order(ivw_summary$P_value), ]
    display_n <- min(20, nrow(ivw_sorted))
    print(ivw_sorted[1:display_n, c("Gene", "OR", "P_value")], row.names = FALSE)
    
    sig_genes <- ivw_sorted[ivw_sorted$P_value < 0.05, ]
    if (nrow(sig_genes) > 0) {
      cat("\n\n显著关联基因 (P < 0.05):\n")
      print(sig_genes[, c("Gene", "OR", "P_value")], row.names = FALSE)
    }
  }
  
  cat("\n\n成功分析基因数:", length(results_list), "/", length(genes), "\n")
  cat("结果保存至:", output_dir, "\n")
}

cat("=" , paste(rep("=", 80), collapse = ""), "\n", sep = "")
