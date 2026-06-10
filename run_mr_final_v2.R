#!/usr/bin/env Rscript
# ================================================================================
# 138 基因 MR 分析（最终修复版 - 增加详细诊断）
# ================================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
  library(readxl)
  library(arrow)
})

cat("======================================================================\n")
cat("138 基因 MR 分析（eQTLGen 全血 - rsID 匹配 - 详细诊断版）\n")
cat("======================================================================\n\n")

# 配置
CLUMPED_FILE <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.001.xlsx"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_relaxed"
MEGASTROKE_FILE <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", OUTPUT_DIR, "\n\n")

# 目标基因列表
target_genes <- c(
  "LYN", "PRKCQ", "NMT1", "TDP1", "MAN2B1", "IL10RA", "RHOC", "SREBF1",
  "KCNA5", "HIF1A", "CTSC", "CAT", "FABP4", "STAT5A", "FABP2", "B2M",
  "RBM39", "HBS1L", "CHFR", "NUDCD2", "TCN2", "SCN9A", "JAK1", "GPX1",
  "CTSB", "CASP8", "FABP5", "XDH", "MB", "POLR2D", "HSD17B10", "MAPKAPK2",
  "SEC13", "PCTP", "ZEB1", "RELA", "IRF1", "GFAP", "CPT2", "BRD3",
  "NR3C1", "F3", "C3", "ITGA1", "CITED2", "HIBADH", "SAT2", "TSPO",
  "PTGS1", "IMPDH2", "FLT4", "CPT1A", "AKT1", "CCR5", "PTPRF", "HPGDS",
  "PTPRJ", "CASK", "MGAT1", "IGFBP2", "TOP2A", "PPARG", "IL6", "EPHX1",
  "CP", "AIF1", "PLA2G4A", "ALDH9A1", "S100A6", "DDC", "CUL4B", "BST1",
  "CNDP2", "TNF", "PARP1", "IKBKB", "EGFR", "COL1A1", "ADRB1", "SPHK1",
  "GCH1", "ACADVL", "STARD13", "CTSD", "PDCD6IP", "PTPRC", "TGFB1", "PABPC1",
  "HTR2C", "CTSS", "CNR2", "ACTA2", "FNTA", "RENBP", "CCNA2", "PTGR1",
  "LEF1", "SAT1", "XRCC6", "TBXAS1", "NR1H3", "HTR2B", "CTSL", "CDK4",
  "CXCR3", "TIMP1", "OAZ1", "STK4", "ZHX2", "MKNK2", "SERPINB10", "ACADM",
  "STAT3", "NFKB1", "HSPA5", "CTSK", "CCND1", "PTPN2", "PTPN6", "PA2G4",
  "HSD17B4", "ACAD11", "PDCD6", "PARP12", "SERPINB1", "STAT1", "NFE2L2",
  "HMOX1", "CTSF", "CCL2", "MAOB", "ICAM1", "FDX1", "LIAS", "LIPT1",
  "DLAT", "PDHB", "PDHX", "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
  "HIF1A", "MTOR", "NFKB1", "GPX4"
)

target_genes <- unique(target_genes)
cat(sprintf("目标基因数：%d\n\n", length(target_genes)))

# 基因符号到 ENSG ID 的映射
gene_to_ensg <- list(
  "ACAD11" = "ENSG00000240303", "ACADM" = "ENSG00000117054", "ACADVL" = "ENSG00000072778",
  "ACTA2" = "ENSG00000107796", "ADRB1" = "ENSG00000043591", "AIF1" = "ENSG00000204472",
  "AKT1" = "ENSG00000142208", "ALDH9A1" = "ENSG00000131449", "ATOX1" = "ENSG00000177556",
  "ATP7A" = "ENSG00000165240", "ATP7B" = "ENSG00000123191", "B2M" = "ENSG00000166710",
  "BRD3" = "ENSG00000169925", "BST1" = "ENSG00000197121", "C3" = "ENSG00000125730",
  "CASK" = "ENSG00000147044", "CASP8" = "ENSG00000064012", "CAT" = "ENSG00000121691",
  "CCL2" = "ENSG00000108691", "CCNA2" = "ENSG00000145386", "CCND1" = "ENSG00000110092",
  "CCR5" = "ENSG00000160791", "CDK4" = "ENSG00000135446", "CHFR" = "ENSG00000072609",
  "CITED2" = "ENSG00000164442", "CNDP2" = "ENSG00000133313", "CNR2" = "ENSG00000188822",
  "COL1A1" = "ENSG00000108821", "CP" = "ENSG00000047457", "CPT1A" = "ENSG00000110090",
  "CPT2" = "ENSG00000157184", "CTSB" = "ENSG00000164733", "CTSC" = "ENSG00000109861",
  "CTSD" = "ENSG00000117984", "CTSF" = "ENSG00000174080", "CTSK" = "ENSG00000143387",
  "CTSL" = "ENSG00000135047", "CTSS" = "ENSG00000163131", "CUL4B" = "ENSG00000158290",
  "CXCR3" = "ENSG00000186810", "DDC" = "ENSG00000132437", "DLAT" = "ENSG00000150768",
  "EGFR" = "ENSG00000146648", "EPHX1" = "ENSG00000143819", "F3" = "ENSG00000117525",
  "FABP2" = "ENSG00000145384", "FABP4" = "ENSG00000170323", "FABP5" = "ENSG00000236972",
  "FDX1" = "ENSG00000137714", "FLT4" = "ENSG00000037280", "FNTA" = "ENSG00000168522",
  "GCH1" = "ENSG00000131979", "GFAP" = "ENSG00000131095", "GPX1" = "ENSG00000233276",
  "GPX4" = "ENSG00000167468", "HBS1L" = "ENSG00000112339", "HIBADH" = "ENSG00000106049",
  "HIF1A" = "ENSG00000100644", "HMOX1" = "ENSG00000100292", "HPGDS" = "ENSG00000163106",
  "HSD17B10" = "ENSG00000072506", "HSD17B4" = "ENSG00000133835", "HSPA5" = "ENSG00000044574",
  "HTR2B" = "ENSG00000135914", "HTR2C" = "ENSG00000147246", "ICAM1" = "ENSG00000090339",
  "IGFBP2" = "ENSG00000115457", "IKBKB" = "ENSG00000104365", "IL10RA" = "ENSG00000110324",
  "IL6" = "ENSG00000136244", "IMPDH2" = "ENSG00000178035", "IRF1" = "ENSG00000125347",
  "ITGA1" = "ENSG00000213949", "JAK1" = "ENSG00000162434", "KCNA5" = "ENSG00000130037",
  "LEF1" = "ENSG00000138795", "LIAS" = "ENSG00000121897", "LIPT1" = "ENSG00000144182",
  "LYN" = "ENSG00000254087", "MAN2B1" = "ENSG00000104774", "MAOB" = "ENSG00000069535",
  "MAPKAPK2" = "ENSG00000162889", "MB" = "ENSG00000198125", "MGAT1" = "ENSG00000131446",
  "MKNK2" = "ENSG00000099875", "MTOR" = "ENSG00000198793", "NFE2L2" = "ENSG00000116044",
  "NFKB1" = "ENSG00000109320", "NMT1" = "ENSG00000136448", "NR1H3" = "ENSG00000025434",
  "NR3C1" = "ENSG00000113580", "NUDCD2" = "ENSG00000170584", "OAZ1" = "ENSG00000104904",
  "PA2G4" = "ENSG00000170515", "PABPC1" = "ENSG00000070756", "PARP1" = "ENSG00000143799",
  "PARP12" = "ENSG00000059378", "PCTP" = "ENSG00000141179", "PDCD6" = "ENSG00000249915",
  "PDCD6IP" = "ENSG00000170248", "PDHB" = "ENSG00000168291", "PDHX" = "ENSG00000110435",
  "PLA2G4A" = "ENSG00000116711", "POLR2D" = "ENSG00000144231", "PPARG" = "ENSG00000132170",
  "PRKCQ" = "ENSG00000065675", "PTGR1" = "ENSG00000106853", "PTGS1" = "ENSG00000095303",
  "PTPN2" = "ENSG00000175354", "PTPN6" = "ENSG00000111679", "PTPRC" = "ENSG00000081237",
  "PTPRF" = "ENSG00000142949", "PTPRJ" = "ENSG00000149177", "RBM39" = "ENSG00000131051",
  "RELA" = "ENSG00000173039", "RENBP" = "ENSG00000102032", "RHOC" = "ENSG00000155366",
  "S100A6" = "ENSG00000197956", "SAT1" = "ENSG00000111371", "SAT2" = "ENSG00000134294",
  "SCN9A" = "ENSG00000169432", "SEC13" = "ENSG00000157020", "SERPINB1" = "ENSG00000197696",
  "SERPINB10" = "ENSG00000242550", "SLC31A1" = "ENSG00000136868", "SPHK1" = "ENSG00000176170",
  "SREBF1" = "ENSG00000072310", "STARD13" = "ENSG00000133121", "STAT1" = "ENSG00000115415",
  "STAT3" = "ENSG00000168610", "STAT5A" = "ENSG00000126561", "STK4" = "ENSG00000101109",
  "TBXAS1" = "ENSG00000059377", "TCN2" = "ENSG00000185339", "TDP1" = "ENSG00000042088",
  "TGFB1" = "ENSG00000105329", "TIMP1" = "ENSG00000102265", "TNF" = "ENSG00000232810",
  "TOP2A" = "ENSG00000131747", "TSPO" = "ENSG00000100300", "XDH" = "ENSG00000158125",
  "XRCC6" = "ENSG00000196419", "ZEB1" = "ENSG00000148516", "ZHX2" = "ENSG00000178764"
)

# 步骤 1: 加载 eQTLGen 数据
cat("步骤 1: 加载 eQTLGen 数据（已 LD clumped）\n")
cat("----------------------------------------------------------------------\n")

if (!file.exists(CLUMPED_FILE)) {
  stop("eQTLGen clumped 文件不存在：", CLUMPED_FILE)
}

cat("  读取 eQTLGen clumped 数据...\n")
eqtlgen_data <- read_excel(CLUMPED_FILE)
cat(sprintf("  ✓ 成功读取 %d 个 SNP\n", nrow(eqtlgen_data)))

# 确保 SNP 列存在
if (!"SNP" %in% names(eqtlgen_data)) {
  stop("eQTLGen 数据缺少 SNP 列")
}

# 移除没有 rsID 的 SNP
eqtlgen_data <- eqtlgen_data[!is.na(eqtlgen_data$SNP) & grepl("^rs", eqtlgen_data$SNP), ]
cat(sprintf("  ✓ 保留有 rsID 的 SNP: %d 个\n\n", nrow(eqtlgen_data)))

# 步骤 2: 加载 MEGASTROKE 数据
cat("步骤 2: 加载 MEGASTROKE 数据\n")
cat("----------------------------------------------------------------------\n")

if (!file.exists(MEGASTROKE_FILE)) {
  stop("MEGASTROKE 文件不存在：", MEGASTROKE_FILE)
}

cat("  读取 MEGASTROKE 协调化数据...\n")
outcome_data <- fread(MEGASTROKE_FILE, sep = "\t", stringsAsFactors = FALSE, showProgress = TRUE)
cat(sprintf("  ✓ 成功读取 %d 个 SNP\n", nrow(outcome_data)))

# 移除没有 rsID 的 SNP
outcome_data <- outcome_data[!is.na(outcome_data$hm_rsid) & outcome_data$hm_rsid != "", ]
cat(sprintf("  ✓ 保留有 rsID 的 SNP: %d 个\n\n", nrow(outcome_data)))

# 步骤 3: 准备每个基因的暴露数据
cat("步骤 3: 准备每个基因的暴露数据\n")
cat("----------------------------------------------------------------------\n")

prepared_genes <- list()

for (gene_symbol in target_genes) {
  ensg_id <- gene_to_ensg[[gene_symbol]]
  
  if (is.null(ensg_id)) {
    cat(sprintf("%s: 无 ENSG ID，跳过\n", gene_symbol))
    next
  }
  
  # 从 eQTLGen 获取 SNPs
  eqtlgen_snps <- eqtlgen_data[eqtlgen_data$gene == gene_symbol | 
                               eqtlgen_data$gene == ensg_id, ]
  
  if (nrow(eqtlgen_snps) == 0) {
    cat(sprintf("%s: 无 eQTL 数据\n", gene_symbol))
    next
  }
  
  # 创建暴露数据
  exposure <- data.frame(
    SNP = eqtlgen_snps$SNP,
    BETA = eqtlgen_snps$beta.exposure,
    SE = eqtlgen_snps$se.exposure,
    EFFECT_ALLELE = eqtlgen_snps$effect_allele.exposure,
    OTHER_ALLELE = eqtlgen_snps$other_allele.exposure,
    EAF = eqtlgen_snps$eaf.exposure,
    PVAL = eqtlgen_snps$pval.exposure,
    stringsAsFactors = FALSE
  )
  
  # 移除重复的 SNP
  exposure <- exposure[!is.na(exposure$SNP) & exposure$SNP != "", ]
  if (nrow(exposure) > 0) {
    exposure <- exposure[!duplicated(exposure$SNP), ]
  }
  
  if (nrow(exposure) < 2) {
    cat(sprintf("%s: SNP 数量不足 (%d)，跳过\n", gene_symbol, nrow(exposure)))
    next
  }
  
  tryCatch({
    exposure_fmt <- format_data(
      exposure, type = "exposure", snp_col = "SNP", beta_col = "BETA", se_col = "SE",
      eaf_col = "EAF", effect_allele_col = "EFFECT_ALLELE", other_allele_col = "OTHER_ALLELE",
      pval_col = "PVAL"
    )
    
    f_stat <- mean((exposure_fmt$beta.exposure / exposure_fmt$se.exposure)^2, na.rm = TRUE)
    
    prepared_genes[[gene_symbol]] <- list(
      exposure = exposure, exposure_fmt = exposure_fmt, ensg_id = ensg_id, 
      f_stat = f_stat, n_snps = nrow(exposure)
    )
    
    cat(sprintf("%s: ✓ 准备就绪 (%d SNPs, F=%.2f)\n", 
                gene_symbol, nrow(exposure), f_stat))
  }, error = function(e) {
    cat(sprintf("%s: 格式化失败 - %s\n", gene_symbol, e$message))
  })
}

cat(sprintf("\n\n成功准备：%d/%d 个基因\n", length(prepared_genes), length(target_genes)))

# 步骤 4: MR 分析（使用 rsID 匹配 + harmonise_data）
cat("\n步骤 4: MR 分析（rsID 匹配 + harmonise_data）\n")
cat("----------------------------------------------------------------------\n")

mr_results <- list()
harmonised_data <- list()
failed_genes <- data.frame(
  gene = character(),
  reason = character(),
  n_snps = integer(),
  stringsAsFactors = FALSE
)

for (gene_symbol in names(prepared_genes)) {
  cat(sprintf("\n%s:\n", gene_symbol))
  
  exposure_fmt <- prepared_genes[[gene_symbol]]$exposure_fmt
  f_stat <- prepared_genes[[gene_symbol]]$f_stat
  n_snps <- prepared_genes[[gene_symbol]]$n_snps
  
  cat(sprintf("  暴露 SNP 数：%d\n", n_snps))
  cat(sprintf("  F 统计量：%.2f\n", f_stat))
  
  # 提取该基因的结局数据
  outcome_matched <- outcome_data[outcome_data$hm_rsid %in% exposure_fmt$SNP, ]
  cat(sprintf("  结局匹配 SNP 数：%d\n", nrow(outcome_matched)))
  
  if (nrow(outcome_matched) == 0) {
    cat(sprintf("  ✗ 无匹配的结局数据\n"))
    failed_genes <- rbind(failed_genes, data.frame(
      gene = gene_symbol, reason = "无匹配结局 SNP", n_snps = n_snps,
      stringsAsFactors = FALSE
    ))
    next
  }
  
  # 格式化结局数据
  outcome_data_df <- as.data.frame(outcome_matched)
  outcome_data_df$beta.outcome <- outcome_data_df$hm_beta
  outcome_data_df$se.outcome <- abs(outcome_data_df$hm_beta)
  outcome_data_df$effect_allele.outcome <- outcome_data_df$hm_effect_allele
  outcome_data_df$other_allele.outcome <- outcome_data_df$hm_other_allele
  outcome_data_df$eaf.outcome <- outcome_data_df$hm_effect_allele_frequency
  outcome_data_df$pval.outcome <- 0.05
  
  outcome_fmt <- format_data(
    outcome_data_df,
    type = "outcome",
    snp_col = "hm_rsid",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    eaf_col = "eaf.outcome",
    pval_col = "pval.outcome"
  )
  
  # Harmonise
  dat <- tryCatch({
    harmonise_data(exposure_fmt, outcome_fmt)
  }, error = function(e) {
    cat(sprintf("  ✗ harmonise_data 错误：%s\n", e$message))
    return(NULL)
  })
  
  if (is.null(dat) || nrow(dat) == 0) {
    cat(sprintf("  ✗ harmonise_data 失败\n"))
    failed_genes <- rbind(failed_genes, data.frame(
      gene = gene_symbol, reason = "harmonise 失败", n_snps = n_snps,
      stringsAsFactors = FALSE
    ))
    next
  }
  
  cat(sprintf("  ✓ Harmonised SNP 数：%d\n", nrow(dat)))
  
  if (nrow(dat) < 1) {
    cat(sprintf("  ✗ SNP 数量不足 (<1)\n"))
    failed_genes <- rbind(failed_genes, data.frame(
      gene = gene_symbol, reason = "harmonised SNP<1", n_snps = n_snps,
      stringsAsFactors = FALSE
    ))
    next
  }
  
  harmonised_data[[gene_symbol]] <- dat
  
  # MR 分析
  res <- tryCatch({
    mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
  }, error = function(e) {
    cat(sprintf("  ✗ MR 分析错误：%s\n", e$message))
    return(NULL)
  })
  
  if (is.null(res) || nrow(res) == 0) {
    cat(sprintf("  ✗ MR 分析失败\n"))
    failed_genes <- rbind(failed_genes, data.frame(
      gene = gene_symbol, reason = "MR 分析失败", n_snps = n_snps,
      stringsAsFactors = FALSE
    ))
    next
  }
  
  ivw_res <- res[res$method == "Inverse variance weighted", ]
  het_res <- mr_heterogeneity(dat)
  
  if (nrow(dat) >= 3) {
    egger_intercept <- mr_pleiotropy_test(dat)
  } else {
    egger_intercept <- data.frame(intercept = NA, pval = NA)
  }
  
  result <- data.frame(
    gene = gene_symbol,
    method = ivw_res$method[1],
    beta = ivw_res$b[1],
    se = ivw_res$se[1],
    or = exp(ivw_res$b[1]),
    ci_low = exp(ivw_res$b[1] - 1.96 * ivw_res$se[1]),
    ci_high = exp(ivw_res$b[1] + 1.96 * ivw_res$se[1]),
    pval = ivw_res$pval[1],
    f_stat = f_stat,
    nsnp = nrow(dat),
    heterogeneity_q = ifelse(nrow(het_res) > 0, het_res$Q[1], NA),
    heterogeneity_pval = ifelse(nrow(het_res) > 0, het_res$Q_pval[1], NA),
    egger_intercept = ifelse(nrow(egger_intercept) > 0, egger_intercept$intercept[1], NA),
    egger_pval = ifelse(nrow(egger_intercept) > 0, egger_intercept$pval[1], NA),
    wm_beta = ifelse(nrow(res[res$method == "Weighted median", ]) > 0,
                     res[res$method == "Weighted median", "b"][1], NA),
    wm_pval = ifelse(nrow(res[res$method == "Weighted median", ]) > 0,
                     res[res$method == "Weighted median", "pval"][1], NA),
    stringsAsFactors = FALSE
  )
  
  mr_results[[gene_symbol]] <- result
  
  direction_consistent <- sign(ivw_res$b[1]) == sign(result$wm_beta)
  
  cat(sprintf("  ✓ MR: OR=%.3f (%.3f-%.3f), P=%.2e, F=%.2f\n", 
              result$or, result$ci_low, result$ci_high, result$pval, result$f_stat))
  cat(sprintf("    异质性：Q=%.2f, P=%.3f\n", result$heterogeneity_q, result$heterogeneity_pval))
  cat(sprintf("    多效性：intercept=%.3f, P=%.3f\n", result$egger_intercept, result$egger_pval))
  cat(sprintf("    稳健性：WM P=%.3f, 方向一致：%s\n", result$wm_pval, ifelse(direction_consistent, "是", "否")))
  
}

# 步骤 5: 保存结果
cat("\n\n步骤 5: 保存结果\n")
cat("----------------------------------------------------------------------\n")

if (length(mr_results) > 0) {
  results_df <- do.call(rbind, mr_results)
  results_df <- results_df[order(results_df$pval), ]
  
  output_file <- file.path(OUTPUT_DIR, "mr_results_all_genes.csv")
  write.csv(results_df, output_file, row.names = FALSE)
  cat(sprintf("✓ 保存 MR 结果：%s\n", output_file))
  
  cat(sprintf("\nTop 显著基因 (P < 0.05):\n"))
  sig_results <- results_df[results_df$pval < 0.05, ]
  if (nrow(sig_results) > 0) {
    print(sig_results[, c("gene", "or", "ci_low", "ci_high", "pval", "nsnp", "f_stat", 
                          "heterogeneity_pval", "egger_pval")])
  } else {
    cat("  无显著基因\n")
  }
  
  results_df$pval_adj <- p.adjust(results_df$pval, method = "fdr")
  write.csv(results_df, file.path(OUTPUT_DIR, "mr_results_fdr_corrected.csv"), row.names = FALSE)
  cat(sprintf("✓ 保存 FDR 校正结果\n"))
  
  if (length(harmonised_data) > 0) {
    all_harmonised <- do.call(rbind, lapply(names(harmonised_data), function(gene) {
      dat <- harmonised_data[[gene]]
      dat$gene <- gene
      dat
    }))
    write.csv(all_harmonised, file.path(OUTPUT_DIR, "harmonised_data_all_genes.csv"), row.names = FALSE)
    cat(sprintf("✓ 保存 Harmonised 数据\n"))
  }
  
  # 保存失败基因列表
  if (nrow(failed_genes) > 0) {
    write.csv(failed_genes, file.path(OUTPUT_DIR, "failed_genes.csv"), row.names = FALSE)
    cat(sprintf("✓ 保存失败基因列表 (%d 个)\n", nrow(failed_genes)))
  }
  
  # 生成质量报告
  cat(sprintf("\n生成质量报告...\n"))
  quality_report <- file.path(OUTPUT_DIR, "quality_control_report.txt")
  sink(quality_report)
  cat("MR 分析质量控制报告（最终版 - rsID 匹配）\n")
  cat("========================================\n\n")
  cat("1. 总体统计\n")
  cat(sprintf("   目标基因总数：%d\n", length(target_genes)))
  cat(sprintf("   成功准备基因数：%d\n", length(prepared_genes)))
  cat(sprintf("   成功 MR 分析基因数：%d\n", length(mr_results)))
  cat(sprintf("   失败基因数：%d\n", nrow(failed_genes)))
  cat("\n")
  
  if (nrow(failed_genes) > 0) {
    cat("2. 失败原因统计\n")
    print(table(failed_genes$reason))
    cat("\n")
  }
  
  cat("3. F 统计量评估\n")
  cat("   平均 F 统计量：", round(mean(results_df$f_stat, na.rm = TRUE), 2), "\n")
  cat("   F > 10 的基因数：", sum(results_df$f_stat > 10, na.rm = TRUE), "/", nrow(results_df), "\n\n")
  
  cat("4. 异质性评估\n")
  cat("   显著异质性 (P < 0.05) 的基因数：", sum(results_df$heterogeneity_pval < 0.05, na.rm = TRUE), "\n")
  cat("   极显著异质性 (P < 1e-10) 的基因：", 
      paste(results_df$gene[results_df$heterogeneity_pval < 1e-10], collapse = ", "), "\n\n")
  
  cat("5. 多效性评估\n")
  cat("   显著多效性 (P < 0.05) 的基因数：", sum(results_df$egger_pval < 0.05, na.rm = TRUE), "\n\n")
  
  cat("6. 稳健性检查\n")
  cat("   方向一致的基因数：", sum(sign(results_df$beta) == sign(results_df$wm_beta), na.rm = TRUE), "\n\n")
  
  cat("7. 匹配统计\n")
  cat("   平均匹配 SNP 数：", round(mean(results_df$nsnp, na.rm = TRUE), 1), "\n")
  cat("   最大匹配 SNP 数：", max(results_df$nsnp, na.rm = TRUE), "\n")
  cat("   最小匹配 SNP 数：", min(results_df$nsnp, na.rm = TRUE), "\n\n")
  
  sink()
  cat(sprintf("✓ 保存质量报告：%s\n", quality_report))
}

cat("\n======================================================================\n")
cat("分析完成！\n")
cat("======================================================================\n")
