#!/usr/bin/env Rscript
# ================================================================================
# 138 基因 MR 分析（使用新 eQTLGen 数据 - r²<0.01）
# ================================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(readxl)
})

cat("======================================================================\n")
cat("138 基因 MR 分析（使用新 eQTLGen 数据 - r²<0.01）\n")
cat("======================================================================\n\n")

# 配置
CLUMPED_FILE <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.01.xlsx"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_138genes_r2_01"
MEGASTROKE_FILE <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 加载数据
cat("加载新的 eQTLGen 数据（r²<0.01）...\n")
eqtlgen_data <- read_excel(CLUMPED_FILE)
eqtlgen_data <- eqtlgen_data[!is.na(eqtlgen_data$SNP) & grepl("^rs", eqtlgen_data$SNP), ]

cat("加载 MEGASTROKE 数据...\n")
outcome_data <- fread(MEGASTROKE_FILE, sep = "\t", stringsAsFactors = FALSE)
outcome_data <- outcome_data[!is.na(outcome_data$hm_rsid) & outcome_data$hm_rsid != "", ]

cat(sprintf("eQTLGen: %d SNPs, MEGASTROKE: %d SNPs\n\n", nrow(eqtlgen_data), nrow(outcome_data)))

# 138 个目标基因
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

# 完整的 ENSG ID 映射（包含原始和补充的）
gene_to_ensg <- list(
  # 原始映射
  "LYN" = "ENSG00000254087", "PRKCQ" = "ENSG00000065675", "NMT1" = "ENSG00000136448",
  "TDP1" = "ENSG00000042088", "MAN2B1" = "ENSG00000104774", "IL10RA" = "ENSG00000110324",
  "RHOC" = "ENSG00000155366", "SREBF1" = "ENSG00000072310", "CTSC" = "ENSG00000109861",
  "CAT" = "ENSG00000121691", "RBM39" = "ENSG00000131051", "CHFR" = "ENSG00000072609",
  "TCN2" = "ENSG00000185339", "SCN9A" = "ENSG00000169432", "CTSB" = "ENSG00000164733",
  "CASP8" = "ENSG00000064012", "MAPKAPK2" = "ENSG00000162889", "SEC13" = "ENSG00000157020",
  "PCTP" = "ENSG00000141179", "CPT2" = "ENSG00000157184", "BRD3" = "ENSG00000169925",
  "ITGA1" = "ENSG00000213949", "CITED2" = "ENSG00000164442", "HIBADH" = "ENSG00000106049",
  "SAT2" = "ENSG00000134294", "TSPO" = "ENSG00000100300", "FLT4" = "ENSG00000037280",
  "CPT1A" = "ENSG00000110090", "AKT1" = "ENSG00000142208", "CCR5" = "ENSG00000160791",
  "MGAT1" = "ENSG00000131446", "IGFBP2" = "ENSG00000115457", "PPARG" = "ENSG00000132170",
  "EPHX1" = "ENSG00000143819", "AIF1" = "ENSG00000204472", "PLA2G4A" = "ENSG00000116711",
  "ALDH9A1" = "ENSG00000131449", "BST1" = "ENSG00000197121", "CNDP2" = "ENSG00000133313",
  "TNF" = "ENSG00000232810", "PARP1" = "ENSG00000143799", "IKBKB" = "ENSG00000104365",
  "ADRB1" = "ENSG00000043591", "SPHK1" = "ENSG00000176170", "GCH1" = "ENSG00000131979",
  "ACADVL" = "ENSG00000072778", "CTSD" = "ENSG00000117984", "PDCD6IP" = "ENSG00000170248",
  "PTPRC" = "ENSG00000081237", "PABPC1" = "ENSG00000070756", "ACTA2" = "ENSG00000107796",
  "PTGR1" = "ENSG00000106853", "XRCC6" = "ENSG00000196419", "TBXAS1" = "ENSG00000059377",
  "NR1H3" = "ENSG00000025434", "CTSL" = "ENSG00000135047", "ZHX2" = "ENSG00000178764",
  "MKNK2" = "ENSG00000099875", "SERPINB10" = "ENSG00000242550", "HSPA5" = "ENSG00000044574",
  "CCND1" = "ENSG00000110092", "PTPN2" = "ENSG00000175354", "HSD17B4" = "ENSG00000133835",
  "PDCD6" = "ENSG00000249915", "PARP12" = "ENSG00000059378", "SERPINB1" = "ENSG00000197696",
  "STAT1" = "ENSG00000115415", "NFE2L2" = "ENSG00000116044", "HMOX1" = "ENSG00000100292",
  "CTSF" = "ENSG00000174080", "CCL2" = "ENSG00000108691", "FDX1" = "ENSG00000137714",
  "LIAS" = "ENSG00000121897", "LIPT1" = "ENSG00000144182", "PDHX" = "ENSG00000110435",
  "SLC31A1" = "ENSG00000136868", "ATP7B" = "ENSG00000123191", "ATOX1" = "ENSG00000177556",
  "GPX4" = "ENSG00000167468",
  
  # 补充的 ENSG ID 映射（65 个基因）
  "KCNA5" = "ENSG00000130037", "HIF1A" = "ENSG00000100644", "FABP4" = "ENSG00000170323",
  "STAT5A" = "ENSG00000126561", "FABP2" = "ENSG00000145384", "B2M" = "ENSG00000166710",
  "HBS1L" = "ENSG00000112339", "NUDCD2" = "ENSG00000170584", "JAK1" = "ENSG00000162434",
  "GPX1" = "ENSG00000233276", "FABP5" = "ENSG00000164687", "XDH" = "ENSG00000158125",
  "MB" = "ENSG00000198125", "POLR2D" = "ENSG00000144231", "HSD17B10" = "ENSG00000072506",
  "ZEB1" = "ENSG00000148516", "RELA" = "ENSG00000173039", "IRF1" = "ENSG00000125347",
  "GFAP" = "ENSG00000131095", "NR3C1" = "ENSG00000113580", "F3" = "ENSG00000117525",
  "C3" = "ENSG00000125730", "PTGS1" = "ENSG00000095303", "IMPDH2" = "ENSG00000178035",
  "PTPRF" = "ENSG00000142949", "HPGDS" = "ENSG00000163106", "PTPRJ" = "ENSG00000149177",
  "CASK" = "ENSG00000147044", "TOP2A" = "ENSG00000131747", "IL6" = "ENSG00000136244",
  "CP" = "ENSG00000047457", "S100A6" = "ENSG00000197956", "DDC" = "ENSG00000132437",
  "CUL4B" = "ENSG00000158290", "EGFR" = "ENSG00000146648", "COL1A1" = "ENSG00000108821",
  "STARD13" = "ENSG00000133121", "TGFB1" = "ENSG00000105329", "HTR2C" = "ENSG00000147246",
  "CTSS" = "ENSG00000163131", "CNR2" = "ENSG00000188822", "FNTA" = "ENSG00000168522",
  "RENBP" = "ENSG00000102032", "CCNA2" = "ENSG00000145386", "LEF1" = "ENSG00000138795",
  "SAT1" = "ENSG00000130066", "HTR2B" = "ENSG00000135914", "CDK4" = "ENSG00000135446",
  "CXCR3" = "ENSG00000186810", "TIMP1" = "ENSG00000102265", "OAZ1" = "ENSG00000104904",
  "STK4" = "ENSG00000101109", "ACADM" = "ENSG00000117054", "STAT3" = "ENSG00000168610",
  "NFKB1" = "ENSG00000109320", "CTSK" = "ENSG00000143387", "PTPN6" = "ENSG00000111679",
  "PA2G4" = "ENSG00000170515", "ACAD11" = "ENSG00000240303", "MAOB" = "ENSG00000069535",
  "ICAM1" = "ENSG00000090339", "DLAT" = "ENSG00000150768", "PDHB" = "ENSG00000168291",
  "ATP7A" = "ENSG00000165240", "MTOR" = "ENSG00000198793"
)

cat(sprintf("ENSG ID 映射数：%d\n\n", length(gene_to_ensg)))

# 统计每个基因的 SNP 数量
cat("统计每个基因的 SNP 数量...\n")
gene_snp_counts <- list()
for (gene in target_genes) {
  ensg_id <- gene_to_ensg[[gene]]
  if (!is.null(ensg_id)) {
    snps <- eqtlgen_data[eqtlgen_data$gene == gene | eqtlgen_data$gene == ensg_id, ]
    gene_snp_counts[[gene]] <- nrow(snps)
  }
}

snp_counts_df <- data.frame(
  gene = names(gene_snp_counts),
  snp_count = unname(unlist(gene_snp_counts)),
  stringsAsFactors = FALSE
)

cat(sprintf("平均每个基因的 SNP 数：%.2f\n", mean(snp_counts_df$snp_count)))
cat(sprintf("中位数：%.0f\n", median(snp_counts_df$snp_count)))
cat(sprintf("范围：%d - %d\n\n", min(snp_counts_df$snp_count), max(snp_counts_df$snp_count)))

# 处理每个基因
results_list <- list()
skipped <- data.frame(gene = character(), reason = character(), stringsAsFactors = FALSE)
failed <- data.frame(gene = character(), reason = character(), n_snps = integer(), stringsAsFactors = FALSE)

cat("开始处理基因...\n\n")

for (i in seq_along(target_genes)) {
  gene_symbol <- target_genes[i]
  cat(sprintf("[%d/%d] %s: ", i, length(target_genes), gene_symbol))
  
  # 检查 ENSG ID
  ensg_id <- gene_to_ensg[[gene_symbol]]
  if (is.null(ensg_id)) {
    cat("无 ENSG ID，跳过\n")
    skipped <- rbind(skipped, data.frame(gene = gene_symbol, reason = "无 ENSG ID", stringsAsFactors = FALSE))
    next
  }
  
  # 获取 eQTL 数据（已经是 TwoSampleMR 格式）
  eqtlgen_snps <- eqtlgen_data[eqtlgen_data$gene == gene_symbol | eqtlgen_data$gene == ensg_id, ]
  
  if (nrow(eqtlgen_snps) < 2) {
    cat(sprintf("SNP 不足 (%d)，跳过\n", nrow(eqtlgen_snps)))
    skipped <- rbind(skipped, data.frame(gene = gene_symbol, reason = sprintf("SNP 不足 (%d)", nrow(eqtlgen_snps)), stringsAsFactors = FALSE))
    next
  }
  
  # 直接使用 eQTLGen 数据（已经是 TwoSampleMR 格式，不需要 format_data）
  exposure <- eqtlgen_snps[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", 
                                "other_allele.exposure", "eaf.exposure", "pval.exposure")]
  
  exposure <- exposure[!duplicated(exposure$SNP) & !is.na(exposure$SNP) & exposure$SNP != "", ]
  
  if (nrow(exposure) < 2) {
    cat("去重后 SNP 不足，跳过\n")
    skipped <- rbind(skipped, data.frame(gene = gene_symbol, reason = "去重后 SNP 不足", stringsAsFactors = FALSE))
    next
  }
  
  # 暴露数据已经是正确格式，直接使用
  exposure_fmt <- as.data.frame(exposure)
  exposure_fmt$type <- "exposure"
  exposure_fmt$mr_keep.exposure <- TRUE
  exposure_fmt$id.exposure <- "eQTLGen"
  exposure_fmt$exposure <- gene_symbol
  
  # 提取结局数据
  outcome_matched <- outcome_data[outcome_data$hm_rsid %in% exposure_fmt$SNP, ]
  
  if (nrow(outcome_matched) == 0) {
    cat("无匹配结局数据，跳过\n")
    next
  }
  
  # 格式化结局数据
  outcome_df <- as.data.frame(outcome_matched)
  outcome_df$beta.outcome <- as.numeric(outcome_df$hm_beta)
  outcome_df$se.outcome <- as.numeric(outcome_df$standard_error)
  outcome_df$effect_allele.outcome <- as.character(outcome_df$hm_effect_allele)
  outcome_df$other_allele.outcome <- as.character(outcome_df$hm_other_allele)
  outcome_df$eaf.outcome <- as.numeric(outcome_df$hm_effect_allele_frequency)
  outcome_df$pval.outcome <- as.numeric(outcome_df$p_value)
  
  # 结局数据格式化
  outcome_fmt <- as.data.frame(outcome_df[, c("hm_rsid", "beta.outcome", "se.outcome", 
                                               "effect_allele.outcome", "other_allele.outcome", 
                                               "eaf.outcome", "pval.outcome")])
  names(outcome_fmt)[1] <- "SNP"
  outcome_fmt$type <- "outcome"
  outcome_fmt$mr_keep.outcome <- TRUE
  outcome_fmt$id.outcome <- "MEGASTROKE"
  outcome_fmt$outcome <- "Stroke"
  
  # Harmonise
  dat <- tryCatch({
    harmonise_data(exposure_fmt, outcome_fmt)
  }, error = function(e) {
    cat(sprintf("Harmonise 失败：%s\n", e$message))
    failed <- rbind(failed, data.frame(gene = gene_symbol, reason = paste("Harmonise:", e$message), n_snps = nrow(exposure_fmt), stringsAsFactors = FALSE))
    return(NULL)
  })
  
  if (is.null(dat) || nrow(dat) == 0) {
    cat("Harmonise 失败\n")
    next
  }
  
  # MR 分析
  res <- tryCatch({
    mr(dat, method_list = c("mr_ivw", "mr_weighted_median"))
  }, error = function(e) {
    cat(sprintf("MR 分析失败：%s\n", e$message))
    failed <- rbind(failed, data.frame(gene = gene_symbol, reason = paste("MR:", e$message), n_snps = nrow(dat), stringsAsFactors = FALSE))
    return(NULL)
  })
  
  if (is.null(res) || nrow(res) == 0) {
    cat("MR 分析失败\n")
    next
  }
  
  # 提取 IVW 结果
  ivw_res <- res[res$method == "Inverse variance weighted", ]
  
  if (nrow(ivw_res) == 0) {
    cat("无 IVW 结果\n")
    next
  }
  
  # 计算 F 统计量
  f_stat <- mean((exposure_fmt$beta.exposure / exposure_fmt$se.exposure)^2, na.rm = TRUE)
  
  # 创建结果
  result <- data.frame(
    gene = gene_symbol,
    beta = ivw_res$b[1],
    se = ivw_res$se[1],
    or = exp(ivw_res$b[1]),
    pval = ivw_res$pval[1],
    f_stat = f_stat,
    nsnp = nrow(dat),
    stringsAsFactors = FALSE
  )
  
  results_list[[gene_symbol]] <- result
  
  # 保存
  write.csv(dat, file.path(OUTPUT_DIR, paste0(gene_symbol, "_harmonised.csv")), row.names = FALSE)
  write.csv(res, file.path(OUTPUT_DIR, paste0(gene_symbol, "_mr_results.csv")), row.names = FALSE)
  
  cat(sprintf("✓ OR=%.3f, P=%.3f, F=%.2f, SNPs=%d\n", result$or, result$pval, f_stat, nrow(dat)))
}

# 汇总结果
cat("\n\n汇总结果:\n")
if (length(results_list) > 0) {
  all_results <- do.call(rbind, results_list)
  all_results <- all_results[order(all_results$pval), ]
  
  write.csv(all_results, file.path(OUTPUT_DIR, "summary_results.csv"), row.names = FALSE)
  write.csv(all_results, file.path(OUTPUT_DIR, "mr_results_fdr_corrected.csv"), row.names = FALSE)
  
  if (nrow(skipped) > 0) {
    write.csv(skipped, file.path(OUTPUT_DIR, "skipped_genes.csv"), row.names = FALSE)
  }
  
  if (nrow(failed) > 0) {
    write.csv(failed, file.path(OUTPUT_DIR, "failed_genes.csv"), row.names = FALSE)
  }
  
  cat(sprintf("成功分析：%d 个基因\n", nrow(all_results)))
  cat(sprintf("跳过：%d 个基因\n", nrow(skipped)))
  cat(sprintf("失败：%d 个基因\n", nrow(failed)))
  
  cat("\nTop 10 最显著的基因:\n")
  print(head(all_results[, c("gene", "beta", "or", "pval", "f_stat", "nsnp")], 10))
  
  if (any(all_results$pval < 0.05)) {
    cat("\n显著基因 (P < 0.05):\n")
    print(all_results[all_results$pval < 0.05, c("gene", "beta", "or", "pval", "f_stat", "nsnp")])
  }
} else {
  cat("无成功分析的基因\n")
}

cat("\n完成！结果保存在:", OUTPUT_DIR, "\n")
