#!/usr/bin/env Rscript
# ================================================================================
# 为 138 个目标基因准备 MR 分析数据并运行
# ================================================================================

library(dplyr)
library(data.table)

cat("======================================================================\n")
cat("为 138 个目标基因准备 MR 分析数据\n")
cat("======================================================================\n\n")

# 目标基因列表（138 个，去重后）
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
  "HSD17B4", "ACAD11", "PDCD6", "PARP12", "SERPINB1A", "STAT1", "NFE2L2",
  "HMOX1", "CTSF", "CCL2", "MAOB", "ICAM1", "FDX1", "LIAS", "LIPT1",
  "DLAT", "PDHB", "PDHX", "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
  "HIF1A", "MTOR", "NFKB1", "GPX4"
)

target_genes <- unique(target_genes)
cat(sprintf("目标基因数：%d\n\n", length(target_genes)))

# 基因符号到 ENSG ID 的映射（通过 Ensembl API 查询，143/144 个基因）
# 仅 SERPINB1A 未找到
gene_to_ensg <- list(
  "ACAD11" = "ENSG00000240303",
  "ACADM" = "ENSG00000117054",
  "ACADVL" = "ENSG00000072778",
  "ACTA2" = "ENSG00000107796",
  "ADRB1" = "ENSG00000043591",
  "AIF1" = "ENSG00000204472",
  "AKT1" = "ENSG00000142208",
  "ALDH9A1" = "ENSG00000143149",
  "ATOX1" = "ENSG00000177556",
  "ATP7A" = "ENSG00000165240",
  "ATP7B" = "ENSG00000123191",
  "B2M" = "ENSG00000166710",
  "BRD3" = "ENSG00000169925",
  "BST1" = "ENSG00000197121",
  "C3" = "ENSG00000125730",
  "CASK" = "ENSG00000147044",
  "CASP8" = "ENSG00000064012",
  "CAT" = "ENSG00000121691",
  "CCL2" = "ENSG00000108691",
  "CCNA2" = "ENSG00000145386",
  "CCND1" = "ENSG00000110092",
  "CCR5" = "ENSG00000160791",
  "CDK4" = "ENSG00000135446",
  "CHFR" = "ENSG00000072609",
  "CITED2" = "ENSG00000164442",
  "CNDP2" = "ENSG00000133313",
  "CNR2" = "ENSG00000188822",
  "COL1A1" = "ENSG00000108821",
  "CP" = "ENSG00000047457",
  "CPT1A" = "ENSG00000110090",
  "CPT2" = "ENSG00000157184",
  "CTSB" = "ENSG00000164733",
  "CTSC" = "ENSG00000109861",
  "CTSD" = "ENSG00000117984",
  "CTSF" = "ENSG00000174080",
  "CTSK" = "ENSG00000143387",
  "CTSL" = "ENSG00000135047",
  "CTSS" = "ENSG00000163131",
  "CUL4B" = "ENSG00000158290",
  "CXCR3" = "ENSG00000186810",
  "DDC" = "ENSG00000132437",
  "DLAT" = "ENSG00000150768",
  "EGFR" = "ENSG00000146648",
  "EPHX1" = "ENSG00000143819",
  "F3" = "ENSG00000117525",
  "FABP2" = "ENSG00000145384",
  "FABP4" = "ENSG00000170323",
  "FABP5" = "ENSG00000236972",
  "FDX1" = "ENSG00000137714",
  "FLT4" = "ENSG00000037280",
  "FNTA" = "ENSG00000168522",
  "GCH1" = "ENSG00000131979",
  "GFAP" = "ENSG00000131095",
  "GPX1" = "ENSG00000233276",
  "GPX4" = "ENSG00000167468",
  "HBS1L" = "ENSG00000112339",
  "HIBADH" = "ENSG00000106049",
  "HIF1A" = "ENSG00000100644",
  "HMOX1" = "ENSG00000100292",
  "HPGDS" = "ENSG00000163106",
  "HSD17B10" = "ENSG00000072506",
  "HSD17B4" = "ENSG00000133835",
  "HSPA5" = "ENSG00000044574",
  "HTR2B" = "ENSG00000135914",
  "HTR2C" = "ENSG00000147246",
  "ICAM1" = "ENSG00000090339",
  "IGFBP2" = "ENSG00000115457",
  "IKBKB" = "ENSG00000104365",
  "IL10RA" = "ENSG00000110324",
  "IL6" = "ENSG00000136244",
  "IMPDH2" = "ENSG00000178035",
  "IRF1" = "ENSG00000125347",
  "ITGA1" = "ENSG00000213949",
  "JAK1" = "ENSG00000162434",
  "KCNA5" = "ENSG00000130037",
  "LEF1" = "ENSG00000138795",
  "LIAS" = "ENSG00000121897",
  "LIPT1" = "ENSG00000144182",
  "LYN" = "ENSG00000254087",
  "MAN2B1" = "ENSG00000104774",
  "MAOB" = "ENSG00000069535",
  "MAPKAPK2" = "ENSG00000162889",
  "MB" = "ENSG00000198125",
  "MGAT1" = "ENSG00000131446",
  "MKNK2" = "ENSG00000099875",
  "MTOR" = "ENSG00000198793",
  "NFE2L2" = "ENSG00000116044",
  "NFKB1" = "ENSG00000109320",
  "NMT1" = "ENSG00000136448",
  "NR1H3" = "ENSG00000025434",
  "NR3C1" = "ENSG00000113580",
  "NUDCD2" = "ENSG00000170584",
  "OAZ1" = "ENSG00000104904",
  "PA2G4" = "ENSG00000170515",
  "PABPC1" = "ENSG00000070756",
  "PARP1" = "ENSG00000143799",
  "PARP12" = "ENSG00000059378",
  "PCTP" = "ENSG00000141179",
  "PDCD6" = "ENSG00000249915",
  "PDCD6IP" = "ENSG00000170248",
  "PDHB" = "ENSG00000168291",
  "PDHX" = "ENSG00000110435",
  "PLA2G4A" = "ENSG00000116711",
  "POLR2D" = "ENSG00000144231",
  "PPARG" = "ENSG00000132170",
  "PRKCQ" = "ENSG00000065675",
  "PTGR1" = "ENSG00000106853",
  "PTGS1" = "ENSG00000095303",
  "PTPN2" = "ENSG00000175354",
  "PTPN6" = "ENSG00000111679",
  "PTPRC" = "ENSG00000081237",
  "PTPRF" = "ENSG00000142949",
  "PTPRJ" = "ENSG00000149177",
  "RBM39" = "ENSG00000131051",
  "RELA" = "ENSG00000173039",
  "RENBP" = "ENSG00000102032",
  "RHOC" = "ENSG00000155366",
  "S100A6" = "ENSG00000197956",
  "SAT1" = "ENSG00000111371",
  "SAT2" = "ENSG00000134294",
  "SCN9A" = "ENSG00000169432",
  "SEC13" = "ENSG00000157020",
  "SERPINB10" = "ENSG00000242550",
  "SLC31A1" = "ENSG00000136868",
  "SPHK1" = "ENSG00000176170",
  "SREBF1" = "ENSG00000072310",
  "STARD13" = "ENSG00000133121",
  "STAT1" = "ENSG00000115415",
  "STAT3" = "ENSG00000168610",
  "STAT5A" = "ENSG00000126561",
  "STK4" = "ENSG00000101109",
  "TBXAS1" = "ENSG00000059377",
  "TCN2" = "ENSG00000185339",
  "TDP1" = "ENSG00000042088",
  "TGFB1" = "ENSG00000105329",
  "TIMP1" = "ENSG00000102265",
  "TNF" = "ENSG00000232810",
  "TOP2A" = "ENSG00000131747",
  "TSPO" = "ENSG00000100300",
  "XDH" = "ENSG00000158125",
  "XRCC6" = "ENSG00000196419",
  "ZEB1" = "ENSG00000148516",
  "ZHX2" = "ENSG00000178764"
)

# 配置
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_data_preparation"
EXPOSURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_matched/matched_data"
OUTCOME_FILE <- "D:/EQTL/mr_results_megastroke/megastroke_outcome_146genes.csv"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", OUTPUT_DIR, "\n\n")

# 步骤 1: 检查已有的暴露数据
cat("步骤 1: 检查已有的暴露数据\n")
cat("----------------------------------------------------------------------\n")

exposure_files <- list.files(EXPOSURE_DIR, pattern = "_exposure\\.csv$", full.names = TRUE)
cat(sprintf("暴露数据文件数：%d\n", length(exposure_files)))

# 建立 ENSG ID 到文件的映射
ensg_to_file <- list()
for (exposure_file in exposure_files) {
  # 提取 ENSG ID（不含版本号）
  ensg_id <- basename(exposure_file)
  ensg_base <- sub("\\..*", "", ensg_id)
  ensg_to_file[[ensg_base]] <- exposure_file
}

cat(sprintf("唯一 ENSG ID 数：%d\n\n", length(ensg_to_file)))

# 步骤 2: 为每个目标基因准备数据
cat("步骤 2: 为每个目标基因准备数据\n")
cat("----------------------------------------------------------------------\n")

prepared_genes <- list()
missing_genes <- character(0)

for (gene in target_genes) {
  ensg_id <- gene_to_ensg[[gene]]
  
  if (is.null(ensg_id)) {
    cat(sprintf("✗ %s: 无 ENSG ID 映射\n", gene))
    missing_genes <- c(missing_genes, gene)
    next
  }
  
  ensg_base <- sub("\\..*", "", ensg_id)
  
  exposure_file <- ensg_to_file[[ensg_base]]
  
  if (is.null(exposure_file)) {
    cat(sprintf("✗ %s (%s): 无暴露数据文件\n", gene, ensg_id))
    missing_genes <- c(missing_genes, gene)
    next
  }
  
  tryCatch({
    data <- fread(exposure_file, stringsAsFactors = FALSE)
    
    if (nrow(data) == 0) {
      cat(sprintf("✗ %s (%s): 暴露数据为空\n", gene, ensg_id))
      missing_genes <- c(missing_genes, gene)
      next
    }
    
    data$gene_symbol <- gene
    
    prepared_genes[[gene]] <- data
    cat(sprintf("✓ %s (%s): %d SNPs\n", gene, ensg_id, nrow(data)))
  }, error = function(e) {
    cat(sprintf("✗ %s (%s): 读取失败 - %s\n", gene, ensg_id, e$message))
    missing_genes <- c(missing_genes, gene)
  })
}

cat(sprintf("\n准备成功：%d/%d 个基因\n", length(prepared_genes), length(target_genes)))
cat(sprintf("缺少数据：%d 个基因\n\n", length(missing_genes)))

if (length(missing_genes) > 0) {
  cat("缺少数据的基因:\n")
  cat(paste(missing_genes, collapse = ", "), "\n\n")
}

# 步骤 3: 保存准备结果
cat("步骤 3: 保存准备结果\n")
cat("----------------------------------------------------------------------\n")

# 保存成功准备的基因列表
if (length(prepared_genes) > 0) {
  prepared_df <- do.call(rbind, lapply(names(prepared_genes), function(gene) {
    data.frame(
      gene = gene,
      ensg_id = gene_to_ensg[[gene]],
      nsnp = nrow(prepared_genes[[gene]]),
      stringsAsFactors = FALSE
    )
  }))
  
  prepared_file <- file.path(OUTPUT_DIR, "prepared_genes.csv")
  write.csv(prepared_df, prepared_file, row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("✓ 保存成功准备的基因：%s\n", prepared_file))
  cat(sprintf("  共 %d 个基因，%d 个 SNP\n\n", nrow(prepared_df), sum(prepared_df$nsnp)))
}

# 保存缺少数据的基因列表
if (length(missing_genes) > 0) {
  missing_df <- data.frame(gene = missing_genes, stringsAsFactors = FALSE)
  missing_file <- file.path(OUTPUT_DIR, "missing_genes.csv")
  write.csv(missing_df, missing_file, row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("✓ 保存缺少数据的基因：%s\n", missing_file))
  cat(sprintf("  共 %d 个基因\n\n", nrow(missing_df)))
}

# 步骤 4: MR 分析
cat("步骤 4: MR 分析\n")
cat("----------------------------------------------------------------------\n")

if (!file.exists(OUTCOME_FILE)) {
  stop("结局数据文件不存在：", OUTCOME_FILE)
}

outcome_data <- fread(OUTCOME_FILE, stringsAsFactors = FALSE)
cat(sprintf("加载结局数据：%d 个 SNP\n\n", nrow(outcome_data)))

# 创建 chr:pos 键
outcome_data$chr_pos <- paste(outcome_data$chr, outcome_data$pos.outcome, sep = ":")
cat(sprintf("创建 %d 个 chr:pos 键\n\n", nrow(outcome_data)))

mr_results <- list()

for (gene_symbol in names(prepared_genes)) {
  exposure <- prepared_genes[[gene_symbol]]
  
  cat(sprintf("\n%s:\n", gene_symbol))
  
  # 创建暴露的 chr:pos 键
  exposure$chr_pos <- paste(exposure$CHR, exposure$BP, sep = ":")
  
  # 按 chr:pos 匹配
  common_chr_pos <- intersect(exposure$chr_pos, outcome_data$chr_pos)
  
  cat(sprintf("  暴露 chr:pos: %d 个\n", length(unique(exposure$chr_pos))))
  cat(sprintf("  匹配的 chr:pos: %d 个\n", length(common_chr_pos)))
  
  if (length(common_chr_pos) < 3) {
    cat(sprintf("  ✗ 匹配的 SNP 太少 (%d)\n", length(common_chr_pos)))
    next
  }
  
  # 提取匹配的数据
  exp_matched <- exposure[exposure$chr_pos %in% common_chr_pos, ]
  out_matched <- outcome_data[outcome_data$chr_pos %in% common_chr_pos, ]
  
  # 排序
  exp_matched <- exp_matched[order(exp_matched$chr_pos), ]
  out_matched <- out_matched[order(out_matched$chr_pos), ]
  
  # 等位基因匹配
  allele_match <- (exp_matched$EFFECT_ALLELE == out_matched$effect_allele.outcome) |
                  (exp_matched$EFFECT_ALLELE == out_matched$other_allele.outcome)
  
  cat(sprintf("  等位基因匹配：%d/%d\n", sum(allele_match), nrow(exp_matched)))
  
  if (sum(allele_match) < 3) {
    cat(sprintf("  ✗ 等位基因匹配的 SNP 太少\n"))
    next
  }
  
  exp_matched <- exp_matched[allele_match, ]
  out_matched <- out_matched[allele_match, ]
  
  # 翻转 beta
  need_flip <- exp_matched$EFFECT_ALLELE != out_matched$effect_allele.outcome
  if (any(need_flip)) {
    exp_matched$BETA[need_flip] <- -exp_matched$BETA[need_flip]
  }
  
  # F 统计量
  f_stat <- mean((exp_matched$BETA / exp_matched$SE)^2, na.rm = TRUE)
  cat(sprintf("  F 统计量：%.2f\n", f_stat))
  
  # IVW 方法
  weights <- 1 / (out_matched$SE^2)
  beta_ivw <- sum(exp_matched$BETA * out_matched$BETA * weights, na.rm = TRUE) / 
              sum(weights, na.rm = TRUE)
  se_ivw <- sqrt(1 / sum(weights, na.rm = TRUE))
  pval_ivw <- 2 * pnorm(-abs(beta_ivw / se_ivw))
  
  # OR 和 CI
  or_ivw <- exp(beta_ivw)
  ci_low <- exp(beta_ivw - 1.96 * se_ivw)
  ci_high <- exp(beta_ivw + 1.96 * se_ivw)
  
  # 组织分布
  tissue_dist <- table(exp_matched$TISSUE)
  
  # 保存结果
  result <- data.frame(
    gene = gene_symbol,
    method = "IVW",
    beta = beta_ivw,
    se = se_ivw,
    or = or_ivw,
    ci_low = ci_low,
    ci_high = ci_high,
    pval = pval_ivw,
    f_stat = f_stat,
    nsnp = nrow(exp_matched),
    n_brain = ifelse("Brain_Cortex" %in% names(tissue_dist), tissue_dist["Brain_Cortex"], 0),
    n_blood = ifelse("Whole_Blood" %in% names(tissue_dist), tissue_dist["Whole_Blood"], 0),
    stringsAsFactors = FALSE
  )
  
  mr_results[[gene_symbol]] <- result
  
  cat(sprintf("  ✓ MR: OR=%.3f (%.3f-%.3f), P=%.2e\n", or_ivw, ci_low, ci_high, pval_ivw))
}

cat(sprintf("\n\n完成 %d 个基因的 MR 分析\n\n", length(mr_results)))

# 步骤 5: 保存 MR 结果
cat("步骤 5: 保存 MR 结果\n")
cat("----------------------------------------------------------------------\n")

if (length(mr_results) > 0) {
  all_results <- do.call(rbind, mr_results)
  all_results$fdr <- p.adjust(all_results$pval, method = "fdr")
  all_results <- all_results[order(all_results$pval), ]
  
  # 保存详细结果
  result_file <- file.path(OUTPUT_DIR, "mr_results_detailed.csv")
  write.csv(all_results, result_file, row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("  ✓ 保存详细结果：%d 个基因\n", nrow(all_results)))
  cat(sprintf("    文件：%s\n", result_file))
  
  # 保存显著结果
  sig_results <- all_results[all_results$fdr < 0.05, ]
  if (nrow(sig_results) > 0) {
    sig_file <- file.path(OUTPUT_DIR, "mr_results_significant.csv")
    write.csv(sig_results, sig_file, row.names = FALSE, fileEncoding = "UTF-8")
    cat(sprintf("  ✓ 保存显著结果 (FDR<0.05)：%d 个基因\n", nrow(sig_results)))
    cat(sprintf("    文件：%s\n", sig_file))
    
    # 打印显著结果
    cat("\n\n显著 MR 结果 (FDR < 0.05):\n")
    cat("----------------------------------------------------------------------\n")
    print(sig_results[, c("gene", "or", "ci_low", "ci_high", "pval", "fdr", "nsnp", "n_brain", "n_blood")])
  }
}

cat("\n======================================================================\n")
cat("完成！\n")
cat("======================================================================\n")
