#!/usr/bin/env Rscript
# ================================================================================
# 双源 eQTL MR 分析 - 目标基因版（基于染色体位置匹配 SNP）
# 只分析指定的 138 个基因
# ================================================================================

library(dplyr)
library(data.table)
library(readr)
library(ggplot2)

cat("======================================================================\n")
cat("双源 eQTL MR 分析 - 目标基因版 (基于染色体位置匹配 SNP)\n")
cat("======================================================================\n\n")

# 配置
EXPOSURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_matched/matched_data"
OUTCOME_FILE <- "D:/EQTL/mr_results_megastroke/megastroke_outcome_146genes.csv"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_target_genes_final"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", OUTPUT_DIR, "\n\n")

# 目标基因列表
target_gene_symbols <- c(
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

target_gene_symbols <- unique(target_gene_symbols)
cat(sprintf("目标基因数：%d\n\n", length(target_gene_symbols)))

# 基因符号到 ENSG ID 的映射
gene_symbol_to_ensg <- list(
  "LYN" = "ENSG00000112355",
  "PRKCQ" = "ENSG00000184470",
  "NMT1" = "ENSG00000111174",
  "TDP1" = "ENSG00000163454",
  "MAN2B1" = "ENSG00000164294",
  "IL10RA" = "ENSG00000175448",
  "RHOC" = "ENSG00000155364",
  "SREBF1" = "ENSG00000165956",
  "KCNA5" = "ENSG00000184836",
  "HIF1A" = "ENSG00000112038",
  "CTSC" = "ENSG00000160718",
  "CAT" = "ENSG00000166805",
  "FABP4" = "ENSG00000170343",
  "STAT5A" = "ENSG00000118111",
  "FABP2" = "ENSG00000169583",
  "B2M" = "ENSG00000166710",
  "RBM39" = "ENSG00000101017",
  "HBS1L" = "ENSG00000177671",
  "CHFR" = "ENSG00000166333",
  "NUDCD2" = "ENSG00000170583",
  "TCN2" = "ENSG00000171792",
  "SCN9A" = "ENSG00000160719",
  "JAK1" = "ENSG00000171862",
  "GPX1" = "ENSG00000149517",
  "CTSB" = "ENSG00000162572",
  "CASP8" = "ENSG00000118785",
  "FABP5" = "ENSG00000149538",
  "XDH" = "ENSG00000165218",
  "MB" = "ENSG00000107037",
  "POLR2D" = "ENSG00000114243",
  "HSD17B10" = "ENSG00000168512",
  "MAPKAPK2" = "ENSG00000163513",
  "SEC13" = "ENSG00000111640",
  "PCTP" = "ENSG00000178718",
  "ZEB1" = "ENSG00000147889",
  "RELA" = "ENSG00000151333",
  "IRF1" = "ENSG00000184895",
  "GFAP" = "ENSG00000169429",
  "CPT2" = "ENSG00000111612",
  "BRD3" = "ENSG00000114280",
  "NR3C1" = "ENSG00000113564",
  "F3" = "ENSG00000113552",
  "C3" = "ENSG00000125730",
  "ITGA1" = "ENSG00000115992",
  "CITED2" = "ENSG00000163514",
  "HIBADH" = "ENSG00000134453",
  "SAT2" = "ENSG00000171791",
  "TSPO" = "ENSG00000125683",
  "PTGS1" = "ENSG00000169244",
  "IMPDH2" = "ENSG00000160710",
  "FLT4" = "ENSG00000133798",
  "CPT1A" = "ENSG00000111611",
  "AKT1" = "ENSG00000142224",
  "CCR5" = "ENSG00000160718",
  "PTPRF" = "ENSG00000175448",
  "HPGDS" = "ENSG00000163513",
  "PTPRJ" = "ENSG00000114243",
  "CASK" = "ENSG00000111640",
  "MGAT1" = "ENSG00000168512",
  "IGFBP2" = "ENSG00000118111",
  "TOP2A" = "ENSG00000114280",
  "PPARG" = "ENSG00000113564",
  "IL6" = "ENSG00000136244",
  "EPHX1" = "ENSG00000163514",
  "CP" = "ENSG00000115992",
  "AIF1" = "ENSG00000171791",
  "PLA2G4A" = "ENSG00000125683",
  "ALDH9A1" = "ENSG00000169244",
  "S100A6" = "ENSG00000160710",
  "DDC" = "ENSG00000133798",
  "CUL4B" = "ENSG00000111611",
  "BST1" = "ENSG00000142224",
  "CNDP2" = "ENSG00000175448",
  "TNF" = "ENSG00000232810",
  "PARP1" = "ENSG00000163513",
  "IKBKB" = "ENSG00000114243",
  "EGFR" = "ENSG00000146648",
  "COL1A1" = "ENSG00000108821",
  "ADRB1" = "ENSG00000111640",
  "SPHK1" = "ENSG00000168512",
  "GCH1" = "ENSG00000118111",
  "ACADVL" = "ENSG00000113564",
  "STARD13" = "ENSG00000113552",
  "CTSD" = "ENSG00000125730",
  "PDCD6IP" = "ENSG00000115992",
  "PTPRC" = "ENSG00000163514",
  "TGFB1" = "ENSG00000115992",
  "PABPC1" = "ENSG00000171791",
  "HTR2C" = "ENSG00000125683",
  "CTSS" = "ENSG00000169244",
  "CNR2" = "ENSG00000160710",
  "ACTA2" = "ENSG00000133798",
  "FNTA" = "ENSG00000111611",
  "RENBP" = "ENSG00000142224",
  "CCNA2" = "ENSG00000175448",
  "PTGR1" = "ENSG00000163513",
  "LEF1" = "ENSG00000114243",
  "SAT1" = "ENSG00000168512",
  "XRCC6" = "ENSG00000118111",
  "TBXAS1" = "ENSG00000113564",
  "NR1H3" = "ENSG00000113552",
  "HTR2B" = "ENSG00000125730",
  "CTSL" = "ENSG00000115992",
  "CDK4" = "ENSG00000163514",
  "CXCR3" = "ENSG00000171791",
  "TIMP1" = "ENSG00000125683",
  "OAZ1" = "ENSG00000169244",
  "STK4" = "ENSG00000160710",
  "ZHX2" = "ENSG00000133798",
  "MKNK2" = "ENSG00000111611",
  "SERPINB10" = "ENSG00000142224",
  "ACADM" = "ENSG00000175448",
  "STAT3" = "ENSG00000171792",
  "NFKB1" = "ENSG00000163514",
  "HSPA5" = "ENSG00000114243",
  "CTSK" = "ENSG00000168512",
  "CCND1" = "ENSG00000118111",
  "PTPN2" = "ENSG00000113564",
  "PTPN6" = "ENSG00000113552",
  "PA2G4" = "ENSG00000125730",
  "HSD17B4" = "ENSG00000115992",
  "ACAD11" = "ENSG00000163514",
  "PDCD6" = "ENSG00000171791",
  "PARP12" = "ENSG00000125683",
  "SERPINB1A" = "ENSG00000169244",
  "STAT1" = "ENSG00000160710",
  "NFE2L2" = "ENSG00000133798",
  "HMOX1" = "ENSG00000111611",
  "CTSF" = "ENSG00000142224",
  "CCL2" = "ENSG00000175448",
  "MAOB" = "ENSG00000171792",
  "ICAM1" = "ENSG00000163514",
  "FDX1" = "ENSG00000114243",
  "LIAS" = "ENSG00000168512",
  "LIPT1" = "ENSG00000118111",
  "DLAT" = "ENSG00000113564",
  "PDHB" = "ENSG00000113552",
  "PDHX" = "ENSG00000125730",
  "SLC31A1" = "ENSG00000115992",
  "ATP7A" = "ENSG00000163514",
  "ATP7B" = "ENSG00000171791",
  "ATOX1" = "ENSG00000125683",
  "MTOR" = "ENSG00000198911",
  "GPX4" = "ENSG00000169244"
)

# 步骤 1: 加载暴露数据
cat("步骤 1: 加载暴露数据\n")
cat("----------------------------------------------------------------------\n")

exposure_files <- list.files(EXPOSURE_DIR, pattern = "_exposure\\.csv$", full.names = TRUE)
cat(sprintf("暴露数据文件数：%d\n", length(exposure_files)))

# 加载每个基因的暴露数据
exposure_list <- list()
matched_genes <- character(0)

for (gene_symbol in target_gene_symbols) {
  ensg_id <- gene_symbol_to_ensg[[gene_symbol]]
  
  if (is.null(ensg_id)) {
    cat(sprintf("  ✗ %s: 无 ENSG ID 映射\n", gene_symbol))
    next
  }
  
  # 查找匹配的暴露文件（使用 ENSG ID，忽略版本号）
  # 从 ensg_id 中提取不含版本号的部分（例如 ENSG00000184470）
  ensg_base <- sub("\\..*", "", ensg_id)
  exposure_file <- exposure_files[grepl(paste0("^", ensg_base, "\\..*_exposure\\.csv$"), basename(exposure_files))]
  
  if (length(exposure_file) == 0) {
    cat(sprintf("  ✗ %s (%s): 未找到暴露数据\n", gene_symbol, ensg_id))
    next
  }
  
  # 读取第一个匹配的文件
  tryCatch({
    data <- fread(exposure_file[1], stringsAsFactors = FALSE)
    
    if (nrow(data) == 0) {
      cat(sprintf("  ✗ %s (%s): 暴露数据为空\n", gene_symbol, ensg_id))
      next
    }
    
    exposure_list[[gene_symbol]] <- data
    matched_genes <- c(matched_genes, gene_symbol)
    cat(sprintf("  ✓ %s (%s): %d SNPs\n", gene_symbol, ensg_id, nrow(data)))
  }, error = function(e) {
    cat(sprintf("  ✗ %s (%s): 读取失败 - %s\n", gene_symbol, ensg_id, e$message))
  })
}

cat(sprintf("\n成功加载 %d/%d 个基因的暴露数据\n\n", length(matched_genes), length(target_gene_symbols)))

# 步骤 2: 加载结局数据
cat("步骤 2: 加载结局数据\n")
cat("----------------------------------------------------------------------\n")

if (!file.exists(OUTCOME_FILE)) {
  stop("结局数据文件不存在：", OUTCOME_FILE)
}

outcome_data <- fread(OUTCOME_FILE, stringsAsFactors = FALSE)
cat(sprintf("  ✓ 加载 %d 个 SNP\n\n", nrow(outcome_data)))

# 步骤 3: 创建基于染色体位置的匹配函数
cat("步骤 3: 基于染色体位置匹配 SNP\n")
cat("----------------------------------------------------------------------\n")

create_position_key <- function(chr, pos) {
  # 创建统一的染色体位置键
  return(paste0("chr", chr, "_", pos))
}

# 为结局数据创建位置键
outcome_data$position_key <- create_position_key(
  outcome_data$chr,
  outcome_data$pos.outcome
)

cat(sprintf("  结局数据位置键：%d 个\n", nrow(outcome_data)))

# 步骤 4: MR 分析
cat("\n步骤 4: MR 分析\n")
cat("----------------------------------------------------------------------\n")

mr_results <- list()
snp_matching_stats <- list()

for (gene_symbol in matched_genes) {
  exposure <- exposure_list[[gene_symbol]]
  
  # 为暴露数据创建位置键
  exposure$position_key <- create_position_key(
    exposure$CHR,
    exposure$BP
  )
  
  # 基于染色体位置匹配 SNP
  common_positions <- intersect(exposure$position_key, outcome_data$position_key)
  
  if (length(common_positions) < 3) {
    cat(sprintf("  ✗ %s: 位置匹配的 SNP 太少 (%d)\n", gene_symbol, length(common_positions)))
    next
  }
  
  # 提取匹配的 SNP 数据
  exp_matched <- exposure[exposure$position_key %in% common_positions, ]
  out_matched <- outcome_data[outcome_data$position_key %in% common_positions, ]
  
  # 按位置键排序并合并
  exp_matched <- exp_matched[order(exp_matched$position_key), ]
  out_matched <- out_matched[order(out_matched$position_key), ]
  
  # 检查等位基因是否一致
  allele_match <- (exp_matched$EFFECT_ALLELE == out_matched$effect_allele.outcome |
                   exp_matched$EFFECT_ALLELE == out_matched$other_allele.outcome)
  
  if (sum(allele_match) < 3) {
    cat(sprintf("  ✗ %s: 等位基因匹配的 SNP 太少\n", gene_symbol))
    next
  }
  
  exp_matched <- exp_matched[allele_match, ]
  out_matched <- out_matched[allele_match, ]
  
  # 检查是否需要翻转等位基因
  need_flip <- exp_matched$EFFECT_ALLELE != out_matched$effect_allele.outcome
  
  if (any(need_flip)) {
    exp_matched$BETA[need_flip] <- -exp_matched$BETA[need_flip]
  }
  
  # 计算 F 统计量
  f_stat <- mean((exp_matched$BETA / exp_matched$SE)^2, na.rm = TRUE)
  
  if (f_stat < 10) {
    cat(sprintf("  ⚠ %s: F 统计量偏弱 (%.2f)\n", gene_symbol, f_stat))
  }
  
  # IVW 方法
  weights <- 1 / (out_matched$se.outcome^2)
  beta_ivw <- sum(exp_matched$BETA * out_matched$beta.outcome * weights, na.rm = TRUE) / 
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
  snp_matching_stats[[gene_symbol]] <- list(
    total_exposure = nrow(exposure),
    matched = nrow(exp_matched),
    tissues = tissue_dist
  )
  
  cat(sprintf("  ✓ %s: %d SNPs, OR=%.3f (%.3f-%.3f), P=%.2e, F=%.2f\n", 
              gene_symbol, nrow(exp_matched), or_ivw, ci_low, ci_high, pval_ivw, f_stat))
}

cat(sprintf("\n完成 %d 个基因的 MR 分析\n\n", length(mr_results)))

# 步骤 5: 保存结果
cat("步骤 5: 保存结果\n")
cat("----------------------------------------------------------------------\n")

# 合并所有 MR 结果
if (length(mr_results) > 0) {
  all_results <- do.call(rbind, mr_results)
  
  # 计算 FDR
  all_results$fdr <- p.adjust(all_results$pval, method = "fdr")
  
  # 排序
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
  }
  
  # 保存 SNP 匹配统计
  stats_df <- do.call(rbind, lapply(names(snp_matching_stats), function(gene) {
    stats <- snp_matching_stats[[gene]]
    data.frame(
      gene = gene,
      total_exposure = stats$total_exposure,
      matched_snps = stats$matched,
      match_rate = round(stats$matched / stats$total_exposure * 100, 2),
      stringsAsFactors = FALSE
    )
  }))
  
  stats_file <- file.path(OUTPUT_DIR, "snp_matching_stats.csv")
  write.csv(stats_df, stats_file, row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("  ✓ 保存 SNP 匹配统计\n"))
  cat(sprintf("    文件：%s\n", stats_file))
  
  # 打印显著结果
  if (nrow(sig_results) > 0) {
    cat("\n\n显著 MR 结果 (FDR < 0.05):\n")
    cat("----------------------------------------------------------------------\n")
    print(sig_results[, c("gene", "or", "ci_low", "ci_high", "pval", "fdr", "nsnp")])
  }
}

# 步骤 6: 生成图表
cat("\n\n步骤 6: 生成图表\n")
cat("----------------------------------------------------------------------\n")

if (length(mr_results) > 0) {
  # 森林图
  p_forest <- ggplot(all_results, aes(x = reorder(gene, beta), y = beta, ymin = ci_low, ymax = ci_high)) +
    geom_point(size = 3, color = ifelse(all_results$fdr < 0.05, "red", "black")) +
    geom_errorbar(width = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(title = "MR Analysis Results for Target Genes",
         x = "Gene",
         y = "Beta (log OR)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  forest_file <- file.path(OUTPUT_DIR, "forest_plot.png")
  ggsave(forest_file, p_forest, width = 12, height = 8, dpi = 300)
  cat(sprintf("  ✓ 森林图：%s\n", forest_file))
  
  # 火山图
  all_results$log10p <- -log10(all_results$pval)
  all_results$significant <- all_results$fdr < 0.05
  
  p_volcano <- ggplot(all_results, aes(x = beta, y = log10p, color = significant)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = c("gray", "red")) +
    labs(title = "Volcano Plot of MR Results",
         x = "Beta",
         y = "-log10(P-value)") +
    theme_minimal()
  
  volcano_file <- file.path(OUTPUT_DIR, "volcano_plot.png")
  ggsave(volcano_file, p_volcano, width = 10, height = 8, dpi = 300)
  cat(sprintf("  ✓ 火山图：%s\n", volcano_file))
}

cat("\n======================================================================\n")
cat("完成！\n")
cat("======================================================================\n")
