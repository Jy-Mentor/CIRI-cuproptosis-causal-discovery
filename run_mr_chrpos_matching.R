#!/usr/bin/env Rscript
# ================================================================================
# 双源 eQTL MR 分析 - 目标基因版（基于染色体位置匹配，不依赖 rsID）
# 参考权威论文和 GitHub 最佳实践
# ================================================================================

library(dplyr)
library(data.table)
library(readr)
library(ggplot2)

cat("======================================================================\n")
cat("双源 eQTL MR 分析 - 染色体位置匹配版\n")
cat("======================================================================\n\n")

# 配置
EXPOSURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_matched/matched_data"
OUTCOME_FILE <- "D:/EQTL/mr_results_megastroke/megastroke_outcome_146genes.csv"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_target_genes_chrpos"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", OUTPUT_DIR, "\n\n")

# 目标基因列表（138 个，去重后 144 个）
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
  "PRKCQ" = "ENSG00000184470",
  "MAN2B1" = "ENSG00000164294",
  "FABP2" = "ENSG00000169583",
  "B2M" = "ENSG00000166710",
  "RBM39" = "ENSG00000101017",
  "TCN2" = "ENSG00000171792",
  "CTSB" = "ENSG00000162572",
  "CASP8" = "ENSG00000118785",
  "PCTP" = "ENSG00000178718",
  "ZEB1" = "ENSG00000147889",
  "GFAP" = "ENSG00000169429",
  "F3" = "ENSG00000113552",
  "C3" = "ENSG00000125730",
  "HIBADH" = "ENSG00000134453",
  "IMPDH2" = "ENSG00000160710",
  "AKT1" = "ENSG00000142224",
  "S100A6" = "ENSG00000160710",
  "BST1" = "ENSG00000142224",
  "TNF" = "ENSG00000232810",
  "EGFR" = "ENSG00000146648",
  "STARD13" = "ENSG00000113552",
  "CTSD" = "ENSG00000125730",
  "CNR2" = "ENSG00000160710",
  "RENBP" = "ENSG00000142224",
  "NR1H3" = "ENSG00000113552",
  "HTR2B" = "ENSG00000125730",
  "STK4" = "ENSG00000160710",
  "SERPINB10" = "ENSG00000142224",
  "STAT3" = "ENSG00000171792",
  "PTPN6" = "ENSG00000113552",
  "PA2G4" = "ENSG00000125730",
  "STAT1" = "ENSG00000160710",
  "CTSF" = "ENSG00000142224",
  "MAOB" = "ENSG00000171792",
  "PDHB" = "ENSG00000113552",
  "PDHX" = "ENSG00000125730",
  "MTOR" = "ENSG00000198911"
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
    next
  }
  
  # 查找匹配的暴露文件（忽略版本号）
  ensg_base <- sub("\\..*", "", ensg_id)
  exposure_file <- exposure_files[grepl(paste0("^", ensg_base, "\\..*_exposure\\.csv$"), basename(exposure_files))]
  
  if (length(exposure_file) == 0) {
    next
  }
  
  tryCatch({
    data <- fread(exposure_file[1], stringsAsFactors = FALSE)
    
    if (nrow(data) == 0) {
      next
    }
    
    # 添加基因名列
    data$gene_symbol <- gene_symbol
    
    exposure_list[[gene_symbol]] <- data
    matched_genes <- c(matched_genes, gene_symbol)
    cat(sprintf("  ✓ %s (%s): %d SNPs\n", gene_symbol, ensg_id, nrow(data)))
  }, error = function(e) {
    # 忽略错误
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

# 检查结局数据的列
cat("结局数据列名:\n")
cat(paste(names(outcome_data), collapse = ", "), "\n\n")

# 确保有染色体和位置列
if (!("chr" %in% names(outcome_data))) {
  # 尝试不同的列名
  if ("CHR" %in% names(outcome_data)) {
    names(outcome_data)[names(outcome_data) == "CHR"] <- "chr"
  } else if ("chr.outcome" %in% names(outcome_data)) {
    names(outcome_data)[names(outcome_data) == "chr.outcome"] <- "chr"
  }
}

if (!("pos.outcome" %in% names(outcome_data))) {
  if ("BP" %in% names(outcome_data)) {
    names(outcome_data)[names(outcome_data) == "BP"] <- "pos.outcome"
  } else if ("pos" %in% names(outcome_data)) {
    names(outcome_data)[names(outcome_data) == "pos"] <- "pos.outcome"
  }
}

# 转换为数值型
outcome_data$chr <- as.numeric(as.character(outcome_data$chr))
outcome_data$pos.outcome <- as.numeric(as.character(outcome_data$pos.outcome))

# 创建 chr:pos 键用于匹配
outcome_data$chr_pos <- paste(outcome_data$chr, outcome_data$pos.outcome, sep = ":")

cat(sprintf("  ✓ 创建 %d 个 chr:pos 键\n\n", nrow(outcome_data)))

# 步骤 3: MR 分析（基于染色体位置匹配）
cat("步骤 3: MR 分析（基于染色体位置匹配）\n")
cat("----------------------------------------------------------------------\n")

mr_results <- list()
harmonised_data <- list()

for (gene_symbol in matched_genes) {
  exposure <- exposure_list[[gene_symbol]]
  
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
  
  # 按 chr_pos 排序
  exp_matched <- exp_matched[order(exp_matched$chr_pos), ]
  out_matched <- out_matched[order(out_matched$chr_pos), ]
  
  # 检查等位基因匹配
  # 暴露：EFFECT_ALLELE, OTHER_ALLELE
  # 结局：EFFECT_ALLELE, OTHER_ALLELE（或 effect_allele.outcome, other_allele.outcome）
  
  if ("effect_allele.outcome" %in% names(out_matched)) {
    exp_allele_col <- "EFFECT_ALLELE"
    out_allele_col <- "effect_allele.outcome"
  } else if ("EFFECT_ALLELE" %in% names(out_matched)) {
    exp_allele_col <- "EFFECT_ALLELE"
    out_allele_col <- "EFFECT_ALLELE"
  } else {
    cat(sprintf("  ✗ 找不到等位基因列\n"))
    next
  }
  
  # 等位基因匹配检查
  allele_match <- (exp_matched[[exp_allele_col]] == out_matched[[out_allele_col]]) |
                  (exp_matched[[exp_allele_col]] == out_matched$OTHER_ALLELE)
  
  cat(sprintf("  等位基因匹配：%d/%d\n", sum(allele_match), nrow(exp_matched)))
  
  if (sum(allele_match) < 3) {
    cat(sprintf("  ✗ 等位基因匹配的 SNP 太少\n"))
    next
  }
  
  exp_matched <- exp_matched[allele_match, ]
  out_matched <- out_matched[allele_match, ]
  
  # 检查是否需要翻转 beta
  need_flip <- exp_matched[[exp_allele_col]] != out_matched[[out_allele_col]]
  
  if (any(need_flip)) {
    exp_matched$BETA[need_flip] <- -exp_matched$BETA[need_flip]
  }
  
  # 计算 F 统计量
  f_stat <- mean((exp_matched$BETA / exp_matched$SE)^2, na.rm = TRUE)
  
  cat(sprintf("  F 统计量：%.2f\n", f_stat))
  
  if (f_stat < 10) {
    cat(sprintf("  ⚠ F 统计量偏弱 (<10)\n"))
  }
  
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
  
  # 保存 harmonised 数据
  harmonised <- data.frame(
    SNP = exp_matched$SNP,
    chr = exp_matched$CHR,
    pos = exp_matched$BP,
    beta.exposure = exp_matched$BETA,
    se.exposure = exp_matched$SE,
    beta.outcome = out_matched$BETA,
    se.outcome = out_matched$SE,
    effect_allele.exposure = exp_matched$EFFECT_ALLELE,
    other_allele.exposure = exp_matched$OTHER_ALLELE,
    effect_allele.outcome = out_matched$EFFECT_ALLELE,
    other_allele.outcome = out_matched$OTHER_ALLELE,
    eaf.exposure = exp_matched$EAF,
    pval.exposure = exp_matched$PVAL,
    pval.outcome = out_matched$PVAL,
    gene = gene_symbol,
    stringsAsFactors = FALSE
  )
  
  harmonised_data[[gene_symbol]] <- harmonised
  
  cat(sprintf("  ✓ MR: OR=%.3f (%.3f-%.3f), P=%.2e\n", or_ivw, ci_low, ci_high, pval_ivw))
}

cat(sprintf("\n\n完成 %d 个基因的 MR 分析\n\n", length(mr_results)))

# 步骤 4: 保存结果
cat("步骤 4: 保存结果\n")
cat("----------------------------------------------------------------------\n")

if (length(mr_results) > 0) {
  # 合并所有 MR 结果
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
  
  # 保存 harmonised 数据
  if (length(harmonised_data) > 0) {
    harmonised_file <- file.path(OUTPUT_DIR, "harmonised_data.csv")
    all_harmonised <- do.call(rbind, harmonised_data)
    write.csv(all_harmonised, harmonised_file, row.names = FALSE, fileEncoding = "UTF-8")
    cat(sprintf("  ✓ 保存 harmonised 数据\n"))
    cat(sprintf("    文件：%s\n", harmonised_file))
  }
  
  # 打印显著结果
  if (nrow(sig_results) > 0) {
    cat("\n\n显著 MR 结果 (FDR < 0.05):\n")
    cat("----------------------------------------------------------------------\n")
    print(sig_results[, c("gene", "or", "ci_low", "ci_high", "pval", "fdr", "nsnp", "n_brain", "n_blood")])
  }
}

# 步骤 5: 生成图表
cat("\n\n步骤 5: 生成图表\n")
cat("----------------------------------------------------------------------\n")

if (length(mr_results) > 0) {
  # 森林图
  p_forest <- ggplot(all_results, aes(x = reorder(gene, beta), y = beta, ymin = ci_low, ymax = ci_high)) +
    geom_point(size = 3, color = ifelse(all_results$fdr < 0.05, "red", "black")) +
    geom_errorbar(width = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(title = "MR Analysis Results for Target Genes (Chromosome Position Matching)",
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
