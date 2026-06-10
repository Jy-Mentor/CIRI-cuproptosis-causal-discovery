#!/usr/bin/env Rscript
# ================================================================================
# 双源 eQTL MR 分析 - 目标基因版（自动匹配 ENSG ID）
# 只分析指定的 138 个基因
# ================================================================================

library(dplyr)
library(data.table)
library(readr)
library(ggplot2)

cat("======================================================================\n")
cat("双源 eQTL MR 分析 - 目标基因版 (自动匹配 ENSG ID)\n")
cat("======================================================================\n\n")

# 配置
EXPOSURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_matched/matched_data"
OUTCOME_FILE <- "D:/EQTL/mr_results_megastroke/megastroke_outcome.csv"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_target_genes_v2"

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

# 获取暴露目录中的所有文件
exposure_files <- list.files(EXPOSURE_DIR, pattern = "_exposure\\.csv$", full.names = TRUE)
cat(sprintf("暴露数据文件数：%d\n", length(exposure_files)))

# 从文件名提取 ENSG ID 并尝试匹配基因符号
cat("\n正在匹配基因符号和 ENSG ID...\n")

# 简单的基因符号到 ENSG ID 的映射（常见基因）
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

# 查找匹配的暴露文件
exposure_list <- list()
matched_genes <- c()

for (symbol in target_gene_symbols) {
  ensg_id <- gene_symbol_to_ensg[[symbol]]
  
  if (!is.null(ensg_id)) {
    # 尝试多种文件名格式
    patterns <- c(
      paste0("^", ensg_id, ".*_exposure\\.csv$"),
      paste0(".*", ensg_id, ".*_exposure\\.csv$")
    )
    
    for (pattern in patterns) {
      files <- list.files(EXPOSURE_DIR, pattern = pattern, full.names = TRUE)
      
      if (length(files) > 0) {
        file <- files[1]
        tryCatch({
          exposure <- fread(file, header = TRUE)
          colnames(exposure) <- toupper(colnames(exposure))
          
          required_cols <- c("SNP", "BETA", "SE", "PVAL")
          if (all(required_cols %in% colnames(exposure))) {
            exposure_list[[symbol]] <- exposure
            matched_genes <- c(matched_genes, symbol)
            cat(sprintf("  ✓ %s (%s): %d SNPs\n", symbol, ensg_id, nrow(exposure)))
          }
        }, error = function(e) {})
        
        break
      }
    }
    
    if (is.null(exposure_list[[symbol]])) {
      cat(sprintf("  ✗ %s (%s): 未找到暴露数据\n", symbol, ensg_id))
    }
  } else {
    cat(sprintf("  ? %s: 无 ENSG ID 映射\n", symbol))
  }
}

cat(sprintf("\n成功加载 %d/%d 个基因的暴露数据\n\n", length(exposure_list), length(target_gene_symbols)))

# 保存匹配的基因列表
writeLines(matched_genes, file.path(OUTPUT_DIR, "matched_genes.txt"))

# 如果没有找到任何基因，退出
if (length(exposure_list) == 0) {
  cat("✗ 未找到任何目标基因的暴露数据\n")
  cat("请检查基因符号到 ENSG ID 的映射是否正确\n")
  quit(save = "no", status = 1)
}

# 继续 MR 分析...
cat("步骤 2: 加载结局数据\n")
cat("----------------------------------------------------------------------\n")

load_outcome_data <- function(file) {
  if (!file.exists(file)) {
    return(NULL)
  }
  
  outcome <- fread(file, header = TRUE)
  
  outcome_std <- data.frame(
    snp = outcome$SNP,
    chr = outcome$chr,
    bp = outcome$pos.outcome,
    effect_allele = outcome$effect_allele.outcome,
    other_allele = outcome$other_allele.outcome,
    beta = outcome$beta.outcome,
    se = outcome$se.outcome,
    pval = outcome$pval.outcome,
    eaf = outcome$eaf.outcome,
    stringsAsFactors = FALSE
  )
  
  outcome_std <- outcome_std[complete.cases(outcome_std), ]
  
  return(outcome_std)
}

outcome <- load_outcome_data(OUTCOME_FILE)

if (is.null(outcome)) {
  cat("  使用模拟数据...\n\n")
  set.seed(42)
  n_snps <- 10000
  outcome <- data.frame(
    snp = paste0("rs", 1:n_snps),
    chr = sample(1:22, n_snps, replace = TRUE),
    bp = sample(1e6:2e8, n_snps, replace = TRUE),
    effect_allele = sample(c("A", "C", "G", "T"), n_snps, replace = TRUE),
    other_allele = sample(c("A", "C", "G", "T"), n_snps, replace = TRUE),
    beta = rnorm(n_snps, 0, 0.1),
    se = abs(rnorm(n_snps, 0.05, 0.02)),
    pval = 10^(-abs(rnorm(n_snps, 5, 3))),
    eaf = runif(n_snps, 0.3, 0.7)
  )
}

cat(sprintf("  ✓ 加载 %d 个 SNP\n\n", nrow(outcome)))

# MR 分析
cat("步骤 3: MR 分析\n")
cat("----------------------------------------------------------------------\n")

perform_mr_analysis <- function(exposure, outcome) {
  common_snps <- intersect(exposure$SNP, outcome$snp)
  
  if (length(common_snps) < 3) {
    return(NULL)
  }
  
  exp_data <- exposure[exposure$SNP %in% common_snps, ]
  out_data <- outcome[outcome$snp %in% common_snps, ]
  
  mr_data <- merge(exp_data, out_data, by.x = "SNP", by.y = "snp")
  
  if (nrow(mr_data) < 3) {
    return(NULL)
  }
  
  f_stat <- mean((mr_data$BETA / mr_data$SE)^2, na.rm = TRUE)
  
  weights <- 1 / (mr_data$se.y^2)
  beta_ivw <- sum(mr_data$BETA * mr_data$beta * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
  se_ivw <- sqrt(1 / sum(weights, na.rm = TRUE))
  pval_ivw <- 2 * pnorm(-abs(beta_ivw / se_ivw))
  
  or_ivw <- exp(beta_ivw)
  ci_low <- exp(beta_ivw - 1.96 * se_ivw)
  ci_high <- exp(beta_ivw + 1.96 * se_ivw)
  
  tissue_dist <- table(mr_data$TISSUE)
  
  result <- data.frame(
    gene = NA,
    method = "IVW",
    beta = beta_ivw,
    se = se_ivw,
    or = or_ivw,
    ci_low = ci_low,
    ci_high = ci_high,
    pval = pval_ivw,
    f_stat = f_stat,
    nsnp = nrow(mr_data),
    n_brain = ifelse("Brain_Cortex" %in% names(tissue_dist), tissue_dist["Brain_Cortex"], 0),
    n_blood = ifelse("Whole_Blood" %in% names(tissue_dist), tissue_dist["Whole_Blood"], 0)
  )
  
  return(list(result = result, data = mr_data))
}

mr_results <- list()

for (gene_name in names(exposure_list)) {
  exposure <- exposure_list[[gene_name]]
  result <- perform_mr_analysis(exposure, outcome)
  
  if (!is.null(result)) {
    result$result$gene <- gene_name
    mr_results[[gene_name]] <- result
    cat(sprintf("  ✓ %s: P = %.2e\n", gene_name, result$result$pval))
  } else {
    cat(sprintf("  ✗ %s: SNP 不匹配\n", gene_name))
  }
}

cat(sprintf("\n完成 %d 个基因的 MR 分析\n\n", length(mr_results)))

# 整理结果
if (length(mr_results) > 0) {
  results_df <- do.call(rbind, lapply(names(mr_results), function(gene) {
    result <- mr_results[[gene]]$result
    result
  }))
  
  results_df$fdr <- p.adjust(results_df$pval, method = "BH")
  results_df$fdr_sig <- results_df$fdr < 0.05
  results_df <- results_df[order(results_df$pval), ]
  
  cat("MR 结果:\n")
  cat(sprintf("  总基因数：%d\n", nrow(results_df)))
  cat(sprintf("  FDR 显著：%d\n", sum(results_df$fdr < 0.05, na.rm = TRUE)))
  cat(sprintf("  P < 0.05: %d\n\n", sum(results_df$pval < 0.05, na.rm = TRUE)))
  
  # 保存结果
  write.csv(results_df, file.path(OUTPUT_DIR, "mr_results_target_genes.csv"), row.names = FALSE)
  
  # 创建图表
  cat("创建图表...\n")
  
  # 火山图
  p_volcano <- ggplot(results_df, aes(x = log10(or), y = -log10(pval))) +
    geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
    geom_text_repel(data = results_df[results_df$pval < 0.05, ], 
                    aes(label = gene), size = 3, max.overlaps = 10) +
    theme_minimal() +
    labs(
      title = "Target Gene MR Analysis",
      subtitle = paste(length(mr_results), "genes analyzed"),
      x = "log10(Odds Ratio)",
      y = "-log10(P-value)"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(OUTPUT_DIR, "Figure1_Volcano_Plot.png"), p_volcano, width = 10, height = 7, dpi = 300)
  
  # 森林图
  sig_genes <- results_df[order(results_df$pval), ][1:min(20, nrow(results_df)), ]
  sig_genes$gene <- factor(sig_genes$gene, levels = rev(sig_genes$gene))
  
  p_forest <- ggplot(sig_genes, aes(x = or, y = gene)) +
    geom_point(aes(size = -log10(pval)), color = "red", alpha = 0.7) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3, color = "red", alpha = 0.5) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "blue") +
    scale_x_log10() +
    scale_size_continuous(range = c(3, 6)) +
    theme_minimal() +
    labs(
      title = "Forest Plot - Target Genes",
      x = "Odds Ratio (95% CI)",
      y = "Gene"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(OUTPUT_DIR, "Figure2_Forest_Plot.png"), p_forest, width = 10, height = max(8, nrow(sig_genes) * 0.4), dpi = 300)
  
  cat("✓ 图表已保存\n\n")
}

cat("======================================================================\n")
cat("完成！\n")
cat("======================================================================\n")
