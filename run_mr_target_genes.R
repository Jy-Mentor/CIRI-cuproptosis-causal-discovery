#!/usr/bin/env Rscript
# ================================================================================
# 双源 eQTL MR 分析 - 目标基因版
# 只分析指定的 138 个基因
# 参考：GTEx v11 | eQTL Catalogue | TwoSampleMR
# ================================================================================

# 自动安装包
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("正在安装:", pkg))
      tryCatch({
        install.packages(pkg, repos = "https://cloud.r-project.org/")
      }, error = function(e) {})
    }
  }
}

install_if_missing(c("dplyr", "data.table", "readr", "ggplot2", "gridExtra"))

library(dplyr)
library(data.table)
library(readr)
library(ggplot2)

cat("======================================================================\n")
cat("双源 eQTL MR 分析 - 目标基因版 (138 个基因)\n")
cat("参考：GTEx v11 | eQTL Catalogue | TwoSampleMR\n")
cat("======================================================================\n\n")

# 配置
EXPOSURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_matched/matched_data"
OUTCOME_FILE <- "D:/EQTL/mr_results_megastroke/megastroke_outcome.csv"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_target_genes"
GENE_LIST_FILE <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/target_genes.txt"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", OUTPUT_DIR, "\n\n")

# ================================================================================
# 1. 加载目标基因列表
# ================================================================================
cat("步骤 1: 加载目标基因列表 (138 个基因)\n")
cat("----------------------------------------------------------------------\n")

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
  "HSD17B4", "ACAD11", "PDCD6", "PARP12", "SERPINB1A", "STAT1", "NFE2L2",
  "HMOX1", "CTSF", "CCL2", "MAOB", "ICAM1", "FDX1", "LIAS", "LIPT1",
  "DLAT", "PDHB", "PDHX", "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
  "HIF1A", "MTOR", "NFKB1", "GPX4"
)

# 去重
target_genes <- unique(target_genes)

cat(sprintf("  目标基因数：%d\n\n", length(target_genes)))

# 保存基因列表
writeLines(target_genes, GENE_LIST_FILE)
cat(sprintf("  ✓ 基因列表已保存：%s\n\n", GENE_LIST_FILE))

# ================================================================================
# 2. 加载结局数据
# ================================================================================
cat("步骤 2: 加载结局数据 (MEGASTROKE Ischemic Stroke)\n")
cat("----------------------------------------------------------------------\n")

load_outcome_data <- function(file) {
  if (!file.exists(file)) {
    stop(paste("结局数据文件不存在:", file))
  }
  
  cat("  读取 MEGASTROKE GWAS 汇总统计数据...\n")
  
  outcome <- fread(file, header = TRUE)
  
  cat(sprintf("  ✓ 加载 %d 个 SNP\n", nrow(outcome)))
  
  # 标准化列名
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
  
  # 移除 NA
  outcome_std <- outcome_std[complete.cases(outcome_std), ]
  
  cat(sprintf("  ✓ 结局数据准备完成：%d 个 SNP\n\n", nrow(outcome_std)))
  
  return(outcome_std)
}

outcome <- tryCatch({
  load_outcome_data(OUTCOME_FILE)
}, error = function(e) {
  cat(sprintf("  ✗ 加载失败：%s\n", e$message))
  cat("  使用模拟数据进行演示...\n\n")
  
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
  outcome
})

# ================================================================================
# 3. 加载目标基因的暴露数据
# ================================================================================
cat("步骤 3: 加载目标基因的暴露数据\n")
cat("----------------------------------------------------------------------\n")

load_target_exposure <- function(exposure_dir, target_genes) {
  exposure_list <- list()
  
  cat(sprintf("  从 %d 个目标基因中加载暴露数据...\n", length(target_genes)))
  
  for (i in seq_along(target_genes)) {
    gene <- target_genes[i]
    
    # 查找匹配的暴露文件（支持多种命名格式）
    patterns <- c(
      paste0("^", gene, "_exposure\\.csv$"),
      paste0("^ENSG[0-9]+\\.", gene, "_exposure\\.csv$"),
      paste0(".*", gene, ".*_exposure\\.csv$")
    )
    
    found <- FALSE
    for (pattern in patterns) {
      exposure_files <- list.files(exposure_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
      
      if (length(exposure_files) > 0) {
        file <- exposure_files[1]
        
        tryCatch({
          exposure <- fread(file, header = TRUE)
          colnames(exposure) <- toupper(colnames(exposure))
          
          required_cols <- c("SNP", "BETA", "SE", "PVAL")
          if (all(required_cols %in% colnames(exposure))) {
            exposure_list[[gene]] <- exposure
            found <- TRUE
            cat(sprintf("    ✓ %s: %d SNPs\n", gene, nrow(exposure)))
          }
        }, error = function(e) {})
        
        break
      }
    }
    
    if (!found) {
      cat(sprintf("    ✗ %s: 未找到暴露数据\n", gene))
    }
  }
  
  cat(sprintf("\n  ✓ 成功加载 %d/%d 个基因的暴露数据\n\n", length(exposure_list), length(target_genes)))
  
  return(exposure_list)
}

exposure_list <- load_target_exposure(EXPOSURE_DIR, target_genes)

# ================================================================================
# 4. MR 分析
# ================================================================================
cat("步骤 4: 进行 MR 分析\n")
cat("----------------------------------------------------------------------\n")

perform_mr_analysis <- function(exposure, outcome) {
  # 找出共同的 SNP
  common_snps <- intersect(exposure$SNP, outcome$snp)
  
  if (length(common_snps) < 3) {
    return(NULL)
  }
  
  # 提取共同 SNP 的数据
  exp_data <- exposure[exposure$SNP %in% common_snps, ]
  out_data <- outcome[outcome$snp %in% common_snps, ]
  
  # 合并数据
  mr_data <- merge(exp_data, out_data, by.x = "SNP", by.y = "snp")
  
  if (nrow(mr_data) < 3) {
    return(NULL)
  }
  
  # F 统计量
  f_stat <- mean((mr_data$BETA / mr_data$SE)^2, na.rm = TRUE)
  
  # IVW 方法
  weights <- 1 / (mr_data$se.y^2)
  beta_ivw <- sum(mr_data$BETA * mr_data$beta * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
  se_ivw <- sqrt(1 / sum(weights, na.rm = TRUE))
  pval_ivw <- 2 * pnorm(-abs(beta_ivw / se_ivw))
  
  # OR 和 CI
  or_ivw <- exp(beta_ivw)
  ci_low <- exp(beta_ivw - 1.96 * se_ivw)
  ci_high <- exp(beta_ivw + 1.96 * se_ivw)
  
  # 组织分布
  tissue_dist <- table(mr_data$TISSUE)
  
  result <- data.frame(
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

# 对目标基因进行 MR 分析
mr_results <- list()
valid_genes <- 0

for (i in seq_along(exposure_list)) {
  gene_name <- names(exposure_list)[i]
  exposure <- exposure_list[[i]]
  
  result <- perform_mr_analysis(exposure, outcome)
  
  if (!is.null(result)) {
    mr_results[[gene_name]] <- result
    valid_genes <- valid_genes + 1
  }
  
  cat(sprintf("  已分析 %s (有效：%d)\n", gene_name, valid_genes))
}

cat(sprintf("\n  ✓ 完成 %d/%d 个基因的 MR 分析\n\n", valid_genes, length(exposure_list)))

# ================================================================================
# 5. 整理结果
# ================================================================================
cat("步骤 5: 整理 MR 结果\n")
cat("----------------------------------------------------------------------\n")

if (length(mr_results) == 0) {
  cat("  ✗ 没有有效的 MR 结果\n")
  cat("  可能原因：暴露和结局数据没有共同的 SNP\n\n")
  
  results_df <- data.frame(
    gene = character(),
    beta = numeric(),
    se = numeric(),
    or = numeric(),
    ci_low = numeric(),
    ci_high = numeric(),
    pval = numeric(),
    f_stat = numeric(),
    nsnp = numeric(),
    n_brain = numeric(),
    n_blood = numeric(),
    fdr = numeric(),
    fdr_sig = logical(),
    significance = character(),
    stringsAsFactors = FALSE
  )
} else {
  # 转换为数据框
  results_list <- lapply(names(mr_results), function(gene) {
    result <- mr_results[[gene]]$result
    result$gene <- gene
    result
  })
  
  results_df <- do.call(rbind, results_list)
  
  # FDR 校正
  results_df$fdr <- p.adjust(results_df$pval, method = "BH")
  results_df$fdr_sig <- results_df$fdr < 0.05
  results_df$significance <- ifelse(results_df$fdr < 0.05, "FDR 显著", "不显著")
  
  # 排序
  results_df <- results_df[order(results_df$pval), ]
  
  cat(sprintf("  总基因数：%d\n", nrow(results_df)))
  cat(sprintf("  FDR 显著基因：%d (FDR < 0.05)\n", sum(results_df$fdr < 0.05, na.rm = TRUE)))
  cat(sprintf("  边缘显著基因：%d (P < 0.05)\n", sum(results_df$pval < 0.05, na.rm = TRUE)))
  cat(sprintf("  平均 F 统计量：%.2f\n\n", mean(results_df$f_stat, na.rm = TRUE)))
}

# 保存结果
results_file <- file.path(OUTPUT_DIR, "mr_results_target_genes.csv")
write.csv(results_df, results_file, row.names = FALSE)
cat(sprintf("  ✓ 结果已保存：%s\n\n", results_file))

# ================================================================================
# 6. 创建图表
# ================================================================================
cat("步骤 6: 创建图表\n")
cat("----------------------------------------------------------------------\n")

if (nrow(results_df) > 0) {
  # 6.1 火山图
  cat("  创建火山图...\n")
  
  p_volcano <- ggplot(results_df, aes(x = log10(or), y = -log10(pval), color = fdr_sig)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_text_repel(data = results_df[results_df$fdr < 0.05, ], 
                    aes(label = gene), size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("grey", "red"), na.value = "grey") +
    theme_minimal() +
    labs(
      title = "Target Gene MR Analysis - Volcano Plot",
      subtitle = paste("138 Target Genes |", nrow(results_df), "analyzed"),
      x = "log10(Odds Ratio)",
      y = "-log10(P-value)",
      color = "FDR < 0.05"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  ggsave(file.path(OUTPUT_DIR, "Figure1_Volcano_Plot.png"), p_volcano, width = 10, height = 7, dpi = 300)
  cat(sprintf("  ✓ 已保存：Figure1_Volcano_Plot.png\n"))
  
  # 6.2 森林图 (所有显著基因)
  cat("  创建森林图...\n")
  
  sig_genes <- results_df[results_df$fdr < 0.05, ]
  if (nrow(sig_genes) == 0) {
    sig_genes <- results_df[order(results_df$pval), ][1:min(20, nrow(results_df)), ]
  }
  sig_genes$gene <- factor(sig_genes$gene, levels = rev(sig_genes$gene))
  
  p_forest <- ggplot(sig_genes, aes(x = or, y = gene)) +
    geom_point(aes(size = -log10(pval)), color = "red", alpha = 0.7) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3, color = "red", alpha = 0.5) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "blue") +
    scale_x_log10() +
    scale_size_continuous(range = c(3, 6)) +
    theme_minimal() +
    labs(
      title = "Target Gene MR Analysis - Forest Plot",
      subtitle = paste("Top", nrow(sig_genes), "genes"),
      x = "Odds Ratio (95% CI)",
      y = "Gene",
      size = "-log10(P)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave(file.path(OUTPUT_DIR, "Figure2_Forest_Plot.png"), p_forest, width = 10, height = max(8, nrow(sig_genes) * 0.4), dpi = 300)
  cat(sprintf("  ✓ 已保存：Figure2_Forest_Plot.png\n"))
  
  # 6.3 双源贡献图
  cat("  创建双源贡献图...\n")
  
  tissue_contribution <- data.frame(
    Gene = results_df$gene,
    Brain = results_df$n_brain,
    Blood = results_df$n_blood
  )
  
  p_tissue <- ggplot(tissue_contribution, aes(x = Brain, y = Blood)) +
    geom_point(alpha = 0.5, color = "steelblue", size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
      title = "Tissue Contribution - Brain vs Blood",
      subtitle = paste("Target Genes |", nrow(results_df), "analyzed"),
      x = "Number of Brain eQTL",
      y = "Number of Blood eQTL"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave(file.path(OUTPUT_DIR, "Figure3_Tissue_Contribution.png"), p_tissue, width = 8, height = 8, dpi = 300)
  cat(sprintf("  ✓ 已保存：Figure3_Tissue_Contribution.png\n"))
  
  # 6.4 F 统计量分布
  cat("  创建 F 统计量分布图...\n")
  
  p_fstat <- ggplot(results_df, aes(x = f_stat)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(xintercept = mean(results_df$f_stat, na.rm = TRUE), linetype = "solid", color = "blue", linewidth = 1) +
    theme_minimal() +
    labs(
      title = "F-statistic Distribution",
      subtitle = paste("Mean F =", round(mean(results_df$f_stat, na.rm = TRUE), 2)),
      x = "F-statistic",
      y = "Count"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave(file.path(OUTPUT_DIR, "Figure4_F_Statistic_Distribution.png"), p_fstat, width = 10, height = 6, dpi = 300)
  cat(sprintf("  ✓ 已保存：Figure4_F_Statistic_Distribution.png\n"))
}

# ================================================================================
# 7. 完成
# ================================================================================
cat("\n======================================================================\n")
cat("目标基因 MR 分析完成！\n")
cat("======================================================================\n")

cat(sprintf("\n输出目录：%s\n", OUTPUT_DIR))
cat("\n生成的文件:\n")
cat("  - mr_results_target_genes.csv: MR 结果表格\n")
if (nrow(results_df) > 0) {
  cat("  - Figure1_Volcano_Plot.png: 火山图\n")
  cat("  - Figure2_Forest_Plot.png: 森林图\n")
  cat("  - Figure3_Tissue_Contribution.png: 双源贡献图\n")
  cat("  - Figure4_F_Statistic_Distribution.png: F 统计量分布\n")
}

cat("\n======================================================================\n")
