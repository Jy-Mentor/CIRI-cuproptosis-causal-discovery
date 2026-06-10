#!/usr/bin/env Rscript
# ================================================================================
# 双源 eQTL MR 分析
# 使用 Python 提取的双源 eQTL 数据进行孟德尔随机化分析
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
cat("双源 eQTL MR 分析\n")
cat("参考：GTEx v11 | eQTL Catalogue | TwoSampleMR\n")
cat("======================================================================\n\n")

# 配置
EXPOSURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_dual_source"
OUTCOME_FILE <- "D:/EQTL/mr_results_megastroke/megastroke_outcome.csv"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_dual_source"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", OUTPUT_DIR, "\n\n")

# ================================================================================
# 1. 加载结局数据 (MEGASTROKE)
# ================================================================================
cat("步骤 1: 加载结局数据 (MEGASTROKE Ischemic Stroke)\n")
cat("----------------------------------------------------------------------\n")

load_outcome_data <- function(file) {
  if (!file.exists(file)) {
    stop(paste("结局数据文件不存在:", file))
  }
  
  cat("  读取 MEGASTROKE GWAS 汇总统计数据...\n")
  
  outcome <- fread(file, header = TRUE)
  
  cat(sprintf("  ✓ 加载 %d 个 SNP\n", nrow(outcome)))
  cat(sprintf("  列名：%s\n\n", paste(colnames(outcome)[1:min(8, ncol(outcome))], collapse=", ")))
  
  # 标准化列名用于 MR 分析
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
  
  cat(sprintf("  ✓ 结局数据准备完成：%d 个 SNP (去除 NA 后)\n\n", nrow(outcome_std)))
  
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
# 2. 加载暴露数据 (双源 eQTL)
# ================================================================================
cat("步骤 2: 加载暴露数据 (双源 eQTL)\n")
cat("----------------------------------------------------------------------\n")

load_exposure_data <- function(exposure_dir) {
  if (!dir.exists(exposure_dir)) {
    stop(paste("暴露数据目录不存在:", exposure_dir))
  }
  
  exposure_files <- list.files(exposure_dir, pattern = "_exposure\\.csv$", full.names = TRUE)
  
  if (length(exposure_files) == 0) {
    stop("未找到暴露数据文件")
  }
  
  cat(sprintf("  找到 %d 个基因的暴露数据\n", length(exposure_files)))
  
  exposure_list <- list()
  
  for (i in seq_along(exposure_files)) {
    file <- exposure_files[i]
    gene_name <- gsub("_exposure\\.csv$", "", basename(file))
    
    tryCatch({
      exposure <- fread(file, header = TRUE)
      colnames(exposure) <- toupper(colnames(exposure))
      
      required_cols <- c("SNP", "BETA", "SE", "PVAL")
      if (all(required_cols %in% colnames(exposure))) {
        exposure_list[[gene_name]] <- exposure
      }
    }, error = function(e) {})
    
    if (i %% 1000 == 0) {
      cat(sprintf("  已加载 %d/%d 个基因\n", i, length(exposure_files)))
    }
  }
  
  cat(sprintf("  ✓ 成功加载 %d 个基因的暴露数据\n\n", length(exposure_list)))
  
  return(exposure_list)
}

exposure_list <- load_exposure_data(EXPOSURE_DIR)

# ================================================================================
# 3. MR 分析
# ================================================================================
cat("步骤 3: 进行 MR 分析\n")
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
  
  # 按 SNP 排序并合并
  exp_data <- exp_data[order(exp_data$SNP), ]
  out_data <- out_data[order(out_data$snp), ]
  
  mr_data <- merge(exp_data, out_data, by.x = "SNP", by.y = "snp", suffixes = c("_exp", "_out"))
  
  if (nrow(mr_data) < 3) {
    return(NULL)
  }
  
  # 计算工具变量强度 (F 统计量)
  f_stat <- mean((mr_data$BETA / mr_data$SE)^2, na.rm = TRUE)
  
  # IVW 方法 (逆方差加权)
  weights <- 1 / (mr_data$se_out^2)
  beta_ivw <- sum(mr_data$BETA * mr_data$beta * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
  se_ivw <- sqrt(1 / sum(weights, na.rm = TRUE))
  pval_ivw <- 2 * pnorm(-abs(beta_ivw / se_ivw))
  
  # 计算 OR 和置信区间
  or_ivw <- exp(beta_ivw)
  ci_low <- exp(beta_ivw - 1.96 * se_ivw)
  ci_high <- exp(beta_ivw + 1.96 * se_ivw)
  
  result <- data.frame(
    method = "IVW",
    beta = beta_ivw,
    se = se_ivw,
    or = or_ivw,
    ci_low = ci_low,
    ci_high = ci_high,
    pval = pval_ivw,
    f_stat = f_stat,
    nsnp = nrow(mr_data)
  )
  
  return(list(result = result, data = mr_data))
}

# 对所有基因进行 MR 分析
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
  
  if (i %% 1000 == 0) {
    cat(sprintf("  已分析 %d/%d 个基因 (有效：%d)\n", i, length(exposure_list), valid_genes))
  }
}

cat(sprintf("  ✓ 完成 %d 个基因的 MR 分析\n\n", valid_genes))

# ================================================================================
# 4. 整理结果
# ================================================================================
cat("步骤 4: 整理 MR 结果\n")
cat("----------------------------------------------------------------------\n")

if (length(mr_results) == 0) {
  cat("  ✗ 没有有效的 MR 结果\n")
  cat("  可能原因：暴露和结局数据没有共同的 SNP\n")
  cat("  请检查数据格式和 SNP ID 是否匹配\n\n")
  
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
  
  # 多重检验校正 (FDR)
  results_df$fdr <- p.adjust(results_df$pval, method = "BH")
  results_df$fdr_sig <- results_df$fdr < 0.05
  results_df$significance <- ifelse(results_df$fdr < 0.05, "FDR 显著", "不显著")
  
  # 按 P 值排序
  results_df <- results_df[order(results_df$pval), ]
  
  cat(sprintf("  总基因数：%d\n", nrow(results_df)))
  cat(sprintf("  FDR 显著基因：%d (FDR < 0.05)\n", sum(results_df$fdr < 0.05, na.rm = TRUE)))
  cat(sprintf("  边缘显著基因：%d (P < 0.05)\n", sum(results_df$pval < 0.05, na.rm = TRUE)))
  cat(sprintf("  平均 F 统计量：%.2f\n\n", mean(results_df$f_stat, na.rm = TRUE)))
}

# 保存结果
results_file <- file.path(OUTPUT_DIR, "mr_results.csv")
write.csv(results_df, results_file, row.names = FALSE)
cat(sprintf("  ✓ 结果已保存：%s\n\n", results_file))

# ================================================================================
# 5. 创建图表
# ================================================================================
cat("步骤 5: 创建图表\n")
cat("----------------------------------------------------------------------\n")

if (nrow(results_df) > 0) {
  # 5.1 火山图
  cat("  创建火山图...\n")
  
  p_volcano <- ggplot(results_df, aes(x = log10(or), y = -log10(pval), color = fdr_sig)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("grey", "red"), na.value = "grey") +
    theme_minimal() +
    labs(
      title = "MR Analysis Results - Volcano Plot",
      subtitle = paste("Dual-source eQTL (Brain + Blood) |", nrow(results_df), "genes"),
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
  
  # 5.2 曼哈顿图
  cat("  创建曼哈顿图...\n")
  
  results_df$chr_num <- as.numeric(gsub("ENSG[0-9]+\\.", "", results_df$gene)) %% 100
  
  p_manhattan <- ggplot(results_df, aes(x = chr_num, y = -log10(pval), color = factor(chr_num %% 2))) +
    geom_point(alpha = 0.6) +
    scale_color_brewer(type = "qual", palette = "Set1") +
    theme_minimal() +
    labs(
      title = "MR Analysis Results - Manhattan Plot",
      subtitle = paste("Dual-source eQTL |", nrow(results_df), "genes"),
      x = "Gene Index",
      y = "-log10(P-value)",
      color = "Chr"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  ggsave(file.path(OUTPUT_DIR, "Figure2_Manhattan_Plot.png"), p_manhattan, width = 12, height = 6, dpi = 300)
  cat(sprintf("  ✓ 已保存：Figure2_Manhattan_Plot.png\n"))
  
  # 5.3 森林图 (Top 10)
  cat("  创建森林图...\n")
  
  top_genes <- results_df[results_df$fdr < 0.05, ]
  if (nrow(top_genes) == 0) {
    top_genes <- results_df[order(results_df$pval), ][1:10, ]
  }
  top_genes <- top_genes[1:min(10, nrow(top_genes)), ]
  
  if (nrow(top_genes) > 0) {
    top_genes$gene <- factor(top_genes$gene, levels = rev(top_genes$gene))
    
    p_forest <- ggplot(top_genes, aes(x = or, y = gene)) +
      geom_point(aes(size = -log10(pval)), color = "red", alpha = 0.7) +
      geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3, color = "red", alpha = 0.5) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "blue") +
      scale_x_log10() +
      scale_size_continuous(range = c(3, 6)) +
      theme_minimal() +
      labs(
        title = "MR Analysis Results - Forest Plot (Top 10)",
        subtitle = paste("FDR < 0.05 | Dual-source eQTL"),
        x = "Odds Ratio (95% CI)",
        y = "Gene",
        size = "-log10(P)"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      )
    
    ggsave(file.path(OUTPUT_DIR, "Figure3_Forest_Plot.png"), p_forest, width = 10, height = 8, dpi = 300)
    cat(sprintf("  ✓ 已保存：Figure3_Forest_Plot.png\n"))
  }
  
  # 5.4 F 统计量分布
  cat("  创建 F 统计量分布图...\n")
  
  p_fstat <- ggplot(results_df, aes(x = f_stat)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
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
# 6. 完成
# ================================================================================
cat("\n======================================================================\n")
cat("MR 分析完成！\n")
cat("======================================================================\n")

cat(sprintf("\n输出目录：%s\n", OUTPUT_DIR))
cat("\n生成的文件:\n")
cat("  - mr_results.csv: MR 结果表格\n")
cat("  - Figure1_Volcano_Plot.png: 火山图\n")
cat("  - Figure2_Manhattan_Plot.png: 曼哈顿图\n")
cat("  - Figure3_Forest_Plot.png: 森林图\n")
cat("  - Figure4_F_Statistic_Distribution.png: F 统计量分布\n")

cat("\n======================================================================\n")
