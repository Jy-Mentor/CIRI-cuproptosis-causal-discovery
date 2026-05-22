# L3 遗传因果锚定：TwoSampleMR 双因果认证
# 输入：L2b 涌现宏节点代表基因
# 主方案：TwoSampleMR（R 包，IEU Open GWAS API）

library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(42)

# ==================== 0. 配置参数 ====================
OUTPUT_DIR <- "../results/L3_genetic_causal"
FIGURE_DIR <- "../figures/L3"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# 宏节点代表基因（示例，需根据 L2b 输出更新）
MACRO_NODE_REPRESENTATIVES <- c("FDX1", "DLAT", "LIAS", "SLC31A1", "ATOX1")

# 结局：CIRI/卒中 GWAS
OUTCOME_ID <- "ebi-a-GCST008526"  # 缺血性卒中 GWAS

# ==================== 1. 提取 eQTL 工具变量 ====================
cat("=== 步骤1: 提取 eQTL 工具变量 ===\n")

extract_eqtl <- function(gene, exposure_id) {
  cat(sprintf("  提取基因 %s 的 eQTL...\n", gene))
  
  # 从 eQTLGen 或 GTEx 提取
  eqtl_data <- extract_instruments(
    outcomes = exposure_id,
    p1 = 5e-8,
    p2 = 5e-8,
    clump = TRUE
  )
  
  if (is.null(eqtl_data) || nrow(eqtl_data) == 0) {
    cat(sprintf("    警告: 基因 %s 无显著 eQTL\n", gene))
    return(NULL)
  }
  
  # 计算 F 统计量
  eqtl_data[["F"]] <- (eqtl_data[["beta.exposure"]] / eqtl_data[["se.exposure"]])^2
  
  cat(sprintf("    找到 %d 个 SNP，F 统计量中位数: %.2f\n", 
              nrow(eqtl_data), median(eqtl_data[["F"]])))
  
  return(eqtl_data)
}

# ==================== 2. 提取结局数据 ====================
cat("\n=== 步骤2: 提取结局数据 ===\n")

extract_outcome <- function(exposure_data, outcome_id) {
  cat(sprintf("  提取结局 %s 的数据...\n", outcome_id))
  
  outcome_data <- extract_outcome_data(
    snps = exposure_data[["SNP"]],
    outcomes = outcome_id,
    proxies = FALSE
  )
  
  if (is.null(outcome_data) || nrow(outcome_data) == 0) {
    cat("    警告: 无匹配结局数据\n")
    return(NULL)
  }
  
  cat(sprintf("    找到 %d 个匹配 SNP\n", nrow(outcome_data)))
  
  return(outcome_data)
}

# ==================== 3. 数据协调 ====================
cat("\n=== 步骤3: 数据协调 ===\n")

harmonise_data <- function(exposure_data, outcome_data) {
  cat("  协调暴露和结局数据...\n")
  
  harmonised_data <- harmonise_data(
    exposure_dat = exposure_data,
    outcome_dat = outcome_data
  )
  
  cat(sprintf("  协调后保留 %d 个 SNP\n", nrow(harmonised_data)))
  
  return(harmonised_data)
}

# ==================== 4. MR 分析 ====================
cat("\n=== 步骤4: MR 分析 ===\n")

run_mr <- function(harmonised_data) {
  cat("  运行 MR 分析...\n")
  
  # 主分析：IVW
  mr_results <- mr(harmonised_data)
  
  # 敏感性分析
  mr_sensitivity <- mr_sensitivity_test(harmonised_data)
  
  # MR-Egger 截距检验（多效性）
  egger_intercept <- mr_pleiotropy_test(harmonised_data)
  
  # 异质性检验
  heterogeneity <- mr_heterogeneity(harmonised_data)
  
  # Steiger 检验（方向性）
  steiger <- directionality_test(harmonised_data)
  
  cat("  MR 结果:\n")
  print(mr_results)
  
  return(list(
    mr = mr_results,
    sensitivity = mr_sensitivity,
    egger = egger_intercept,
    heterogeneity = heterogeneity,
    steiger = steiger
  ))
}

# ==================== 5. 共定位分析 ====================
cat("\n=== 步骤5: 共定位分析 ===\n")

run_coloc <- function(gene, exposure_id, outcome_id) {
  cat(sprintf("  运行基因 %s 的共定位分析...\n", gene))
  
  # 使用 coloc 包
  # 注意：需要区域级别的汇总统计数据
  # 此处为示例框架，实际运行需提供区域数据
  
  cat("    共定位分析需要区域级别的 eQTL 和 GWAS 数据\n")
  cat("    请提供 GTEx v8 区域数据后运行\n")
  
  return(NULL)
}

# ==================== 6. 结果汇总 ====================
cat("\n=== 步骤6: 结果汇总 ===\n")

summarize_mr_results <- function(all_results, macro_representatives) {
  cat("  汇总 MR 结果...\n")
  
  summary_df <- data.frame()
  
  for (gene in macro_representatives) {
    if (gene %in% names(all_results)) {
      res <- all_results[[gene]]
      
      if (!is.null(res$mr) && nrow(res$mr) > 0) {
        ivw_result <- res$mr[res$mr[["method"]] == "Inverse variance weighted", ]
        
        if (nrow(ivw_result) > 0) {
          summary_df <- rbind(summary_df, data.frame(
            gene = gene,
            beta = ivw_result[["b"]],
            se = ivw_result[["se"]],
            pval = ivw_result[["pval"]],
            or = exp(ivw_result[["b"]]),
            or_ci_lower = exp(ivw_result[["b"]] - 1.96 * ivw_result[["se"]]),
            or_ci_upper = exp(ivw_result[["b"]] + 1.96 * ivw_result[["se"]]),
            direction = ifelse(ivw_result[["b"]] > 0, "risk", "protective")
          ))
        }
      }
    }
  }
  
  return(summary_df)
}

# ==================== 7. 可视化 ====================
cat("\n=== 步骤7: 可视化 ===\n")

plot_mr_results <- function(summary_df, figure_dir) {
  cat("  绘制 MR 结果森林图...\n")
  
  if (nrow(summary_df) == 0) {
    cat("    无结果可绘制\n")
    return(NULL)
  }
  
  p <- ggplot(summary_df, aes(x = gene, y = beta, ymin = beta - 1.96 * se, ymax = beta + 1.96 * se)) +
    geom_pointrange(aes(color = direction), size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    coord_flip() +
    labs(
      title = "Mendelian Randomization: Cuproptosis Genes → CIRI",
      x = "Gene",
      y = "Beta (MR-IVW)",
      color = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 10)
    )
  
  ggsave(file.path(figure_dir, "MR_forest_plot.pdf"), plot = p, width = 8, height = 6)
  
  cat("  森林图已保存\n")
}

# ==================== 主流程 ====================
cat("============================================================\n")
cat("L3 遗传因果锚定：TwoSampleMR 双因果认证\n")
cat("============================================================\n")

# 注意：此处为分析框架，实际运行需要：
# 1. 配置 IEU Open GWAS API 访问权限
# 2. 提供 GTEx v8 脑组织 eQTL 数据
# 3. 提供缺血性卒中 GWAS 汇总统计数据

cat("\n分析框架已建立，实际运行需要配置以下数据源：\n")
cat("  ✓ eQTLGen (血液 eQTL, n≈31,000)\n")
cat("  ✓ GTEx v8 (脑组织 eQTL: frontal cortex/hippocampus/cerebellum)\n")
cat("  ✓ IEU Open GWAS API (缺血性卒中 GWAS)\n")
cat("  ✓ GenAsia (东亚人群 eQTL)\n")

cat("\n自检标准:\n")
cat("  ✓ F 统计量 > 10\n")
cat("  ✓ Cochran Q p > 0.05（随机效应 IVW 可放宽）\n")
cat("  ✓ Steiger p < 0.05\n")
cat("  ✓ coloc PP.H4 > 0.75\n")

cat("\nL3 分析框架完成！\n")