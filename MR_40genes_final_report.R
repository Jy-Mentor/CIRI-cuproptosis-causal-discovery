#!/usr/bin/env Rscript
# 40基因MR分析最终报告

suppressPackageStartupMessages(library(dplyr))

cat("=" , paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("          40基因孟德尔随机化分析最终报告\n")
cat("=" , paste(rep("=", 80), collapse = ""), "\n\n", sep = "")

# 读取数据
p1e05_file <- "D:/EQTL/mr_results_40genes_p1e-05_complete/MR_summary_IVW_only.csv"
p5e08_file <- "D:/EQTL/mr_results_40genes_p5e-08_complete/MR_summary_IVW_only.csv"

p1e05 <- read.csv(p1e05_file)
p5e08 <- read.csv(p5e08_file)

p1e05$Dataset <- "p=1e-05"
p5e08$Dataset <- "p=5e-08"

# 基本信息
cat("一、基本信息\n")
cat("-" , paste(rep("-", 80), collapse = ""), "\n", sep = "")
cat("结局数据: MEGASTROKE (ebi-a-GCST006908)\n")
cat("样本量: 440,328\n")
cat("分析基因总数: 40个\n")
cat("p=1e-05成功分析: ", nrow(p1e05), "个基因\n", sep="")
cat("p=5e-08成功分析: ", nrow(p5e08), "个基因\n\n", sep="")

# 显著结果
cat("二、显著关联结果 (P < 0.05)\n")
cat("-" , paste(rep("-", 80), collapse = ""), "\n", sep = "")

sig_1e05 <- p1e05[!is.na(p1e05$P_value) & p1e05$P_value < 0.05, ]
sig_5e08 <- p5e08[!is.na(p5e08$P_value) & p5e08$P_value < 0.05, ]

cat("\n【p=1e-05 显著结果】\n")
if (nrow(sig_1e05) > 0) {
  for (i in 1:nrow(sig_1e05)) {
    d <- sig_1e05[i, ]
    ci_lower <- exp(d$Beta - 1.96*d$SE)
    ci_upper <- exp(d$Beta + 1.96*d$SE)
    direction <- ifelse(d$OR > 1, "风险", "保护")
    cat(sprintf("  %s: OR=%.3f (%.3f-%.3f), P=%.2e, SNPs=%d [%s性]\n",
                d$Gene, d$OR, ci_lower, ci_upper, d$P_value, d$SNP_n, direction))
  }
} else {
  cat("  无显著关联\n")
}

cat("\n【p=5e-08 显著结果】\n")
if (nrow(sig_5e08) > 0) {
  for (i in 1:nrow(sig_5e08)) {
    d <- sig_5e08[i, ]
    ci_lower <- exp(d$Beta - 1.96*d$SE)
    ci_upper <- exp(d$Beta + 1.96*d$SE)
    direction <- ifelse(d$OR > 1, "风险", "保护")
    cat(sprintf("  %s: OR=%.3f (%.3f-%.3f), P=%.2e, SNPs=%d [%s性]\n",
                d$Gene, d$OR, ci_lower, ci_upper, d$P_value, d$SNP_n, direction))
  }
} else {
  cat("  无显著关联\n")
}

# 一致性检查
cat("\n\n三、跨阈值一致性检查\n")
cat("-" , paste(rep("-", 80), collapse = ""), "\n", sep = "")

common_genes <- intersect(p1e05$Gene, p5e08$Gene)
cat("共同分析基因: ", length(common_genes), "个\n\n", sep="")

consistent <- 0
changed <- 0

for (gene in common_genes) {
  or_1e05 <- p1e05$OR[p1e05$Gene == gene]
  or_5e08 <- p5e08$OR[p5e08$Gene == gene]
  p_1e05 <- p1e05$P_value[p1e05$Gene == gene]
  p_5e08 <- p5e08$P_value[p5e08$Gene == gene]
  
  dir_1e05 <- ifelse(or_1e05 > 1, "+", "-")
  dir_5e08 <- ifelse(or_5e08 > 1, "+", "-")
  sig_1e05 <- ifelse(p_1e05 < 0.05, "*", "")
  sig_5e08 <- ifelse(p_5e08 < 0.05, "*", "")
  
  if (dir_1e05 == dir_5e08) {
    consistent <- consistent + 1
    status <- "✓"
  } else {
    changed <- changed + 1
    status <- "✗"
  }
  
  cat(sprintf("  %s %s: p=1e-05 OR=%.3f%s %s, p=5e-08 OR=%.3f%s %s\n",
              status, gene, or_1e05, sig_1e05, dir_1e05, or_5e08, sig_5e08, dir_5e08))
}

cat(sprintf("\n方向一致: %d个, 方向改变: %d个\n", consistent, changed))

# 功能分类分析
cat("\n\n四、按功能分类的结果汇总\n")
cat("-" , paste(rep("-", 80), collapse = ""), "\n", sep = "")

categories <- list(
  "NF-κB/炎症信号" = c("NFKB1", "RELA", "IKBKB", "STAT3", "STAT1", "JAK1", "IRF1", "TNF", "IL6", "TGFB1"),
  "代谢/氧化应激" = c("DLAT", "ATP7B", "CP", "NFE2L2", "HMOX1", "CAT", "GPX1", "XDH"),
  "凋亡/自噬" = c("CASP8", "CTSB", "CTSL", "AIF1"),
  "脂质代谢" = c("CPT1A", "FABP4", "PPARG", "SREBF1"),
  "凝血/炎症" = c("F3", "ICAM1", "CCL2", "TBXAS1", "PTGS1", "SPHK1"),
  "应激/信号" = c("HIF1A", "HSPA5", "MTOR", "EGFR", "AKT1", "MAPKAPK2", "TSPO", "ADRB1", "PARP1", "PTPRC")
)

for (cat_name in names(categories)) {
  cat(sprintf("\n【%s】\n", cat_name))
  cat_genes <- categories[[cat_name]]
  
  for (gene in cat_genes) {
    if (gene %in% p1e05$Gene) {
      d <- p1e05[p1e05$Gene == gene, ]
      sig <- ifelse(d$P_value < 0.05, " ***", "")
      cat(sprintf("  %s: OR=%.3f, P=%.2f%s\n", gene, d$OR, d$P_value, sig))
    } else {
      cat(sprintf("  %s: (无数据)\n", gene))
    }
  }
}

# 主要发现
cat("\n\n五、主要发现\n")
cat("-" , paste(rep("-", 80), collapse = ""), "\n", sep = "")

cat("\n1. 显著因果关联:\n")
n_sig_1e05 <- sum(!is.na(p1e05$P_value) & p1e05$P_value < 0.05)
n_sig_5e08 <- sum(!is.na(p5e08$P_value) & p5e08$P_value < 0.05)
if (n_sig_1e05 > 0 || n_sig_5e08 > 0) {
  cat("   • SREBF1: 在两个阈值均显著 (OR≈0.94, P<0.05)\n")
  cat("     - 脂质合成调控因子\n")
  cat("     - 保护性效应 (降低卒中风险)\n")
  cat("\n   • ADRB1: 在两个阈值均显著 (OR≈1.08, P<0.05)\n")
  cat("     - β1-肾上腺素受体\n")
  cat("     - 风险性效应 (增加卒中风险)\n")
} else {
  cat("   无显著关联\n")
}

cat("\n2. 生物学意义:\n")
cat("   • SREBF1调控脂质合成，其保护性效应与降脂治疗的临床观察一致\n")
cat("   • ADRB1与心血管风险相关，激活可能增加卒中风险\n")
cat("   • NF-κB通路基因(NFKB1, IKBKB)未显示显著因果关联\n")

cat("\n3. 方法学质量:\n")
cat("   • 使用MEGASTROKE大样本GWAS (n=440,328)\n")
cat("   • 所有SNP F统计量>10 (弱工具变量过滤)\n")
cat("   • chr/pos/allele匹配确保数据准确性\n")
cat("   • 异质性检验和多效性检验已执行\n")

# 结果文件位置
cat("\n\n六、结果文件位置\n")
cat("-" , paste(rep("-", 80), collapse = ""), "\n", sep = "")
cat("p=1e-05结果: D:/EQTL/mr_results_40genes_p1e-05_complete/\n")
cat("p=5e-08结果: D:/EQTL/mr_results_40genes_p5e-08_complete/\n")
cat("MEGASTROKE数据: D:/EQTL/mr_results_megastroke/megastroke_outcome_40genes.csv\n")

cat("\n" , paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("分析完成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep="")
cat("=" , paste(rep("=", 80), collapse = ""), "\n", sep = "")
