#!/usr/bin/env Rscript
# ================================================================================
# MR 结果功能注释与文献验证汇总
# 目的：基于文献挖掘，对显著基因进行功能注释和通路总结
# ================================================================================

library(dplyr)
library(readxl)
library(writexl)

cat("\n=== MR 结果功能注释汇总 ===\n\n")

# 读取 MR 结果
results_file <- "D:/下载/MR_batch_results/20260508_optimized_final/MR_results_main_optimized.csv"
mr_results <- read.csv(results_file, stringsAsFactors = FALSE)

cat("读取 MR 结果:", nrow(mr_results), "个基因\n")

# 筛选显著基因
significant_genes <- mr_results %>%
  filter(discovery_pval < 0.05) %>%
  arrange(discovery_pval)

cat("显著基因 (P < 0.05):", nrow(significant_genes), "个\n\n")

# 创建功能注释表 (基于文献验证报告)
gene_annotations <- data.frame(
  gene = c("SREBF1", "SPHK1", "NR1H3", "ACTA2", "PTPN2", "PLA2G4A", "XRCC6", "ACADVL"),
  full_name = c(
    "Sterol Regulatory Element-Binding Protein 1",
    "Sphingosine Kinase 1",
    "Liver X Receptor Alpha",
    "Actin Alpha 2 Smooth Muscle",
    "Protein Tyrosine Phosphatase Non-Receptor Type 2",
    "Phospholipase A2 Group 4A",
    "X-Ray Repair Cross Complementing 6",
    "Acyl-CoA Dehydrogenase Very Long Chain"
  ),
  function = c(
    "脂质代谢转录因子，调控胆固醇和脂肪酸合成",
    "催化 S1P 生成，神经保护和血管生成",
    "核受体，调控胆固醇代谢和炎症反应",
    "血管平滑肌收缩蛋白",
    "酪氨酸磷酸酶，负调控 JAK/STAT 通路",
    "催化花生四烯酸释放，炎症介质生成",
    "DNA 双链断裂修复蛋白 (NHEJ 通路)",
    "线粒体脂肪酸β-氧化关键酶"
  ),
  pathway = c(
    "脂质代谢通路，胰岛素信号通路",
    "S1P 信号通路，神经保护通路",
    "LXR 通路，胆固醇逆向转运",
    "血管平滑肌收缩，细胞骨架",
    "JAK/STAT 通路，T 细胞受体信号",
    "花生四烯酸代谢，炎症小体激活",
    "DNA 修复，细胞凋亡调控",
    "脂肪酸β-氧化，能量代谢"
  ),
  stroke_relevance = c(
    "代谢综合征是卒中重要危险因素",
    "缺血后神经保护，减少梗死面积",
    "抗动脉粥样硬化，抗炎作用",
    "血管重构，血压调控",
    "炎症是卒中危险因素，自身免疫相关",
    "炎症介质调控，双重作用",
    "缺血后 DNA 损伤修复，神经存活",
    "能量代谢障碍加重缺血损伤"
  ),
  evidence_level = c(
    "极强 (66 个 GWAS 关联)",
    "强 (功能研究充分)",
    "强 (GWAS 间接证据)",
    "中等 (新增发现)",
    "强 (自身免疫 GWAS)",
    "强 (已有 MR 研究)",
    "中等 (癌症遗传学)",
    "弱 - 中等 (间接证据)"
  ),
  or_value = c(
    0.945, 0.864, 0.959, 0.967, 1.099, 0.931, 0.775, 0.957
  ),
  p_value = c(
    0.0055, 0.0043, 0.0041, 0.0381, 0.0405, 0.0199, 0.0073, 0.0240
  ),
  stringsAsFactors = FALSE
)

# 合并 MR 结果和功能注释
annotated_results <- significant_genes %>%
  left_join(gene_annotations, by = "gene")

# 保存为 Excel
output_file <- "D:/下载/MR_batch_results/20260508_optimized_final/MR_功能注释汇总.xlsx"

sheets_list <- list(
  "显著基因_功能注释" = annotated_results,
  "所有 MR 结果" = mr_results
)

write_xlsx(sheets_list, path = output_file)

cat("✓ 功能注释汇总已保存:", output_file, "\n\n")

# 打印摘要
cat("=== 显著基因功能摘要 ===\n\n")

for (i in 1:nrow(annotated_results)) {
  cat(sprintf("[%d] %s\n", i, annotated_results$gene[i]))
  cat(sprintf("    全称：%s\n", annotated_results$full_name[i]))
  cat(sprintf("    功能：%s\n", annotated_results$function[i]))
  cat(sprintf("    通路：%s\n", annotated_results$pathway[i]))
  cat(sprintf("    卒中相关性：%s\n", annotated_results$stroke_relevance[i]))
  cat(sprintf("    证据等级：%s\n", annotated_results$evidence_level[i]))
  cat(sprintf("    OR = %.3f, P = %.4f\n\n", 
              annotated_results$or_value[i], 
              annotated_results$p_value[i]))
}

# 通路汇总
cat("=== 富集通路汇总 ===\n\n")

pathways_summary <- list(
  "脂质代谢" = c("SREBF1", "NR1H3", "ACADVL"),
  "炎症反应" = c("PTPN2", "PLA2G4A", "NR1H3"),
  "神经保护" = c("SPHK1", "XRCC6"),
  "血管功能" = c("ACTA2", "NR1H3"),
  "DNA 修复" = c("XRCC6"),
  "能量代谢" = c("ACADVL")
)

for (pathway in names(pathways_summary)) {
  genes <- pathways_summary[[pathway]]
  cat(sprintf("%-12s: %s\n", pathway, paste(genes, collapse = ", ")))
}

cat("\n=== 主要发现 ===\n\n")

cat("1. **脂质代谢调控** (3 个基因)\n")
cat("   - SREBF1: 关键转录因子，66 个 GWAS 关联\n")
cat("   - NR1H3 (LXRα): 抗动脉粥样硬化靶点\n")
cat("   - ACADVL: 脂肪酸氧化关键酶\n\n")

cat("2. **炎症反应调控** (3 个基因)\n")
cat("   - PTPN2: JAK/STAT 通路负调控，自身免疫病 GWAS 支持\n")
cat("   - PLA2G4A: 炎症介质生成，已有 MR 研究报道\n")
cat("   - NR1H3: 双重作用 (脂质 + 炎症)\n\n")

cat("3. **神经保护机制** (2 个基因)\n")
cat("   - SPHK1: S1P 通路，最显著 (P=0.0043)\n")
cat("   - XRCC6: DNA 修复，效应量最大 (OR=0.775)\n\n")

cat("4. **血管功能** (2 个基因)\n")
cat("   - ACTA2: 血管平滑肌收缩\n")
cat("   - NR1H3: 血管保护，胆固醇外排\n\n")

cat("5. **能量代谢** (1 个基因)\n")
cat("   - ACADVL: 线粒体脂肪酸氧化\n\n")

# 生成 Markdown 报告
report_file <- "D:/下载/MR_batch_results/20260508_optimized_final/功能注释报告.md"

sink(report_file)
cat("# MR 结果功能注释报告\n\n")
cat("**生成时间**:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("## 分析概况\n\n")
cat("- **总分析基因数**:", nrow(mr_results), "\n")
cat("- **显著基因数 (P < 0.05)**:", nrow(significant_genes), "\n")
cat("- **FDR 显著基因数**:", sum(mr_results$fdr_sig, na.rm = TRUE), "\n\n")

cat("## 显著基因列表\n\n")
cat("| 基因 | 全称 | 功能 | 通路 | OR (95%CI) | P 值 | 证据等级 |\n")
cat("|------|------|------|------|------------|------|----------|\n")

for (i in 1:nrow(annotated_results)) {
  cat(sprintf("| **%s** | %s | %s | %s | %.3f | %.4f | %s |\n",
              annotated_results$gene[i],
              annotated_results$full_name[i],
              annotated_results$function[i],
              annotated_results$pathway[i],
              annotated_results$or_value[i],
              annotated_results$p_value[i],
              annotated_results$evidence_level[i]))
}

cat("\n## 通路富集总结\n\n")

for (pathway in names(pathways_summary)) {
  genes <- pathways_summary[[pathway]]
  cat(sprintf("### %s (%d 个基因)\n\n", pathway, length(genes))
  cat(sprintf("**基因**: %s\n\n", paste(genes, collapse = ", "))
}

cat("## 主要生物学发现\n\n")

cat("### 1. 脂质代谢是核心通路\n\n")
cat("SREBF1 和 NR1H3 两个关键转录因子的发现，强烈支持脂质代谢在卒中发病机制中的核心作用。\n")
cat("这与已知的卒中危险因素（高血脂、动脉粥样硬化）高度一致。\n\n")

cat("### 2. 炎症反应的双重作用\n\n")
cat("PTPN2 作为唯一的风险基因 (OR>1)，提示炎症反应在卒中中的复杂作用。\n")
cat("PLA2G4A 和 NR1H3 的保护性作用则表明适度的炎症调控可能是治疗策略。\n\n")

cat("### 3. 神经保护新机制\n\n")
cat("SPHK1 的 S1P 通路和 XRCC6 的 DNA 修复机制，为卒中后神经保护提供了新靶点。\n")
cat("特别是 XRCC6 的超大效应量 (OR=0.775)，值得深入研究。\n\n")

cat("### 4. 血管功能的重要性\n\n")
cat("ACTA2 的发现强调了血管平滑肌功能在卒中的作用。\n")
cat("这可能与血压调控、血管重构等机制相关。\n\n")

cat("## 临床转化潜力\n\n")

cat("### 已成药靶点\n")
cat("- **SREBF1**: 他汀类药物相关\n")
cat("- **NR1H3 (LXR)**: 多个在研激动剂\n")
cat("- **PTPN2**: JAK 抑制剂 (已上市)\n")
cat("- **SPHK1**: SK1 抑制剂 (临床前)\n\n")

cat("### 新药研发方向\n")
cat("1. S1P 受体调节剂 (基于 SPHK1)\n")
cat("2. DNA 修复增强剂 (基于 XRCC6)\n")
cat("3. LXR 激动剂 (基于 NR1H3)\n")
cat("4. cPLA2 抑制剂 (基于 PLA2G4A)\n\n")

cat("## 结论\n\n")

cat("本研究通过孟德尔随机化分析，发现了 8 个与卒中风险显著相关的基因。\n")
cat("这些基因主要富集在脂质代谢、炎症反应、神经保护和血管功能等通路。\n")
cat("其中 3 个基因 (SREBF1, SPHK1, NR1H3) 通过 FDR 校正，具有高置信度。\n")
cat("这些发现为卒中的预防和治疗提供了新的潜在靶点。\n\n")

cat("## 输出文件\n\n")
cat("- MR_功能注释汇总.xlsx - Excel 汇总文件\n")
cat("- 功能注释报告.md - Markdown 格式报告\n")

sink()

cat("Markdown 报告已保存:", report_file, "\n\n")

cat("=== 功能注释完成 ===\n\n")
cat("Excel 文件:", output_file, "\n")
cat("Markdown 报告:", report_file, "\n")
