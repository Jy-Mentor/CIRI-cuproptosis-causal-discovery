#!/usr/bin/env Rscript
# 调试脚本：检查匹配问题

library(readxl)
library(data.table)

# 读取暴露数据
exposure_file <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.01.xlsx"
exposure_data <- read_excel(exposure_file)

genes <- c("NFKB1", "STAT3", "HIF1A", "HSPA5", "HMOX1",
           "RELA", "NFE2L2", "CP", "LIAS", "IKBKB",
           "JAK1", "PARP1", "CASP8", "MTOR", "PTPRC")

# 筛选目标基因
exposure_snps <- exposure_data[toupper(exposure_data$gene) %in% toupper(genes), ]
cat("目标基因SNP数:", nrow(exposure_snps), "\n")

# 创建chr_pos
exposure_snps$chr_pos <- paste(exposure_snps$chr.exposure, exposure_snps$pos.exposure, sep = ":")
cat("暴露数据chr_pos示例:\n")
print(head(exposure_snps$chr_pos))

# 读取结局数据
outcome_file <- "D:/EQTL/mr_results_megastroke/megastroke_outcome.csv"
outcome_raw <- fread(outcome_file)

cat("\n结局数据列名:\n")
print(names(outcome_raw))

# 检查结局数据的chr和pos列
cat("\n结局数据chr列名检查:\n")
print("chr" %in% names(outcome_raw))
print("chr.outcome" %in% names(outcome_raw))

# 创建chr_pos
cat("\n创建结局数据chr_pos...\n")
outcome_raw$chr_pos <- paste(outcome_raw$chr, outcome_raw$pos.outcome, sep = ":")
cat("结局数据chr_pos示例:\n")
print(head(outcome_raw$chr_pos))

# 检查匹配
cat("\n=== 匹配检查 ===\n")
exp_chr_pos <- unique(exposure_snps$chr_pos)
out_chr_pos <- unique(outcome_raw$chr_pos)

common <- intersect(exp_chr_pos, out_chr_pos)
cat("暴露唯一chr_pos数:", length(exp_chr_pos), "\n")
cat("结局唯一chr_pos数:", length(out_chr_pos), "\n")
cat("共同chr_pos数:", length(common), "\n")

if (length(common) > 0) {
  cat("\n匹配的chr_pos示例:\n")
  print(head(common))
}
