#!/usr/bin/env Rscript
# 详细调试匹配问题

library(readxl)
library(data.table)
library(dplyr)

# 读取数据
exposure_file <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.01.xlsx"
exposure_data <- read_excel(exposure_file)

genes <- c("NFKB1", "STAT3", "HIF1A", "HSPA5", "HMOX1",
           "RELA", "NFE2L2", "CP", "LIAS", "IKBKB",
           "JAK1", "PARP1", "CASP8", "MTOR", "PTPRC")

# 筛选目标基因
exposure_snps <- exposure_data[toupper(exposure_data$gene) %in% toupper(genes), ]
exposure_snps$chr <- as.integer(exposure_snps$chr.exposure)
exposure_snps$pos <- as.numeric(exposure_snps$pos.exposure)
exposure_snps$gene_name <- toupper(exposure_snps$gene)
exposure_snps$chr_pos <- paste(exposure_snps$chr, exposure_snps$pos, sep = ":")

# 读取结局数据
outcome_file <- "D:/EQTL/mr_results_megastroke/megastroke_outcome.csv"
outcome_raw <- fread(outcome_file)
outcome_raw$chr.outcome <- as.integer(outcome_raw$chr)
outcome_raw$pos.outcome <- as.numeric(outcome_raw$pos.outcome)
outcome_raw$chr_pos <- paste(outcome_raw$chr.outcome, outcome_raw$pos.outcome, sep = ":")

cat("=== 整体匹配检查 ===\n")
cat("暴露chr_pos数:", length(unique(exposure_snps$chr_pos)), "\n")
cat("结局chr_pos数:", length(unique(outcome_raw$chr_pos)), "\n")
cat("共同chr_pos数:", length(intersect(unique(exposure_snps$chr_pos), unique(outcome_raw$chr_pos))), "\n\n")

# 测试单个基因
gene <- "CASP8"
exp_gene <- exposure_snps[exposure_snps$gene_name == gene, ]

cat("=== 基因:", gene, "===\n")
cat("暴露SNP数:", nrow(exp_gene), "\n")
cat("暴露chr:pos:\n")
print(exp_gene$chr_pos)

cat("\n暴露chr类型:", class(exp_gene$chr), "\n")
cat("暴露pos类型:", class(exp_gene$pos), "\n")
cat("暴露chr_pos类型:", class(exp_gene$chr_pos), "\n")

exp_chr_pos <- exp_gene$chr_pos
cat("\nexp_chr_pos内容:\n")
print(exp_chr_pos)

cat("\n结局chr_pos前10个:\n")
print(head(outcome_raw$chr_pos, 10))

cat("\n结局chr_pos类型:", class(outcome_raw$chr_pos), "\n")

# 检查是否有任何匹配
matches <- outcome_raw$chr_pos %in% exp_chr_pos
cat("\n匹配结果:", sum(matches), "个匹配\n")

if (sum(matches) > 0) {
  cat("匹配的chr_pos:\n")
  print(outcome_raw$chr_pos[matches])
} else {
  # 检查具体问题
  cat("\n检查第一个暴露chr_pos是否在结局中:\n")
  test_pos <- exp_chr_pos[1]
  cat("测试位置:", test_pos, "\n")
  cat("是否在结局中:", test_pos %in% outcome_raw$chr_pos, "\n")
  
  # 手动检查
  parts <- strsplit(test_pos, ":")[[1]]
  test_chr <- as.integer(parts[1])
  test_pos_num <- as.numeric(parts[2])
  
  cat("分解: chr=", test_chr, ", pos=", test_pos_num, "\n")
  cat("结局中是否有chr=", test_chr, ":", any(outcome_raw$chr.outcome == test_chr), "\n")
  cat("结局中是否有pos=", test_pos_num, ":", any(outcome_raw$pos.outcome == test_pos_num), "\n")
  
  # 找到这个位置的结局数据
  matching_rows <- outcome_raw$chr.outcome == test_chr & outcome_raw$pos.outcome == test_pos_num
  cat("同时匹配chr和pos的行数:", sum(matching_rows), "\n")
}
