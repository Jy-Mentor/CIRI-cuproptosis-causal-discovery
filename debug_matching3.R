#!/usr/bin/env Rscript
# 调试主脚本中的具体问题

library(readxl)
library(data.table)
library(dplyr)

# 完全复制主脚本的逻辑
exposure_file <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.01.xlsx"
exposure_data <- read_excel(exposure_file)

genes <- c("NFKB1", "STAT3", "HIF1A", "HSPA5", "HMOX1",
           "RELA", "NFE2L2", "CP", "LIAS", "IKBKB",
           "JAK1", "PARP1", "CASP8", "MTOR", "PTPRC")

gene_col <- "gene"

# 筛选目标基因（如主脚本第65行）
exposure_snps <- exposure_data[toupper(exposure_data[[gene_col]]) %in% toupper(genes), ]
cat("筛选后SNP数:", nrow(exposure_snps), "\n")

# 重命名列（如主脚本第125-149行）
chr_col <- "chr.exposure"
pos_col <- "pos.exposure"
beta_col <- "beta.exposure"
se_col <- "se.exposure"
ea_col <- "effect_allele.exposure"
oa_col <- "other_allele.exposure"
pval_col <- "pval.exposure"
eaf_col <- "eaf.exposure"
snp_col <- "SNP"

exposure_snps$chr <- as.integer(exposure_snps[[chr_col]])
exposure_snps$pos <- as.numeric(exposure_snps[[pos_col]])
exposure_snps$gene_name <- toupper(as.character(exposure_snps[[gene_col]]))
exposure_snps$beta <- as.numeric(exposure_snps[[beta_col]])
exposure_snps$se <- as.numeric(exposure_snps[[se_col]])
exposure_snps$effect_allele <- exposure_snps[[ea_col]]
exposure_snps$other_allele <- exposure_snps[[oa_col]]
exposure_snps$pval <- as.numeric(exposure_snps[[pval_col]])
exposure_snps$eaf <- as.numeric(exposure_snps[[eaf_col]])
exposure_snps$SNP <- exposure_snps[[snp_col]]

cat("转换后SNP数:", nrow(exposure_snps), "\n")
cat("chr列NA数:", sum(is.na(exposure_snps$chr)), "\n")
cat("pos列NA数:", sum(is.na(exposure_snps$pos)), "\n")

# 移除NA值（如主脚本第152-157行）
valid <- !is.na(exposure_snps$chr) & !is.na(exposure_snps$pos) & 
         !is.na(exposure_snps$beta) & !is.na(exposure_snps$se) &
         !is.na(exposure_snps$effect_allele) & !is.na(exposure_snps$other_allele)
exposure_snps <- exposure_snps[valid, ]
cat("移除NA后SNP数:", nrow(exposure_snps), "\n")

# 创建chr_pos（如主脚本第160行）
exposure_snps$chr_pos <- paste(exposure_snps$chr, exposure_snps$pos, sep = ":")
cat("创建chr_pos后SNP数:", nrow(exposure_snps), "\n")

# 读取结局数据
outcome_file <- "D:/EQTL/mr_results_megastroke/megastroke_outcome.csv"
outcome_raw <- fread(outcome_file)
outcome_raw$chr.outcome <- as.integer(outcome_raw$chr)
outcome_raw$pos.outcome <- as.numeric(outcome_raw$pos.outcome)
outcome_raw$chr_pos <- paste(outcome_raw$chr.outcome, outcome_raw$pos.outcome, sep = ":")

# 测试CASP8（如主脚本第203-221行）
gene <- "CASP8"
cat("\n========== 测试基因:", gene, "==========\n")

exp_gene <- exposure_snps[exposure_snps$gene_name == toupper(gene), ]
cat("exp_gene行数:", nrow(exp_gene), "\n")

if (nrow(exp_gene) > 0) {
  cat("exp_gene$chr:\n")
  print(exp_gene$chr)
  cat("exp_gene$pos:\n")
  print(exp_gene$pos)
  
  exp_chr_pos <- paste(exp_gene$chr, exp_gene$pos, sep = ":")
  cat("exp_chr_pos:\n")
  print(exp_chr_pos)
  cat("唯一exp_chr_pos数:", length(unique(exp_chr_pos)), "\n")
  
  # 检查匹配
  cat("\n结局数据chr_pos前10个:\n")
  print(head(outcome_raw$chr_pos, 10))
  
  matches <- outcome_raw$chr_pos %in% exp_chr_pos
  cat("匹配数:", sum(matches), "\n")
  
  if (sum(matches) > 0) {
    cat("匹配的chr_pos:\n")
    print(outcome_raw$chr_pos[matches])
  }
}
