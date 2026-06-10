#!/usr/bin/env Rscript
# 测试本地GTEx v11数据读取

set.seed(42)

# 加载必要的包
suppressPackageStartupMessages({
  if (!require("TwoSampleMR", quietly = TRUE)) install.packages("TwoSampleMR")
  library(TwoSampleMR)
})

# 测试函数：从本地GTEx v11数据获取基因eQTL
test_get_gtex_eqtl <- function(gene_symbol) {
  cat("Testing", gene_symbol, "...\n")
  
  tryCatch({
    # 本地GTEx v11数据文件路径
    egenes_file <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Brain_Cortex.v11.eGenes.txt"
    
    if (!file.exists(egenes_file)) {
      cat("ERROR: GTEx eGenes file not found\n")
      return(NULL)
    }
    
    # 读取eGenes文件
    cat("Reading eGenes file...\n")
    eqtl_data <- read.table(egenes_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # 按基因名称筛选
    gene_data <- eqtl_data[eqtl_data$gene_name == gene_symbol, ]
    
    if (nrow(gene_data) == 0) {
      cat("WARNING: No eQTLs found for", gene_symbol, "\n")
      return(NULL)
    }
    
    cat("SUCCESS: Found", nrow(gene_data), "eQTLs for", gene_symbol, "\n")
    
    # 显示前5行数据
    cat("First 5 eQTLs:\n")
    print(head(gene_data[, c("rs_id_dbSNP157_GRCh38p14", "slope", "slope_se", "pval_nominal", "af", "alt", "ref")], 5))
    
    # 转换为TwoSampleMR格式
    exp_dat <- data.frame(
      SNP = gene_data$rs_id_dbSNP157_GRCh38p14,
      beta.exposure = gene_data$slope,
      se.exposure = gene_data$slope_se,
      pval.exposure = gene_data$pval_nominal,
      eaf.exposure = gene_data$af,
      effect_allele.exposure = gene_data$alt,
      other_allele.exposure = gene_data$ref,
      exposure = paste0(gene_symbol, " (Brain Cortex GTEx v11)"),
      stringsAsFactors = FALSE
    )
    
    # 过滤掉没有rsID的SNP
    exp_dat <- exp_dat[exp_dat$SNP != "", ]
    
    if (nrow(exp_dat) == 0) {
      cat("WARNING: No SNPs with rsID found\n")
      return(NULL)
    }
    
    cat("After filtering SNPs with rsID:", nrow(exp_dat), "SNPs\n")
    
    # 格式转换
    exp_fmt <- format_data(
      exp_dat,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "beta.exposure",
      se_col = "se.exposure",
      pval_col = "pval.exposure",
      eaf_col = "eaf.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure"
    )
    
    cat("SUCCESS: Data formatted for TwoSampleMR\n")
    return(exp_fmt)
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# 测试几个基因
cat("==================================================\n")
cat("Testing GTEx v11 Local Data Reader\n")
cat("==================================================\n\n")

test_genes <- c("FDX1", "NFKB1", "STAT3")

for (gene in test_genes) {
  cat("\n", rep("=", 50), "\n", sep = "")
  result <- test_get_gtex_eqtl(gene)
  if (!is.null(result)) {
    cat("\nFormatted data sample:\n")
    print(head(result, 3))
  }
  cat(rep("=", 50), "\n", sep = "")
}

cat("\n==================================================\n")
cat("Test Complete!\n")
cat("==================================================\n")
