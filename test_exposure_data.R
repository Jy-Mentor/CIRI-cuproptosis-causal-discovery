#!/usr/bin/env Rscript
# =============================================================================
# 测试暴露数据问题 - 参考 GitHub 权威实践
# =============================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(readxl)
  library(arrow)
})

cat("======================================================================\n")
cat("暴露数据问题诊断测试\n")
cat("======================================================================\n\n")

# 设置 JWT Token
token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3ODMwNDA4MywiZXhwIjoxNzc5NTEzNjgzfQ.ZtcIUEx_xYtrVD_EE-UboKyLlC-lZBq2pjn-iYhJzxocqHdA-02K9n_Qbw-5ngQ07GHjjIYqVtmkZfJ3OJl1yI-tOMBBFzVKe0nkwDcB6-yBgjgBaxVm8vq_pbNrMwy_ZezY5ys9jx7I8T4bZYg9KeUbSwj04OfNP82kGcKXIOErXXVy-Ie3dbUogDRSjnCT-_32yNQxuWpiyYnPWSrWQbQ2HlUQiiDTdFGzWeJJKfSRvjQzdp5g3nccxht0m5A0UsPCdvkyHFEpvPVZ-NpjCjkgy8GbZBv4cmDMSc5JJL6HLO0eV508SRKxMdp-gL6qVdhGiJ2i9XmZQE27-aq_7g"
Sys.setenv(OPENGWAS_JWT = token)

# 验证 Token
cat("1. 验证 OpenGWAS API Token...\n")
tryCatch({
  user_info <- ieugwasr::user()
  cat("   ✓ Token 有效\n")
  cat("   用户:", user_info$user[[1]], "\n")
  cat("   Tier:", user_info$tier[[1]], "\n\n")
}, error = function(e) {
  cat("   ✗ Token 无效:", conditionMessage(e), "\n\n")
})

# =============================================================================
# 测试 1: 检查本地 eQTLGen 数据格式
# =============================================================================
cat("2. 检查本地 eQTLGen 数据格式...\n")
eqtlgen_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-8_kb_10000_r2_0.001.xlsx"

if (file.exists(eqtlgen_file)) {
  eqtlgen_data <- read_excel(eqtlgen_file)
  cat("   文件:", basename(eqtlgen_file), "\n")
  cat("   SNP 数量:", nrow(eqtlgen_data), "\n")
  cat("   列名:", paste(names(eqtlgen_data), collapse = ", "), "\n")
  
  # 检查关键列
  required_cols <- c("SNP", "chr.exposure", "pos.exposure", "beta.exposure", "se.exposure", 
                     "effect_allele.exposure", "other_allele.exposure", 
                     "eaf.exposure", "pval.exposure")
  missing_cols <- setdiff(required_cols, names(eqtlgen_data))
  
  if (length(missing_cols) > 0) {
    cat("   ✗ 缺少列:", paste(missing_cols, collapse = ", "), "\n")
    cat("   提示：需要将 chr.exposure 重命名为 CHR, pos.exposure 重命名为 BP\n")
  } else {
    cat("   ✓ 所有必需列都存在\n")
  }
  
  # 显示前几行
  cat("\n   前 3 行数据:\n")
  print(head(eqtlgen_data[, 1:min(10, ncol(eqtlgen_data))], 3))
  cat("\n")
  
  # 检查 SNP ID 格式
  cat("   SNP ID 示例:\n")
  print(eqtlgen_data$SNP[1:5])
  cat("\n")
  
  # 检查基因分布
  if ("gene" %in% names(eqtlgen_data)) {
    cat("   基因分布:\n")
    gene_counts <- table(eqtlgen_data$gene)
    print(head(sort(gene_counts, decreasing = TRUE), 10))
    cat("\n")
  }
  
} else {
  cat("   ✗ 文件不存在:", eqtlgen_file, "\n\n")
}

# =============================================================================
# 测试 2: 检查 GTEx 本地数据格式
# =============================================================================
cat("3. 检查 GTEx 本地数据格式...\n")

gtex_egenes_file <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Brain_Cortex.v11.eGenes.txt"
gtex_pairs_file <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Brain_Cortex.v11.eQTLs.signif_pairs.parquet"

if (file.exists(gtex_egenes_file)) {
  cat("   读取 eGenes 文件...\n")
  gtex_egenes <- fread(gtex_egenes_file, stringsAsFactors = FALSE)
  cat("   ✓ eGenes 数量:", nrow(gtex_egenes), "\n")
  cat("   列名:", paste(names(gtex_egenes), collapse = ", "), "\n")
  
  # 显示前几行
  cat("\n   eGenes 前 3 行:\n")
  print(head(gtex_egenes[, 1:min(10, ncol(gtex_egenes))], 3))
  cat("\n")
  
  # 检查 gene_id 格式
  cat("   gene_id 示例:\n")
  print(head(gtex_egenes$gene_id, 5))
  cat("\n")
  
} else {
  cat("   ✗ eGenes 文件不存在\n\n")
}

if (file.exists(gtex_pairs_file)) {
  cat("   读取 parquet 文件...\n")
  gtex_pairs <- read_parquet(gtex_pairs_file)
  cat("   ✓ eQTL 对数量:", nrow(gtex_pairs), "\n")
  cat("   列名:", paste(names(gtex_pairs), collapse = ", "), "\n")
  
  # 显示前几行
  cat("\n   eQTL 对前 3 行:\n")
  print(head(gtex_pairs[, 1:min(10, ncol(gtex_pairs))], 3))
  cat("\n")
  
  # 检查 variant_id 格式
  cat("   variant_id 示例:\n")
  print(head(gtex_pairs$variant_id, 5))
  cat("\n")
  
  # 检查 phenotype_id 与 eGenes 的 gene_id 是否匹配
  if (file.exists(gtex_egenes_file)) {
    unique_phenotypes <- unique(gtex_pairs$phenotype_id)
    unique_genes <- unique(gtex_egenes$gene_id)
    
    cat("   匹配检查:\n")
    cat("   - eGenes 中 unique gene_id 数量:", length(unique_genes), "\n")
    cat("   - parquet 中 unique phenotype_id 数量:", length(unique_phenotypes), "\n")
    
    # 检查是否有重叠
    common_ids <- intersect(unique_phenotypes, unique_genes)
    cat("   - 共同的 ID 数量:", length(common_ids), "\n")
    
    if (length(common_ids) == 0) {
      cat("   ✗ 警告：eGenes 和 parquet 文件的 ID 完全不匹配！\n")
      cat("      这可能是 GTEx 数据无法匹配的关键原因\n\n")
    } else {
      cat("   ✓ ID 匹配正常\n\n")
    }
  }
  
} else {
  cat("   ✗ parquet 文件不存在\n\n")
}

# =============================================================================
# 测试 3: 从 IEU OpenGWAS 获取 GTEx 数据（权威方法）
# =============================================================================
cat("4. 测试 IEU OpenGWAS API 获取 GTEx 数据（参考 GitHub 最佳实践）...\n")

# GTEx v8 脑组织数据集
gtex_datasets <- list(
  "Brain_Frontal_Cortex" = "ieu-b-4171",
  "Brain_Cortex" = "ieu-b-4178",
  "Whole_Blood" = "ieu-b-4180"
)

test_gene <- "FDX1"

for (tissue_name in names(gtex_datasets)) {
  dataset_id <- gtex_datasets[[tissue_name]]
  cat("\n   测试", tissue_name, "(", dataset_id, ")...\n")
  
  tryCatch({
    # 提取该基因的所有 eQTL
    eqtl_dat <- extract_instruments(
      outcomes = dataset_id,
      p1 = 5e-8,
      clump = TRUE,
      r2 = 0.001,
      kb = 1000
    )
    
    if (!is.null(eqtl_dat) && nrow(eqtl_dat) > 0) {
      cat("   ✓ 成功获取", nrow(eqtl_dat), "个 SNP\n")
      
      # 筛选目标基因
      if ("exposure" %in% names(eqtl_dat)) {
        eqtl_dat$gene <- gsub(" .*", "", eqtl_dat$exposure)
        gene_snps <- eqtl_dat[grepl(paste0("^", test_gene, "$"), eqtl_dat$gene, ignore.case = TRUE), ]
        
        if (nrow(gene_snps) > 0) {
          cat("   ✓", test_gene, "有", nrow(gene_snps), "个 SNP\n")
          cat("   SNP 示例:", gene_snps$SNP[1], "\n")
          cat("   染色体:", gene_snps$CHR[1], "\n")
          cat("   位置:", gene_snps$BP[1], "\n")
        } else {
          cat("   ✗", test_gene, "无数据\n")
        }
      }
    } else {
      cat("   ✗ 无数据返回\n")
    }
    
  }, error = function(e) {
    cat("   ✗ 错误:", conditionMessage(e), "\n")
  })
}

cat("\n")

# =============================================================================
# 测试 4: 检查 MEGASTROKE 结局数据格式
# =============================================================================
cat("5. 检查 MEGASTROKE 结局数据格式...\n")

megastroke_file <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"

if (file.exists(megastroke_file)) {
  outcome_data <- fread(megastroke_file, sep = "\t", stringsAsFactors = FALSE)
  cat("   ✓ 成功读取", nrow(outcome_data), "个 SNP\n")
  cat("   列名:", paste(names(outcome_data), collapse = ", "), "\n")
  
  # 显示前几行
  cat("\n   前 3 行:\n")
  print(head(outcome_data[, 1:min(10, ncol(outcome_data))], 3))
  cat("\n")
  
  # 检查染色体分布
  if ("chromosome" %in% names(outcome_data)) {
    cat("   染色体分布:\n")
    chr_counts <- table(outcome_data$chromosome)
    print(chr_counts)
    cat("\n")
  }
  
  # 检查 SNP ID 格式
  if ("variant_id" %in% names(outcome_data)) {
    cat("   variant_id 示例:\n")
    print(head(outcome_data$variant_id, 10))
    cat("\n")
  }
  
} else {
  cat("   ✗ 文件不存在\n\n")
}

# =============================================================================
# 总结
# =============================================================================
cat("======================================================================\n")
cat("诊断总结\n")
cat("======================================================================\n\n")

cat("关键问题排查:\n")
cat("1. eQTLGen 数据格式是否正确？\n")
cat("2. GTEx eGenes 和 parquet 文件的 ID 是否匹配？\n")
cat("3. IEU OpenGWAS API 是否能获取 GTEx 数据？\n")
cat("4. MEGASTROKE 数据的染色体位置是否与 eQTL 重叠？\n\n")

cat("建议的解决方案:\n")
cat("- 如果 GTEx 本地数据 ID 不匹配，改用 IEU OpenGWAS API 获取\n")
cat("- 如果 SNP ID 格式不匹配，使用染色体位置匹配而非 rsID\n")
cat("- 如果 MEGASTROKE 区域不重叠，检查坐标系统是否一致 (GRCh37 vs GRCh38)\n\n")

cat("测试完成！\n")
