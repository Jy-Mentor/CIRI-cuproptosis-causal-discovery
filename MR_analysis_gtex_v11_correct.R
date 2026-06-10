#!/usr/bin/env Rscript
# MR分析脚本 - GTEx v11 Whole Blood数据集（修正版）
# 使用chr/pos/allele匹配，不依赖rsID
# 直接分块读取结局数据避免内存问题

library(readxl)
library(data.table)
library(readr)
library(TwoSampleMR)
library(dplyr)
library(purrr)

# ============================================
# 参数设置
# ============================================
signif_pairs_file <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eQTLs.signif_pairs.parquet"
outcome_file  <- "D:/EQTL/eqtlgen_ieu_outcome.csv"
genes         <- c("NFKB1", "FDX1", "STAT3", "HIF1A", "HMOX1",
                   "GPX4", "HSPA5", "AGER", "DLAT")
output_dir    <- "D:/EQTL/mr_results_gtex_v11"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== GTEx v11 Whole Blood MR 分析（chr/pos/allele匹配版） ==========\n")
cat("signif_pairs文件:", signif_pairs_file, "\n")
cat("结局文件:", outcome_file, "\n")
cat("分析基因:", paste(genes, collapse = ", "), "\n\n")

# ============================================
# 1. 读取signif_pairs.parquet文件
# ============================================
cat("步骤1: 读取signif_pairs数据...\n")

if (!requireNamespace("arrow", quietly = TRUE)) {
  cat("安装arrow包...\n")
  install.packages("arrow", repos = "https://cran.r-project.org")
}
library(arrow)

signif_pairs <- read_parquet(signif_pairs_file)
cat("signif_pairs总行数:", nrow(signif_pairs), "\n")
cat("signif_pairs列名:", paste(names(signif_pairs), collapse = ", "), "\n\n")

# ============================================
# 2. 创建基因名到phenotype_id的映射
# ============================================
cat("步骤2: 创建基因映射...\n")

# 从phenotype_id提取基因名（去除版本号）
gene_mapping <- data.frame(
  phenotype_id = unique(signif_pairs$phenotype_id),
  stringsAsFactors = FALSE
)
gene_mapping$gene_base <- sub("\\.[0-9]+$", "", gene_mapping$phenotype_id)

# ENSG ID到目标基因名的映射
ensg_to_gene <- c(
  "ENSG00000109320" = "NFKB1",
  "ENSG00000204305" = "AGER", 
  "ENSG00000044574" = "HSPA5",
  "ENSG00000137714" = "FDX1",
  "ENSG00000150768" = "HIF1A",
  "ENSG00000100644" = "STAT3",
  "ENSG00000168610" = "GPX4",
  "ENSG00000167468" = "HMOX1",
  "ENSG00000100292" = "DLAT"
)

gene_mapping$gene_name <- ensg_to_gene[gene_mapping$gene_base]
gene_mapping <- gene_mapping[!is.na(gene_mapping$gene_name), ]

cat("找到", nrow(gene_mapping), "个目标基因\n")
cat("目标基因:", paste(gene_mapping$gene_name, collapse = ", "), "\n\n")

# ============================================
# 3. 筛选目标基因的SNP并解析variant_id
# ============================================
cat("步骤3: 筛选并解析SNP（提取chr/pos/ref/alt）...\n")

# 筛选目标基因的SNP
target_phenotype_ids <- gene_mapping$phenotype_id
exposure_snps <- signif_pairs[signif_pairs$phenotype_id %in% target_phenotype_ids, ]
cat("目标基因SNP数:", nrow(exposure_snps), "\n")

# 添加基因名
phenotype_to_gene <- setNames(gene_mapping$gene_name, gene_mapping$phenotype_id)
exposure_snps$gene_name <- phenotype_to_gene[exposure_snps$phenotype_id]

# 从variant_id解析chr/pos/ref/alt
# variant_id格式: chr1_665098_G_A_b38 (chr_pos_ref_alt_build)
parse_variant <- function(variant_id) {
  parts <- strsplit(as.character(variant_id), "_")[[1]]
  if (length(parts) >= 4) {
    chr <- as.integer(gsub("chr", "", parts[1]))
    pos <- as.numeric(parts[2])
    ref <- parts[3]
    alt <- parts[4]
    return(c(chr = chr, pos = pos, ref = ref, alt = alt))
  } else {
    return(c(chr = NA, pos = NA, ref = NA, alt = NA))
  }
}

# 应用解析（向量化以提高效率）
variant_parts <- t(sapply(exposure_snps$variant_id, parse_variant))
exposure_snps$chr <- as.integer(variant_parts[, "chr"])
exposure_snps$pos <- as.numeric(variant_parts[, "pos"])
exposure_snps$ref <- variant_parts[, "ref"]
exposure_snps$alt <- variant_parts[, "alt"]

# 移除无法解析的SNP
valid_snps <- !is.na(exposure_snps$chr) & !is.na(exposure_snps$pos) & 
              !is.na(exposure_snps$ref) & !is.na(exposure_snps$alt)
exposure_snps <- exposure_snps[valid_snps, ]
cat("有效SNP数（成功解析chr/pos/ref/alt）:", nrow(exposure_snps), "\n")

# 检查每个基因的SNP数
gene_snp_counts <- table(exposure_snps$gene_name)
cat("各基因SNP数:\n")
print(gene_snp_counts)
cat("\n")

if (nrow(exposure_snps) == 0) {
  cat("错误: 没有有效SNP\n")
  quit(status = 1)
}

# 创建chr:pos用于匹配
exposure_snps$chr_pos <- paste(exposure_snps$chr, exposure_snps$pos, sep = ":")
unique_chr_pos <- unique(exposure_snps$chr_pos)
cat("唯一chr:pos数:", length(unique_chr_pos), "\n")

# ============================================
# 4. 读取Python预过滤的结局数据
# ============================================
cat("\n步骤4: 读取预过滤的结局数据...\n")

# 使用Python预过滤的小文件
filtered_file <- file.path(output_dir, "outcome_filtered_gtex.csv")

if (!file.exists(filtered_file)) {
  cat("错误: 未找到预过滤文件:", filtered_file, "\n")
  cat("请先运行Python过滤脚本\n")
  quit(status = 1)
}

outcome_raw <- fread(filtered_file)
cat("结局数据行数:", nrow(outcome_raw), "\n")
cat("结局数据列名:", paste(names(outcome_raw), collapse = ", "), "\n\n")

# 确保chr和pos是数值型
outcome_raw$chr.outcome <- as.integer(outcome_raw$chr.outcome)
outcome_raw$pos.outcome <- as.numeric(outcome_raw$pos.outcome)

# ============================================
# 5. 循环每个基因跑MR（使用chr/pos/allele匹配）
# ============================================
cat("\n步骤5: 逐个基因跑MR（chr/pos/allele匹配）...\n\n")

results_list <- list()

for (gene in genes) {
  cat("==========", gene, "==========\n")
  
  # 提取该基因暴露数据
  exp_gene <- exposure_snps[exposure_snps$gene_name == gene, ]
  cat("  暴露SNP:", nrow(exp_gene), "个\n")
  
  if (nrow(exp_gene) == 0) {
    cat("  跳过: 暴露数据中无该基因\n\n")
    next
  }
  
  # 提取暴露的chr:pos用于匹配
  exp_chr_pos <- paste(exp_gene$chr, exp_gene$pos, sep = ":")
  cat("  暴露chr:pos数:", length(unique(exp_chr_pos)), "个\n")
  
  # 在结局数据中匹配chr:pos
  out_matched <- outcome_raw[paste(outcome_raw$chr.outcome, outcome_raw$pos.outcome, sep = ":") %in% exp_chr_pos, ]
  cat("  按chr:pos匹配到结局SNP:", nrow(out_matched), "个\n")
  
  if (nrow(out_matched) == 0) {
    cat("  跳过: 未匹配到结局SNP\n\n")
    next
  }
  
  # 准备暴露数据框（用于合并）
  exp_df <- data.frame(
    chr = exp_gene$chr,
    pos = exp_gene$pos,
    ref = exp_gene$ref,
    alt = exp_gene$alt,
    beta.exp = exp_gene$slope,
    se.exp = exp_gene$slope_se,
    pval.exp = exp_gene$pval_nominal,
    eaf.exp = exp_gene$af,
    gene = exp_gene$gene_name,
    variant_id = exp_gene$variant_id,
    stringsAsFactors = FALSE
  )
  
  # 准备结局数据框（用于合并）
  out_df <- data.frame(
    chr = out_matched$chr.outcome,
    pos = out_matched$pos.outcome,
    beta.out = out_matched$beta.outcome,
    se.out = out_matched$se.outcome,
    pval.out = out_matched$pval.outcome,
    eaf.out = out_matched$eaf.outcome,
    effect_allele.out = out_matched$effect_allele.outcome,
    other_allele.out = out_matched$other_allele.outcome,
    SNP.out = out_matched$SNP,
    samplesize.out = out_matched$samplesize.outcome,
    outcome = out_matched$outcome,
    stringsAsFactors = FALSE
  )
  
  # 按chr+pos合并
  merged <- merge(exp_df, out_df, by = c("chr", "pos"), suffixes = c(".exp", ".out"))
  cat("  合并后行数:", nrow(merged), "\n")
  
  if (nrow(merged) == 0) {
    cat("  跳过: 合并后无数据\n\n")
    next
  }
  
  # 检查等位基因是否匹配（允许翻转）
  merged$allele_match <- with(merged,
    (alt == effect_allele.out & ref == other_allele.out) |
    (alt == other_allele.out & ref == effect_allele.out)
  )
  
  cat("  等位基因匹配/翻转:", sum(merged$allele_match), "/", nrow(merged), "\n")
  
  # 只保留等位基因匹配的行
  merged_clean <- merged[merged$allele_match, ]
  
  if (nrow(merged_clean) == 0) {
    cat("  跳过: 无等位基因匹配的SNP\n\n")
    next
  }
  
  # 处理翻转（如果暴露和结局的效应等位基因相反，beta取反）
  merged_clean$beta.out.flip <- with(merged_clean,
    ifelse(alt == effect_allele.out, beta.out, -beta.out)
  )
  
  merged_clean$eaf.out.flip <- with(merged_clean,
    ifelse(alt == effect_allele.out, eaf.out, 1 - eaf.out)
  )
  
  # 构建harmonised数据框
  dat <- data.frame(
    SNP = merged_clean$SNP.out,
    beta.exposure = merged_clean$beta.exp,
    beta.outcome = merged_clean$beta.out.flip,
    se.exposure = merged_clean$se.exp,
    se.outcome = merged_clean$se.out,
    effect_allele.exposure = merged_clean$alt,
    other_allele.exposure = merged_clean$ref,
    effect_allele.outcome = merged_clean$effect_allele.out,
    other_allele.outcome = merged_clean$other_allele.out,
    eaf.exposure = merged_clean$eaf.exp,
    eaf.outcome = merged_clean$eaf.out.flip,
    pval.exposure = merged_clean$pval.exp,
    pval.outcome = merged_clean$pval.out,
    samplesize.outcome = merged_clean$samplesize.out,
    gene.exposure = merged_clean$gene,
    id.exposure = merged_clean$gene,
    exposure = merged_clean$gene,
    outcome = merged_clean$outcome,
    id.outcome = "Stroke",
    mr_keep = TRUE,
    stringsAsFactors = FALSE
  )
  
  cat("  Harmonised SNP:", nrow(dat), "个\n")
  
  # 检查SNP数量
  if (nrow(dat) < 1) {
    cat("  跳过: SNP不足\n\n")
    next
  }
  
  # 计算F统计量并过滤弱工具变量（F>10）
  dat$F_stat <- (dat$beta.exposure / dat$se.exposure)^2
  dat_filtered <- dat[dat$F_stat > 10, ]
  
  if (nrow(dat_filtered) == 0) {
    cat("  跳过: F统计量均<10（弱工具变量）\n\n")
    next
  }
  
  cat("  F>10的SNP:", nrow(dat_filtered), "个\n")
  
  # 跑MR
  res <- tryCatch({
    mr(dat_filtered, method_list = c("mr_ivw", "mr_egger_regression",
                                     "mr_weighted_median", "mr_weighted_mode"))
  }, error = function(e) {
    cat("  MR分析失败:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(res) || nrow(res) == 0) {
    cat("  MR分析无结果\n\n")
    next
  }
  
  # 敏感性分析
  hetero <- tryCatch(mr_heterogeneity(dat_filtered), error = function(e) NULL)
  pleio  <- tryCatch(mr_pleiotropy_test(dat_filtered), error = function(e) NULL)
  
  # 保存结果
  results_list[[gene]] <- list(
    gene    = gene,
    dat     = dat_filtered,
    res     = res,
    hetero  = hetero,
    pleio   = pleio
  )
  
  # 打印关键结果
  ivw <- res[res$method == "Inverse variance weighted", ]
  if (nrow(ivw) > 0 && !is.na(ivw$b[1])) {
    cat("  IVW Beta:", round(ivw$b[1], 4),
        "| OR:", round(exp(ivw$b[1]), 3),
        "| P:", format(ivw$pval[1], digits = 3, scientific = TRUE), "\n")
    
    # 计算统计功效
    if (!is.na(ivw$se[1]) && ivw$se[1] > 0) {
      power <- pnorm(-1.96 + abs(ivw$b[1])/ivw$se[1]) + pnorm(-1.96 - abs(ivw$b[1])/ivw$se[1])
      cat("  统计功效:", round(power * 100, 1), "%\n")
    }
    cat("\n")
  } else {
    cat("  IVW结果无效或缺失\n\n")
  }
  
  # 保存该基因结果
  write.csv(res, file.path(output_dir, paste0(gene, "_MR_results.csv")), row.names = FALSE)
  write.csv(dat_filtered, file.path(output_dir, paste0(gene, "_harmonised_data.csv")), row.names = FALSE)
}

# ============================================
# 6. 汇总结果表
# ============================================
cat("\n步骤6: 生成汇总结果...\n")

summary_df <- data.frame(
  Gene = character(),
  Method = character(),
  Beta = numeric(),
  SE = numeric(),
  OR = numeric(),
  P_value = numeric(),
  Power = numeric(),
  stringsAsFactors = FALSE
)

for (gene in names(results_list)) {
  res <- results_list[[gene]]$res
  for (i in 1:nrow(res)) {
    beta_val <- res$b[i]
    se_val <- res$se[i]
    if (!is.na(beta_val) && is.numeric(beta_val)) {
      # 计算功效
      power <- NA
      if (!is.na(se_val) && se_val > 0) {
        power <- pnorm(-1.96 + abs(beta_val)/se_val) + pnorm(-1.96 - abs(beta_val)/se_val)
      }
      summary_df <- rbind(summary_df, data.frame(
        Gene = gene,
        Method = res$method[i],
        Beta = beta_val,
        SE = se_val,
        OR = exp(beta_val),
        P_value = res$pval[i],
        Power = round(power * 100, 1)
      ))
    }
  }
}

write.csv(summary_df, file.path(output_dir, "MR_summary_all_genes.csv"), row.names = FALSE)

# 只保留IVW结果
ivw_summary <- summary_df[summary_df$Method == "Inverse variance weighted", ]
if (nrow(ivw_summary) > 0) {
  write.csv(ivw_summary, file.path(output_dir, "MR_summary_IVW_only.csv"), row.names = FALSE)
  
  cat("\n========== IVW汇总结果 ==========\n")
  print(ivw_summary[, c("Gene", "Beta", "OR", "P_value", "Power")])
}

cat("\n========== 分析完成 ==========\n")
cat("输出目录:", output_dir, "\n")
cat("成功分析基因数:", length(results_list), "\n")
if (length(results_list) > 0) {
  cat("成功分析的基因:", paste(names(results_list), collapse = ", "), "\n")
}
cat("汇总结果已保存\n")
