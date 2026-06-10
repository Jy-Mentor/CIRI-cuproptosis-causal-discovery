#!/usr/bin/env Rscript
# MR分析脚本 - eQTLGen p=5e-08数据集
# 使用chr/pos/allele匹配，不依赖rsID

library(readxl)
library(data.table)
library(readr)
library(TwoSampleMR)
library(dplyr)
library(purrr)

# ============================================
# 参数设置
# ============================================
exposure_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-08_kb_1000_r2_0.01.xlsx"
outcome_file  <- "D:/EQTL/eqtlgen_ieu_outcome.csv"
genes         <- c("NFKB1", "FDX1", "STAT3", "HIF1A", "HMOX1",
                   "GPX4", "HSPA5", "AGER", "DLAT")
output_dir    <- "D:/EQTL/mr_results_p5e-08"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== eQTLGen p=5e-08 MR 分析（chr/pos/allele匹配版） ==========\n")
cat("暴露文件:", exposure_file, "\n")
cat("结局文件:", outcome_file, "\n")
cat("分析基因:", paste(genes, collapse = ", "), "\n\n")

# ============================================
# 1. 读取暴露数据（Excel文件）
# ============================================
cat("步骤1: 读取暴露数据...\n")

exposure_data <- read_excel(exposure_file)
cat("暴露数据总行数:", nrow(exposure_data), "\n")
cat("暴露数据列名:", paste(names(exposure_data), collapse = ", "), "\n\n")

# 查看前几行了解数据结构
cat("暴露数据前5行:\n")
print(head(exposure_data, 5))
cat("\n")

# ============================================
# 2. 筛选目标基因
# ============================================
cat("步骤2: 筛选目标基因...\n")

# 确定基因名列
gene_col <- NULL
possible_gene_cols <- c("Gene", "gene", "GENE", "gene_name", "GeneName", "gene_symbol")
for (col in possible_gene_cols) {
  if (col %in% names(exposure_data)) {
    gene_col <- col
    break
  }
}

if (is.null(gene_col)) {
  cat("错误: 无法找到基因名列。可用列名:", paste(names(exposure_data), collapse = ", "), "\n")
  quit(status = 1)
}

cat("基因名列:", gene_col, "\n")

# 筛选目标基因
exposure_snps <- exposure_data[exposure_data[[gene_col]] %in% genes, ]
cat("目标基因SNP数:", nrow(exposure_snps), "\n")

# 检查每个基因的SNP数
gene_snp_counts <- table(exposure_snps[[gene_col]])
cat("各基因SNP数:\n")
print(gene_snp_counts)
cat("\n")

if (nrow(exposure_snps) == 0) {
  cat("错误: 没有目标基因的SNP\n")
  quit(status = 1)
}

# ============================================
# 3. 确定列名并解析数据
# ============================================
cat("步骤3: 解析暴露数据列...\n")

cols <- names(exposure_snps)

# 查找各列
find_col <- function(patterns, cols) {
  for (p in patterns) {
    if (p %in% cols) return(p)
  }
  return(NULL)
}

chr_col <- find_col(c("chr.exposure", "Chr", "chr", "CHR", "Chromosome", "chromosome"), cols)
pos_col <- find_col(c("pos.exposure", "Pos", "pos", "POS", "Position", "position", "BP", "bp"), cols)
snp_col <- find_col(c("SNP", "snp", "RSID", "rsid", "rsId"), cols)
beta_col <- find_col(c("beta.exposure", "beta", "BETA", "Beta", "slope", "SLOPE", "effect"), cols)
se_col <- find_col(c("se.exposure", "se", "SE", "std_err", "StdErr", "slope_se"), cols)
pval_col <- find_col(c("pval.exposure", "pval", "Pval", "PVAL", "pvalue", "Pvalue", "P", "p"), cols)
ea_col <- find_col(c("effect_allele.exposure", "A1", "a1", "effect_allele", "EffectAllele", "alt", "ALT"), cols)
oa_col <- find_col(c("other_allele.exposure", "A2", "a2", "other_allele", "OtherAllele", "ref", "REF"), cols)
eaf_col <- find_col(c("eaf.exposure", "eaf", "EAF", "freq", "FREQ", "Freq", "AF", "af", "maf", "MAF"), cols)

cat("列映射:\n")
cat("  基因:", gene_col, "\n")
cat("  chr:", chr_col, "\n")
cat("  pos:", pos_col, "\n")
cat("  SNP:", snp_col, "\n")
cat("  beta:", beta_col, "\n")
cat("  se:", se_col, "\n")
cat("  pval:", pval_col, "\n")
cat("  effect allele:", ea_col, "\n")
cat("  other allele:", oa_col, "\n")
cat("  EAF:", eaf_col, "\n\n")

# 检查必需列
required_cols <- c(chr_col, pos_col, beta_col, se_col, ea_col, oa_col)
missing <- required_cols[sapply(required_cols, is.null)]
if (length(missing) > 0) {
  cat("错误: 缺失必需列:", paste(missing, collapse = ", "), "\n")
  quit(status = 1)
}

# 重命名列为标准名称
exposure_snps$chr <- as.integer(exposure_snps[[chr_col]])
exposure_snps$pos <- as.numeric(exposure_snps[[pos_col]])
exposure_snps$gene_name <- exposure_snps[[gene_col]]
exposure_snps$beta <- as.numeric(exposure_snps[[beta_col]])
exposure_snps$se <- as.numeric(exposure_snps[[se_col]])
exposure_snps$effect_allele <- exposure_snps[[ea_col]]
exposure_snps$other_allele <- exposure_snps[[oa_col]]

if (!is.null(pval_col)) {
  exposure_snps$pval <- as.numeric(exposure_snps[[pval_col]])
} else {
  exposure_snps$pval <- 1e-10  # 默认值
}

if (!is.null(eaf_col)) {
  exposure_snps$eaf <- as.numeric(exposure_snps[[eaf_col]])
} else {
  exposure_snps$eaf <- 0.5  # 默认值
}

if (!is.null(snp_col)) {
  exposure_snps$SNP <- exposure_snps[[snp_col]]
} else {
  exposure_snps$SNP <- paste0("chr", exposure_snps$chr, "_", exposure_snps$pos)
}

# 移除NA值
valid <- !is.na(exposure_snps$chr) & !is.na(exposure_snps$pos) & 
         !is.na(exposure_snps$beta) & !is.na(exposure_snps$se) &
         !is.na(exposure_snps$effect_allele) & !is.na(exposure_snps$other_allele)
exposure_snps <- exposure_snps[valid, ]

cat("有效SNP数:", nrow(exposure_snps), "\n\n")

# 创建chr:pos用于匹配
exposure_snps$chr_pos <- paste(exposure_snps$chr, exposure_snps$pos, sep = ":")
unique_chr_pos <- unique(exposure_snps$chr_pos)
cat("唯一chr:pos数:", length(unique_chr_pos), "\n")

# 保存chr:pos列表
chr_pos_file <- file.path(output_dir, "exposure_chr_pos_list.txt")
writeLines(unique_chr_pos, chr_pos_file)
cat("chr:pos列表已保存到:", chr_pos_file, "\n\n")

# ============================================
# 4. 读取Python预过滤的结局数据或使用原始结局文件
# ============================================
cat("步骤4: 读取结局数据...\n")

# 首先尝试使用已有的过滤文件
filtered_file <- file.path(output_dir, "outcome_filtered.csv")

if (file.exists(filtered_file)) {
  cat("使用预过滤的结局文件...\n")
  outcome_raw <- fread(filtered_file)
} else {
  # 需要先用Python过滤
  cat("错误: 未找到预过滤文件:", filtered_file, "\n")
  cat("请先运行Python过滤脚本过滤结局数据\n")
  cat("chr:pos列表已保存，可用于Python过滤\n")
  quit(status = 1)
}

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
    beta.exp = exp_gene$beta,
    se.exp = exp_gene$se,
    pval.exp = exp_gene$pval,
    eaf.exp = exp_gene$eaf,
    effect_allele.exp = exp_gene$effect_allele,
    other_allele.exp = exp_gene$other_allele,
    SNP.exp = exp_gene$SNP,
    gene = exp_gene$gene_name,
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
    (effect_allele.exp == effect_allele.out & other_allele.exp == other_allele.out) |
    (effect_allele.exp == other_allele.out & other_allele.exp == effect_allele.out)
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
    ifelse(effect_allele.exp == effect_allele.out, beta.out, -beta.out)
  )
  
  merged_clean$eaf.out.flip <- with(merged_clean,
    ifelse(effect_allele.exp == effect_allele.out, eaf.out, 1 - eaf.out)
  )
  
  # 构建harmonised数据框
  dat <- data.frame(
    SNP = merged_clean$SNP.out,
    beta.exposure = merged_clean$beta.exp,
    beta.outcome = merged_clean$beta.out.flip,
    se.exposure = merged_clean$se.exp,
    se.outcome = merged_clean$se.out,
    effect_allele.exposure = merged_clean$effect_allele.exp,
    other_allele.exposure = merged_clean$other_allele.exp,
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
