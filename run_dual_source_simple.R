#!/usr/bin/env Rscript
# ================================================================================
# 双源 eQTL 整合 MR 分析
# 整合 eQTLGen 全血与 GTEx 脑组织数据
# 参考实践：
#   - GTEx Portal v11: https://gtexportal.org/
#   - eQTL Catalogue: https://eqtlcatalogue.org/
#   - TwoSampleMR: https://mrcieu.github.io/TwoSampleMR/
#   - MR-Base: https://www.mrbase.org/
# ================================================================================

# 包管理
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("正在安装:", pkg))
      tryCatch({
        install.packages(pkg, repos = "https://cloud.r-project.org/")
      }, error = function(e) {})
    }
  }
}

install_if_missing(c("arrow", "dplyr", "data.table", "readr"))

library(arrow)
library(dplyr)
library(data.table)
library(readr)

cat("======================================================================\n")
cat("双源 eQTL 整合 MR 分析\n")
cat("参考：GTEx v11 | eQTL Catalogue | TwoSampleMR\n")
cat("======================================================================\n\n")

# 配置
GTEx_BRAIN_FILE <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Brain_Cortex.v11.eQTLs.signif_pairs.parquet"
GTEx_BLOOD_FILE <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eQTLs.signif_pairs.parquet"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_dual_source"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", OUTPUT_DIR, "\n\n")

# ================================================================================
# 1. 加载 eGenes 列表
# ================================================================================
cat("步骤 1: 加载 eGenes 列表\n")
cat("----------------------------------------------------------------------\n")

brain_genes_df <- fread("C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Brain_Cortex.v11.eGenes.txt")
blood_genes_df <- fread("C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eGenes.txt")

brain_genes <- unique(brain_genes_df$gene_id)
blood_genes <- unique(blood_genes_df$gene_id)

cat(sprintf("  ✓ 脑皮层 eGenes: %d 个基因\n", length(brain_genes)))
cat(sprintf("  ✓ 全血 eGenes: %d 个基因\n", length(blood_genes)))
cat(sprintf("  ✓ 共同基因：%d 个\n\n", length(intersect(brain_genes, blood_genes))))

# ================================================================================
# 2. 加载 eQTL 数据
# ================================================================================
cat("步骤 2: 加载 eQTL 数据\n")
cat("----------------------------------------------------------------------\n")

cat("  加载脑皮层 eQTL 数据...\n")
brain_eqtl <- read_parquet(GTEx_BRAIN_FILE)
cat(sprintf("  ✓ 脑皮层：%d 个 eQTL 对\n", nrow(brain_eqtl)))

cat("  加载全血 eQTL 数据...\n")
blood_eqtl <- read_parquet(GTEx_BLOOD_FILE)
cat(sprintf("  ✓ 全血：%d 个 eQTL 对\n\n", nrow(blood_eqtl)))

# ================================================================================
# 3. 格式化 eQTL 数据为 MR 输入
# 参考 GTEx v11 和 TwoSampleMR 标准格式
# ================================================================================
cat("步骤 3: 格式化 eQTL 数据为 MR 输入\n")
cat("----------------------------------------------------------------------\n")

format_eqtl_for_mr <- function(eqtl_data, gene_name, tissue_type) {
  # 格式化 eQTL 数据为 MR 分析输入格式
  # 参考：GTEx v11 标准 | TwoSampleMR 格式要求
  
  if (is.null(eqtl_data) || nrow(eqtl_data) == 0) {
    return(NULL)
  }
  
  # GTEx v11 使用 phenotype_id 格式：geneID_position_strand
  # 需要从中提取 gene_id
  # 参考：GTEx Portal v11 Documentation
  gene_data <- eqtl_data %>%
    mutate(gene_id_temp = sapply(strsplit(as.character(phenotype_id), "_"), function(x) x[1])) %>%
    filter(gene_id_temp == gene_name)
  
  if (nrow(gene_data) == 0) {
    return(NULL)
  }
  
  # 选择最强的 eQTL (最低 P 值) - 限制 Top 50
  # 参考 GTEx 推荐做法和 MR-Base 标准
  gene_data <- gene_data %>% 
    arrange(pval_nominal) %>% 
    slice_head(n = 50)
  
  # 解析 variant_id 获取染色体和位置
  # variant_id 格式：chr_pos_ref_alt_b38 (GTEx v11 标准)
  parse_variant_id <- function(vid) {
    tryCatch({
      parts <- strsplit(as.character(vid), "_")[[1]]
      if (length(parts) >= 4) {
        return(list(
          chr = gsub("chr", "", parts[1]),
          pos = as.integer(parts[2]),
          ref = parts[3],
          alt = parts[4]
        ))
      }
    }, error = function(e) {})
    return(list(chr="NA", pos=0, ref="NA", alt="NA"))
  }
  
  # 解析变异信息
  variant_info <- lapply(gene_data$variant_id, parse_variant_id)
  gene_data$CHR <- sapply(variant_info, function(x) x$chr)
  gene_data$BP <- sapply(variant_info, function(x) x$pos)
  gene_data$REF <- sapply(variant_info, function(x) x$ref)
  gene_data$ALT <- sapply(variant_info, function(x) x$alt)
  
  # 格式化输出 (TwoSampleMR 标准格式)
  # 参考：https://mrcieu.github.io/TwoSampleMR/reference/format_data.html
  mr_format <- gene_data %>%
    transmute(
      SNP = variant_id,
      CHR = CHR,
      BP = BP,
      EFFECT_ALLELE = ALT,
      OTHER_ALLELE = REF,
      BETA = slope,
      SE = slope_se,
      PVAL = pval_nominal,
      EAF = af,
      GENE = gene_name,
      TISSUE = tissue_type
    )
  
  return(mr_format)
}

# ================================================================================
# 4. 处理所有基因
# ================================================================================
cat("步骤 4: 处理所有基因\n")
cat("----------------------------------------------------------------------\n")

# 目标基因列表（所有在任一组织中表达的基因）
target_genes <- unique(c(brain_genes, blood_genes))
cat(sprintf("目标基因总数：%d\n\n", length(target_genes)))

stats <- list(both=0, brain_only=0, blood_only=0, neither=0)
processed <- 0

for (i in seq_along(target_genes)) {
  gene <- target_genes[i]
  output_file <- file.path(OUTPUT_DIR, paste0(gene, "_exposure.csv"))
  
  in_brain <- gene %in% brain_genes
  in_blood <- gene %in% blood_genes
  
  all_data <- list()
  
  # 添加脑组织数据
  if (in_brain) {
    brain_data <- format_eqtl_for_mr(brain_eqtl, gene, "Brain_Cortex")
    if (!is.null(brain_data) && nrow(brain_data) > 0) {
      all_data[[length(all_data) + 1]] <- brain_data
    }
  }
  
  # 添加全血数据
  if (in_blood) {
    blood_data <- format_eqtl_for_mr(blood_eqtl, gene, "Whole_Blood")
    if (!is.null(blood_data) && nrow(blood_data) > 0) {
      all_data[[length(all_data) + 1]] <- blood_data
    }
  }
  
  # 保存
  if (length(all_data) > 0) {
    combined <- bind_rows(all_data)
    fwrite(combined, output_file)
    processed <- processed + 1
    
    # 更新统计
    if (in_brain && in_blood) {
      stats$both <- stats$both + 1
    } else if (in_brain) {
      stats$brain_only <- stats$brain_only + 1
    } else if (in_blood) {
      stats$blood_only <- stats$blood_only + 1
    }
  } else {
    stats$neither <- stats$neither + 1
  }
  
  # 进度显示
  if (i %% 100 == 0) {
    cat(sprintf("  进度：%d/%d (已处理：%d)\n", i, length(target_genes), processed))
  }
}

# ================================================================================
# 5. 输出统计信息
# ================================================================================
cat("\n======================================================================\n")
cat("处理完成\n")
cat("======================================================================\n")

cat(sprintf("\n总基因数：%d\n", length(target_genes)))
cat(sprintf("成功处理：%d 个基因\n\n", processed))

cat("组织分布:\n")
cat(sprintf("  - 双组织都有：%d 个基因\n", stats$both))
cat(sprintf("  - 仅脑组织：  %d 个基因\n", stats$brain_only))
cat(sprintf("  - 仅全血：    %d 个基因\n", stats$blood_only))
cat(sprintf("  - 都无数据：  %d 个基因\n\n", stats$neither))

cat(sprintf("输出目录：%s\n", OUTPUT_DIR))

# 保存统计文件
stats_text <- sprintf(
  "双源 eQTL 整合统计\n",
  "参考：GTEx v11 | eQTL Catalogue\n\n",
  "总基因数：%d\n",
  "成功处理：%d\n\n",
  "组织分布:\n",
  "  - 双组织都有：%d\n",
  "  - 仅脑组织：%d\n",
  "  - 仅全血：%d\n",
  "  - 都无数据：%d\n",
  "输出目录：%s\n",
  length(target_genes), processed, 
  stats$both, stats$brain_only, stats$blood_only, stats$neither,
  OUTPUT_DIR
)
writeLines(stats_text, file.path(OUTPUT_DIR, "stats.txt"))

cat("\n✓ 统计文件已保存\n")

# 显示样本文件预览
if (processed > 0) {
  sample_files <- list.files(OUTPUT_DIR, pattern = "_exposure.csv$", full.names = TRUE)
  if (length(sample_files) > 0) {
    cat("\n样本文件预览:\n")
    sample_data <- fread(sample_files[1])
    cat(sprintf("  文件：%s\n", basename(sample_files[1])))
    cat(sprintf("  维度：%d 行 × %d 列\n", nrow(sample_data), ncol(sample_data)))
    cat(sprintf("  列名：%s\n\n", paste(colnames(sample_data), collapse=", ")))
  }
}

cat("======================================================================\n")
cat("下一步\n")
cat("======================================================================\n")
cat("
1. 检查暴露数据
   目录:", OUTPUT_DIR, "

2. 运行双源 MR 分析
   使用 TwoSampleMR 包进行 MR 分析
   命令：Rscript run_mr_analysis.R

3. 比较单源 vs 双源结果
   - 查看新增的显著基因
   - 分析组织特异性效应

参考资源:
- GTEx Portal v11: https://gtexportal.org/
- eQTL Catalogue: https://eqtlcatalogue.org/
- TwoSampleMR: https://mrcieu.github.io/TwoSampleMR/
- MR-Base: https://www.mrbase.org/
")

cat("\n======================================================================\n")
