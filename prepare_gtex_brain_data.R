#!/usr/bin/env Rscript
# ================================================================================
# 准备 GTEx v11 脑组织 eQTL 数据
# 参考 GTEx Portal 和 eQTL Catalogue 的标准处理流程
# 来源：https://gtexportal.org/ | https://eqtlcatalogue.org/
# ================================================================================

# 自动安装包
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("正在安装:", pkg))
      tryCatch({
        install.packages(pkg, repos = "https://cloud.r-project.org/")
      }, error = function(e) {
        message(paste("安装失败:", pkg, "-", e$message))
      })
    }
  }
}

# 安装必要的包
install_if_missing(c("arrow", "dplyr", "data.table", "readr"))

library(arrow)
library(dplyr)
library(data.table)
library(readr)

# 配置 - 使用用户实际数据路径
GTEx_BRAIN_FILE <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Brain_Cortex.v11.eQTLs.signif_pairs.parquet"
GTEx_BLOOD_FILE <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eQTLs.signif_pairs.parquet"
GTEx_BRAIN_EGENES <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Brain_Cortex.v11.eGenes.txt"
GTEx_BLOOD_EGENES <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eGenes.txt"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_gtex_brain"
GENE_LIST_FILE <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/gene_list_optimized.txt"

cat("=" , rep("=", 69), "\n", sep="")
cat("准备 GTEx v11 脑组织 eQTL 数据\n")
cat("参考：GTEx Portal v11 | eQTL Catalogue\n")
cat("=", rep("=", 69), "\n", sep="")

# 创建输出目录
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("\n输出目录:", OUTPUT_DIR, "\n")

# ================================================================================
# 加载基因列表
# ================================================================================

gene_list <- character(0)
if (file.exists(GENE_LIST_FILE)) {
  gene_list <- readLines(GENE_LIST_FILE, warn = FALSE)
  gene_list <- gene_list[gene_list != ""]
  cat(sprintf("✓ 加载基因列表：%d 个基因\n", length(gene_list)))
} else {
  cat(sprintf("✗ 基因列表文件不存在：%s\n", GENE_LIST_FILE))
  cat("  将使用 eQTL 数据中的所有基因\n")
}

# ================================================================================
# 加载 GTEx 数据 (参考 GTEx v11 标准格式)
# ================================================================================

cat("\n", "=", rep("=", 69), "\n", sep="")
cat("加载 GTEx v11 数据\n")
cat("=", rep("=", 69), "\n", sep="")

# 加载脑组织 eGenes (基因列表)
brain_genes_df <- NULL
if (file.exists(GTEx_BRAIN_EGENES)) {
  cat("\n加载脑皮层 eGenes 列表...\n")
  tryCatch({
    brain_genes_df <- fread(GTEx_BRAIN_EGENES, header=TRUE, sep="\t")
    cat(sprintf("  ✓ 脑皮层 eGenes: %d 个基因\n", nrow(brain_genes_df)))
  }, error = function(e) {
    cat(sprintf("  ✗ 加载失败：%s\n", e$message))
    brain_genes_df <- NULL
  })
}

# 加载全血 eGenes (基因列表)
blood_genes_df <- NULL
if (file.exists(GTEx_BLOOD_EGENES)) {
  cat("\n加载全血 eGenes 列表...\n")
  tryCatch({
    blood_genes_df <- fread(GTEx_BLOOD_EGENES, header=TRUE, sep="\t")
    cat(sprintf("  ✓ 全血 eGenes: %d 个基因\n", nrow(blood_genes_df)))
  }, error = function(e) {
    cat(sprintf("  ✗ 加载失败：%s\n", e$message))
    blood_genes_df <- NULL
  })
}

# 加载脑组织数据 (Brain Cortex)
brain_eqtl <- NULL
if (file.exists(GTEx_BRAIN_FILE)) {
  cat("\n加载脑皮层 eQTL 数据...\n")
  cat("  文件:", GTEx_BRAIN_FILE, "\n")
  tryCatch({
    brain_eqtl <- read_parquet(GTEx_BRAIN_FILE)
    cat(sprintf("  ✓ 脑皮层数据：%d 个 eQTL 对\n", nrow(brain_eqtl)))
    cat(sprintf("    列：%s\n", paste(colnames(brain_eqtl), collapse=", ")))
  }, error = function(e) {
    cat(sprintf("  ✗ 加载失败：%s\n", e$message))
    brain_eqtl <- NULL
  })
} else {
  cat(sprintf("  ✗ 文件不存在：%s\n", GTEx_BRAIN_FILE))
}

# 加载全血数据 (Whole Blood)
blood_eqtl <- NULL
if (file.exists(GTEx_BLOOD_FILE)) {
  cat("\n加载全血 eQTL 数据...\n")
  cat("  文件:", GTEx_BLOOD_FILE, "\n")
  tryCatch({
    blood_eqtl <- read_parquet(GTEx_BLOOD_FILE)
    cat(sprintf("  ✓ 全血数据：%d 个 eQTL 对\n", nrow(blood_eqtl)))
    cat(sprintf("    列：%s\n", paste(colnames(blood_eqtl), collapse=", ")))
  }, error = function(e) {
    cat(sprintf("  ✗ 加载失败：%s\n", e$message))
    blood_eqtl <- NULL
  })
} else {
  cat(sprintf("  ✗ 文件不存在：%s\n", GTEx_BLOOD_FILE))
}

# ================================================================================
# 格式化数据为 MR 输入
# 参考 TwoSampleMR 和 MR-Base 标准格式
# ================================================================================

format_eqtl_for_mr <- function(eqtl_data, gene_name, tissue_type) {
  # 格式化 eQTL 数据为 MR 分析输入格式
  # 参考：https://mrcieu.github.io/TwoSampleMR/
  
  if (is.null(eqtl_data) || nrow(eqtl_data) == 0) {
    return(NULL)
  }
  
  # 筛选该基因的 eQTL
  # GTEx v11 使用 phenotype_id 格式：gene_id|phenotype_id
  # 我们需要从 phenotype_id 中提取 gene_id
  # phenotype_id 格式：ENSG00000000003.14_136872842_+
  gene_data <- eqtl_data %>%
    mutate(gene_id = sapply(strsplit(phenotype_id, "_"), function(x) x[1])) %>%
    filter(gene_id == gene_name)
  
  if (nrow(gene_data) == 0) {
    return(NULL)
  }
  
  # 选择最强的 eQTL (最低 P 值) - 限制 Top 50
  # 参考 GTEx 推荐做法
  gene_data <- gene_data %>% 
    arrange(pval_nominal) %>% 
    slice_head(n = 50)
  
  # 解析 variant_id 获取染色体和位置
  # variant_id 格式：chr_pos_ref_alt_b38 (GTEx v11 标准)
  parse_variant_id <- function(variant_id) {
    tryCatch({
      parts <- strsplit(as.character(variant_id), "_")[[1]]
      if (length(parts) >= 4) {
        chr_ <- gsub("chr", "", parts[1])
        pos <- as.integer(parts[2])
        ref <- parts[3]
        alt <- parts[4]
        return(list(chr=chr_, pos=pos, ref=ref, alt=alt))
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
      TISSUE = tissue_type,
      TSS_DISTANCE = start_distance
    )
  
  return(mr_format)
}

# ================================================================================
# 为每个基因创建暴露文件
# ================================================================================

cat("\n", "=", rep("=", 69), "\n", sep="")
cat("格式化 eQTL 数据为 MR 输入\n")
cat("参考：TwoSampleMR / MR-Base 标准\n")
cat("=", rep("=", 69), "\n", sep="")

# 获取所有可用基因
brain_genes <- character(0)
if (!is.null(brain_genes_df)) {
  brain_genes <- unique(brain_genes_df$gene_id)
  cat(sprintf("\n脑皮层 eQTL 数据包含：%d 个基因\n", length(brain_genes)))
}

blood_genes <- character(0)
if (!is.null(blood_genes_df)) {
  blood_genes <- unique(blood_genes_df$gene_id)
  cat(sprintf("全血 eQTL 数据包含：%d 个基因\n", length(blood_genes)))
}

# 交集和并集
common_genes <- intersect(brain_genes, blood_genes)
cat(sprintf("\n双组织共同基因：%d 个\n", length(common_genes)))

# 如果有基因列表，使用基因列表；否则使用所有 eQTL 基因
if (length(gene_list) > 0) {
  target_genes <- gene_list
  cat(sprintf("\n目标基因：%d 个（来自基因列表）\n", length(target_genes)))
} else {
  # 使用所有在任一组织中表达的基因
  target_genes <- unique(c(brain_genes, blood_genes))
  cat(sprintf("\n目标基因：%d 个（来自 eQTL 数据）\n", length(target_genes)))
}

# 统计
stats <- list(
  brain_only = 0,
  blood_only = 0,
  both = 0,
  neither = 0
)

# 处理每个基因
cat(sprintf("\n开始处理 %d 个基因...\n", length(target_genes)))
processed_count <- 0

for (i in seq_along(target_genes)) {
  gene <- target_genes[i]
  output_file <- file.path(OUTPUT_DIR, paste0(gene, "_exposure.csv"))
  
  # 检查基因是否在数据中
  in_brain <- gene %in% brain_genes
  in_blood <- gene %in% blood_genes
  
  all_data <- list()
  
  # 添加脑组织数据
  if (in_brain && !is.null(brain_eqtl)) {
    brain_data <- format_eqtl_for_mr(brain_eqtl, gene, "Brain_Cortex")
    if (!is.null(brain_data) && nrow(brain_data) > 0) {
      all_data[[1]] <- brain_data
    }
  }
  
  # 添加全血数据
  if (in_blood && !is.null(blood_eqtl)) {
    blood_data <- format_eqtl_for_mr(blood_eqtl, gene, "Whole_Blood")
    if (!is.null(blood_data) && nrow(blood_data) > 0) {
      all_data[[length(all_data) + 1]] <- blood_data
    }
  }
  
  # 保存
  if (length(all_data) > 0) {
    # 合并所有组织的数据
    combined <- bind_rows(all_data)
    fwrite(combined, output_file)
    processed_count <- processed_count + 1
    
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
  if (i %% 50 == 0) {
    cat(sprintf("  处理进度：%d/%d (已处理：%d)\n", i, length(target_genes), processed_count))
  }
}

# ================================================================================
# 输出统计信息
# ================================================================================

cat("\n", "=", rep("=", 69), "\n", sep="")
cat("处理完成\n")
cat("=", rep("=", 69), "\n", sep="")

cat(sprintf("\n总基因数：%d\n", length(target_genes)))
cat(sprintf("成功处理：%d 个基因\n", processed_count))
cat("\n组织分布:\n")
cat(sprintf("  - 双组织都有： %d 个基因\n", stats$both))
cat(sprintf("  - 仅脑组织：   %d 个基因\n", stats$brain_only))
cat(sprintf("  - 仅全血：     %d 个基因\n", stats$blood_only))
cat(sprintf("  - 都无数据：   %d 个基因\n", stats$neither))
cat(sprintf("\n输出目录：%s\n", OUTPUT_DIR))

# 创建统计文件
stats_file <- file.path(OUTPUT_DIR, "processing_stats.txt")
sink(stats_file)
cat("GTEx v11 eQTL 数据处理统计\n")
cat("参考：GTEx Portal v11 | eQTL Catalogue\n")
cat(paste(rep("=", 50), collapse=""), "\n\n")
cat(sprintf("总基因数：%d\n", length(target_genes)))
cat(sprintf("成功处理：%d\n\n", processed_count))
cat("组织分布:\n")
cat(sprintf("  - 双组织都有：%d\n", stats$both))
cat(sprintf("  - 仅脑组织：  %d\n", stats$brain_only))
cat(sprintf("  - 仅全血：    %d\n", stats$blood_only))
cat(sprintf("  - 都无数据：  %d\n\n", stats$neither))
cat(sprintf("输出目录：%s\n", OUTPUT_DIR))
sink()

cat(sprintf("\n✓ 已保存统计信息：%s\n", stats_file))

# 创建基因列表文件
gene_list_output <- file.path(OUTPUT_DIR, "gene_list_processed.txt")
writeLines(target_genes, gene_list_output)
cat(sprintf("✓ 已保存基因列表：%s\n", gene_list_output))

# 创建样本文件预览
if (processed_count > 0) {
  sample_gene <- target_genes[1]
  sample_file <- file.path(OUTPUT_DIR, paste0(sample_gene, "_exposure.csv"))
  if (file.exists(sample_file)) {
    cat("\n\n样本文件预览 (", sample_gene, "_exposure.csv):\n", sep="")
    sample_data <- fread(sample_file)
    cat(sprintf("维度：%d 行 × %d 列\n", nrow(sample_data), ncol(sample_data)))
    cat("列名:", paste(colnames(sample_data), collapse=", "), "\n")
    cat("\n前 6 行:\n")
    print(head(sample_data))
  }
}

cat("\n", "=", rep("=", 69), "\n", sep="")
cat("下一步\n")
cat("=", rep("=", 69), "\n", sep="")
cat("
1. 检查输出文件
   目录:", OUTPUT_DIR, "

2. 运行双源 MR 分析
   命令：& \"C:\\R\\R-4.5.2\\bin\\Rscript.exe\" run_dual_source_mr.R

3. 比较单源 vs 双源结果
   - 查看新增的显著基因
   - 分析组织特异性效应

参考资源:
- GTEx Portal: https://gtexportal.org/
- eQTL Catalogue: https://eqtlcatalogue.org/
- TwoSampleMR: https://mrcieu.github.io/TwoSampleMR/
- MR-Base: https://www.mrbase.org/
")

cat("\n", "=", rep("=", 69), "\n", sep="")
