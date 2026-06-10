#!/usr/bin/env Rscript
# ================================================================================
# GTEx 数据整合 - 使用简单独立 SNP 选择方法
# ================================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(dplyr)
})

cat("======================================================================\n")
cat("GTEx 数据整合 - 使用简单独立 SNP 选择方法\n")
cat("======================================================================\n\n")

# 配置
GTEx_DIR <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL"
OUTPUT_DIR <- "D:/EQTL/clump"
P_THRESHOLD <- 1e-5
KB_WINDOW <- 1000

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 138 个目标基因
target_genes <- c(
  "LYN", "PRKCQ", "NMT1", "TDP1", "MAN2B1", "IL10RA", "RHOC", "SREBF1",
  "KCNA5", "HIF1A", "CTSC", "CAT", "FABP4", "STAT5A", "FABP2", "B2M",
  "RBM39", "HBS1L", "CHFR", "NUDCD2", "TCN2", "SCN9A", "JAK1", "GPX1",
  "CTSB", "CASP8", "FABP5", "XDH", "MB", "POLR2D", "HSD17B10", "MAPKAPK2",
  "SEC13", "PCTP", "ZEB1", "RELA", "IRF1", "GFAP", "CPT2", "BRD3",
  "NR3C1", "F3", "C3", "ITGA1", "CITED2", "HIBADH", "SAT2", "TSPO",
  "PTGS1", "IMPDH2", "FLT4", "CPT1A", "AKT1", "CCR5", "PTPRF", "HPGDS",
  "PTPRJ", "CASK", "MGAT1", "IGFBP2", "TOP2A", "PPARG", "IL6", "EPHX1",
  "CP", "AIF1", "PLA2G4A", "ALDH9A1", "S100A6", "DDC", "CUL4B", "BST1",
  "CNDP2", "TNF", "PARP1", "IKBKB", "EGFR", "COL1A1", "ADRB1", "SPHK1",
  "GCH1", "ACADVL", "STARD13", "CTSD", "PDCD6IP", "PTPRC", "TGFB1", "PABPC1",
  "HTR2C", "CTSS", "CNR2", "ACTA2", "FNTA", "RENBP", "CCNA2", "PTGR1",
  "LEF1", "SAT1", "XRCC6", "TBXAS1", "NR1H3", "HTR2B", "CTSL", "CDK4",
  "CXCR3", "TIMP1", "OAZ1", "STK4", "ZHX2", "MKNK2", "SERPINB10", "ACADM",
  "STAT3", "NFKB1", "HSPA5", "CTSK", "CCND1", "PTPN2", "PTPN6", "PA2G4",
  "HSD17B4", "ACAD11", "PDCD6", "PARP12", "SERPINB1", "STAT1", "NFE2L2",
  "HMOX1", "CTSF", "CCL2", "MAOB", "ICAM1", "FDX1", "LIAS", "LIPT1",
  "DLAT", "PDHB", "PDHX", "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
  "HIF1A", "MTOR", "NFKB1", "GPX4"
)

target_genes <- unique(target_genes)
cat(sprintf("目标基因数：%d\n\n", length(target_genes)))

# 从 eGenes 文件中提取基因名到 ENSG 的映射
cat("从 GTEx eGenes 文件提取基因映射...\n")
blood_genes_file <- file.path(GTEx_DIR, "Whole_Blood.v11.eGenes.txt")
blood_genes <- fread(blood_genes_file, sep = "\t", header = TRUE)

blood_genes_mapping <- blood_genes[, c("gene_id", "gene_name")]
names(blood_genes_mapping) <- c("ensg_id", "gene_symbol")
blood_genes_mapping <- unique(blood_genes_mapping)
blood_genes_mapping <- blood_genes_mapping[!is.na(blood_genes_mapping$gene_symbol) & blood_genes_mapping$gene_symbol != "", ]

cat(sprintf("  全血 eGenes 映射数：%d\n\n", nrow(blood_genes_mapping)))

# 过滤目标基因的 ENSG ID
target_ensg_blood <- blood_genes_mapping[blood_genes_mapping$gene_symbol %in% target_genes, ]
cat(sprintf("  目标基因在全血中的 ENSG ID 数：%d\n\n", nrow(target_ensg_blood)))

# ================================================================================
# 处理全血 eQTL - 使用简单独立 SNP 选择
# ================================================================================
cat("处理全血 eQTL 数据（简单独立 SNP 选择）...\n")

blood_eqtl_file <- file.path(GTEx_DIR, "blood_eqtl.csv")

cat("  加载 eQTL 数据...\n")
blood_eqtl <- fread(blood_eqtl_file, header = TRUE)
cat(sprintf("  原始数据：%d 行\n", nrow(blood_eqtl)))

cat("  过滤显著 eQTL (P < 1e-5)...\n")
blood_sig <- blood_eqtl[blood_eqtl$pval_nominal < P_THRESHOLD, ]
cat(sprintf("  显著 eQTL 数：%d\n", nrow(blood_sig)))

cat("  过滤目标基因...\n")
target_ensg_ids <- target_ensg_blood$ensg_id
blood_target <- blood_sig[blood_sig$phenotype_id %in% target_ensg_ids, ]
cat(sprintf("  目标基因 eQTL 数：%d\n", nrow(blood_target)))

if (nrow(blood_target) > 0) {
  # 提取 rsID 和染色体位置信息
  cat("  提取 rsID 和染色体位置信息...\n")
  
  blood_target$rs_id <- NA
  blood_target$chr <- NA
  blood_target$pos <- NA
  
  for (i in 1:nrow(blood_target)) {
    gene_variant <- blood_genes[blood_genes$gene_id == blood_target$phenotype_id[i] & 
                                   blood_genes$variant_id == blood_target$variant_id[i], ]
    if (nrow(gene_variant) > 0 && !is.na(gene_variant$rs_id_dbSNP157_GRCh38p14) && 
        gene_variant$rs_id_dbSNP157_GRCh38p14 != "") {
      blood_target$rs_id[i] <- gene_variant$rs_id_dbSNP157_GRCh38p14
    }
    
    # 提取染色体和位置
    variant_parts <- strsplit(as.character(blood_target$variant_id[i]), "_")[[1]]
    if (length(variant_parts) >= 4) {
      chr <- gsub("chr", "", variant_parts[1])
      pos <- as.numeric(variant_parts[2])
      blood_target$chr[i] <- chr
      blood_target$pos[i] <- pos
      
      if (is.na(blood_target$rs_id[i])) {
        blood_target$rs_id[i] <- paste0("chr", chr, "_", pos, "_", variant_parts[3], "_", variant_parts[4])
      }
    }
  }
  
  cat(sprintf("  有 rsID 的 eQTL 数：%d\n", sum(!is.na(blood_target$rs_id))))
  cat(sprintf("  有染色体位置信息的 eQTL 数：%d\n", sum(!is.na(blood_target$chr))))
  
  # 转换为 TwoSampleMR 格式
  cat("  转换为 TwoSampleMR 格式...\n")
  blood_formatted <- data.frame(
    SNP = blood_target$rs_id,
    beta.exposure = blood_target$slope,
    se.exposure = blood_target$slope_se,
    effect_allele.exposure = sapply(strsplit(as.character(blood_target$variant_id), "_"), function(x) x[3]),
    other_allele.exposure = sapply(strsplit(as.character(blood_target$variant_id), "_"), function(x) x[4]),
    eaf.exposure = blood_target$af,
    pval.exposure = blood_target$pval_nominal,
    gene = blood_target$phenotype_id,
    chr = blood_target$chr,
    pos = blood_target$pos,
    exposure = "GTEx_Whole_Blood",
    stringsAsFactors = FALSE
  )
  
  blood_formatted <- blood_formatted[!is.na(blood_formatted$SNP) & 
                                      !is.na(blood_formatted$chr) & 
                                      !is.na(blood_formatted$pos), ]
  cat(sprintf("  格式化后：%d 行\n", nrow(blood_formatted)))
  
  # 简单独立 SNP 选择：按基因分组，每 KB_WINDOW kb 选择一个最显著的 SNP
  cat(sprintf("  进行简单独立 SNP 选择（每 %d kb 选择一个 SNP）...\n", KB_WINDOW))
  
  blood_formatted <- blood_formatted %>%
    arrange(gene, chr, pos, pval.exposure) %>%
    group_by(gene, chr, bin = floor(pos / KB_WINDOW)) %>%
    filter(pval.exposure == min(pval.exposure)) %>%
    ungroup() %>%
    select(-bin)
  
  cat(sprintf("  独立 SNP 选择后：%d 行\n\n", nrow(blood_formatted)))
  
  # 保存
  blood_output <- file.path(OUTPUT_DIR, "GTEx_Whole_Blood_p_1e-05_kb_1000_simple.tsv")
  write.table(blood_formatted, blood_output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  保存到：%s\n\n", blood_output))
} else {
  cat("  无目标基因 eQTL，跳过\n\n")
  blood_formatted <- NULL
}

# ================================================================================
# 处理脑组织 eQTL
# ================================================================================
cat("处理脑组织 eQTL 数据（简单独立 SNP 选择）...\n")

brain_genes_file <- file.path(GTEx_DIR, "Brain_Cortex.v11.eGenes.txt")
brain_genes <- fread(brain_genes_file, sep = "\t", header = TRUE)

brain_genes_mapping <- brain_genes[, c("gene_id", "gene_name")]
names(brain_genes_mapping) <- c("ensg_id", "gene_symbol")
brain_genes_mapping <- unique(brain_genes_mapping)
brain_genes_mapping <- brain_genes_mapping[!is.na(brain_genes_mapping$gene_symbol) & brain_genes_mapping$gene_symbol != "", ]

target_ensg_brain <- brain_genes_mapping[brain_genes_mapping$gene_symbol %in% target_genes, ]
cat(sprintf("  目标基因在脑组织中的 ENSG ID 数：%d\n\n", nrow(target_ensg_brain)))

brain_eqtl_file <- file.path(GTEx_DIR, "brain_eqtl.csv")

cat("  加载 eQTL 数据...\n")
brain_eqtl <- fread(brain_eqtl_file, header = TRUE)
cat(sprintf("  原始数据：%d 行\n", nrow(brain_eqtl)))

cat("  过滤显著 eQTL (P < 1e-5)...\n")
brain_sig <- brain_eqtl[brain_eqtl$pval_nominal < P_THRESHOLD, ]
cat(sprintf("  显著 eQTL 数：%d\n", nrow(brain_sig)))

cat("  过滤目标基因...\n")
target_ensg_ids_brain <- target_ensg_brain$ensg_id
brain_target <- brain_sig[brain_sig$phenotype_id %in% target_ensg_ids_brain, ]
cat(sprintf("  目标基因 eQTL 数：%d\n", nrow(brain_target)))

if (nrow(brain_target) > 0) {
  cat("  提取 rsID 和染色体位置信息...\n")
  
  brain_target$rs_id <- NA
  brain_target$chr <- NA
  brain_target$pos <- NA
  
  for (i in 1:nrow(brain_target)) {
    gene_variant <- brain_genes[brain_genes$gene_id == brain_target$phenotype_id[i] & 
                                 brain_genes$variant_id == brain_target$variant_id[i], ]
    if (nrow(gene_variant) > 0 && !is.na(gene_variant$rs_id_dbSNP157_GRCh38p14) && 
        gene_variant$rs_id_dbSNP157_GRCh38p14 != "") {
      brain_target$rs_id[i] <- gene_variant$rs_id_dbSNP157_GRCh38p14
    }
    
    variant_parts <- strsplit(as.character(brain_target$variant_id[i]), "_")[[1]]
    if (length(variant_parts) >= 4) {
      chr <- gsub("chr", "", variant_parts[1])
      pos <- as.numeric(variant_parts[2])
      brain_target$chr[i] <- chr
      brain_target$pos[i] <- pos
      
      if (is.na(brain_target$rs_id[i])) {
        brain_target$rs_id[i] <- paste0("chr", chr, "_", pos, "_", variant_parts[3], "_", variant_parts[4])
      }
    }
  }
  
  cat(sprintf("  有 rsID 的 eQTL 数：%d\n", sum(!is.na(brain_target$rs_id))))
  cat(sprintf("  有染色体位置信息的 eQTL 数：%d\n", sum(!is.na(brain_target$chr))))
  
  cat("  转换为 TwoSampleMR 格式...\n")
  brain_formatted <- data.frame(
    SNP = brain_target$rs_id,
    beta.exposure = brain_target$slope,
    se.exposure = brain_target$slope_se,
    effect_allele.exposure = sapply(strsplit(as.character(brain_target$variant_id), "_"), function(x) x[3]),
    other_allele.exposure = sapply(strsplit(as.character(brain_target$variant_id), "_"), function(x) x[4]),
    eaf.exposure = brain_target$af,
    pval.exposure = brain_target$pval_nominal,
    gene = brain_target$phenotype_id,
    chr = brain_target$chr,
    pos = brain_target$pos,
    exposure = "GTEx_Brain_Cortex",
    stringsAsFactors = FALSE
  )
  
  brain_formatted <- brain_formatted[!is.na(brain_formatted$SNP) & 
                                      !is.na(brain_formatted$chr) & 
                                      !is.na(brain_formatted$pos), ]
  cat(sprintf("  格式化后：%d 行\n", nrow(brain_formatted)))
  
  cat(sprintf("  进行简单独立 SNP 选择（每 %d kb 选择一个 SNP）...\n", KB_WINDOW))
  
  brain_formatted <- brain_formatted %>%
    arrange(gene, chr, pos, pval.exposure) %>%
    group_by(gene, chr, bin = floor(pos / KB_WINDOW)) %>%
    filter(pval.exposure == min(pval.exposure)) %>%
    ungroup() %>%
    select(-bin)
  
  cat(sprintf("  独立 SNP 选择后：%d 行\n\n", nrow(brain_formatted)))
  
  brain_output <- file.path(OUTPUT_DIR, "GTEx_Brain_Cortex_p_1e-05_kb_1000_simple.tsv")
  write.table(brain_formatted, brain_output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  保存到：%s\n\n", brain_output))
} else {
  cat("  无目标基因 eQTL，跳过\n\n")
  brain_formatted <- NULL
}

# ================================================================================
# 整合所有数据源
# ================================================================================
cat("整合所有 eQTL 数据源...\n")

eqtlgen_file <- file.path(OUTPUT_DIR, "eQTLgen_allgene_p_1e-05_kb_1000_r2_0.01.xlsx")
eqtlgen_data <- read_excel(eqtlgen_file)
cat(sprintf("  eQTLGen: %d SNPs\n", nrow(eqtlgen_data)))

gtex_sources <- list()
if (!is.null(blood_formatted)) {
  blood_clumped_df <- fread(blood_output, sep = "\t")
  cat(sprintf("  GTEx 全血：%d SNPs\n", nrow(blood_clumped_df)))
  gtex_sources[["blood"]] <- blood_clumped_df
}

if (!is.null(brain_formatted)) {
  brain_clumped_df <- fread(brain_output, sep = "\t")
  cat(sprintf("  GTEx 脑组织：%d SNPs\n", nrow(brain_clumped_df)))
  gtex_sources[["brain"]] <- brain_clumped_df
}

# 合并所有数据
all_data <- list(eqtlgen_data)
if (length(gtex_sources) > 0) {
  all_data <- c(all_data, gtex_sources)
}

merged_eqtl <- do.call(rbind, lapply(all_data, function(df) {
  df[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", 
         "other_allele.exposure", "eaf.exposure", "pval.exposure", "gene")]
}))

# 去重（优先保留 P 值最小的）
merged_eqtl <- merged_eqtl[order(merged_eqtl$pval.exposure), ]
merged_eqtl <- merged_eqtl[!duplicated(merged_eqtl$SNP), ]

cat(sprintf("  整合后：%d SNPs\n\n", nrow(merged_eqtl)))

# 保存整合数据
merged_output <- file.path(OUTPUT_DIR, "merged_eqtl_all_sources_p_1e-05_simple.tsv")
write.table(merged_eqtl, merged_output, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  整合数据保存到：%s\n\n", merged_output))

# 统计每个数据源的贡献
cat("各数据源贡献:\n")

# 通过检查 SNP ID 格式来判断数据源
# eQTLGen 的 SNP 是 rsID 格式，GTEx 的 SNP 包含 chr 前缀
eqtlgen_count <- sum(grepl("^rs", merged_eqtl$SNP))
blood_count <- sum(grepl("Whole_Blood", merged_eqtl$gene))
brain_count <- sum(grepl("Brain_Cortex", merged_eqtl$gene))

cat(sprintf("  eQTLGen: %d SNPs (%.1f%%)\n", eqtlgen_count, eqtlgen_count / nrow(merged_eqtl) * 100))
cat(sprintf("  GTEx 全血：%d SNPs (%.1f%%)\n", blood_count, blood_count / nrow(merged_eqtl) * 100))
cat(sprintf("  GTEx 脑组织：%d SNPs (%.1f%%)\n", brain_count, brain_count / nrow(merged_eqtl) * 100))

# 统计每个基因的 SNP 数
cat("\n每个基因的 SNP 统计:\n")
gene_counts <- table(merged_eqtl$gene)
cat(sprintf("  有 eQTL 的基因数：%d\n", length(gene_counts)))
cat(sprintf("  平均每个基因：%.0f 个 SNP\n", mean(gene_counts)))
cat(sprintf("  中位数：%.0f 个 SNP\n", median(gene_counts)))
cat(sprintf("  范围：%d - %d 个 SNP\n\n", min(gene_counts), max(gene_counts)))

cat("======================================================================\n")
cat("完成！\n")
cat("======================================================================\n")
cat("\n下一步:\n")
cat("1. 使用整合后的数据运行 MR 分析\n")
cat("2. 脚本位置：run_mr_138genes_complete.R\n")
cat("3. 修改 CLUMPED_FILE 为：%s\n", merged_output)
cat("\n")
