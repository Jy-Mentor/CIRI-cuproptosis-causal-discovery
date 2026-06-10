#!/usr/bin/env Rscript
# ================================================================================
# 准备 GTEx 数据用于 MR 分析（修复版）
# ================================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(readxl)
})

cat("======================================================================\n")
cat("准备 GTEx 数据用于 MR 分析（修复版）\n")
cat("======================================================================\n\n")

# 配置
GTEx_DIR <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL"
OUTPUT_DIR <- "D:/EQTL/clump"
P_THRESHOLD <- 1e-5
R2_THRESHOLD <- 0.01
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

# 提取 gene_name（HGNC 符号）
blood_genes_mapping <- blood_genes[, c("gene_id", "gene_name")]
names(blood_genes_mapping) <- c("ensg_id", "gene_symbol")
blood_genes_mapping <- unique(blood_genes_mapping)
blood_genes_mapping <- blood_genes_mapping[!is.na(blood_genes_mapping$gene_symbol) & blood_genes_mapping$gene_symbol != "", ]

cat(sprintf("  全血 eGenes 映射数：%d\n\n", nrow(blood_genes_mapping)))

# 过滤目标基因的 ENSG ID
target_ensg_blood <- blood_genes_mapping[blood_genes_mapping$gene_symbol %in% target_genes, ]
cat(sprintf("  目标基因在全血中的 ENSG ID 数：%d\n", nrow(target_ensg_blood)))
cat("  前 10 个:\n")
print(head(target_ensg_blood, 10))
cat("\n")

# ================================================================================
# 处理全血 eQTL
# ================================================================================
cat("处理全血 eQTL 数据...\n")

blood_eqtl_file <- file.path(GTEx_DIR, "blood_eqtl.csv")

cat("  加载 eQTL 数据...\n")
blood_eqtl <- fread(blood_eqtl_file, header = TRUE)
cat(sprintf("  原始数据：%d 行\n", nrow(blood_eqtl)))

# 过滤显著 eQTL
cat("  过滤显著 eQTL (P < 1e-5)...\n")
blood_sig <- blood_eqtl[blood_eqtl$pval_nominal < P_THRESHOLD, ]
cat(sprintf("  显著 eQTL 数：%d\n", nrow(blood_sig)))

# 过滤目标基因
cat("  过滤目标基因...\n")
target_ensg_ids <- target_ensg_blood$ensg_id
blood_target <- blood_sig[blood_sig$phenotype_id %in% target_ensg_ids, ]
cat(sprintf("  目标基因 eQTL 数：%d\n", nrow(blood_target)))

if (nrow(blood_target) > 0) {
  # 提取 rsID（从 variant_id 中提取或使用 eGenes 文件的映射）
  cat("  提取 rsID...\n")
  
  # 从 variant_id 提取染色体位置信息
  blood_target$rs_id <- NA
  for (i in 1:nrow(blood_target)) {
    variant_parts <- strsplit(as.character(blood_target$variant_id[i]), "_")[[1]]
    if (length(variant_parts) >= 4) {
      # 格式：chr1_64764_C_T_b38
      chr <- variant_parts[1]
      pos <- variant_parts[2]
      ref <- variant_parts[3]
      alt <- variant_parts[4]
      
      # 在 eGenes 映射中查找
      gene_variant <- blood_genes[blood_genes$gene_id == blood_target$phenotype_id[i] & 
                                   blood_genes$variant_id == blood_target$variant_id[i], ]
      if (nrow(gene_variant) > 0 && !is.na(gene_variant$rs_id_dbSNP157_GRCh38p14) && 
          gene_variant$rs_id_dbSNP157_GRCh38p14 != "") {
        blood_target$rs_id[i] <- gene_variant$rs_id_dbSNP157_GRCh38p14
      } else {
        # 如果没有 rsID，使用 chr_pos_ref_alt 格式
        blood_target$rs_id[i] <- paste0(chr, "_", pos, "_", ref, "_", alt)
      }
    }
  }
  
  cat(sprintf("  有 rsID 的 eQTL 数：%d\n", sum(!is.na(blood_target$rs_id))))
  
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
    exposure = "GTEx_Whole_Blood",
    stringsAsFactors = FALSE
  )
  
  # 过滤掉没有 SNP 标识的
  blood_formatted <- blood_formatted[!is.na(blood_formatted$SNP) & blood_formatted$SNP != "", ]
  cat(sprintf("  格式化后：%d 行\n", nrow(blood_formatted)))
  
  # LD clumping
  if (nrow(blood_formatted) > 1) {
    cat("  进行 LD clumping (r²<0.01, kb=1000)...\n")
    blood_clumped <- clump_data(
      blood_formatted,
      clump_r2 = R2_THRESHOLD,
      clump_kb = KB_WINDOW,
      clump_p = P_THRESHOLD
    )
    cat(sprintf("  Clumping 后：%d 行\n\n", nrow(blood_clumped)))
    
    # 保存
    blood_output <- file.path(OUTPUT_DIR, "GTEx_Whole_Blood_p_1e-05_kb_1000_r2_0.01.tsv")
    write.table(blood_clumped, blood_output, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  保存到：%s\n\n", blood_output))
  } else {
    cat("  SNP 不足，跳过 clumping\n\n")
    blood_clumped <- NULL
  }
} else {
  cat("  无目标基因 eQTL，跳过\n\n")
  blood_clumped <- NULL
}

# ================================================================================
# 处理脑组织 eQTL
# ================================================================================
cat("处理脑组织 eQTL 数据...\n")

brain_genes_file <- file.path(GTEx_DIR, "Brain_Cortex.v11.eGenes.txt")
brain_genes <- fread(brain_genes_file, sep = "\t", header = TRUE)

brain_genes_mapping <- brain_genes[, c("gene_id", "gene_name")]
names(brain_genes_mapping) <- c("ensg_id", "gene_symbol")
brain_genes_mapping <- unique(brain_genes_mapping)
brain_genes_mapping <- brain_genes_mapping[!is.na(brain_genes_mapping$gene_symbol) & brain_genes_mapping$gene_symbol != "", ]

cat(sprintf("  脑组织 eGenes 映射数：%d\n", nrow(brain_genes_mapping)))

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
  cat("  提取 rsID...\n")
  
  brain_target$rs_id <- NA
  for (i in 1:nrow(brain_target)) {
    gene_variant <- brain_genes[brain_genes$gene_id == brain_target$phenotype_id[i] & 
                                 brain_genes$variant_id == brain_target$variant_id[i], ]
    if (nrow(gene_variant) > 0 && !is.na(gene_variant$rs_id_dbSNP157_GRCh38p14) && 
        gene_variant$rs_id_dbSNP157_GRCh38p14 != "") {
      brain_target$rs_id[i] <- gene_variant$rs_id_dbSNP157_GRCh38p14
    } else {
      variant_parts <- strsplit(as.character(brain_target$variant_id[i]), "_")[[1]]
      if (length(variant_parts) >= 4) {
        chr <- variant_parts[1]
        pos <- variant_parts[2]
        ref <- variant_parts[3]
        alt <- variant_parts[4]
        brain_target$rs_id[i] <- paste0(chr, "_", pos, "_", ref, "_", alt)
      }
    }
  }
  
  cat(sprintf("  有 rsID 的 eQTL 数：%d\n", sum(!is.na(brain_target$rs_id))))
  
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
    exposure = "GTEx_Brain_Cortex",
    stringsAsFactors = FALSE
  )
  
  brain_formatted <- brain_formatted[!is.na(brain_formatted$SNP) & brain_formatted$SNP != "", ]
  cat(sprintf("  格式化后：%d 行\n", nrow(brain_formatted)))
  
  if (nrow(brain_formatted) > 1) {
    cat("  进行 LD clumping (r²<0.01, kb=1000)...\n")
    brain_clumped <- clump_data(
      brain_formatted,
      clump_r2 = R2_THRESHOLD,
      clump_kb = KB_WINDOW,
      clump_p = P_THRESHOLD
    )
    cat(sprintf("  Clumping 后：%d 行\n\n", nrow(brain_clumped)))
    
    brain_output <- file.path(OUTPUT_DIR, "GTEx_Brain_Cortex_p_1e-05_kb_1000_r2_0.01.tsv")
    write.table(brain_clumped, brain_output, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  保存到：%s\n\n", brain_output))
  } else {
    cat("  SNP 不足，跳过 clumping\n\n")
    brain_clumped <- NULL
  }
} else {
  cat("  无目标基因 eQTL，跳过\n\n")
  brain_clumped <- NULL
}

# ================================================================================
# 整合所有数据源
# ================================================================================
cat("整合所有 eQTL 数据源...\n")

eqtlgen_file <- file.path(OUTPUT_DIR, "eQTLgen_allgene_p_1e-05_kb_1000_r2_0.01.xlsx")
eqtlgen_data <- read_excel(eqtlgen_file)
cat(sprintf("  eQTLGen: %d SNPs\n", nrow(eqtlgen_data)))

# 读取 GTEx 数据
gtex_sources <- list()
if (!is.null(blood_clumped)) {
  blood_clumped_df <- fread(blood_output, sep = "\t")
  cat(sprintf("  GTEx 全血：%d SNPs\n", nrow(blood_clumped_df)))
  gtex_sources[["blood"]] <- blood_clumped_df
}

if (!is.null(brain_clumped)) {
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
merged_output <- file.path(OUTPUT_DIR, "merged_eqtl_all_sources_p_1e-05_kb_1000_r2_0.01.tsv")
write.table(merged_eqtl, merged_output, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  整合数据保存到：%s\n\n", merged_output))

# 统计每个数据源的贡献
cat("各数据源贡献:\n")
eqtlgen_count <- sum(!is.na(merged_eqtl$exposure) | !grepl("GTEx", merged_eqtl$exposure))
blood_count <- sum(grepl("Whole_Blood", merged_eqtl$exposure), na.rm = TRUE)
brain_count <- sum(grepl("Brain_Cortex", merged_eqtl$exposure), na.rm = TRUE)

cat(sprintf("  eQTLGen: %d SNPs (%.1f%%)\n", eqtlgen_count, eqtlgen_count / nrow(merged_eqtl) * 100))
cat(sprintf("  GTEx 全血：%d SNPs (%.1f%%)\n", blood_count, blood_count / nrow(merged_eqtl) * 100))
cat(sprintf("  GTEx 脑组织：%d SNPs (%.1f%%)\n", brain_count, brain_count / nrow(merged_eqtl) * 100))

cat("\n完成！\n")
