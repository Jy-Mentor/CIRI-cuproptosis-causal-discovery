#!/usr/bin/env Rscript
# ================================================================================
# 138 基因 MR 分析（双源 eQTL 整合）
# 参考权威论文和 GitHub 最佳实践
# - eQTLGen: 已 LD clumped (r² < 0.001, kb = 10000, p < 5e-8) - 直接用
# - GTEx 脑组织：使用所有显著 eQTL（无需 clumping）
# - 完整的敏感性分析（异质性、多效性）
# ================================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
  library(readxl)
  library(arrow)
})

cat("======================================================================\n")
cat("138 基因 MR 分析（双源 eQTL 整合 - eQTLGen + GTEx 脑组织）\n")
cat("======================================================================\n\n")

# 配置（已确认文件路径正确）
CLUMPED_FILE <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-8_kb_10000_r2_0.001.xlsx"  # ✓ 已验证
gtex_base <- "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL"
GTEx_EGENES_FILE <- file.path(gtex_base, "Brain_Cortex.v11.eGenes.txt")
GTEx_PAIRS_FILE <- file.path(gtex_base, "Brain_Cortex.v11.eQTLs.signif_pairs.parquet")
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_dual_source"
MEGASTROKE_FILE <- "D:/下载/29531354-GCST006906-EFO_0000712.h.tsv.gz"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", OUTPUT_DIR, "\n\n")

# 目标基因列表
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

# 基因符号到 ENSG ID 的映射
gene_to_ensg <- list(
  "ACAD11" = "ENSG00000240303", "ACADM" = "ENSG00000117054", "ACADVL" = "ENSG00000072778",
  "ACTA2" = "ENSG00000107796", "ADRB1" = "ENSG00000043591", "AIF1" = "ENSG00000204472",
  "AKT1" = "ENSG00000142208", "ALDH9A1" = "ENSG00000143149", "ATOX1" = "ENSG00000177556",
  "ATP7A" = "ENSG00000165240", "ATP7B" = "ENSG00000123191", "B2M" = "ENSG00000166710",
  "BRD3" = "ENSG00000169925", "BST1" = "ENSG00000197121", "C3" = "ENSG00000125730",
  "CASK" = "ENSG00000147044", "CASP8" = "ENSG00000064012", "CAT" = "ENSG00000121691",
  "CCL2" = "ENSG00000108691", "CCNA2" = "ENSG00000145386", "CCND1" = "ENSG00000110092",
  "CCR5" = "ENSG00000160791", "CDK4" = "ENSG00000135446", "CHFR" = "ENSG00000072609",
  "CITED2" = "ENSG00000164442", "CNDP2" = "ENSG00000133313", "CNR2" = "ENSG00000188822",
  "COL1A1" = "ENSG00000108821", "CP" = "ENSG00000047457", "CPT1A" = "ENSG00000110090",
  "CPT2" = "ENSG00000157184", "CTSB" = "ENSG00000164733", "CTSC" = "ENSG00000109861",
  "CTSD" = "ENSG00000117984", "CTSF" = "ENSG00000174080", "CTSK" = "ENSG00000143387",
  "CTSL" = "ENSG00000135047", "CTSS" = "ENSG00000163131", "CUL4B" = "ENSG00000158290",
  "CXCR3" = "ENSG00000186810", "DDC" = "ENSG00000132437", "DLAT" = "ENSG00000150768",
  "EGFR" = "ENSG00000146648", "EPHX1" = "ENSG00000143819", "F3" = "ENSG00000117525",
  "FABP2" = "ENSG00000145384", "FABP4" = "ENSG00000170323", "FABP5" = "ENSG00000236972",
  "FDX1" = "ENSG00000137714", "FLT4" = "ENSG00000037280", "FNTA" = "ENSG00000168522",
  "GCH1" = "ENSG00000131979", "GFAP" = "ENSG00000131095", "GPX1" = "ENSG00000233276",
  "GPX4" = "ENSG00000167468", "HBS1L" = "ENSG00000112339", "HIBADH" = "ENSG00000106049",
  "HIF1A" = "ENSG00000100644", "HMOX1" = "ENSG00000100292", "HPGDS" = "ENSG00000163106",
  "HSD17B10" = "ENSG00000072506", "HSD17B4" = "ENSG00000133835", "HSPA5" = "ENSG00000044574",
  "HTR2B" = "ENSG00000135914", "HTR2C" = "ENSG00000147246", "ICAM1" = "ENSG00000090339",
  "IGFBP2" = "ENSG00000115457", "IKBKB" = "ENSG00000104365", "IL10RA" = "ENSG00000110324",
  "IL6" = "ENSG00000136244", "IMPDH2" = "ENSG00000178035", "IRF1" = "ENSG00000125347",
  "ITGA1" = "ENSG00000213949", "JAK1" = "ENSG00000162434", "KCNA5" = "ENSG00000130037",
  "LEF1" = "ENSG00000138795", "LIAS" = "ENSG00000121897", "LIPT1" = "ENSG00000144182",
  "LYN" = "ENSG00000254087", "MAN2B1" = "ENSG00000104774", "MAOB" = "ENSG00000069535",
  "MAPKAPK2" = "ENSG00000162889", "MB" = "ENSG00000198125", "MGAT1" = "ENSG00000131446",
  "MKNK2" = "ENSG00000099875", "MTOR" = "ENSG00000198793", "NFE2L2" = "ENSG00000116044",
  "NFKB1" = "ENSG00000109320", "NMT1" = "ENSG00000136448", "NR1H3" = "ENSG00000025434",
  "NR3C1" = "ENSG00000113580", "NUDCD2" = "ENSG00000170584", "OAZ1" = "ENSG00000104904",
  "PA2G4" = "ENSG00000170515", "PABPC1" = "ENSG00000070756", "PARP1" = "ENSG00000143799",
  "PARP12" = "ENSG00000059378", "PCTP" = "ENSG00000141179", "PDCD6" = "ENSG00000249915",
  "PDCD6IP" = "ENSG00000170248", "PDHB" = "ENSG00000168291", "PDHX" = "ENSG00000110435",
  "PLA2G4A" = "ENSG00000116711", "POLR2D" = "ENSG00000144231", "PPARG" = "ENSG00000132170",
  "PRKCQ" = "ENSG00000065675", "PTGR1" = "ENSG00000106853", "PTGS1" = "ENSG00000095303",
  "PTPN2" = "ENSG00000175354", "PTPN6" = "ENSG00000111679", "PTPRC" = "ENSG00000081237",
  "PTPRF" = "ENSG00000142949", "PTPRJ" = "ENSG00000149177", "RBM39" = "ENSG00000131051",
  "RELA" = "ENSG00000173039", "RENBP" = "ENSG00000102032", "RHOC" = "ENSG00000155366",
  "S100A6" = "ENSG00000197956", "SAT1" = "ENSG00000111371", "SAT2" = "ENSG00000134294",
  "SCN9A" = "ENSG00000169432", "SEC13" = "ENSG00000157020", "SERPINB1" = "ENSG00000197696",
  "SERPINB10" = "ENSG00000242550", "SLC31A1" = "ENSG00000136868", "SPHK1" = "ENSG00000176170",
  "SREBF1" = "ENSG00000072310", "STARD13" = "ENSG00000133121", "STAT1" = "ENSG00000115415",
  "STAT3" = "ENSG00000168610", "STAT5A" = "ENSG00000126561", "STK4" = "ENSG00000101109",
  "TBXAS1" = "ENSG00000059377", "TCN2" = "ENSG00000185339", "TDP1" = "ENSG00000042088",
  "TGFB1" = "ENSG00000105329", "TIMP1" = "ENSG00000102265", "TNF" = "ENSG00000232810",
  "TOP2A" = "ENSG00000131747", "TSPO" = "ENSG00000100300", "XDH" = "ENSG00000158125",
  "XRCC6" = "ENSG00000196419", "ZEB1" = "ENSG00000148516", "ZHX2" = "ENSG00000178764"
)

# 步骤 1: 加载 eQTLGen 数据（已 LD clumped）
cat("步骤 1: 加载 eQTLGen 数据（已 LD clumped）\n")
cat("----------------------------------------------------------------------\n")

if (!file.exists(CLUMPED_FILE)) {
  stop("eQTLGen clumped 文件不存在：", CLUMPED_FILE)
}

cat("  读取 eQTLGen clumped 数据...\n")
cat("  文件：", basename(CLUMPED_FILE), "\n")
cat("  参数：P < 5e-08, KB = 10000, r² < 0.001（严格去 LD）\n\n")

eqtlgen_data <- read_excel(CLUMPED_FILE)
cat(sprintf("  ✓ 成功读取 %d 个 SNP\n\n", nrow(eqtlgen_data)))

# 步骤 2: 加载 GTEx 脑组织数据（使用所有显著 eQTL，无需 clumping）
cat("步骤 2: 加载 GTEx 脑组织 eQTL 数据（所有显著 eQTL）\n")
cat("----------------------------------------------------------------------\n")

cat("  读取 GTEx eGenes 数据...\n")
gtex_egenes <- fread(GTEx_EGENES_FILE, stringsAsFactors = FALSE)
cat(sprintf("  ✓ 成功读取 %d 个 eGenes\n", nrow(gtex_egenes)))

# 读取所有显著 eQTL 对（无需 clumping，直接用所有 p<5e-8 的 SNP）
cat("  读取 GTEx 显著 eQTL 对...\n")
gtex_pairs <- read_parquet(GTEx_PAIRS_FILE)
cat(sprintf("  ✓ 成功读取 %d 个 eQTL 对\n", nrow(gtex_pairs)))

# 筛选显著的 eQTL（p < 5e-8）
gtex_signif <- gtex_pairs[gtex_pairs$pval_nominal < 5e-08, ]
cat(sprintf("  ✓ 筛选显著 eQTL (P < 5e-08): %d 个\n\n", nrow(gtex_signif)))

# 步骤 3: 加载 MEGASTROKE 数据
cat("步骤 3: 加载 MEGASTROKE 数据\n")
cat("----------------------------------------------------------------------\n")

if (!file.exists(MEGASTROKE_FILE)) {
  stop("MEGASTROKE 文件不存在：", MEGASTROKE_FILE)
}

cat("  读取 MEGASTROKE 协调化数据...\n")
outcome_data <- fread(MEGASTROKE_FILE, sep = "\t", stringsAsFactors = FALSE, showProgress = TRUE)
cat(sprintf("  ✓ 成功读取 %d 个 SNP\n\n", nrow(outcome_data)))

# 重命名以匹配 TwoSampleMR 格式
col_mapping <- list(
  variant_id = "SNP", chromosome = "chr", base_pair_location = "pos.outcome",
  beta = "beta.outcome", standard_error = "se.outcome",
  effect_allele = "effect_allele.outcome", other_allele = "other_allele.outcome",
  effect_allele_frequency = "eaf.outcome", `p-value` = "pval.outcome"
)

for (old_name in names(col_mapping)) {
  if (old_name %in% names(outcome_data)) {
    setnames(outcome_data, old_name, col_mapping[[old_name]])
  }
}

outcome_data$phenotype <- "MEGASTROKE_Stroke"
outcome_data$id.outcome <- "MEGASTROKE_Stroke"
outcome_data$outcome <- "MEGASTROKE_Stroke"
outcome_data$samplesize.outcome <- 520000
outcome_data$ncase.outcome <- 67162
outcome_data$ncontrol.outcome <- 454450

# 步骤 4: 准备每个基因的双源暴露数据
cat("步骤 4: 准备每个基因的双源暴露数据\n")
cat("----------------------------------------------------------------------\n")

prepared_genes <- list()

for (gene_symbol in target_genes) {
  ensg_id <- gene_to_ensg[[gene_symbol]]
  
  if (is.null(ensg_id)) {
    cat(sprintf("%s: 无 ENSG ID，跳过\n", gene_symbol))
    next
  }
  
  # 从 eQTLGen 获取 SNPs（已 clumped）
  eqtlgen_snps <- eqtlgen_data[eqtlgen_data$gene == gene_symbol | 
                               eqtlgen_data$gene == ensg_id, ]
  
  # 从 GTEx 获取所有显著 eQTL（无需 clumping）
  # 注意：parquet 文件中 phenotype_id 就是 gene_id
  gtex_snps <- gtex_signif[gtex_signif$phenotype_id == ensg_id, ]
  
  # 合并双源数据
  if (nrow(eqtlgen_snps) > 0) {
    # eQTLGen 数据（已经是 TwoSampleMR 格式）
    eqtlgen_df <- data.frame(
      SNP = eqtlgen_snps$SNP,
      BETA = eqtlgen_snps$beta.exposure,
      SE = eqtlgen_snps$se.exposure,
      EFFECT_ALLELE = eqtlgen_snps$effect_allele.exposure,
      OTHER_ALLELE = eqtlgen_snps$other_allele.exposure,
      EAF = eqtlgen_snps$eaf.exposure,
      PVAL = eqtlgen_snps$pval.exposure,
      CHR = eqtlgen_snps$chr.exposure,
      BP = eqtlgen_snps$pos.exposure,
      TISSUE = "Whole_Blood",
      stringsAsFactors = FALSE
    )
    
    if (nrow(gtex_snps) > 0) {
      # GTEx 数据（需要转换格式）
      # variant_id 格式：chr_pos_ref_alt_b38
      gtex_df <- data.frame(
        SNP = gtex_snps$variant_id,
        BETA = gtex_snps$slope,
        SE = gtex_snps$slope_se,
        EFFECT_ALLELE = sapply(strsplit(gtex_snps$variant_id, "_"), function(x) x[4]),
        OTHER_ALLELE = sapply(strsplit(gtex_snps$variant_id, "_"), function(x) x[3]),
        EAF = gtex_snps$af,
        PVAL = gtex_snps$pval_nominal,
        CHR = as.numeric(gsub("chr", "", sapply(strsplit(gtex_snps$variant_id, "_"), function(x) x[1]))),
        BP = as.numeric(sapply(strsplit(gtex_snps$variant_id, "_"), function(x) x[2])),
        TISSUE = "Brain_Cortex",
        stringsAsFactors = FALSE
      )
      
      # 合并双源
      exposure <- rbind(eqtlgen_df, gtex_df)
    } else {
      exposure <- eqtlgen_df
    }
  } else if (nrow(gtex_snps) > 0) {
    # 只有 GTEx 数据
    gtex_df <- data.frame(
      SNP = gtex_snps$variant_id,
      BETA = gtex_snps$slope,
      SE = gtex_snps$slope_se,
      EFFECT_ALLELE = sapply(strsplit(gtex_snps$variant_id, "_"), function(x) x[4]),
      OTHER_ALLELE = sapply(strsplit(gtex_snps$variant_id, "_"), function(x) x[3]),
      EAF = gtex_snps$af,
      PVAL = gtex_snps$pval_nominal,
      CHR = as.numeric(gsub("chr", "", sapply(strsplit(gtex_snps$variant_id, "_"), function(x) x[1]))),
      BP = as.numeric(sapply(strsplit(gtex_snps$variant_id, "_"), function(x) x[2])),
      TISSUE = "Brain_Cortex",
      stringsAsFactors = FALSE
    )
    exposure <- gtex_df
  } else {
    cat(sprintf("%s: 无 eQTL 数据\n", gene_symbol))
    next
  }
  
  if (nrow(exposure) < 2) {
    cat(sprintf("%s: SNP 数量不足 (%d)\n", gene_symbol, nrow(exposure)))
    next
  }
  
  # 格式化暴露数据
  tryCatch({
    exposure_fmt <- format_data(
      exposure, type = "exposure", snp_col = "SNP", beta_col = "BETA", se_col = "SE",
      eaf_col = "EAF", effect_allele_col = "EFFECT_ALLELE", other_allele_col = "OTHER_ALLELE",
      pval_col = "PVAL"
    )
    
    # 计算 F 统计量（每个 SNP 的 F 统计量平均值）
    f_stat <- mean((exposure_fmt$beta.exposure / exposure_fmt$se.exposure)^2, na.rm = TRUE)
    
    # 组织分布
    tissue_dist <- table(exposure$TISSUE)
    n_blood <- ifelse("Whole_Blood" %in% names(tissue_dist), sum(tissue_dist["Whole_Blood"]), 0)
    n_brain <- ifelse("Brain_Cortex" %in% names(tissue_dist), sum(tissue_dist["Brain_Cortex"]), 0)
    
    prepared_genes[[gene_symbol]] <- list(
      exposure = exposure, exposure_fmt = exposure_fmt, ensg_id = ensg_id, 
      f_stat = f_stat, n_blood = n_blood, n_brain = n_brain
    )
    
    cat(sprintf("%s: ✓ 准备就绪 (%d SNPs: %d 血，%d 脑，F=%.2f)\n", 
                gene_symbol, nrow(exposure), n_blood, n_brain, f_stat))
  }, error = function(e) {
    cat(sprintf("%s: 格式化失败 - %s\n", gene_symbol, e$message))
  })
}

cat(sprintf("\n\n成功准备：%d/%d 个基因\n", length(prepared_genes), length(target_genes)))

# 步骤 5: MR 分析（包含敏感性分析）
cat("\n步骤 5: MR 分析（包含敏感性分析）\n")
cat("----------------------------------------------------------------------\n")

mr_results <- list()
harmonised_data <- list()

for (gene_symbol in names(prepared_genes)) {
  cat(sprintf("\n%s:\n", gene_symbol))
  
  exposure <- prepared_genes[[gene_symbol]]$exposure
  exposure_fmt <- prepared_genes[[gene_symbol]]$exposure_fmt
  f_stat <- prepared_genes[[gene_symbol]]$f_stat
  n_blood <- prepared_genes[[gene_symbol]]$n_blood
  n_brain <- prepared_genes[[gene_symbol]]$n_brain
  
  chr <- as.character(unique(exposure$CHR)[1])
  min_bp <- min(exposure$BP, na.rm = TRUE) - 5e5
  max_bp <- max(exposure$BP, na.rm = TRUE) + 5e5
  
  cat(sprintf("  染色体：%s\n", chr))
  cat(sprintf("  位置范围：%d-%d\n", min_bp, max_bp))
  cat(sprintf("  F 统计量：%.2f\n", f_stat))
  cat(sprintf("  组织来源：%d 血，%d 脑\n", n_blood, n_brain))
  
  tryCatch({
    # 筛选区域数据
    region_data <- outcome_data[outcome_data$chr == chr & 
                                outcome_data$pos.outcome >= min_bp & 
                                outcome_data$pos.outcome <= max_bp, ]
    
    if (nrow(region_data) == 0) {
      cat(sprintf("  ✗ 该区域无 MEGASTROKE 数据\n"))
      next
    }
    
    # 按 SNP ID 匹配
    matched_snps <- intersect(exposure_fmt$SNP, region_data$SNP)
    
    if (length(matched_snps) < 2) {
      cat(sprintf("  ✗ 匹配的 SNP 太少 (%d)\n", length(matched_snps)))
      next
    }
    
    cat(sprintf("  匹配的 SNP：%d\n", length(matched_snps)))
    
    # 提取匹配的数据
    exp_matched <- exposure_fmt[exposure_fmt$SNP %in% matched_snps, ]
    out_matched <- region_data[region_data$SNP %in% matched_snps, ]
    
    # 排序
    exp_matched <- exp_matched[order(exp_matched$SNP), ]
    out_matched <- out_matched[order(out_matched$SNP), ]
    
    # Harmonise
    dat <- harmonise_data(exp_matched, out_matched)
    
    if (nrow(dat) < 2) {
      cat(sprintf("  ✗ Harmonise 后 SNP 不足 (%d)\n", nrow(dat)))
      next
    }
    
    harmonised_data[[gene_symbol]] <- dat
    
    # MR 分析（多种方法）
    res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
    
    if (nrow(res) == 0) {
      cat(sprintf("  ✗ MR 分析失败\n"))
      next
    }
    
    # 提取 IVW 结果
    ivw_res <- res[res$method == "Inverse variance weighted", ]
    
    # 异质性检查
    het_res <- mr_heterogeneity(dat)
    
    # 多效性检查
    if (nrow(dat) >= 3) {
      egger_intercept <- mr_pleiotropy_test(dat)
    } else {
      egger_intercept <- data.frame(intercept = NA, pval = NA)
    }
    
    # 保存结果
    result <- data.frame(
      gene = gene_symbol,
      method = ivw_res$method[1],
      beta = ivw_res$b[1],
      se = ivw_res$se[1],
      or = exp(ivw_res$b[1]),
      ci_low = exp(ivw_res$b[1] - 1.96 * ivw_res$se[1]),
      ci_high = exp(ivw_res$b[1] + 1.96 * ivw_res$se[1]),
      pval = ivw_res$pval[1],
      f_stat = f_stat,
      nsnp = nrow(dat),
      n_blood = n_blood,
      n_brain = n_brain,
      # 敏感性分析结果
      heterogeneity_q = ifelse(nrow(het_res) > 0, het_res$Q[1], NA),
      heterogeneity_pval = ifelse(nrow(het_res) > 0, het_res$Q_pval[1], NA),
      egger_intercept = ifelse(nrow(egger_intercept) > 0, egger_intercept$intercept[1], NA),
      egger_pval = ifelse(nrow(egger_intercept) > 0, egger_intercept$pval[1], NA),
      # Weighted Median
      wm_beta = ifelse(nrow(res[res$method == "Weighted median", ]) > 0,
                       res[res$method == "Weighted median", "b"][1], NA),
      wm_pval = ifelse(nrow(res[res$method == "Weighted median", ]) > 0,
                       res[res$method == "Weighted median", "pval"][1], NA),
      stringsAsFactors = FALSE
    )
    
    mr_results[[gene_symbol]] <- result
    
    direction_consistent <- sign(ivw_res$b[1]) == sign(result$wm_beta)
    
    cat(sprintf("  ✓ MR: OR=%.3f (%.3f-%.3f), P=%.2e, F=%.2f\n", 
                result$or, result$ci_low, result$ci_high, result$pval, result$f_stat))
    cat(sprintf("    异质性：Q=%.2f, P=%.3f\n", result$heterogeneity_q, result$heterogeneity_pval))
    cat(sprintf("    多效性：intercept=%.3f, P=%.3f\n", result$egger_intercept, result$egger_pval))
    cat(sprintf("    稳健性：WM P=%.3f, 方向一致：%s\n", result$wm_pval, ifelse(direction_consistent, "是", "否")))
    
  }, error = function(e) {
    cat(sprintf("  ✗ 错误：%s\n", e$message))
  })
}

# 步骤 6: 保存结果
cat("\n\n步骤 6: 保存结果\n")
cat("----------------------------------------------------------------------\n")

if (length(mr_results) > 0) {
  results_df <- do.call(rbind, mr_results)
  results_df <- results_df[order(results_df$pval), ]
  
  output_file <- file.path(OUTPUT_DIR, "mr_results_all_genes.csv")
  write.csv(results_df, output_file, row.names = FALSE)
  cat(sprintf("✓ 保存 MR 结果：%s\n", output_file))
  
  cat(sprintf("\nTop 显著基因 (P < 0.05):\n"))
  sig_results <- results_df[results_df$pval < 0.05, ]
  if (nrow(sig_results) > 0) {
    print(sig_results[, c("gene", "or", "ci_low", "ci_high", "pval", "nsnp", "f_stat", 
                          "heterogeneity_pval", "egger_pval", "n_blood", "n_brain")])
  } else {
    cat("  无显著基因\n")
  }
  
  results_df$pval_adj <- p.adjust(results_df$pval, method = "fdr")
  write.csv(results_df, file.path(OUTPUT_DIR, "mr_results_fdr_corrected.csv"), row.names = FALSE)
  cat(sprintf("✓ 保存 FDR 校正结果\n"))
  
  if (length(harmonised_data) > 0) {
    all_harmonised <- do.call(rbind, lapply(names(harmonised_data), function(gene) {
      dat <- harmonised_data[[gene]]
      dat$gene <- gene
      dat
    }))
    write.csv(all_harmonised, file.path(OUTPUT_DIR, "harmonised_data_all_genes.csv"), row.names = FALSE)
    cat(sprintf("✓ 保存 Harmonised 数据\n"))
  }
  
  # 生成质量报告
  cat(sprintf("\n生成质量报告...\n"))
  quality_report <- file.path(OUTPUT_DIR, "quality_control_report.txt")
  sink(quality_report)
  cat("MR 分析质量控制报告（双源 eQTL）\n")
  cat("=============================\n\n")
  cat("1. F 统计量评估\n")
  cat("   平均 F 统计量：", round(mean(results_df$f_stat, na.rm = TRUE), 2), "\n")
  cat("   F > 10 的基因数：", sum(results_df$f_stat > 10, na.rm = TRUE), "/", nrow(results_df), "\n\n")
  
  cat("2. 异质性评估\n")
  cat("   显著异质性 (P < 0.05) 的基因数：", sum(results_df$heterogeneity_pval < 0.05, na.rm = TRUE), "\n")
  cat("   极显著异质性 (P < 1e-10) 的基因：", 
      paste(results_df$gene[results_df$heterogeneity_pval < 1e-10], collapse = ", "), "\n\n")
  
  cat("3. 多效性评估\n")
  cat("   显著多效性 (P < 0.05) 的基因数：", sum(results_df$egger_pval < 0.05, na.rm = TRUE), "\n\n")
  
  cat("4. 稳健性检查\n")
  cat("   方向一致的基因数：", sum(sign(results_df$beta) == sign(results_df$wm_beta), na.rm = TRUE), "\n\n")
  
  cat("5. 双源整合统计\n")
  cat("   同时有血和脑数据的基因数：", sum(results_df$n_blood > 0 & results_df$n_brain > 0), "\n")
  cat("   只有血数据的基因数：", sum(results_df$n_blood > 0 & results_df$n_brain == 0), "\n")
  cat("   只有脑数据的基因数：", sum(results_df$n_blood == 0 & results_df$n_brain > 0), "\n\n")
  
  sink()
  cat(sprintf("✓ 保存质量报告：%s\n", quality_report))
}

cat("\n======================================================================\n")
cat("分析完成！\n")
cat("======================================================================\n")
