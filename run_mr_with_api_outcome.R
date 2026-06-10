#!/usr/bin/env Rscript
# ================================================================================
# 为 138 个目标基因运行 MR 分析
# 使用 IEU OpenGWAS API 动态提取 MEGASTROKE 结局数据
# ================================================================================

library(dplyr)
library(data.table)
library(ieugwasr)
library(TwoSampleMR)

cat("======================================================================\n")
cat("138 基因 MR 分析（使用 IEU API 提取结局数据）\n")
cat("======================================================================\n\n")

# 目标基因列表（138 个，去重后）
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

# 基因符号到 ENSG ID 的映射（通过 Ensembl API 查询，143/144 个基因）
# 仅 SERPINB1A 未找到
gene_to_ensg <- list(
  "ACAD11" = "ENSG00000240303",
  "ACADM" = "ENSG00000117054",
  "ACADVL" = "ENSG00000072778",
  "ACTA2" = "ENSG00000107796",
  "ADRB1" = "ENSG00000043591",
  "AIF1" = "ENSG00000204472",
  "AKT1" = "ENSG00000142208",
  "ALDH9A1" = "ENSG00000143149",
  "ATOX1" = "ENSG00000177556",
  "ATP7A" = "ENSG00000165240",
  "ATP7B" = "ENSG00000123191",
  "B2M" = "ENSG00000166710",
  "BRD3" = "ENSG00000169925",
  "BST1" = "ENSG00000197121",
  "C3" = "ENSG00000125730",
  "CASK" = "ENSG00000147044",
  "CASP8" = "ENSG00000064012",
  "CAT" = "ENSG00000121691",
  "CCL2" = "ENSG00000108691",
  "CCNA2" = "ENSG00000145386",
  "CCND1" = "ENSG00000110092",
  "CCR5" = "ENSG00000160791",
  "CDK4" = "ENSG00000135446",
  "CHFR" = "ENSG00000072609",
  "CITED2" = "ENSG00000164442",
  "CNDP2" = "ENSG00000133313",
  "CNR2" = "ENSG00000188822",
  "COL1A1" = "ENSG00000108821",
  "CP" = "ENSG00000047457",
  "CPT1A" = "ENSG00000110090",
  "CPT2" = "ENSG00000157184",
  "CTSB" = "ENSG00000164733",
  "CTSC" = "ENSG00000109861",
  "CTSD" = "ENSG00000117984",
  "CTSF" = "ENSG00000174080",
  "CTSK" = "ENSG00000143387",
  "CTSL" = "ENSG00000135047",
  "CTSS" = "ENSG00000163131",
  "CUL4B" = "ENSG00000158290",
  "CXCR3" = "ENSG00000186810",
  "DDC" = "ENSG00000132437",
  "DLAT" = "ENSG00000150768",
  "EGFR" = "ENSG00000146648",
  "EPHX1" = "ENSG00000143819",
  "F3" = "ENSG00000117525",
  "FABP2" = "ENSG00000145384",
  "FABP4" = "ENSG00000170323",
  "FABP5" = "ENSG00000236972",
  "FDX1" = "ENSG00000137714",
  "FLT4" = "ENSG00000037280",
  "FNTA" = "ENSG00000168522",
  "GCH1" = "ENSG00000131979",
  "GFAP" = "ENSG00000131095",
  "GPX1" = "ENSG00000233276",
  "GPX4" = "ENSG00000167468",
  "HBS1L" = "ENSG00000112339",
  "HIBADH" = "ENSG00000106049",
  "HIF1A" = "ENSG00000100644",
  "HMOX1" = "ENSG00000100292",
  "HPGDS" = "ENSG00000163106",
  "HSD17B10" = "ENSG00000072506",
  "HSD17B4" = "ENSG00000133835",
  "HSPA5" = "ENSG00000044574",
  "HTR2B" = "ENSG00000135914",
  "HTR2C" = "ENSG00000147246",
  "ICAM1" = "ENSG00000090339",
  "IGFBP2" = "ENSG00000115457",
  "IKBKB" = "ENSG00000104365",
  "IL10RA" = "ENSG00000110324",
  "IL6" = "ENSG00000136244",
  "IMPDH2" = "ENSG00000178035",
  "IRF1" = "ENSG00000125347",
  "ITGA1" = "ENSG00000213949",
  "JAK1" = "ENSG00000162434",
  "KCNA5" = "ENSG00000130037",
  "LEF1" = "ENSG00000138795",
  "LIAS" = "ENSG00000121897",
  "LIPT1" = "ENSG00000144182",
  "LYN" = "ENSG00000254087",
  "MAN2B1" = "ENSG00000104774",
  "MAOB" = "ENSG00000069535",
  "MAPKAPK2" = "ENSG00000162889",
  "MB" = "ENSG00000198125",
  "MGAT1" = "ENSG00000131446",
  "MKNK2" = "ENSG00000099875",
  "MTOR" = "ENSG00000198793",
  "NFE2L2" = "ENSG00000116044",
  "NFKB1" = "ENSG00000109320",
  "NMT1" = "ENSG00000136448",
  "NR1H3" = "ENSG00000025434",
  "NR3C1" = "ENSG00000113580",
  "NUDCD2" = "ENSG00000170584",
  "OAZ1" = "ENSG00000104904",
  "PA2G4" = "ENSG00000170515",
  "PABPC1" = "ENSG00000070756",
  "PARP1" = "ENSG00000143799",
  "PARP12" = "ENSG00000059378",
  "PCTP" = "ENSG00000141179",
  "PDCD6" = "ENSG00000249915",
  "PDCD6IP" = "ENSG00000170248",
  "PDHB" = "ENSG00000168291",
  "PDHX" = "ENSG00000110435",
  "PLA2G4A" = "ENSG00000116711",
  "POLR2D" = "ENSG00000144231",
  "PPARG" = "ENSG00000132170",
  "PRKCQ" = "ENSG00000065675",
  "PTGR1" = "ENSG00000106853",
  "PTGS1" = "ENSG00000095303",
  "PTPN2" = "ENSG00000175354",
  "PTPN6" = "ENSG00000111679",
  "PTPRC" = "ENSG00000081237",
  "PTPRF" = "ENSG00000142949",
  "PTPRJ" = "ENSG00000149177",
  "RBM39" = "ENSG00000131051",
  "RELA" = "ENSG00000173039",
  "RENBP" = "ENSG00000102032",
  "RHOC" = "ENSG00000155366",
  "S100A6" = "ENSG00000197956",
  "SAT1" = "ENSG00000111371",
  "SAT2" = "ENSG00000134294",
  "SCN9A" = "ENSG00000169432",
  "SEC13" = "ENSG00000157020",
  "SERPINB1" = "ENSG00000197696",
  "SERPINB10" = "ENSG00000242550",
  "SLC31A1" = "ENSG00000136868",
  "SPHK1" = "ENSG00000176170",
  "SREBF1" = "ENSG00000072310",
  "STARD13" = "ENSG00000133121",
  "STAT1" = "ENSG00000115415",
  "STAT3" = "ENSG00000168610",
  "STAT5A" = "ENSG00000126561",
  "STK4" = "ENSG00000101109",
  "TBXAS1" = "ENSG00000059377",
  "TCN2" = "ENSG00000185339",
  "TDP1" = "ENSG00000042088",
  "TGFB1" = "ENSG00000105329",
  "TIMP1" = "ENSG00000102265",
  "TNF" = "ENSG00000232810",
  "TOP2A" = "ENSG00000131747",
  "TSPO" = "ENSG00000100300",
  "XDH" = "ENSG00000158125",
  "XRCC6" = "ENSG00000196419",
  "ZEB1" = "ENSG00000148516",
  "ZHX2" = "ENSG00000178764"
)

# 配置
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_api"
EXPOSURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_matched/matched_data"
MEGASTROKE_ID <- "ieu-a-798"  # MEGASTROKE 在 IEU OpenGWAS 的 ID
JWT_TOKEN <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3ODMwNDA4MywiZXhwIjoxNzc5NTEzNjgzfQ.ZtcIUEx_xYtrVD_EE-UboKyLlC-lZBq2pjn-iYhJzxocqHdA-02K9n_Qbw-5ngQ07GHjjIYqVtmkZfJ3OJl1yI-tOMBBFzVKe0nkwDcB6-yBgjgBaxVm8vq_pbNrMwy_ZezY5ys9jx7I8T4bZYg9KeUbSwj04OfNP82kGcKXIOErXXVy-Ie3dbUogDRSjnCT-_32yNQxuWpiyYnPWSrWQbQ2HlUQiiDTdFGzWeJJKfSRvjQzdp5g3nccxht0m5A0UsPCdvkyHFEpvPVZ-NpjCjkgy8GbZBv4cmDMSc5JJL6HLO0eV508SRKxMdp-gL6qVdhGiJ2i9XmZQE27-aq_7g"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", OUTPUT_DIR, "\n\n")

# 步骤 1: 读取暴露数据
cat("步骤 1: 读取暴露数据\n")
cat("----------------------------------------------------------------------\n")

exposure_files <- list.files(EXPOSURE_DIR, pattern = "_exposure\\.csv$", full.names = TRUE)
cat(sprintf("暴露数据文件数：%d\n", length(exposure_files)))

# 建立 ENSG ID 到文件的映射
ensg_to_file <- list()
for (exposure_file in exposure_files) {
  ensg_id <- basename(exposure_file)
  ensg_base <- sub("\\..*", "", ensg_id)
  ensg_to_file[[ensg_base]] <- exposure_file
}

cat(sprintf("建立 ENSG ID 映射：%d 个文件\n\n", length(ensg_to_file)))

# 步骤 2: 准备 MR 分析
cat("步骤 2: 准备 MR 分析\n")
cat("----------------------------------------------------------------------\n")

prepared_genes <- list()
missing_exposure <- character()

for (gene_symbol in target_genes) {
  ensg_id <- gene_to_ensg[[gene_symbol]]
  
  if (is.null(ensg_id)) {
    cat(sprintf("%s: 无 ENSG ID，跳过\n", gene_symbol))
    next
  }
  
  # 查找暴露文件
  file <- ensg_to_file[[ensg_id]]
  
  if (is.null(file) || !file.exists(file)) {
    cat(sprintf("%s (%s): 无暴露数据\n", gene_symbol, ensg_id))
    missing_exposure <- c(missing_exposure, gene_symbol)
    next
  }
  
  # 读取暴露数据
  exposure <- fread(file, stringsAsFactors = FALSE)
  
  if (nrow(exposure) < 3) {
    cat(sprintf("%s: SNP 数量不足 (%d)\n", gene_symbol, nrow(exposure)))
    next
  }
  
  # 转换为 data.frame（format_data 要求）
  exposure <- as.data.frame(exposure)
  
  # 格式化暴露数据
  tryCatch({
    exposure_fmt <- format_data(
      exposure,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "EAF",
      effect_allele_col = "EFFECT_ALLELE",
      other_allele_col = "OTHER_ALLELE",
      pval_col = "PVAL"
    )
    
    prepared_genes[[gene_symbol]] <- list(
      exposure = exposure,
      exposure_fmt = exposure_fmt,
      ensg_id = ensg_id
    )
    
    cat(sprintf("%s: ✓ 准备就绪 (%d SNPs)\n", gene_symbol, nrow(exposure)))
  }, error = function(e) {
    cat(sprintf("%s: 格式化失败 - %s\n", gene_symbol, e$message))
  })
}

cat(sprintf("\n\n成功准备：%d/%d 个基因\n", length(prepared_genes), length(target_genes)))
cat(sprintf("缺失暴露：%d 个基因\n\n", length(missing_exposure)))

# 保存缺失暴露的基因列表
if (length(missing_exposure) > 0) {
  writeLines(missing_exposure, file.path(OUTPUT_DIR, "missing_exposure_genes.txt"))
}

# 步骤 3: MR 分析（使用 IEU API 提取结局数据）
cat("步骤 3: MR 分析（使用 IEU API 提取结局数据）\n")
cat("----------------------------------------------------------------------\n")

mr_results <- list()
harmonised_data <- list()

for (gene_symbol in names(prepared_genes)) {
  cat(sprintf("\n%s:\n", gene_symbol))
  
  exposure <- prepared_genes[[gene_symbol]]$exposure
  exposure_fmt <- prepared_genes[[gene_symbol]]$exposure_fmt
  
  # 获取该基因的染色体和位置范围
  chr <- as.character(unique(exposure$CHR)[1])
  min_bp <- min(exposure$BP, na.rm = TRUE) - 5e5  # 扩展 500kb
  max_bp <- max(exposure$BP, na.rm = TRUE) + 5e5
  
  cat(sprintf("  染色体：%s\n", chr))
  cat(sprintf("  位置范围：%d-%d\n", min_bp, max_bp))
  
  # 使用 IEU API 提取该区域的 MEGASTROKE SNP
  tryCatch({
    # 方法 1: 优先按 RSID 提取（使用 JWT token）
    outcome_fmt <- extract_outcome_data(
      snps = exposure_fmt$SNP,
      outcomes = MEGASTROKE_ID,
      token = JWT_TOKEN
    )
    
    # 如果按 RSID 提取失败或数量太少，用区域提取
    if (is.null(outcome_fmt) || nrow(outcome_fmt) < 3) {
      cat(sprintf("  RSID 提取失败 (%d SNPs)，尝试区域提取...\n", 
                  ifelse(is.null(outcome_fmt), 0, nrow(outcome_fmt))))
      
      # 用染色体：起始 - 终止 提取
      variants <- paste0("chr", chr, ":", min_bp, "-", max_bp)
      region_snps <- associations(
        id = MEGASTROKE_ID,
        variants = variants,
        token = JWT_TOKEN
      )
      
      if (!is.null(region_snps) && nrow(region_snps) > 0) {
        outcome_fmt <- format_data(region_snps, type = "outcome")
        cat(sprintf("  区域提取成功：%d SNPs\n", nrow(outcome_fmt)))
      } else {
        cat(sprintf("  区域提取失败\n"))
        next
      }
    }
    
    if (is.null(outcome_fmt) || nrow(outcome_fmt) < 3) {
      cat(sprintf("  ✗ 结局数据不足 (%d SNPs)\n", nrow(outcome_fmt)))
      next
    }
    
    # Harmonise 数据
    dat <- harmonise_data(exposure_fmt, outcome_fmt)
    
    if (nrow(dat) < 3) {
      cat(sprintf("  ✗ Harmonise 后 SNP 不足 (%d)\n", nrow(dat)))
      next
    }
    
    harmonised_data[[gene_symbol]] <- dat
    
    # MR 分析（IVW 方法）
    res <- mr(dat, method_list = c("mr_ivw"))
    
    if (nrow(res) == 0) {
      cat(sprintf("  ✗ MR 分析失败\n"))
      next
    }
    
    # 提取 IVW 结果
    ivw_res <- res[res$method == "Inverse variance weighted", ]
    
    # 计算 F 统计量
    f_stat <- mean((exposure_fmt$beta.outcome / exposure_fmt$se.outcome)^2, na.rm = TRUE)
    
    # 组织分布
    tissue_dist <- table(exposure$TISSUE)
    n_brain <- ifelse("Brain_Cortex" %in% names(tissue_dist), tissue_dist["Brain_Cortex"], 0)
    n_blood <- ifelse("Whole_Blood" %in% names(tissue_dist), tissue_dist["Whole_Blood"], 0)
    
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
      n_brain = n_brain,
      n_blood = n_blood,
      stringsAsFactors = FALSE
    )
    
    mr_results[[gene_symbol]] <- result
    
    cat(sprintf("  ✓ MR: OR=%.3f (%.3f-%.3f), P=%.2e, F=%.2f\n", 
                result$or, result$ci_low, result$ci_high, result$pval, result$f_stat))
    
  }, error = function(e) {
    cat(sprintf("  ✗ 错误：%s\n", e$message))
  })
}

# 步骤 4: 保存结果
cat("\n\n步骤 4: 保存结果\n")
cat("----------------------------------------------------------------------\n")

if (length(mr_results) > 0) {
  results_df <- do.call(rbind, mr_results)
  results_df <- results_df[order(results_df$pval), ]
  
  # 保存为 CSV
  output_file <- file.path(OUTPUT_DIR, "mr_results_all_genes.csv")
  write.csv(results_df, output_file, row.names = FALSE)
  cat(sprintf("✓ 保存 MR 结果：%s\n", output_file))
  
  # 显示 top 结果
  cat(sprintf("\nTop 显著基因 (P < 0.05):\n"))
  sig_results <- results_df[results_df$pval < 0.05, ]
  if (nrow(sig_results) > 0) {
    print(sig_results[, c("gene", "or", "ci_low", "ci_high", "pval", "nsnp", "f_stat")])
  } else {
    cat("  无显著基因\n")
  }
  
  # 保存 FDR 校正后的 P 值
  results_df$pval_adj <- p.adjust(results_df$pval, method = "fdr")
  write.csv(results_df, file.path(OUTPUT_DIR, "mr_results_fdr_corrected.csv"), row.names = FALSE)
  cat(sprintf("✓ 保存 FDR 校正结果\n"))
}

# 保存 harmonised 数据
if (length(harmonised_data) > 0) {
  all_harmonised <- do.call(rbind, lapply(names(harmonised_data), function(gene) {
    dat <- harmonised_data[[gene]]
    dat$gene <- gene
    dat
  }))
  write.csv(all_harmonised, file.path(OUTPUT_DIR, "harmonised_data_all_genes.csv"), row.names = FALSE)
  cat(sprintf("✓ 保存 Harmonised 数据\n"))
}

cat("\n======================================================================\n")
cat("分析完成！\n")
cat("======================================================================\n")
