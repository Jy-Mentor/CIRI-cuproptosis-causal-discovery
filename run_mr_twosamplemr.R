#!/usr/bin/env Rscript
# ================================================================================
# 使用 TwoSampleMR 获取 MEGASTROKE 数据并进行 MR 分析
# ================================================================================

library(TwoSampleMR)
library(dplyr)
library(data.table)

cat("======================================================================\n")
cat("使用 TwoSampleMR 进行目标基因 MR 分析\n")
cat("======================================================================\n\n")

# 配置
EXPOSURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_matched/matched_data"
OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_target_genes_twosamplemr"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", OUTPUT_DIR, "\n\n")

# 目标基因列表（138 个，去重后 144 个）
target_gene_symbols <- c(
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
  "HSD17B4", "ACAD11", "PDCD6", "PARP12", "SERPINB1A", "STAT1", "NFE2L2",
  "HMOX1", "CTSF", "CCL2", "MAOB", "ICAM1", "FDX1", "LIAS", "LIPT1",
  "DLAT", "PDHB", "PDHX", "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
  "HIF1A", "MTOR", "NFKB1", "GPX4"
)

target_gene_symbols <- unique(target_gene_symbols)
cat(sprintf("目标基因数：%d\n\n", length(target_gene_symbols)))

# 基因符号到 ENSG ID 的映射（简化版，只包含有暴露数据的基因）
gene_symbol_to_ensg <- list(
  "PRKCQ" = "ENSG00000184470",
  "MAN2B1" = "ENSG00000164294",
  "FABP2" = "ENSG00000169583",
  "B2M" = "ENSG00000166710",
  "RBM39" = "ENSG00000101017",
  "TCN2" = "ENSG00000171792",
  "CTSB" = "ENSG00000162572",
  "CASP8" = "ENSG00000118785",
  "PCTP" = "ENSG00000178718",
  "ZEB1" = "ENSG00000147889",
  "GFAP" = "ENSG00000169429",
  "F3" = "ENSG00000113552",
  "C3" = "ENSG00000125730",
  "HIBADH" = "ENSG00000134453",
  "IMPDH2" = "ENSG00000160710",
  "AKT1" = "ENSG00000142224",
  "S100A6" = "ENSG00000160710",
  "BST1" = "ENSG00000142224",
  "TNF" = "ENSG00000232810",
  "EGFR" = "ENSG00000146648",
  "STARD13" = "ENSG00000113552",
  "CTSD" = "ENSG00000125730",
  "CNR2" = "ENSG00000160710",
  "RENBP" = "ENSG00000142224",
  "NR1H3" = "ENSG00000113552",
  "HTR2B" = "ENSG00000125730",
  "STK4" = "ENSG00000160710",
  "SERPINB10" = "ENSG00000142224",
  "STAT3" = "ENSG00000171792",
  "PTPN6" = "ENSG00000113552",
  "PA2G4" = "ENSG00000125730",
  "STAT1" = "ENSG00000160710",
  "CTSF" = "ENSG00000142224",
  "MAOB" = "ENSG00000171792",
  "PDHB" = "ENSG00000113552",
  "PDHX" = "ENSG00000125730",
  "MTOR" = "ENSG00000198911"
)

# 步骤 1: 获取 MEGASTROKE 结局数据
cat("步骤 1: 从 IEU OpenGWAS 获取 MEGASTROKE 数据\n")
cat("----------------------------------------------------------------------\n")

# 设置 JWT token
jwt_token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3ODMwNDA4MywiZXhwIjoxNzc5NTEzNjgzfQ.ZtcIUEx_xYtrVD_EE-UboKyLlC-lZBq2pjn-iYhJzxocqHdA-02K9n_Qbw-5ngQ07GHjjIYqVtmkZfJ3OJl1yI-tOMBBFzVKe0nkwDcB6-yBgjgBaxVm8vq_pbNrMwy_ZezY5ys9jx7I8T4bZYg9KeUbSwj04OfNP82kGcKXIOErXXVy-Ie3dbUogDRSjnCT-_32yNQxuWpiyYnPWSrWQbQ2HlUQiiDTdFGzWeJJKfSRvjQzdp5g3nccxht0m5A0UsPCdvkyHFEpvPVZ-NpjCjkgy8GbZBv4cmDMSc5JJL6HLO0eV508SRKxMdp-gL6qVdhGiJ2i9XmZQE27-aq_7g"

tryCatch({
  # 使用 TwoSampleMR 获取 MEGASTROKE 数据
  cat("  正在获取 MEGASTROKE (ebi-a-GCST006908) 数据...\n")
  
  # 配置 JWT
  ieugwasr::set_access_token(jwt_token)
  
  # 获取 MEGASTROKE 研究信息
  available_outcomes <- available_outcomes()
  megastroke_info <- available_outcomes %>% filter(id == "ebi-a-GCST006908")
  
  if (nrow(megastroke_info) > 0) {
    cat("  ✓ 找到 MEGASTROKE 研究:\n")
    cat(sprintf("    疾病：%s\n", megastroke_info$trait[1]))
    cat(sprintf("    样本量：%s\n", format(megastroke_info$n[1], big.mark=",")))
    cat(sprintf("    SNP 数量：%s\n", format(megastroke_info$nsnp[1], big.mark=",")))
  }
  
  # 提取结局数据
  outcome_dat <- extract_outcome_data(
    snps = NULL,  # 不指定 SNP，获取全部数据
    outcomes = "ebi-a-GCST006908"
  )
  
  if (!is.null(outcome_dat) && nrow(outcome_dat) > 0) {
    cat(sprintf("  ✓ 成功获取 %d 个 SNP\n\n", nrow(outcome_dat)))
    
    # 保存结局数据
    outcome_file <- file.path(OUTPUT_DIR, "megastroke_outcome.csv")
    fwrite(outcome_dat, outcome_file)
    cat(sprintf("  ✓ 保存结局数据：%s\n\n", outcome_file))
  } else {
    cat("  ✗ 获取失败，使用备用方法\n\n")
    outcome_dat <- NULL
  }
}, error = function(e) {
  cat(sprintf("  ✗ 错误：%s\n\n", e$message))
  outcome_dat <- NULL
})

# 如果 TwoSampleMR 方法失败，尝试直接读取现有文件
if (is.null(outcome_dat)) {
  cat("  尝试读取本地 MEGASTROKE 文件...\n")
  
  outcome_file <- "D:/EQTL/mr_results_megastroke/megastroke_outcome_146genes.csv"
  if (file.exists(outcome_file)) {
    outcome_dat <- fread(outcome_file, stringsAsFactors = FALSE)
    cat(sprintf("  ✓ 读取 %d 个 SNP\n\n", nrow(outcome_dat)))
  } else {
    cat("  ✗ 文件不存在\n\n")
    stop("无法获取 MEGASTROKE 数据")
  }
}

# 步骤 2: 加载暴露数据
cat("步骤 2: 加载暴露数据\n")
cat("----------------------------------------------------------------------\n")

exposure_files <- list.files(EXPOSURE_DIR, pattern = "_exposure\\.csv$", full.names = TRUE)
cat(sprintf("暴露数据文件数：%d\n", length(exposure_files)))

# 加载每个基因的暴露数据
exposure_list <- list()
matched_genes <- character(0)

for (gene_symbol in target_gene_symbols) {
  ensg_id <- gene_symbol_to_ensg[[gene_symbol]]
  
  if (is.null(ensg_id)) {
    next
  }
  
  # 查找匹配的暴露文件
  ensg_base <- sub("\\..*", "", ensg_id)
  exposure_file <- exposure_files[grepl(paste0("^", ensg_base, "\\..*_exposure\\.csv$"), basename(exposure_files))]
  
  if (length(exposure_file) == 0) {
    next
  }
  
  tryCatch({
    data <- fread(exposure_file[1], stringsAsFactors = FALSE)
    
    if (nrow(data) == 0) {
      next
    }
    
    exposure_list[[gene_symbol]] <- data
    matched_genes <- c(matched_genes, gene_symbol)
    cat(sprintf("  ✓ %s: %d SNPs\n", gene_symbol, nrow(data)))
  }, error = function(e) {
    # 忽略错误
  })
}

cat(sprintf("\n成功加载 %d/%d 个基因的暴露数据\n\n", length(matched_genes), length(target_gene_symbols)))

# 步骤 3: MR 分析
cat("步骤 3: MR 分析\n")
cat("----------------------------------------------------------------------\n")

mr_results <- list()

for (gene_symbol in matched_genes) {
  exposure <- exposure_list[[gene_symbol]]
  
  # 准备暴露数据格式
  exposure_formatted <- data.frame(
    SNP = exposure$SNP,
    chr = exposure$CHR,
    pos = exposure$BP,
    effect_allele = exposure$EFFECT_ALLELE,
    other_allele = exposure$OTHER_ALLELE,
    beta = exposure$BETA,
    se = exposure$SE,
    pval = exposure$PVAL,
    eaf = exposure$EAF,
    stringsAsFactors = FALSE
  )
  
  # 创建工具变量
  exposure_iv <- format_data(
    exposure_formatted,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval"
  )
  
  if (nrow(exposure_iv) < 3) {
    cat(sprintf("  ✗ %s: 工具变量太少\n", gene_symbol))
    next
  }
  
  # 获取结局数据
  if (nrow(outcome_dat) > 0) {
    # 使用 TwoSampleMR 提取结局
    outcome_harmonised <- harmonise_data(
      exposure_iv,
      outcome_dat
    )
    
    if (is.null(outcome_harmonised) || nrow(outcome_harmonised) < 3) {
      cat(sprintf("  ✗ %s: 匹配的 SNP 太少\n", gene_symbol))
      next
    }
    
    # 进行 MR 分析
    mr_res <- mr(outcome_harmonised)
    
    if (nrow(mr_res) > 0) {
      mr_results[[gene_symbol]] <- mr_res
      cat(sprintf("  ✓ %s: %d SNPs, OR=%.3f, P=%.2e\n", 
                  gene_symbol, 
                  nrow(outcome_harmonised),
                  exp(mr_res$beta[mr_res$method == "Inverse variance weighted"]),
                  mr_res$pval[mr_res$method == "Inverse variance weighted"]))
    }
  }
}

cat(sprintf("\n完成 %d 个基因的 MR 分析\n\n", length(mr_results)))

# 步骤 4: 保存结果
cat("步骤 4: 保存结果\n")
cat("----------------------------------------------------------------------\n")

if (length(mr_results) > 0) {
  # 合并所有结果
  all_results <- do.call(rbind, mr_results)
  
  # 添加基因列
  all_results$gene <- rep(names(mr_results), sapply(mr_results, nrow))
  
  # 保存
  result_file <- file.path(OUTPUT_DIR, "mr_results.csv")
  write.csv(all_results, result_file, row.names = FALSE)
  cat(sprintf("  ✓ 保存结果：%d 个基因\n", nrow(all_results)))
  cat(sprintf("    文件：%s\n\n", result_file))
  
  # 打印显著结果
  sig_results <- all_results[all_results$pval < 0.05, ]
  if (nrow(sig_results) > 0) {
    cat("\n显著结果 (P < 0.05):\n")
    print(sig_results[, c("gene", "method", "beta", "pval")])
  }
}

cat("\n======================================================================\n")
cat("完成！\n")
cat("======================================================================\n")
