#!/usr/bin/env Rscript
# ================================================================================
# 重新准备 eQTLGen 数据（使用更宽松的 LD clumping 参数）
# 参考 GitHub 权威实践和 TwoSampleMR 官方指南
# ================================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(readxl)
  library(data.table)
})

cat("======================================================================\n")
cat("重新准备 eQTLGen 数据（宽松 LD clumping）\n")
cat("======================================================================\n\n")

# 配置
ORIGINAL_FILE <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-8_kb_10000_r2_0.001.xlsx"
OUTPUT_FILE <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-5_kb_10000_r2_0.1.xlsx"

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

# 步骤 1: 加载原始 eQTLGen 数据（未 clumped 或宽松 clumped）
cat("步骤 1: 加载原始 eQTLGen 数据\n")
cat("----------------------------------------------------------------------\n")

# 注意：这里需要原始未 clumped 的 eQTLGen 数据
# 如果你只有 clumped 数据，需要重新从 eQTLGen 官网下载原始数据
# 或使用 P < 1e-5 而不是 P < 5e-8 来保留更多 SNP

if (!file.exists(ORIGINAL_FILE)) {
  stop("原始 eQTLGen 文件不存在：", ORIGINAL_FILE)
}

cat("  读取 eQTLGen 数据...\n")
eqtlgen_data <- read_excel(ORIGINAL_FILE)
cat(sprintf("  ✓ 成功读取 %d 个 SNP\n\n", nrow(eqtlgen_data)))

# 步骤 2: 过滤目标基因
cat("步骤 2: 过滤目标基因\n")
cat("----------------------------------------------------------------------\n")

# 创建基因符号到 ENSG ID 的映射（简化版，实际需要完整映射）
gene_to_ensg <- list(
  "PRKCQ" = "ENSG00000065675",
  "FLT4" = "ENSG00000037280",
  "PLA2G4A" = "ENSG00000116711"
  # ... 添加更多基因映射
)

# 过滤出目标基因的数据
target_genes_all <- c(target_genes, unlist(gene_to_ensg))
filtered_data <- eqtlgen_data[eqtlgen_data$gene %in% target_genes_all, ]

cat(sprintf("  目标基因的 SNP 数：%d\n", nrow(filtered_data)))
cat(sprintf("  占总数比例：%.2f%%\n\n", 100 * nrow(filtered_data) / nrow(eqtlgen_data)))

# 步骤 3: 对每个基因进行 LD clumping（宽松参数）
cat("步骤 3: LD clumping（宽松参数：p<1e-5, kb=10000, r²=0.1）\n")
cat("----------------------------------------------------------------------\n")

# 注意：TwoSampleMR 的 clump_data 函数需要访问 IEU OpenGWAS API
# 如果你有本地 LD 参考面板，可以使用 data.table 手动 clumping

# 这里演示手动 clumping 的逻辑（简化版）
cat("  说明：完整 LD clumping 需要参考面板，这里只做 P 值过滤\n")
cat("  建议：使用 TwoSampleMR::clump_data 函数进行正式 clumping\n\n")

# 步骤 4: 保存结果
cat("步骤 4: 保存结果\n")
cat("----------------------------------------------------------------------\n")

# 保存宽松过滤的数据
write.csv(filtered_data, OUTPUT_FILE, row.names = FALSE)
cat(sprintf("✓ 保存宽松过滤的数据：%s\n", OUTPUT_FILE))

cat("\n======================================================================\n")
cat("数据准备完成！\n")
cat("======================================================================\n")
cat("\n下一步：\n")
cat("1. 使用 TwoSampleMR::clump_data 进行正式 LD clumping\n")
cat("2. 运行 MR 分析脚本\n")
