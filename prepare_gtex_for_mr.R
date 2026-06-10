#!/usr/bin/env Rscript
# ================================================================================
# 准备 GTEx 数据用于 MR 分析
# ================================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(readxl)
  library(arrow)
})

cat("======================================================================\n")
cat("准备 GTEx 数据用于 MR 分析\n")
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

# 完整的 ENSG ID 映射
gene_to_ensg <- list(
  "LYN" = "ENSG00000254087", "PRKCQ" = "ENSG00000065675", "NMT1" = "ENSG00000136448",
  "TDP1" = "ENSG00000042088", "MAN2B1" = "ENSG00000104774", "IL10RA" = "ENSG00000110324",
  "RHOC" = "ENSG00000155366", "SREBF1" = "ENSG00000072310", "CTSC" = "ENSG00000109861",
  "CAT" = "ENSG00000121691", "RBM39" = "ENSG00000131051", "CHFR" = "ENSG00000072609",
  "TCN2" = "ENSG00000185339", "SCN9A" = "ENSG00000169432", "CTSB" = "ENSG00000164733",
  "CASP8" = "ENSG00000064012", "MAPKAPK2" = "ENSG00000162889", "SEC13" = "ENSG00000157020",
  "PCTP" = "ENSG00000141179", "CPT2" = "ENSG00000157184", "BRD3" = "ENSG00000169925",
  "ITGA1" = "ENSG00000213949", "CITED2" = "ENSG00000164442", "HIBADH" = "ENSG00000106049",
  "SAT2" = "ENSG00000134294", "TSPO" = "ENSG00000100300", "FLT4" = "ENSG00000037280",
  "CPT1A" = "ENSG00000110090", "AKT1" = "ENSG00000142208", "CCR5" = "ENSG00000160791",
  "MGAT1" = "ENSG00000131446", "IGFBP2" = "ENSG00000115457", "PPARG" = "ENSG00000132170",
  "EPHX1" = "ENSG00000143819", "AIF1" = "ENSG00000204472", "PLA2G4A" = "ENSG00000116711",
  "ALDH9A1" = "ENSG00000131449", "BST1" = "ENSG00000197121", "CNDP2" = "ENSG00000133313",
  "TNF" = "ENSG00000232810", "PARP1" = "ENSG00000143799", "IKBKB" = "ENSG00000104365",
  "ADRB1" = "ENSG00000043591", "SPHK1" = "ENSG00000176170", "GCH1" = "ENSG00000131979",
  "ACADVL" = "ENSG00000072778", "CTSD" = "ENSG00000117984", "PDCD6IP" = "ENSG00000170248",
  "PTPRC" = "ENSG00000081237", "PABPC1" = "ENSG00000070756", "ACTA2" = "ENSG00000107796",
  "PTGR1" = "ENSG00000106853", "XRCC6" = "ENSG00000196419", "TBXAS1" = "ENSG00000059377",
  "NR1H3" = "ENSG00000025434", "CTSL" = "ENSG00000135047", "ZHX2" = "ENSG00000178764",
  "MKNK2" = "ENSG00000099875", "SERPINB10" = "ENSG00000242550", "HSPA5" = "ENSG00000044574",
  "CCND1" = "ENSG00000110092", "PTPN2" = "ENSG00000175354", "HSD17B4" = "ENSG00000133835",
  "PDCD6" = "ENSG00000249915", "PARP12" = "ENSG00000059378", "SERPINB1" = "ENSG00000197696",
  "STAT1" = "ENSG00000115415", "NFE2L2" = "ENSG00000116044", "HMOX1" = "ENSG00000100292",
  "CTSF" = "ENSG00000174080", "CCL2" = "ENSG00000108691", "FDX1" = "ENSG00000137714",
  "LIAS" = "ENSG00000121897", "LIPT1" = "ENSG00000144182", "PDHX" = "ENSG00000110435",
  "SLC31A1" = "ENSG00000136868", "ATP7B" = "ENSG00000123191", "ATOX1" = "ENSG00000177556",
  "GPX4" = "ENSG00000167468",
  "KCNA5" = "ENSG00000130037", "HIF1A" = "ENSG00000100644", "FABP4" = "ENSG00000170323",
  "STAT5A" = "ENSG00000126561", "FABP2" = "ENSG00000145384", "B2M" = "ENSG00000166710",
  "HBS1L" = "ENSG00000112339", "NUDCD2" = "ENSG00000170584", "JAK1" = "ENSG00000162434",
  "GPX1" = "ENSG00000233276", "FABP5" = "ENSG00000164687", "XDH" = "ENSG00000158125",
  "MB" = "ENSG00000198125", "POLR2D" = "ENSG00000144231", "HSD17B10" = "ENSG00000072506",
  "ZEB1" = "ENSG00000148516", "RELA" = "ENSG00000173039", "IRF1" = "ENSG00000125347",
  "GFAP" = "ENSG00000131095", "NR3C1" = "ENSG00000113580", "F3" = "ENSG00000117525",
  "C3" = "ENSG00000125730", "PTGS1" = "ENSG00000095303", "IMPDH2" = "ENSG00000178035",
  "PTPRF" = "ENSG00000142949", "HPGDS" = "ENSG00000163106", "PTPRJ" = "ENSG00000149177",
  "CASK" = "ENSG00000147044", "TOP2A" = "ENSG00000131747", "IL6" = "ENSG00000136244",
  "CP" = "ENSG00000047457", "S100A6" = "ENSG00000197956", "DDC" = "ENSG00000132437",
  "CUL4B" = "ENSG00000158290", "EGFR" = "ENSG00000146648", "COL1A1" = "ENSG00000108821",
  "STARD13" = "ENSG00000133121", "TGFB1" = "ENSG00000105329", "HTR2C" = "ENSG00000147246",
  "CTSS" = "ENSG00000163131", "CNR2" = "ENSG00000188822", "FNTA" = "ENSG00000168522",
  "RENBP" = "ENSG00000102032", "CCNA2" = "ENSG00000145386", "LEF1" = "ENSG00000138795",
  "SAT1" = "ENSG00000130066", "HTR2B" = "ENSG00000135914", "CDK4" = "ENSG00000135446",
  "CXCR3" = "ENSG00000186810", "TIMP1" = "ENSG00000102265", "OAZ1" = "ENSG00000104904",
  "STK4" = "ENSG00000101109", "ACADM" = "ENSG00000117054", "STAT3" = "ENSG00000168610",
  "NFKB1" = "ENSG00000109320", "CTSK" = "ENSG00000143387", "PTPN6" = "ENSG00000111679",
  "PA2G4" = "ENSG00000170515", "ACAD11" = "ENSG00000240303", "MAOB" = "ENSG00000069535",
  "ICAM1" = "ENSG00000090339", "DLAT" = "ENSG00000150768", "PDHB" = "ENSG00000168291",
  "ATP7A" = "ENSG00000165240", "MTOR" = "ENSG00000198793"
)

# 从 rsID 到 variant_id 的映射函数
extract_rs_id <- function(variant_id) {
  # variant_id 格式：chr1_64764_C_T_b38 或 ENSG00000310526.1
  # 从 eGenes 文件中提取 rs_id_dbSNP157_GRCh38p14
  if (grepl("^rs", variant_id)) {
    return(variant_id)
  } else {
    return(NA)
  }
}

# ================================================================================
# 处理全血 eQTL
# ================================================================================
cat("处理全血 eQTL 数据...\n")

# 读取 CSV 格式（更快）
blood_eqtl_file <- file.path(GTEx_DIR, "blood_eqtl.csv")
blood_genes_file <- file.path(GTEx_DIR, "Whole_Blood.v11.eGenes.txt")

cat("  加载 eQTL 数据...\n")
blood_eqtl <- fread(blood_eqtl_file, header = TRUE)
cat(sprintf("  原始数据：%d 行\n", nrow(blood_eqtl)))

# 读取 eGenes 获取 rsID 映射
cat("  加载 eGenes 获取 rsID 映射...\n")
blood_genes <- fread(blood_genes_file, sep = "\t", header = TRUE)

# 创建 variant_id 到 rsID 的映射
cat("  创建 rsID 映射...\n")
rsid_map <- blood_genes[, c("variant_id", "rs_id_dbSNP157_GRCh38p14")]
names(rsid_map)[2] <- "rs_id"
rsid_map <- rsid_map[!is.na(rsid_map$rs_id) & rsid_map$rs_id != "", ]
rsid_map <- unique(rsid_map)
cat(sprintf("  rsID 映射数：%d\n", nrow(rsid_map)))

# 过滤显著 eQTL
cat("  过滤显著 eQTL (P < 1e-5)...\n")
blood_sig <- blood_eqtl[blood_eqtl$pval_nominal < P_THRESHOLD, ]
cat(sprintf("  显著 eQTL 数：%d\n", nrow(blood_sig)))

# 添加 rsID
cat("  添加 rsID...\n")
blood_sig <- merge(blood_sig, rsid_map, by = "variant_id", all.x = TRUE)
cat(sprintf("  有 rsID 的 eQTL 数：%d\n", sum(!is.na(blood_sig$rs_id))))

# 转换为 TwoSampleMR 格式
cat("  转换为 TwoSampleMR 格式...\n")
blood_formatted <- data.frame(
  SNP = blood_sig$rs_id,
  beta.exposure = blood_sig$slope,
  se.exposure = blood_sig$slope_se,
  effect_allele.exposure = sapply(strsplit(as.character(blood_sig$variant_id), "_"), function(x) x[3]),
  other_allele.exposure = sapply(strsplit(as.character(blood_sig$variant_id), "_"), function(x) x[4]),
  eaf.exposure = blood_sig$af,
  pval.exposure = blood_sig$pval_nominal,
  gene = blood_sig$phenotype_id,
  exposure = "GTEx_Whole_Blood",
  stringsAsFactors = FALSE
)

# 过滤掉没有 rsID 的
blood_formatted <- blood_formatted[!is.na(blood_formatted$SNP) & blood_formatted$SNP != "", ]
cat(sprintf("  格式化后：%d 行\n", nrow(blood_formatted)))

# 按基因过滤（只保留目标基因）
cat("  按目标基因过滤...\n")
target_ensg <- unname(unlist(gene_to_ensg))
blood_formatted <- blood_formatted[blood_formatted$gene %in% target_ensg, ]
cat(sprintf("  目标基因相关 eQTL 数：%d\n", nrow(blood_formatted)))

# LD clumping
cat("  进行 LD clumping (r²<0.01, kb=1000)...\n")
blood_clumped <- clump_data(
  blood_formatted,
  clump_r2 = R2_THRESHOLD,
  clump_kb = KB_WINDOW,
  clump_p = P_THRESHOLD
)
cat(sprintf("  Clumping 后：%d 行\n\n", nrow(blood_clumped)))

# 保存
blood_output <- file.path(OUTPUT_DIR, "GTEx_Whole_Blood_p_1e-05_kb_1000_r2_0.01.xlsx")
write.table(blood_clumped, blood_output, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  保存到：%s\n\n", blood_output))

# ================================================================================
# 处理脑组织 eQTL
# ================================================================================
cat("处理脑组织 eQTL 数据...\n")

brain_eqtl_file <- file.path(GTEx_DIR, "brain_eqtl.csv")
brain_genes_file <- file.path(GTEx_DIR, "Brain_Cortex.v11.eGenes.txt")

cat("  加载 eQTL 数据...\n")
brain_eqtl <- fread(brain_eqtl_file, header = TRUE)
cat(sprintf("  原始数据：%d 行\n", nrow(brain_eqtl)))

cat("  加载 eGenes 获取 rsID 映射...\n")
brain_genes <- fread(brain_genes_file, sep = "\t", header = TRUE)

cat("  创建 rsID 映射...\n")
brain_rsid_map <- brain_genes[, c("variant_id", "rs_id_dbSNP157_GRCh38p14")]
names(brain_rsid_map)[2] <- "rs_id"
brain_rsid_map <- brain_rsid_map[!is.na(brain_rsid_map$rs_id) & brain_rsid_map$rs_id != "", ]
brain_rsid_map <- unique(brain_rsid_map)
cat(sprintf("  rsID 映射数：%d\n", nrow(brain_rsid_map)))

cat("  过滤显著 eQTL (P < 1e-5)...\n")
brain_sig <- brain_eqtl[brain_eqtl$pval_nominal < P_THRESHOLD, ]
cat(sprintf("  显著 eQTL 数：%d\n", nrow(brain_sig)))

cat("  添加 rsID...\n")
brain_sig <- merge(brain_sig, brain_rsid_map, by = "variant_id", all.x = TRUE)
cat(sprintf("  有 rsID 的 eQTL 数：%d\n", sum(!is.na(brain_sig$rs_id))))

cat("  转换为 TwoSampleMR 格式...\n")
brain_formatted <- data.frame(
  SNP = brain_sig$rs_id,
  beta.exposure = brain_sig$slope,
  se.exposure = brain_sig$slope_se,
  effect_allele.exposure = sapply(strsplit(as.character(brain_sig$variant_id), "_"), function(x) x[3]),
  other_allele.exposure = sapply(strsplit(as.character(brain_sig$variant_id), "_"), function(x) x[4]),
  eaf.exposure = brain_sig$af,
  pval.exposure = brain_sig$pval_nominal,
  gene = brain_sig$phenotype_id,
  exposure = "GTEx_Brain_Cortex",
  stringsAsFactors = FALSE
)

brain_formatted <- brain_formatted[!is.na(brain_formatted$SNP) & brain_formatted$SNP != "", ]
cat(sprintf("  格式化后：%d 行\n", nrow(brain_formatted)))

target_ensg <- unname(unlist(gene_to_ensg))
brain_formatted <- brain_formatted[brain_formatted$gene %in% target_ensg, ]
cat(sprintf("  目标基因相关 eQTL 数：%d\n", nrow(brain_formatted)))

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

# ================================================================================
# 整合 eQTLGen + GTEx 全血 + GTEx 脑组织
# ================================================================================
cat("整合所有 eQTL 数据源...\n")

eqtlgen_file <- file.path(OUTPUT_DIR, "eQTLgen_allgene_p_1e-05_kb_1000_r2_0.01.xlsx")
eqtlgen_data <- read_excel(eqtlgen_file)
cat(sprintf("  eQTLGen: %d SNPs\n", nrow(eqtlgen_data)))

# 读取 GTEx 数据
blood_clumped_df <- fread(blood_output, sep = "\t")
cat(sprintf("  GTEx 全血：%d SNPs\n", nrow(blood_clumped_df)))

brain_clumped_df <- fread(brain_output, sep = "\t")
cat(sprintf("  GTEx 脑组织：%d SNPs\n", nrow(brain_clumped_df)))

# 合并（去重）
merged_eqtl <- rbind(
  eqtlgen_data[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", 
                   "other_allele.exposure", "eaf.exposure", "pval.exposure", "gene")],
  blood_clumped_df[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", 
                       "other_allele.exposure", "eaf.exposure", "pval.exposure", "gene")],
  brain_clumped_df[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", 
                       "other_allele.exposure", "eaf.exposure", "pval.exposure", "gene")]
)

# 去重（优先保留 eQTLGen，其次 GTEx 全血，最后脑组织）
merged_eqtl <- merged_eqtl[order(merged_eqtl$pval.exposure), ]
merged_eqtl <- merged_eqtl[!duplicated(merged_eqtl$SNP), ]

cat(sprintf("  整合后：%d SNPs\n\n", nrow(merged_eqtl)))

# 保存整合数据
merged_output <- file.path(OUTPUT_DIR, "merged_eqtl_all_sources_p_1e-05_kb_1000_r2_0.01.xlsx")
write.table(merged_eqtl, merged_output, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  整合数据保存到：%s\n\n", merged_output))

# 统计每个数据源的贡献
cat("各数据源贡献:\n")
cat(sprintf("  eQTLGen: %d SNPs (%.1f%%)\n", 
            sum(grepl("eQTLGen", merged_eqtl$exposure) | is.na(merged_eqtl$exposure)),
            sum(grepl("eQTLGen", merged_eqtl$exposure) | is.na(merged_eqtl$exposure)) / nrow(merged_eqtl) * 100))
cat(sprintf("  GTEx 全血：%d SNPs (%.1f%%)\n", 
            sum(grepl("Whole_Blood", merged_eqtl$exposure), na.rm = TRUE),
            sum(grepl("Whole_Blood", merged_eqtl$exposure), na.rm = TRUE) / nrow(merged_eqtl) * 100))
cat(sprintf("  GTEx 脑组织：%d SNPs (%.1f%%)\n", 
            sum(grepl("Brain_Cortex", merged_eqtl$exposure), na.rm = TRUE),
            sum(grepl("Brain_Cortex", merged_eqtl$exposure), na.rm = TRUE) / nrow(merged_eqtl) * 100))

cat("\n完成！\n")
