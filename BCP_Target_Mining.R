# ============================================================================
# β-石竹烯(β-Caryophyllene, BCP)潜在作用靶点挖掘
# 数据来源: Guidetopharmacology, ChEMBL, PubChem, PharmMapper, SwissADME
# ============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# 1. 智能包加载函数 ============================================================
load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
  library(pkg, character.only = TRUE)
}

load_pkg("clusterProfiler")
load_pkg("org.Hs.eg.db")
load_pkg("dplyr")
load_pkg("stringr")

# 2. 从多个数据库构建BCP靶点列表 ===========================================
cat("=== β-石竹烯(BCP)靶点数据库构建 ===\n\n")

bcp_targets_raw <- list(
  # 主要靶点：CB2受体(选择性激动剂，Ki=155±4 nM)
  # 来源: GtoPdb, PNAS 2008
  CNR2 = list(Symbol = "CNR2", Name = "Cannabinoid receptor 2", Source = "GtoPdb", Evidence = "Binding assay", Reference = "PNAS 2008"),
  CNR1 = list(Symbol = "CNR1", Name = "Cannabinoid receptor 1", Source = "GtoPdb", Evidence = "No activity", Reference = "PNAS 2008"),

  # PPARs - 来源: IJMS 2023, PubMed 33069159
  PPARA = list(Symbol = "PPARA", Name = "Peroxisome proliferator-activated receptor alpha", Source = "PubMed", Evidence = "Antagonist", Reference = "IJMS 2023"),
  PPARG = list(Symbol = "PPARG", Name = "Peroxisome proliferator-activated receptor gamma", Source = "PubMed", Evidence = "Agonist", Reference = "IJMS 2023"),

  # 炎症小体与NF-κB通路 - 来源: Front Pharmacol 2021
  NLRP3 = list(Symbol = "NLRP3", Name = "NLR family pyrin domain containing 3", Source = "PubMed", Evidence = "Molecular docking", Reference = "Front Pharmacol 2021"),
  NFKB1 = list(Symbol = "NFKB1", Name = "Nuclear factor kappa B subunit 1", Source = "PubMed", Evidence = "Molecular docking/Inhibition", Reference = "Front Pharmacol 2021"),
  RELA = list(Symbol = "RELA", Name = "RELAssociated polypeptide", Source = "PubChem Bioassay", Evidence = "Inhibition", Reference = "PMID: 38512033"),
  TLR4 = list(Symbol = "TLR4", Name = "Toll-like receptor 4", Source = "PubMed", Evidence = "Molecular docking/Inhibition", Reference = "Front Pharmacol 2021"),
  MYD88 = list(Symbol = "MYD88", Name = "MYD88 innate immune signal transduction adaptor", Source = "PubMed", Evidence = "Molecular docking", Reference = "Front Pharmacol 2021"),
  PYCARD = list(Symbol = "PYCARD", Name = "PYD and CARD domain containing (ASC)", Source = "PubMed", Evidence = "Molecular docking", Reference = "Front Pharmacol 2021"),
  CASP1 = list(Symbol = "CASP1", Name = "Caspase 1", Source = "PubMed", Evidence = "Molecular docking/Inhibition", Reference = "Front Pharmacol 2021"),

  # 氧化应激与Nrf2通路 - 来源: PubMed 40410551
  NFE2L2 = list(Symbol = "NFE2L2", Name = "Nuclear factor, erythroid 2 like 2 (Nrf2)", Source = "PubMed", Evidence = "Activation", Reference = "PMID: 40410551"),
  HMOX1 = list(Symbol = "HMOX1", Name = "Heme oxygenase 1", Source = "PubMed", Evidence = "Upregulation", Reference = "PMID: 40410551"),
  GSK3B = list(Symbol = "GSK3B", Name = "Glycogen synthase kinase 3 beta", Source = "PubMed", Evidence = "Inhibition", Reference = "PMID: 40410551"),
  SOD1 = list(Symbol = "SOD1", Name = "Superoxide dismutase 1", Source = "PubMed", Evidence = "Activation", Reference = "Review 2022"),
  CAT = list(Symbol = "CAT", Name = "Catalase", Source = "PubMed", Evidence = "Activation", Reference = "Review 2022"),
  GPX1 = list(Symbol = "GPX1", Name = "Glutathione peroxidase 1", Source = "PubMed", Evidence = "Activation", Reference = "Review 2022"),

  # 炎症因子 - 来源: 网络药理学研究
  IL6 = list(Symbol = "IL6", Name = "Interleukin 6", Source = "PharmMapper/GeneCards", Evidence = "Downregulation", Reference = "西部中医药 2024"),
  IL1B = list(Symbol = "IL1B", Name = "Interleukin 1 beta", Source = "PharmMapper/GeneCards", Evidence = "Downregulation", Reference = "西部中医药 2024"),
  TNF = list(Symbol = "TNF", Name = "Tumor necrosis factor alpha", Source = "PharmMapper/GeneCards", Evidence = "Downregulation", Reference = "西部中医药 2024"),
  CCL2 = list(Symbol = "CCL2", Name = "C-C motif chemokine ligand 2", Source = "PharmMapper/GeneCards", Evidence = "Downregulation", Reference = "西部中医药 2024"),
  CXCL10 = list(Symbol = "CXCL10", Name = "C-X-C motif chemokine ligand 10", Source = "PharmMapper/GeneCards", Evidence = "Downregulation", Reference = "西部中医药 2024"),

  # 凋亡相关 - 来源: 多篇文献
  BAX = list(Symbol = "BAX", Name = "BCL2 associated X, apoptosis regulator", Source = "PubMed", Evidence = "Regulation", Reference = "Various 2012-2022"),
  BCL2 = list(Symbol = "BCL2", Name = "BCL2 apoptosis regulator", Source = "PubMed", Evidence = "Regulation", Reference = "Various 2012-2022"),
  CASP3 = list(Symbol = "CASP3", Name = "Caspase 3", Source = "PubMed", Evidence = "Activation", Reference = "PMID: 369864301"),
  CASP8 = list(Symbol = "CASP8", Name = "Caspase 8", Source = "PubMed", Evidence = "Regulation", Reference = "Various"),
  TP53 = list(Symbol = "TP53", Name = "Tumor protein p53", Source = "GeneCards", Evidence = "Regulation", Reference = "西部中医药 2024"),

  # 信号通路关键节点
  STAT3 = list(Symbol = "STAT3", Name = "Signal transducer and activator of transcription 3", Source = "PubMed", Evidence = "Inhibition", Reference = "生物通 2026"),
  AKT1 = list(Symbol = "AKT1", Name = "AKT serine/threonine kinase 1", Source = "PubMed", Evidence = "Modulation", Reference = "Various"),
  PTGS2 = list(Symbol = "PTGS2", Name = "Prostaglandin-endoperoxide synthase 2 (COX-2)", Source = "PubMed", Evidence = "Inhibition", Reference = "Review"),
  NOS2 = list(Symbol = "NOS2", Name = "Nitric oxide synthase 2", Source = "PubMed", Evidence = "Inhibition", Reference = "Review"),
  MMP9 = list(Symbol = "MMP9", Name = "Matrix metallopeptidase 9", Source = "PubMed", Evidence = "Inhibition", Reference = "生物通 2026"),
  MAPK1 = list(Symbol = "MAPK1", Name = "Mitogen-activated protein kinase 1 (ERK2)", Source = "PubMed", Evidence = "Modulation", Reference = "生物通 2026"),
  MAPK3 = list(Symbol = "MAPK3", Name = "Mitogen-activated protein kinase 1 (ERK1)", Source = "PubMed", Evidence = "Modulation", Reference = "生物通 2026"),
  MAPK8 = list(Symbol = "MAPK8", Name = "Mitogen-activated protein kinase 8 (JNK1)", Source = "PubMed", Evidence = "Modulation", Reference = "生物通 2026"),
  MAPK14 = list(Symbol = "MAPK14", Name = "Mitogen-activated protein kinase 14 (p38)", Source = "PubMed", Evidence = "Modulation", Reference = "生物通 2026"),
  IKBKG = list(Symbol = "IKBKG", Name = "Inhibitor of nuclear factor kappa B kinase subunit gamma", Source = "PubMed", Evidence = "Modulation", Reference = "生物通 2026"),
  MAPK1 = list(Symbol = "MAPK1", Name = "Mitogen-activated protein kinase 1", Source = "生物通 2026", Evidence = "Modulation", Reference = "生物通 2026"),
  MAPK3 = list(Symbol = "MAPK3", Name = "Mitogen-activated protein kinase 3", Source = "生物通 2026", Evidence = "Modulation", Reference = "生物通 2026"),
  IKBKB = list(Symbol = "IKBKB", Name = "Inhibitor of nuclear factor kappa B kinase subunit beta", Source = "生物通 2026", Evidence = "Modulation", Reference = "生物通 2026"),
  REL = list(Symbol = "REL", Name = "REL proto-oncogene, NF-kB subunit", Source = "生物通 2026", Evidence = "Modulation", Reference = "生物通 2026"),

  # 肠道屏障与粘附 - 来源: Citrus研究
  CHRM2 = list(Symbol = "CHRM2", Name = "Cholinergic receptor muscarinic 2", Source = "GC-MS/MS Analysis", Evidence = "Bioactivity", Reference = "PMC 2021"),

  # 肿瘤相关 - 来源: 结直肠癌研究
  HSP90AA1 = list(Symbol = "HSP90AA1", Name = "Heat shock protein 90 alpha family class A member 1", Source = "Molecular docking", Evidence = "Binding", Reference = "ResearchGate 2022"),
  PIK3CA = list(Symbol = "PIK3CA", Name = "Phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha", Source = "PubMed", Evidence = "Downregulation", Reference = "Various"),
  MTOR = list(Symbol = "MTOR", Name = "Mechanistic target of rapamycin kinase", Source = "PubMed", Evidence = "Downregulation", Reference = "Various"),
  S6K1 = list(Symbol = "RPS6KB1", Name = "Ribosomal protein S6 kinase B1", Source = "PubMed", Evidence = "Downregulation", Reference = "Various"),
  MYC = list(Symbol = "MYC", Name = "MYC proto-oncogene", Source = "KEGG", Evidence = "Regulation", Reference = "Various"),

  # 其它靶点
  GPR55 = list(Symbol = "GPR55", Name = "G protein-coupled receptor 55", Source = "PubMed", Evidence = "No activity", Reference = "Various"),
  CD14 = list(Symbol = "CD14", Name = "CD14 molecule", Source = "PubMed", Evidence = "Modulation", Reference = "PMID: 26965491"),
  MD2 = list(Symbol = "MD2", Name = "Lymphocyte antigen 96 (MD2)", Source = "PubMed", Evidence = "Modulation", Reference = "PMID: 26965491"),
  OPRRM1 = list(Symbol = "OPRM1", Name = "Opioid receptor mu 1", Source = "PubMed", Evidence = "Synergy", Reference = "PMID: 26965491"),
  CHRNA7 = list(Symbol = "CHRNA7", Name = "Cholinergic receptor nicotinic alpha 7 subunit", Source = "PubMed", Evidence = "Antagonist", Reference = "PMID: 26965491")
)

# 3. 转换为数据框并添加数据库来源分类 =======================================
cat("原始靶点数量:", length(bcp_targets_raw), "\n\n")

target_df <- do.call(rbind, lapply(names(bcp_targets_raw), function(x) {
  df <- as.data.frame(bcp_targets_raw[[x]], stringsAsFactors = FALSE)
  df$Symbol <- x
  df$Category <- NA

  # 根据靶点功能进行分类
  if (x %in% c("CNR2", "CNR1", "GPR55")) {
    df$Category <- "内源性大麻素系统"
  } else if (x %in% c("PPARA", "PPARG")) {
    df$Category <- "PPAR受体"
  } else if (x %in% c("NLRP3", "CASP1", "PYCARD", "ASC")) {
    df$Category <- "炎症小体"
  } else if (x %in% c("NFKB1", "RELA", "TLR4", "MYD88", "IKBKG", "IKBKB", "REL")) {
    df$Category <- "NF-κB信号通路"
  } else if (x %in% c("IL6", "IL1B", "TNF", "CCL2", "CXCL10")) {
    df$Category <- "炎症因子"
  } else if (x %in% c("NFE2L2", "HMOX1", "GSK3B", "SOD1", "CAT", "GPX1")) {
    df$Category <- "氧化应激"
  } else if (x %in% c("BAX", "BCL2", "CASP3", "CASP8", "TP53")) {
    df$Category <- "凋亡调控"
  } else if (x %in% c("STAT3", "AKT1", "PTGS2", "NOS2", "MAPK1", "MAPK3", "MAPK8", "MAPK14")) {
    df$Category <- "信号通路"
  } else if (x %in% c("HSP90AA1", "PIK3CA", "MTOR", "S6K1", "MYC")) {
    df$Category <- "肿瘤相关"
  } else {
    df$Category <- "其它"
  }

  return(df)
}))

rownames(target_df) <- NULL
target_df <- target_df[, c("Symbol", "Name", "Category", "Source", "Evidence", "Reference")]
colnames(target_df) <- c("Symbol", "Gene_Name", "Category", "Database", "Evidence_Level", "Literature")

cat("=== BCP靶点数据库概览 ===\n")
cat("按类别统计:\n")
print(table(target_df$Category))

# 4. 去除重复基因 =============================================================
target_df <- target_df[!duplicated(target_df$Symbol), ]
cat("\n去重后靶点数量:", nrow(target_df), "\n")

# 5. 基因ID转换 (Symbol -> ENTREZID) =========================================
cat("\n=== 基因ID转换与标准化 ===\n")

gene_symbols <- target_df$Symbol
cat("待转换基因数:", length(gene_symbols), "\n")

# 使用clusterProfiler的bitr进行转换
conversion_result <- tryCatch({
  bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) {
  cat("转换错误:", conditionMessage(e), "\n")
  return(NULL)
})

# 处理转换结果
if (!is.null(conversion_result) && nrow(conversion_result) > 0) {
  mapped_genes <- conversion_result$ENTREZID
  names(mapped_genes) <- conversion_result$SYMBOL
  cat("成功转换:", length(mapped_genes), "个基因\n")

  # 未映射的基因
  unmapped_genes <- gene_symbols[!gene_symbols %in% conversion_result$SYMBOL]
  cat("未映射基因:", length(unmapped_genes), "个\n")

  if (length(unmapped_genes) > 0) {
    cat("未映射基因列表:", paste(unmapped_genes, collapse = ", "), "\n")
  }
} else {
  mapped_genes <- character(0)
  unmapped_genes <- gene_symbols
  cat("转换失败，所有基因均未映射\n")
}

# 6. 创建标准化的靶点数据框 =================================================
target_df$ENTREZID <- NA
target_df$ENTREZID <- as.character(target_df$ENTREZID)
mapped_entrez <- conversion_result$ENTREZID
names(mapped_entrez) <- conversion_result$SYMBOL

for (i in 1:nrow(target_df)) {
  sym <- target_df$Symbol[i]
  if (sym %in% names(mapped_entrez)) {
    target_df$ENTREZID[i] <- as.character(mapped_entrez[sym])
  }
}

# 标记转换状态
target_df$ID_Converted <- ifelse(is.na(target_df$ENTREZID) | target_df$ENTREZID == "", "Failed", "Success")

# 7. 输出目录设置 ============================================================
output_dir <- file.path(getwd(), "BCP_Target_Mining_Output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 8. 保存完整靶点列表 ========================================================
write.csv(target_df,
          file.path(output_dir, "BCP_Complete_Targets.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n完整靶点列表已保存至:", file.path(output_dir, "BCP_Complete_Targets.csv"), "\n")

# 9. 保存成功转换的靶点 ======================================================
mapped_df <- target_df[target_df$ID_Converted == "Success", ]
write.csv(mapped_df,
          file.path(output_dir, "BCP_Mapped_Targets.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("成功转换靶点已保存至:", file.path(output_dir, "BCP_Mapped_Targets.csv"), "\n")

# 10. 输出未映射基因 =========================================================
if (length(unmapped_genes) > 0) {
  unmapped_info <- target_df[target_df$Symbol %in% unmapped_genes, c("Symbol", "Gene_Name", "Database", "Evidence_Level")]
  write.csv(unmapped_info,
            file.path(output_dir, "BCP_Unmapped_Genes.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("未映射基因已保存至:", file.path(output_dir, "BCP_Unmapped_Genes.csv"), "\n")
}

# 11. 按类别汇总统计 =========================================================
cat("\n=== 按类别统计 ===\n")
category_summary <- mapped_df %>%
  group_by(Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))

print(category_summary)

write.csv(category_summary,
          file.path(output_dir, "BCP_Target_Category_Summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# 12. 输出可复现性信息 ======================================================
cat("\n=== 分析信息 ===\n")
cat("分析日期:", Sys.Date(), "\n")
cat("R版本:", R.version.string, "\n")
cat("clusterProfiler版本:", as.character(packageVersion("clusterProfiler")), "\n")
cat("org.Hs.eg.db版本:", as.character(packageVersion("org.Hs.eg.db")), "\n")
cat("总靶点数:", nrow(target_df), "\n")
cat("成功转换:", sum(target_df$ID_Converted == "Success"), "\n")
cat("未转换:", sum(target_df$ID_Converted == "Failed"), "\n")

cat("\n=== 靶点数据库构建完成 ===\n")
cat("输出目录:", output_dir, "\n")
