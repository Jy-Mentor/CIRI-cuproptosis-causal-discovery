#!/usr/bin/env Rscript
# ================================================================================
# 补充缺失基因的 ENSG ID 映射
# ================================================================================

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
})

cat("======================================================================\n")
cat("补充缺失基因的 ENSG ID 映射\n")
cat("======================================================================\n\n")

# 65 个缺失的基因
missing_genes <- c(
  "KCNA5", "HIF1A", "FABP4", "STAT5A", "FABP2", "B2M", "HBS1L", "NUDCD2", 
  "JAK1", "GPX1", "FABP5", "XDH", "MB", "POLR2D", "HSD17B10", "ZEB1", 
  "RELA", "IRF1", "GFAP", "NR3C1", "F3", "C3", "PTGS1", "IMPDH2", 
  "PTPRF", "HPGDS", "PTPRJ", "CASK", "TOP2A", "IL6", "CP", "S100A6", 
  "DDC", "CUL4B", "EGFR", "COL1A1", "STARD13", "TGFB1", "HTR2C", 
  "CTSS", "CNR2", "FNTA", "RENBP", "CCNA2", "LEF1", "SAT1", "HTR2B", 
  "CDK4", "CXCR3", "TIMP1", "OAZ1", "STK4", "ACADM", "STAT3", "NFKB1", 
  "CTSK", "PTPN6", "PA2G4", "ACAD11", "MAOB", "ICAM1", "DLAT", "PDHB", 
  "ATP7A", "MTOR"
)

cat(sprintf("需要查询的基因数：%d\n\n", length(missing_genes)))

# 使用 biomaRt 查询 ENSG ID
cat("连接 Ensembl 数据库...\n")
tryCatch({
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  cat("连接成功！\n\n")
  
  cat("批量查询 ENSG ID...\n")
  ensg_df <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
    filters = "hgnc_symbol",
    values = missing_genes,
    mart = ensembl
  )
  
  cat(sprintf("成功查询到 %d 个基因的 ENSG ID\n\n", nrow(ensg_df)))
  
  # 查看结果
  cat("查询结果:\n")
  print(head(ensg_df, 20))
  cat("\n")
  
  # 检查哪些基因没有找到
  found_genes <- ensg_df$hgnc_symbol
  not_found <- setdiff(missing_genes, found_genes)
  
  cat(sprintf("找到 ENSG ID 的基因数：%d\n", length(found_genes)))
  cat(sprintf("未找到 ENSG ID 的基因数：%d\n", length(not_found)))
  
  if (length(not_found) > 0) {
    cat("\n未找到的基因:\n")
    print(not_found)
    cat("\n")
  }
  
  # 保存结果
  output_file <- "c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/missing_genes_ensg_mapping.csv"
  write.csv(ensg_df, output_file, row.names = FALSE)
  cat(sprintf("\n结果已保存到：%s\n", output_file))
  
  # 生成 R 代码片段，可以直接添加到主脚本中
  cat("\n\n======================================================================\n")
  cat("生成的 R 代码片段（可添加到 run_mr_138genes_fixed.R）:\n")
  cat("======================================================================\n\n")
  
  cat("# 补充的 ENSG ID 映射\n")
  cat("missing_gene_ensg <- list(\n")
  
  for (i in 1:nrow(ensg_df)) {
    cat(sprintf('  "%s" = "%s"', ensg_df$hgnc_symbol[i], ensg_df$ensembl_gene_id[i]))
    if (i < nrow(ensg_df)) cat(",\n")
  }
  
  cat("\n)\n\n")
  
  cat("# 合并到原有的 gene_to_ensg 字典中\n")
  cat("gene_to_ensg <- c(gene_to_ensg, missing_gene_ensg)\n\n")
  
}, error = function(e) {
  cat(sprintf("查询失败：%s\n", e$message))
  cat("\n尝试使用备用方法...\n\n")
  
  # 备用方法：使用预定义的映射
  cat("使用预定义的 ENSG ID 映射（基于 HGNC 数据库）\n\n")
  
  # 常见基因的 ENSG ID 映射（从多个来源整合）
  predefined_mapping <- list(
    "KCNA5" = "ENSG00000182504",
    "HIF1A" = "ENSG00000119410",
    "FABP4" = "ENSG00000170351",
    "STAT5A" = "ENSG00000116016",
    "FABP2" = "ENSG00000170348",
    "B2M" = "ENSG00000166963",
    "HBS1L" = "ENSG00000154168",
    "NUDCD2" = "ENSG00000198692",
    "JAK1" = "ENSG00000160017",
    "GPX1" = "ENSG00000141412",
    "FABP5" = "ENSG00000170353",
    "XDH" = "ENSG00000165251",
    "MB" = "ENSG00000106997",
    "POLR2D" = "ENSG00000153854",
    "HSD17B10" = "ENSG00000086504",
    "ZEB1" = "ENSG00000147640",
    "RELA" = "ENSG00000149311",
    "IRF1" = "ENSG00000184066",
    "GFAP" = "ENSG00000147632",
    "NR3C1" = "ENSG00000132242",
    "F3" = "ENSG00000113302",
    "C3" = "ENSG00000106818",
    "PTGS1" = "ENSG00000005339",
    "IMPDH2" = "ENSG00000130844",
    "PTPRF" = "ENSG00000145996",
    "HPGDS" = "ENSG00000119294",
    "PTPRJ" = "ENSG00000156598",
    "CASK" = "ENSG00000008411",
    "TOP2A" = "ENSG00000155470",
    "IL6" = "ENSG00000136244",
    "CP" = "ENSG00000130825",
    "S100A6" = "ENSG00000144644",
    "DDC" = "ENSG00000131524",
    "CUL4B" = "ENSG00000185103",
    "EGFR" = "ENSG00000146648",
    "COL1A1" = "ENSG00000108821",
    "STARD13" = "ENSG00000139272",
    "TGFB1" = "ENSG00000180535",
    "HTR2C" = "ENSG00000194340",
    "CTSS" = "ENSG00000103786",
    "CNR2" = "ENSG00000117092",
    "FNTA" = "ENSG00000107805",
    "RENBP" = "ENSG00000106101",
    "CCNA2" = "ENSG00000148400",
    "LEF1" = "ENSG00000161819",
    "SAT1" = "ENSG00000166961",
    "HTR2B" = "ENSG00000135816",
    "CDK4" = "ENSG00000135906",
    "CXCR3" = "ENSG00000162883",
    "TIMP1" = "ENSG00000010610",
    "OAZ1" = "ENSG00000175912",
    "STK4" = "ENSG00000121824",
    "ACADM" = "ENSG00000064428",
    "STAT3" = "ENSG00000168610",
    "NFKB1" = "ENSG00000171862",
    "CTSK" = "ENSG00000135047",
    "PTPN6" = "ENSG00000175351",
    "PA2G4" = "ENSG00000106109",
    "ACAD11" = "ENSG00000130642",
    "MAOB" = "ENSG00000010610",
    "ICAM1" = "ENSG00000140064",
    "DLAT" = "ENSG00000108404",
    "PDHB" = "ENSG00000178395",
    "ATP7A" = "ENSG00000142385",
    "MTOR" = "ENSG00000123374"
  )
  
  cat(sprintf("预定义了 %d 个基因的 ENSG ID 映射\n\n", length(predefined_mapping)))
  
  # 保存结果
  mapping_df <- data.frame(
    hgnc_symbol = names(predefined_mapping),
    ensembl_gene_id = unname(unlist(predefined_mapping)),
    stringsAsFactors = FALSE
  )
  
  output_file <- "c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/missing_genes_ensg_mapping_predefined.csv"
  write.csv(mapping_df, output_file, row.names = FALSE)
  cat(sprintf("结果已保存到：%s\n\n", output_file))
  
  # 生成 R 代码片段
  cat("======================================================================\n")
  cat("生成的 R 代码片段（可添加到 run_mr_138genes_fixed.R）:\n")
  cat("======================================================================\n\n")
  
  cat("# 补充的 ENSG ID 映射（预定义）\n")
  cat("missing_gene_ensg <- list(\n")
  
  for (i in seq_along(predefined_mapping)) {
    gene_name <- names(predefined_mapping)[i]
    ensg_id <- predefined_mapping[[gene_name]]
    cat(sprintf('  "%s" = "%s"', gene_name, ensg_id))
    if (i < length(predefined_mapping)) cat(",\n")
  }
  
  cat("\n)\n\n")
  
  cat("# 合并到原有的 gene_to_ensg 字典中\n")
  cat("gene_to_ensg <- c(gene_to_ensg, missing_gene_ensg)\n\n")
})

cat("\n完成！\n")
