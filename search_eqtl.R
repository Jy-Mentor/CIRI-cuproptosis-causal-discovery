#!/usr/bin/env Rscript
# 搜索eQTL相关数据集

suppressPackageStartupMessages(library(ieugwasr))

token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3NzAxMjQ5NiwiZXhwIjoxNzc4MjIyMDk2fQ.sz-zieAniOZLfPGDiQRo3W6z5PQ6HfEvsfSjvKCVLOlUJumnpM-jE9Se6QAqXezfH-ffEEbM0heC3pQu6MjBF1ttqaRpIdOv8S7_pq-s-HLVQUJZ3Ge5Ort66N89riaj9nRS-FdGidJAft58CEmjk-Q1lpJpJbOSGJiKCRQ_Man5vbCBiTH3p_e490wbapqqpxaXWs4ki-xYUDlQXuc6u2lZIh7_6TSei6Chmew0EWy8wcVgSTJGAmAiT8jrFTq5ydxsPvzZMnMZgEezAfrHIi9i3zJqbTzggs8Okqy2owfzLzcE7MKt_EcUSHYvMzRsL3YcbILipkodoy-sbzNmpw"
Sys.setenv(OPENGWAS_JWT = token)

cat("Fetching all GWAS info (this may take a minute)...\n")
all_info <- gwasinfo()

cat("\nTotal datasets:", nrow(all_info), "\n")

# 搜索eQTL相关
cat("\n=== eQTL datasets ===\n")
eqtl_datasets <- all_info[grepl("eqtl|eQTL", all_info$trait, ignore.case = TRUE), ]
cat("Found", nrow(eqtl_datasets), "eQTL datasets\n")

if (nrow(eqtl_datasets) > 0) {
  # 显示前20个
  for (i in 1:min(nrow(eqtl_datasets), 20)) {
    cat("\n", i, ". ID:", eqtl_datasets$id[i], "\n")
    cat("   Trait:", substr(eqtl_datasets$trait[i], 1, 80), "\n")
    cat("   Sample:", ifelse(!is.na(eqtl_datasets$sample_size[i]), eqtl_datasets$sample_size[i], "N/A"), "\n")
  }
}

# 搜索基因表达相关
cat("\n\n=== Gene expression datasets ===\n")
expr_datasets <- all_info[grepl("expression|transcriptome|RNA", all_info$trait, ignore.case = TRUE), ]
cat("Found", nrow(expr_datasets), "gene expression datasets\n")

if (nrow(expr_datasets) > 0) {
  for (i in 1:min(nrow(expr_datasets), 10)) {
    cat("\n", i, ". ID:", expr_datasets$id[i], "\n")
    cat("   Trait:", substr(expr_datasets$trait[i], 1, 80), "\n")
  }
}

# 搜索特定基因
cat("\n\n=== Searching for our 9 hub genes ===\n")
genes <- c("NFKB1", "FDX1", "STAT3", "HIF1A", "HMOX1", "GPX4", "TNF", "IL6", "AGER")
for (gene in genes) {
  gene_matches <- all_info[grepl(gene, all_info$trait, ignore.case = TRUE), ]
  cat("\n", gene, ": Found", nrow(gene_matches), "datasets\n")
  if (nrow(gene_matches) > 0) {
    for (i in 1:min(nrow(gene_matches), 3)) {
      cat("   -", gene_matches$id[i], ":", substr(gene_matches$trait[i], 1, 60), "\n")
    }
  }
}

# 搜索组织相关
cat("\n\n=== Tissue-specific datasets ===\n")
tissues <- c("blood", "brain", "liver", "lung", "muscle", "adipose", "skin")
for (tissue in tissues) {
  tissue_matches <- all_info[grepl(tissue, all_info$trait, ignore.case = TRUE), ]
  cat(tissue, ":", nrow(tissue_matches), "datasets\n")
}

cat("\n\nDone!\n")
