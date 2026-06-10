#!/usr/bin/env Rscript
# 检查IEU OpenGWAS中可用的GTEx数据集

suppressPackageStartupMessages({
  library(ieugwasr)
  library(TwoSampleMR)
})

# 设置token
token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3NzAxMjQ5NiwiZXhwIjoxNzc4MjIyMDk2fQ.sz-zieAniOZLfPGDiQRo3W6z5PQ6HfEvsfSjvKCVLOlUJumnpM-jE9Se6QAqXezfH-ffEEbM0heC3pQu6MjBF1ttqaRpIdOv8S7_pq-s-HLVQUJZ3Ge5Ort66N89riaj9nRS-FdGidJAft58CEmjk-Q1lpJpJbOSGJiKCRQ_Man5vbCBiTH3p_e490wbapqqpxaXWs4ki-xYUDlQXuc6u2lZIh7_6TSei6Chmew0EWy8wcVgSTJGAmAiT8jrFTq5ydxsPvzZMnMZgEezAfrHIi9i3zJqbTzggs8Okqy2owfzLzcE7MKt_EcUSHYvMzRsL3YcbILipkodoy-sbzNmpw"
Sys.setenv(OPENGWAS_JWT = token)

cat("Checking OpenGWAS API access...\n")
tryCatch({
  user_info <- ieugwasr::user()
  cat("User:", user_info$user[[1]], "\n")
  cat("Token valid!\n\n")
}, error = function(e) {
  cat("Token error:", conditionMessage(e), "\n")
  quit(status = 1)
})

# 获取所有可用的GWAS数据集信息
cat("Fetching GWAS info...\n")
all_info <- gwasinfo()
cat("Total datasets:", nrow(all_info), "\n\n")

# 搜索GTEx相关数据集
cat("Searching for GTEx datasets...\n")
gtex_datasets <- all_info[grepl("gtex|GTEx", all_info$trait, ignore.case = TRUE), ]
cat("Found", nrow(gtex_datasets), "GTEx-related datasets\n\n")

# 显示脑组织相关的数据集
cat("Brain tissue datasets:\n")
brain_keywords <- c("brain", "cortex", "hippocampus", "cerebellum", "putamen", "amygdala")
for (keyword in brain_keywords) {
  matches <- gtex_datasets[grepl(keyword, gtex_datasets$trait, ignore.case = TRUE), ]
  if (nrow(matches) > 0) {
    cat("\n", toupper(keyword), "\n")
    for (i in 1:min(nrow(matches), 5)) {
      cat("  ID:", matches$id[i], "\n")
      cat("  Trait:", matches$trait[i], "\n")
      cat("  Sample:", ifelse(!is.na(matches$sample_size[i]), matches$sample_size[i], "N/A"), "\n\n")
    }
  }
}

# 测试几个特定的数据集ID
test_ids <- c("ieu-b-4171", "ieu-b-4178", "ieu-b-4166", "eqtl-a-ENSG00000137730")
cat("\n\nTesting specific dataset IDs...\n")
for (id in test_ids) {
  cat("\nTesting:", id, "\n")
  tryCatch({
    # 尝试获取数据集信息
    info <- all_info[all_info$id == id, ]
    if (nrow(info) > 0) {
      cat("  Found:", info$trait[1], "\n")
      cat("  Year:", ifelse(!is.na(info$year), info$year, "N/A"), "\n")
      
      # 尝试提取instruments
      test_dat <- extract_instruments(outcomes = id, p1 = 5e-8)
      if (!is.null(test_dat) && nrow(test_dat) > 0) {
        cat("  SUCCESS: Extracted", nrow(test_dat), "instruments\n")
        cat("  Columns:", paste(names(test_dat), collapse = ", "), "\n")
        if ("exposure" %in% names(test_dat)) {
          cat("  Example exposure:", test_dat$exposure[1], "\n")
        }
      } else {
        cat("  No instruments found (may need different p-value threshold)\n")
      }
    } else {
      cat("  Dataset ID not found in catalog\n")
    }
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })
}

cat("\n\nDone!\n")
