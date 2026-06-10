#!/usr/bin/env Rscript
# 搜索脑组织相关数据集

suppressPackageStartupMessages(library(ieugwasr))

token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3NzAxMjQ5NiwiZXhwIjoxNzc4MjIyMDk2fQ.sz-zieAniOZLfPGDiQRo3W6z5PQ6HfEvsfSjvKCVLOlUJumnpM-jE9Se6QAqXezfH-ffEEbM0heC3pQu6MjBF1ttqaRpIdOv8S7_pq-s-HLVQUJZ3Ge5Ort66N89riaj9nRS-FdGidJAft58CEmjk-Q1lpJpJbOSGJiKCRQ_Man5vbCBiTH3p_e490wbapqqpxaXWs4ki-xYUDlQXuc6u2lZIh7_6TSei6Chmew0EWy8wcVgSTJGAmAiT8jrFTq5ydxsPvzZMnMZgEezAfrHIi9i3zJqbTzggs8Okqy2owfzLzcE7MKt_EcUSHYvMzRsL3YcbILipkodoy-sbzNmpw"
Sys.setenv(OPENGWAS_JWT = token)

cat("Fetching GWAS info...\n")
all_info <- gwasinfo()

cat("\n=== Brain-related datasets ===\n")
brain_datasets <- all_info[grepl("brain|cortex|cerebellum|hippocampus|putamen", 
                                  all_info$trait, ignore.case = TRUE), ]

cat("Found", nrow(brain_datasets), "datasets\n\n")

for (i in 1:min(nrow(brain_datasets), 30)) {
  cat(i, ". ID:", brain_datasets$id[i], "\n")
  cat("   Trait:", substr(brain_datasets$trait[i], 1, 80), "\n")
  cat("   Year:", ifelse(!is.na(brain_datasets$year[i]), brain_datasets$year[i], "N/A"), "\n")
  cat("   Sample:", ifelse(!is.na(brain_datasets$sample_size[i]), brain_datasets$sample_size[i], "N/A"), "\n\n")
}

cat("\n\nDone!\n")
