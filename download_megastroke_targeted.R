#!/usr/bin/env Rscript
# ================================================================================
# 根据目标基因 eQTL 位置下载 MEGASTROKE SNP
# 只下载需要的 SNP，避免下载完整 GWAS 数据
# ================================================================================

library(dplyr)
library(data.table)
library(httr)

cat("======================================================================\n")
cat("批量下载目标基因的 MEGASTROKE SNP\n")
cat("======================================================================\n\n")

# JWT Token
JWT_TOKEN <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3ODMwNDA4MywiZXhwIjoxNzc5NTEzNjgzfQ.ZtcIUEx_xYtrVD_EE-UboKyLlC-lZBq2pjn-iYhJzxocqHdA-02K9n_Qbw-5ngQ07GHjjIYqVtmkZfJ3OJl1yI-tOMBBFzVKe0nkwDcB6-yBgjgBaxVm8vq_pbNrMwy_ZezY5ys9jx7I8T4bZYg9KeUbSwj04OfNP82kGcKXIOErXXVy-Ie3dbUogDRSjnCT-_32yNQxuWpiyYnPWSrWQbQ2HlUQiiDTdFGzWeJJKfSRvjQzdp5g3nccxht0m5A0UsPCdvkyHFEpvPVZ-NpjCjkgy8GbZBv4cmDMSc5JJL6HLO0eV508SRKxMdp-gL6qVdhGiJ2i9XmZQE27-aq_7g"

# 配置
EXPOSURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/exposure_matched/matched_data"
OUTPUT_FILE <- "D:/EQTL/mr_results_megastroke/megastroke_target_genes.csv"

# 步骤 1: 收集所有 SNP 位置
cat("步骤 1: 从暴露数据收集 SNP 位置...\n")
cat("----------------------------------------------------------------------\n")

exposure_files <- list.files(EXPOSURE_DIR, pattern = "_exposure\\.csv$", full.names = TRUE)
cat(sprintf("找到 %d 个暴露数据文件\n", length(exposure_files)))

# 收集所有位置
all_positions <- data.frame(CHR = character(), BP = integer(), stringsAsFactors = FALSE)

for (exposure_file in exposure_files) {
  tryCatch({
    data <- fread(exposure_file, stringsAsFactors = FALSE)
    
    positions <- data.frame(
      CHR = as.character(data$CHR),
      BP = as.integer(data$BP),
      stringsAsFactors = FALSE
    )
    
    all_positions <- bind_rows(all_positions, positions)
    
    if (nrow(all_positions) %% 500 == 0) {
      cat(sprintf("  已收集 %d 个位置...\n", nrow(all_positions)))
    }
    
  }, error = function(e) {
    cat(sprintf("  跳过 %s: %s\n", basename(exposure_file), e$message))
  })
}

# 去重
all_positions <- distinct(all_positions)
cat(sprintf("\n✓ 共 %d 个唯一位置\n\n", nrow(all_positions)))

# 按染色体分组
chr_list <- unique(all_positions$CHR)
cat(sprintf("覆盖 %d 条染色体: %s\n\n", length(chr_list), paste(sort(chr_list), collapse = ", ")))

# 步骤 2: 按染色体区域批量下载
cat("步骤 2: 从 MEGASTROKE API 下载 SNP...\n")
cat("----------------------------------------------------------------------\n")

API_BASE <- "https://api.openGWAS.io"
DATASET_ID <- "ebi-a-GCST006908"

headers <- add_headers(
  Authorization = paste("Bearer", JWT_TOKEN)
)

# 存储结果
all_snps <- list()
request_count <- 0

# 合并相邻位置的函数
merge_regions <- function(positions, window_size = 500000) {
  if (length(positions) == 0) return(list())
  
  positions <- sort(positions)
  regions <- list()
  region_start <- positions[1]
  region_end <- positions[1]
  
  for (pos in positions[-1]) {
    if (pos - region_end < window_size) {
      region_end <- pos
    } else {
      regions[[length(regions) + 1]] <- c(region_start, region_end)
      region_start <- pos
      region_end <- pos
    }
  }
  
  regions[[length(regions) + 1]] <- c(region_start, region_end)
  return(regions)
}

# 对每条染色体下载
for (chr_ in sort(chr_list, method = "radix")) {
  positions <- all_positions$BP[all_positions$CHR == chr_]
  
  # 合并区域
  regions <- merge_regions(positions, window_size = 500000)  # 500kb 窗口
  
  cat(sprintf("\n染色体 %s: %d 个区域\n", chr_, length(regions)))
  
  for (i in seq_along(regions)) {
    region <- regions[[i]]
    region_start <- region[1]
    region_end <- region[2]
    
    tryCatch({
      # 查询该区域的 SNP
      response <- GET(
        url = paste0(API_BASE, "/api/v1/associations"),
        query = list(
          dataset_id = DATASET_ID,
          chromosome = paste0("chr", chr_),
          start = region_start,
          end = region_end,
          page_size = 10000
        ),
        headers = headers,
        timeout = 60
      )
      
      request_count <<- request_count + 1
      
      if (status_code(response) == 200) {
        snps <- content(response, "parsed")
        
        if (length(snps) > 0) {
          all_snps[[length(all_snps) + 1]] <- snps
          cat(sprintf("  区域 %d: %d 个 SNP\n", i, length(snps)))
        }
      } else {
        cat(sprintf("  区域 %d: HTTP %d\n", i, status_code(response)))
      }
      
      # 避免请求过快
      if (request_count %% 10 == 0) {
        Sys.sleep(1)
      } else {
        Sys.sleep(0.3)
      }
      
    }, error = function(e) {
      cat(sprintf("  区域 %d: 错误 - %s\n", i, e$message))
    })
  }
}

# 合并所有 SNP
if (length(all_snps) > 0) {
  all_snps_df <- do.call(rbind, lapply(all_snps, function(snps) {
    do.call(rbind, lapply(snps, function(snp) {
      data.frame(
        SNP = snp$snp,
        CHR = snp$chromosome,
        BP = snp$position,
        EFFECT_ALLELE = snp$effect_allele,
        OTHER_ALLELE = snp$other_allele,
        EAF = snp$eaf,
        BETA = snp$beta,
        SE = snp$se,
        PVAL = snp$pval,
        N = snp$n,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  cat(sprintf("\n✓ 共下载 %d 个 SNP\n", nrow(all_snps_df)))
  cat(sprintf("✓ 发送 %d 个请求\n\n", request_count))
  
  # 步骤 3: 保存数据
  cat("步骤 3: 保存数据...\n")
  cat("----------------------------------------------------------------------\n")
  
  # 去重
  original_count <- nrow(all_snps_df)
  all_snps_df <- distinct(all_snps_df, SNP, .keep_all = TRUE)
  cat(sprintf("去重：%d -> %d 个 SNP\n", original_count, nrow(all_snps_df)))
  
  # 添加结局信息
  all_snps_df$outcome <- "Ischemic Stroke"
  all_snps_df$id.outcome <- DATASET_ID
  
  # 保存
  fwrite(all_snps_df, OUTPUT_FILE)
  
  cat(sprintf("✓ 保存到：%s\n", OUTPUT_FILE))
  cat(sprintf("✓ 文件大小：%.2f MB\n\n", file.size(OUTPUT_FILE) / 1e6))
  
  # 显示覆盖情况
  cat("覆盖的染色体:\n")
  for (chr_ in sort(unique(all_snps_df$CHR), method = "radix")) {
    count <- sum(all_snps_df$CHR == chr_)
    cat(sprintf("  染色体 %s: %d 个 SNP\n", chr_, count))
  }
  
} else {
  cat("\n✗ 未下载任何 SNP\n")
  cat("请检查:\n")
  cat("  1. JWT Token 是否有效\n")
  cat("  2. 网络连接\n")
  cat("  3. API 限制\n")
}

cat("\n======================================================================\n")
cat("完成！\n")
cat("======================================================================\n")
