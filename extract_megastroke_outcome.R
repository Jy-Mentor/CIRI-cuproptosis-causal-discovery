#!/usr/bin/env Rscript
# 从IEU OpenGWAS提取MEGASTROKE结局数据
# 正确的GWAS ID: ebi-a-GCST006908

library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(data.table)

# 设置IEU API token
# 方法1: 从.Renviron读取
if (file.exists(".Renviron")) {
  lines <- readLines(".Renviron")
  for (line in lines) {
    if (grepl("^IEUGWASR_TOKEN=", line)) {
      token <- sub("^IEUGWASR_TOKEN=", "", line)
      Sys.setenv(IEUGWASR_TOKEN = token)
      cat("IEU API token已设置\n")
      break
    }
  }
}

# 验证token
if (Sys.getenv("IEUGWASR_TOKEN") == "") {
  stop("错误: IEUGWASR_TOKEN环境变量未设置")
}

# ============================================
# 步骤 1：读取暴露数据，获取 SNP 列表
# ============================================
exposure_file <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.01.xlsx"

if (!file.exists(exposure_file)) {
  stop("错误：暴露文件不存在: ", exposure_file)
}

exposure_data <- readxl::read_excel(exposure_file)
cat("暴露数据总行数:", nrow(exposure_data), "\n")

# 统一基因名大写（列名是'gene'不是'Gene'）
exposure_data$gene_upper <- toupper(exposure_data$gene)

# 15个基因
genes <- c("NFKB1", "STAT3", "HIF1A", "HSPA5", "HMOX1",
           "RELA", "NFE2L2", "CP", "LIAS", "IKBKB",
           "JAK1", "PARP1", "CASP8", "MTOR", "PTPRC")

# 筛选目标基因
exp_snps <- exposure_data[exposure_data$gene_upper %in% genes, ]
cat("目标基因SNP数:", nrow(exp_snps), "\n")

# 获取唯一 rsID（用于查询 GWAS）
rsid_list <- unique(exp_snps$SNP)
rsid_list <- rsid_list[!is.na(rsid_list) & rsid_list != ""]
cat("暴露 SNP 总数:", length(rsid_list), "\n")

if (length(rsid_list) == 0) {
  stop("错误：没有有效的rsID")
}

# ============================================
# 步骤 2：分批提取 MEGASTROKE 结局数据
# ============================================
batch_size <- 800
n_batches <- ceiling(length(rsid_list) / batch_size)

outcome_list <- list()

for(i in 1:n_batches) {
  start <- (i-1) * batch_size + 1
  end <- min(i * batch_size, length(rsid_list))
  batch_snps <- rsid_list[start:end]
  
  cat("\n正在提取批次", i, "/", n_batches, "(", length(batch_snps), "个 SNP)...\n")
  
  out <- tryCatch({
    extract_outcome_data(
      snps = batch_snps,
      outcomes = "ebi-a-GCST006908",  # MEGASTROKE 正确 ID
      proxies = FALSE                  # 不需要代理 SNP
    )
  }, error = function(e) {
    cat("  批次", i, "失败:", conditionMessage(e), "\n")
    NULL
  })
  
  if(!is.null(out) && nrow(out) > 0) {
    outcome_list[[i]] <- out
    cat("  ✓ 成功提取", nrow(out), "行\n")
  } else {
    cat("  该批次无匹配 SNP\n")
  }
  
  # 避免 API 限流，暂停 1 秒
  if(i < n_batches) Sys.sleep(1)
}

# ============================================
# 步骤 3：合并并保存
# ============================================
if(length(outcome_list) == 0) {
  stop("错误：未提取到任何结局数据！请检查 SNP 格式或网络连接。")
}

outcome_data <- bind_rows(outcome_list)
cat("\n========== 提取完成 ==========\n")
cat("总匹配 SNP 数:", nrow(outcome_data), "\n")
cat("唯一 SNP 数:", length(unique(outcome_data$SNP)), "\n")
cat("结局列名:", paste(names(outcome_data), collapse = ", "), "\n")

# 保存为 CSV
output_dir <- "D:/EQTL/mr_results_megastroke"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fwrite(outcome_data, file.path(output_dir, "megastroke_outcome.csv"))
cat("\n已保存到:", file.path(output_dir, "megastroke_outcome.csv"), "\n")

# ============================================
# 步骤 4：快速验证（确认是 GWAS 不是 eQTL）
# ============================================
cat("\n========== 数据验证 ==========\n")
cat("id.outcome 示例:", paste(unique(outcome_data$id.outcome)[1:3], collapse = ", "), "\n")
cat("outcome 列示例:", paste(unique(outcome_data$outcome)[1:3], collapse = ", "), "\n")

if("samplesize.outcome" %in% names(outcome_data)) {
  ss_range <- range(outcome_data$samplesize.outcome, na.rm = TRUE)
  cat("samplesize.outcome 范围:", ss_range[1], "-", ss_range[2], "\n")
  
  # MEGASTROKE 总样本量 ~40万
  if(ss_range[2] > 100000) {
    cat("✅ 验证通过：这是 GWAS 数据（大样本，符合 MEGASTROKE）\n")
  } else {
    cat("⚠️ 警告：样本量过小，可能仍是 eQTL 数据\n")
  }
} else {
  cat("ℹ️ 无 samplesize 信息，无法验证\n")
}

cat("\n下一步：修改 MR 脚本，使用新的结局文件\n")
cat("outcome_file <- '", file.path(output_dir, "megastroke_outcome.csv"), "'\n")
