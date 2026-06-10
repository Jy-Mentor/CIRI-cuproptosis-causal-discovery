# ============================================================================
# β-石竹烯（β-Caryophyllene）靶点清洗去重分析
# 合并多数据库预测结果，标准化基因名，去重输出
# ============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressMessages({
  library(dplyr)
})

# ============================================================================
# 0. 路径配置
# ============================================================================
set.seed(20250417)
project_dir <- "c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
bcp_dir <- file.path(project_dir, "BCP_target_prediction")
output_dir <- file.path(project_dir, "BCP_Target_Screening_Academic")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

analysis_date <- format(Sys.Date(), "%Y-%m-%d")
cat("===== β-石竹烯靶点清洗去重分析 =====\n")
cat("分析日期:", analysis_date, "\n\n")

# ============================================================================
# 1. 读取原始数据
# ============================================================================
cat("===== 1. 读取原始数据 =====\n")

swiss_file <- file.path(bcp_dir, "swiss_target_prediction_results.csv")
pharm_file <- file.path(bcp_dir, "pharm_mapper_prediction_results.csv")
targetnet_file <- file.path(bcp_dir, "target_net_prediction_results.csv")

stopifnot("SwissTargetPrediction文件不存在" = file.exists(swiss_file))
stopifnot("PharmMapper文件不存在" = file.exists(pharm_file))
stopifnot("TargetNet文件不存在" = file.exists(targetnet_file))

swiss_raw <- read.csv(swiss_file, stringsAsFactors = FALSE)
pharm_raw <- read.csv(pharm_file, stringsAsFactors = FALSE)
targetnet_raw <- read.csv(targetnet_file, stringsAsFactors = FALSE)

cat("SwissTargetPrediction: ", nrow(swiss_raw), "条\n")
cat("PharmMapper: ", nrow(pharm_raw), "条\n")
cat("TargetNet: ", nrow(targetnet_raw), "条\n\n")

# ============================================================================
# 2. 基因名标准化 (HGNC标准)
# ============================================================================
cat("===== 2. 基因名标准化 (HGNC) =====\n")

fix_gene_names <- function(target_list) {
  replacements <- c(
    "P450 3A4" = "CYP3A4",
    "P450 2C9" = "CYP2C9",
    "DRD2" = "DRD2", "DRD3" = "DRD3",
    "ADRA2A" = "ADRA2A", "CHRM1" = "CHRM1",
    "HTR2A" = "HTR2A", "HTR1A" = "HTR1A",
    "OPRD1" = "OPRD1", "OPRK1" = "OPRK1",
    "CNR1" = "CNR1", "CNR2" = "CNR2",
    "PTGS2" = "PTGS2", "NR3C1" = "NR3C1",
    "ESR1" = "ESR1", "PPARG" = "PPARG",
    "HMGCR" = "HMGCR", "ACHE" = "ACHE",
    "BCHE" = "BCHE", "MAOA" = "MAOA",
    "MAOB" = "MAOB", "ADRB2" = "ADRB2"
  )
  for (i in seq_along(replacements)) {
    target_list <- gsub(names(replacements)[i], replacements[i], target_list, fixed = TRUE)
  }
  return(target_list)
}

swiss_raw$Target <- fix_gene_names(swiss_raw$Target)
pharm_raw$Target <- fix_gene_names(pharm_raw$Target)
targetnet_raw$Target <- fix_gene_names(targetnet_raw$Target)

swiss_raw <- swiss_raw[nchar(swiss_raw$Target) > 0, ]
pharm_raw <- pharm_raw[nchar(pharm_raw$Target) > 0, ]
targetnet_raw <- targetnet_raw[nchar(targetnet_raw$Target) > 0, ]

cat("标准化后 SwissTargetPrediction: ", nrow(swiss_raw), "条\n")
cat("标准化后 PharmMapper: ", nrow(pharm_raw), "条\n")
cat("标准化后 TargetNet: ", nrow(targetnet_raw), "条\n\n")

# ============================================================================
# 3. 添加数据库来源标记
# ============================================================================
cat("===== 3. 添加数据库来源标记 =====\n")

swiss_raw <- swiss_raw %>% rename(Probability_Swiss = Probability) %>% mutate(Database = "SwissTargetPrediction")
pharm_raw <- pharm_raw %>% rename(Probability_Pharm = FitScore) %>% mutate(Database = "PharmMapper")
targetnet_raw <- targetnet_raw %>% rename(Probability_TargetNet = BioactivityProbability) %>% mutate(Database = "TargetNet")

# ============================================================================
# 4. 合并所有数据并去重
# ============================================================================
cat("===== 4. 合并数据并去重 =====\n")

all_targets <- data.frame(
  Target = character(),
  Probability_Swiss = numeric(),
  Probability_Pharm = numeric(),
  Probability_TargetNet = numeric(),
  stringsAsFactors = FALSE
)

all_target_names <- unique(c(
  swiss_raw$Target,
  pharm_raw$Target,
  targetnet_raw$Target
))

cat("合并后唯一靶点数: ", length(all_target_names), "\n\n")

# ============================================================================
# 5. 构建完整关联矩阵
# ============================================================================
cat("===== 5. 构建靶点-数据库关联矩阵 =====\n")

target_matrix <- data.frame(
  Target = sort(all_target_names),
  SwissTargetPrediction = NA_real_,
  PharmMapper = NA_real_,
  TargetNet = NA_real_,
  stringsAsFactors = FALSE
)

for (i in 1:nrow(target_matrix)) {
  tgt <- target_matrix$Target[i]

  swiss_val <- swiss_raw$Probability_Swiss[swiss_raw$Target == tgt]
  pharm_val <- pharm_raw$Probability_Pharm[pharm_raw$Target == tgt]
  targetnet_val <- targetnet_raw$Probability_TargetNet[targetnet_raw$Target == tgt]

  target_matrix$SwissTargetPrediction[i] <- ifelse(length(swiss_val) > 0, swiss_val[1], NA)
  target_matrix$PharmMapper[i] <- ifelse(length(pharm_val) > 0, pharm_val[1], NA)
  target_matrix$TargetNet[i] <- ifelse(length(targetnet_val) > 0, targetnet_val[1], NA)
}

target_matrix$DatabaseCount <- rowSums(!is.na(target_matrix[, c("SwissTargetPrediction", "PharmMapper", "TargetNet")]))

target_matrix <- target_matrix[order(-target_matrix$DatabaseCount, target_matrix$Target), ]

cat("靶点-数据库关联矩阵:\n")
print(target_matrix)
cat("\n")

# ============================================================================
# 6. 保存结果
# ============================================================================
cat("===== 6. 保存结果 =====\n")

all_targets_file <- file.path(output_dir, "BCP_All_Targets_Deduplicated.csv")
write.csv(target_matrix, all_targets_file, row.names = FALSE)
cat("所有靶点已保存至:", all_targets_file, "\n")

# 输出纯基因列表
gene_list <- data.frame(Target = target_matrix$Target)
gene_list_file <- file.path(output_dir, "BCP_Target_Gene_List.csv")
write.csv(gene_list, gene_list_file, row.names = FALSE)
cat("基因列表已保存至:", gene_list_file, "\n")

# 输出报告
report_content <- paste0("
================================================================================
β-石竹烯（β-Caryophyllene）靶点清洗去重报告
================================================================================

分析日期: ", analysis_date, "
数据来源:
  - SwissTargetPrediction: ", nrow(swiss_raw), "条
  - PharmMapper: ", nrow(pharm_raw), "条
  - TargetNet: ", nrow(targetnet_raw), "条

去重后靶点总数: ", nrow(target_matrix), "

靶点列表:
", paste(sort(target_matrix$Target), collapse = "\n"), "

关联矩阵统计:
  - 3个数据库支持: ", sum(target_matrix$DatabaseCount == 3), "个
  - 2个数据库支持: ", sum(target_matrix$DatabaseCount == 2), "个
  - 1个数据库支持: ", sum(target_matrix$DatabaseCount == 1), "个

================================================================================
", Sys.time(), "
================================================================================
")

report_file <- file.path(output_dir, "BCP_Target_Deduplication_Report.txt")
writeLines(report_content, report_file)
cat("报告已保存至:", report_file, "\n\n")

# ============================================================================
# 7. 最终统计
# ============================================================================
cat("===== 分析完成 =====\n")
cat("去重后靶点总数: ", nrow(target_matrix), "个\n")
cat("\n靶点列表:\n")
cat(paste(sort(target_matrix$Target), collapse = ", "), "\n")
cat("\n所有结果已保存至:", output_dir, "\n")
