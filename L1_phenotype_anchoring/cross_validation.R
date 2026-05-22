#!/usr/bin/env Rscript

# L1 交叉验证：横向（同时间点）+ 纵向（时间动态）
# 修复：HTML实体解码 + 根据log2FC推导方向

suppressPackageStartupMessages({
  library(openxlsx)
})

set.seed(42)

# HTML实体解码函数
decode_html_entities <- function(x) {
  if (is.null(x) || all(is.na(x))) return(x)
  x <- as.character(x)
  result <- character(length(x))
  for (i in seq_along(x)) {
    val <- x[i]
    if (grepl("&#\\d+;", val, perl = TRUE)) {
      matches <- gregexpr("&#(\\d+);", val, perl = TRUE)
      starts <- matches[[1]]
      if (starts[1] != -1) {
        lengths <- attr(matches[[1]], "match.length")
        for (j in rev(seq_along(starts))) {
          matched <- substr(val, starts[j], starts[j] + lengths[j] - 1)
          num <- as.integer(gsub("[^0-9]", "", matched))
          char_val <- intToUtf8(num)
          val <- paste0(
            substr(val, 1, starts[j] - 1),
            char_val,
            substr(val, starts[j] + lengths[j], nchar(val))
          )
        }
      }
    }
    result[i] <- val
  }
  return(result)
}

# 根据log2FC推导方向
logfc_to_dir <- function(logfc) {
  ifelse(!is.na(logfc), ifelse(logfc > 0, "上调", "下调"), NA)
}

base_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery"
results_dir <- file.path(base_dir, "results/L1_phenotype_anchoring")
l1_dir <- file.path(base_dir, "L1_phenotype_anchoring")

# scRNA-seq 铜死亡基因（列: A=类别, B=Gene, C=log2FC, D=P.Value, E=adj.P.Val, F=显著性, G=上调/下调）
scrna_excel <- file.path(results_dir, "L1_scRNA_GSE174574_Summary.xlsx")
scrna_raw <- read.xlsx(scrna_excel, sheet = "Cuproptosis_Genes", startRow = 3, colNames = TRUE)

scrna <- data.frame(
  gene = toupper(as.character(scrna_raw[["Gene"]])),
  log2FC = as.numeric(scrna_raw[["log2FC"]]),
  # 第7列=上调/下调（解码HTML实体），若解码失败则用log2FC推导
  direction = decode_html_entities(as.character(scrna_raw[[7]])),
  stringsAsFactors = FALSE
)
# 对解码后仍无效的，用log2FC推导
scrna$direction <- ifelse(is.na(scrna$direction) | scrna$direction == "NA" | scrna$direction == "",
  logfc_to_dir(scrna$log2FC), scrna$direction)
scrna <- scrna[!is.na(scrna$gene) & scrna$gene != "", ]

cat("scRNA-seq 方向样本:", paste(unique(scrna$direction), collapse=", "), "\n")

# GSE97537 (24h Rat) - 使用log2FC推导方向（避免编码问题）
gse97537 <- read.csv(file.path(l1_dir, "GSE97537_cuproptosis_DEGs.csv"), stringsAsFactors = FALSE)
gse97537$Human_Gene <- toupper(gse97537$Human_Gene)
gse97537$Direction <- logfc_to_dir(gse97537$log2FC)

# GSE61616 (7d Mouse) - 使用log2FC推导方向
gse61616_excel <- file.path(results_dir, "L1_Bulk_GSE61616_Summary.xlsx")
gse61616_raw <- read.xlsx(gse61616_excel, sheet = "Cuproptosis_Genes", startRow = 3, colNames = TRUE)

gse61616 <- data.frame(
  gene = toupper(as.character(gse61616_raw[[2]])),
  log2FC = as.numeric(gse61616_raw[[4]]),
  direction = logfc_to_dir(as.numeric(gse61616_raw[[4]])),
  stringsAsFactors = FALSE
)
gse61616 <- gse61616[!is.na(gse61616$gene) & gse61616$gene != "", ]

cat("GSE61616 方向样本:", paste(unique(gse61616$direction), collapse=", "), "\n")

cat("\n=== 数据加载 ===\n")
cat("  scRNA-seq:", nrow(scrna), "基因\n")
cat("  GSE97537 (24h Rat):", nrow(gse97537), "基因\n")
cat("  GSE61616 (7d Mouse):", nrow(gse61616), "基因\n")

scrna_map_dir <- setNames(scrna$direction, scrna$gene)
gse61616_map_dir <- setNames(gse61616$direction, gse61616$gene)
gse61616_map_logfc <- setNames(gse61616$log2FC, gse61616$gene)

# 横向对比
cat("\n=== 横向对比：24h 急性期 ===\n")
cat("  GSE97537 (24h Rat) vs scRNA-seq (24h Mouse)\n\n")

horizontal <- data.frame(
  Gene = gse97537$Human_Gene,
  GSE97537_Status = gse97537$Status,
  GSE97537_log2FC = gse97537$log2FC,
  GSE97537_adjP = gse97537$adj.P.Val,
  GSE97537_Dir = gse97537$Direction,
  GSE97537_Sig = gse97537$Significant,
  scRNA_Dir = scrna_map_dir[gse97537$Human_Gene],
  stringsAsFactors = FALSE
)

horizontal$Dir_Consistent <- ifelse(
  !is.na(horizontal$GSE97537_Dir) & !is.na(horizontal$scRNA_Dir) &
  horizontal$GSE97537_Dir == horizontal$scRNA_Dir,
  "一致",
  ifelse(!is.na(horizontal$GSE97537_Dir) & !is.na(horizontal$scRNA_Dir), "不一致", "NA")
)

cat("  铜死亡基因总数:", nrow(horizontal), "\n")
cat("  GSE97537检出:", sum(horizontal$GSE97537_Status == "检出"), "\n")
cat("  双平台检出:", sum(horizontal$GSE97537_Status == "检出" & !is.na(horizontal$scRNA_Dir)), "\n")
cat("  方向一致:", sum(horizontal$Dir_Consistent == "一致"), "/", sum(!is.na(horizontal$Dir_Consistent)), "\n\n")

for (i in seq_len(nrow(horizontal))) {
  if (horizontal$GSE97537_Status[i] == "检出") {
    cat(sprintf("  %-10s | 24hRat: %+7.4f (%s, %s) | scRNA: %s | %s\n",
      horizontal$Gene[i], horizontal$GSE97537_log2FC[i],
      horizontal$GSE97537_Dir[i], horizontal$GSE97537_Sig[i],
      horizontal$scRNA_Dir[i], horizontal$Dir_Consistent[i]))
  }
}

# 纵向对比
cat("\n=== 纵向对比：时间动态 ===\n")
cat("  GSE97537 (24h Rat) vs GSE61616 (7d Mouse)\n\n")

vertical <- data.frame(
  Gene = gse97537$Human_Gene,
  GSE97537_log2FC = gse97537$log2FC,
  GSE97537_adjP = gse97537$adj.P.Val,
  GSE97537_Dir = gse97537$Direction,
  GSE97537_Sig = gse97537$Significant,
  GSE61616_log2FC = gse61616_map_logfc[gse97537$Human_Gene],
  GSE61616_Dir = gse61616_map_dir[gse97537$Human_Gene],
  stringsAsFactors = FALSE
)

vertical$Dir_Consistent <- ifelse(
  !is.na(vertical$GSE97537_Dir) & !is.na(vertical$GSE61616_Dir) &
  vertical$GSE97537_Dir == vertical$GSE61616_Dir,
  "一致",
  ifelse(!is.na(vertical$GSE97537_Dir) & !is.na(vertical$GSE61616_Dir), "不一致", "NA")
)

vertical$Time_Dynamic <- ifelse(
  !is.na(vertical$GSE97537_Dir) & !is.na(vertical$GSE61616_Dir),
  ifelse(vertical$GSE97537_Dir == vertical$GSE61616_Dir, "持续响应", "方向反转"),
  ifelse(!is.na(vertical$GSE97537_Dir), "仅24h响应",
    ifelse(!is.na(vertical$GSE61616_Dir), "仅7d响应", "无响应")))

cat("  双时间点检出:", sum(!is.na(vertical$GSE97537_Dir) & !is.na(vertical$GSE61616_Dir)), "\n")
cat("  持续响应:", sum(vertical$Time_Dynamic == "持续响应"), "\n")
cat("  方向反转:", sum(vertical$Time_Dynamic == "方向反转"), "\n")
cat("  仅24h响应:", sum(vertical$Time_Dynamic == "仅24h响应"), "\n")
cat("  仅7d响应:", sum(vertical$Time_Dynamic == "仅7d响应"), "\n\n")

for (i in seq_len(nrow(vertical))) {
  cat(sprintf("  %-10s | 24h: %+7.4f (%s, %s) | 7d: %+7.4f (%s) | %s\n",
    vertical$Gene[i],
    ifelse(is.na(vertical$GSE97537_log2FC[i]), 0, vertical$GSE97537_log2FC[i]),
    ifelse(is.na(vertical$GSE97537_Dir[i]), "未检出", vertical$GSE97537_Dir[i]),
    ifelse(is.na(vertical$GSE97537_Sig[i]), "?", vertical$GSE97537_Sig[i]),
    ifelse(is.na(vertical$GSE61616_log2FC[i]), 0, vertical$GSE61616_log2FC[i]),
    ifelse(is.na(vertical$GSE61616_Dir[i]), "未检出", vertical$GSE61616_Dir[i]),
    vertical$Time_Dynamic[i]))
}

# 输出 Excel
out_file <- file.path(results_dir, "L1_Cross_Validation.xlsx")
wb <- createWorkbook()

addWorksheet(wb, "Horizontal_24h")
writeData(wb, "Horizontal_24h", horizontal)

addWorksheet(wb, "Vertical_TimeCourse")
writeData(wb, "Vertical_TimeCourse", vertical)

addWorksheet(wb, "Summary")
summary_stats <- data.frame(
  Metric = c("横向: GSE97537检出基因数", "横向: 双平台检出基因数", "横向: 方向一致性",
    "纵向: 双时间点检出基因数", "纵向: 持续响应", "纵向: 方向反转",
    "纵向: 仅24h响应", "纵向: 仅7d响应"),
  Value = c(sum(horizontal$GSE97537_Status == "检出"),
    sum(horizontal$GSE97537_Status == "检出" & !is.na(horizontal$scRNA_Dir)),
    paste0(sum(horizontal$Dir_Consistent == "一致"), "/", sum(!is.na(horizontal$Dir_Consistent))),
    sum(!is.na(vertical$GSE97537_Dir) & !is.na(vertical$GSE61616_Dir)),
    sum(vertical$Time_Dynamic == "持续响应"),
    sum(vertical$Time_Dynamic == "方向反转"),
    sum(vertical$Time_Dynamic == "仅24h响应"),
    sum(vertical$Time_Dynamic == "仅7d响应")),
  stringsAsFactors = FALSE
)
writeData(wb, "Summary", summary_stats)

saveWorkbook(wb, out_file, overwrite = TRUE)

cat("\n✅ 交叉验证完成！\n")
cat("  输出:", out_file, "\n")
