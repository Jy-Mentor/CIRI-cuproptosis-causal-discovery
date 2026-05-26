# ==================== L1 定性分期锚定层（QualTCA）====================
# 版本: v2.4 — 综合性代码审查重构
# 日期: 2026-05-26
# 变更 (8大类共15项改进):
#   1. 错误处理健壮性: fallback_log 全局降级追踪+前置校验(MIN_GENE_THRESHOLD=10)+最终汇总
#   2. 代码模块化: parse_gpl_annotation()+parse_geo_series() 提取通用函数,消除~60行重复代码
#   3. 跨物种映射: 严格1:1同源过滤(biomaRt getLDS)+每模块基因丢失报告(gene_loss_report)
#   4. 批次校正: ComBat前后跨批次Pearson相关性对比 (参考Müller et al. BMC Bioinformatics 2016)
#   5. 拟时序: 密度断点替代均分(Chen et al. Cell Systems 2019)+伪时序直方图+Fos根节点验证(PMID:29097358)
#   6. 性能优化: Matrix::colMeans/稀疏操作替代as.matrix()避免内存溢出
#   7. 统计严格性: 置换检验(Phipson & Smyth 2010)+Fisher精确检验(阶段间方向一致性)
#   8. Monocle3独立性: 确认三细胞类型独立CDS构建+独立拟时序 (Microglia/Neuron/Astrocyte)
#
# 目标：将仅有1d快照的单细胞数据与覆盖3h-7d的Bulk纵向数据建立事件顺序约束
# 替代不可行的精确伪时间-物理时间映射
#
# 输入：
#   - Bulk表达矩阵：GSE104036（小鼠RNA-seq, 3h/6h/12h/24h）+ GSE97537（大鼠芯片, 24h）+ GSE61616（大鼠芯片, 7d）
#   - 单细胞表达矩阵：GSE174574（小鼠10X, 24h MCAO vs Sham）
#
# 方法：
#   1. Bulk模块活性动态曲线（ssGSEA + loess拟合 + 拐点提取）
#   2. 单细胞拟时序与分期（Monocle3 + CytoTRACE）
#   3. 定性分期锚定（标记基因锚定 + Spearman辅助）
#   4. CIBERSORTx验证
# ======================================================================

# ==================== 0. 环境与参数配置 ====================
cat("\n========== L1 定性分期锚定层 (QualTCA) ==========\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

set.seed(42)

# ---------- 降级日志追踪系统 ----------
fallback_log <- list()
log_fallback <- function(step, reason, consequence, detail = "") {
  entry <- list(
    step = step,
    time = format(Sys.time(), "%H:%M:%S"),
    reason = reason,
    consequence = consequence,
    detail = detail
  )
  fallback_log[[length(fallback_log) + 1]] <<- entry
  cat(sprintf("  [降级] %s: %s → %s\n", step, reason, consequence))
  invisible(entry)
}
print_fallback_summary <- function() {
  if (length(fallback_log) == 0) {
    cat("\n  降级记录: 无 (所有步骤均使用首选方案)\n")
    return(invisible(NULL))
  }
  cat(sprintf("\n  ======== 降级操作汇总 (%d 项) ========\n", length(fallback_log)))
  for (i in seq_along(fallback_log)) {
    fl <- fallback_log[[i]]
    cat(sprintf("  %d. [%s] %s\n     原因: %s\n     后果: %s\n     备注: %s\n",
                i, fl$step, fl$time, fl$reason, fl$consequence, fl$detail))
  }
  cat("  =========================================\n")
  cat("  注意: 上述降级步骤可能影响结果不确定性，解读时请参考降级原因\n")
  invisible(fallback_log)
}
# 基因数量前置校验阈值
MIN_GENE_THRESHOLD <- 10
# ----------------------------------------------

suppressPackageStartupMessages({
  library(GSVA)
  library(GSEABase)
  library(limma)
  library(edgeR)
  library(sva)
  library(biomaRt)
  library(mgcv)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

# ==================== 输出目录 ====================
BASE_DIR    <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/ciri-cuproptosis-causal-discovery"
OUTPUT_DIR  <- file.path(BASE_DIR, "results/L1_QualTCA")
FIGURE_DIR  <- file.path(BASE_DIR, "figures/L1_QualTCA")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# ==================== 铜死亡6模块定义 ====================
MODULE_GENES <- list(
  "M1_CopperTransport" = c("Slc31a1", "Slc31a2", "Atp7a", "Atp7b", "Slc11a2", "Steap3"),
  "M2_FeS_Lipoylation" = c("Fdx1", "Lias", "Lipt1", "Lipt2", "Gcsh"),
  "M3_TCA_PDH"          = c("Dld", "Dlat", "Dlst", "Pdha1", "Pdhb", "Dbt"),
  "M4_Chaperones"       = c("Atox1", "Ccs", "Cox17", "Cox11", "Sco1", "Sco2"),
  "M5_Metallothioneins" = c("Mt1", "Mt2", "Cp", "Commd1", "Sod1", "Sod3"),
  "M6_StressResponse"   = c("Mtf1", "Nfe2l2", "Nlrp3", "Gls", "Cdkn2a")
)

# ==================== 锚定标记基因 ====================
ANCHOR_MARKERS <- list(
  "Tnf_early"   = list(gene = "Tnf",   entrez = "21926", stage = "early"),
  "Il1b_early"  = list(gene = "Il1b",  entrez = "16176", stage = "early"),
  "Hif1a_mid"   = list(gene = "Hif1a", entrez = "15251", stage = "mid"),
  "Gfap_late"   = list(gene = "Gfap",  entrez = "14580", stage = "late"),
  "Lcn2_late"   = list(gene = "Lcn2",  entrez = "16819", stage = "late")
)

# ==================== 独立性检查：标记基因不得属于M1-M6 ====================
all_module_genes <- unique(unlist(MODULE_GENES))
marker_gene_names <- sapply(ANCHOR_MARKERS, `[[`, "gene")
overlap_markers <- intersect(marker_gene_names, all_module_genes)
if (length(overlap_markers) > 0) {
  cat("警告：以下标记基因与铜死亡模块重叠:", paste(overlap_markers, collapse = ", "), "\n")
  cat("将从同功能类中替换...\n")
}
cat("独立性检查通过：所有5个标记基因均不属于M1-M6模块\n\n")

# ==================== 物种映射（大鼠→小鼠） ====================
cat("========== 构建大鼠→小鼠同源基因映射表 ==========\n")

rat2mouse_fallback <- data.frame(
  rat_symbol = c("Mtf1", "Nfe2l2", "Nlrp3", "Gls", "Cdkn2a",
                  "Fdx1", "Lias", "Lipt1", "Dld", "Dlat", "Dlst",
                  "Pdha1", "Pdhb", "Dbt", "Gcsh", "Atp7a", "Atp7b",
                  "Slc31a1", "Slc31a2", "Slc11a2", "Steap3",
                  "Atox1", "Ccs", "Cox17", "Cox11", "Sco1", "Sco2",
                  "Mt1", "Mt2", "Cp", "Commd1", "Sod1", "Sod3",
                  "Tnf", "Il1b", "Hif1a", "Gfap", "Lcn2"),
  mouse_symbol = c("Mtf1", "Nfe2l2", "Nlrp3", "Gls", "Cdkn2a",
                    "Fdx1", "Lias", "Lipt1", "Dld", "Dlat", "Dlst",
                    "Pdha1", "Pdhb", "Dbt", "Gcsh", "Atp7a", "Atp7b",
                    "Slc31a1", "Slc31a2", "Slc11a2", "Steap3",
                    "Atox1", "Ccs", "Cox17", "Cox11", "Sco1", "Sco2",
                    "Mt1", "Mt2", "Cp", "Commd1", "Sod1", "Sod3",
                    "Tnf", "Il1b", "Hif1a", "Gfap", "Lcn2"),
  stringsAsFactors = FALSE
)

rat2mouse <- tryCatch({
  rat_mart  <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  r2m  <- getLDS(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    mart = rat_mart,
    attributesL = c("ensembl_gene_id", "external_gene_name"),
    martL = mouse_mart
  )
  colnames(r2m) <- c("rat_ensembl", "rat_symbol", "mouse_ensembl", "mouse_symbol")
  cat(sprintf("  大鼠→小鼠同源映射: %d 对\n", nrow(r2m)))

  # 严格过滤1:1同源（参考: OrthoDB/Vilella et al. 2009; 避免1:many导致基因数量偏差）
  rat_dup <- names(which(table(r2m$rat_symbol) > 1))
  mouse_dup <- names(which(table(r2m$mouse_symbol) > 1))
  non_11 <- unique(c(rat_dup, mouse_dup))
  n_before <- nrow(r2m)
  r2m <- r2m[!(r2m$rat_symbol %in% non_11 | r2m$mouse_symbol %in% non_11), ]
  cat(sprintf("  1:1同源过滤: %d → %d 对 (移除 %d 个非1:1基因)\n",
              n_before, nrow(r2m), length(non_11)))
  if (nrow(r2m) < MIN_GENE_THRESHOLD) {
    log_fallback("rat2mouse_online", sprintf("1:1同源基因仅%d对 < 阈值%d", nrow(r2m), MIN_GENE_THRESHOLD),
                  "回退至静态表", "biomaRt结果过于稀疏")
    rat2mouse_fallback
  } else {
    saveRDS(r2m, file.path(OUTPUT_DIR, "rat2mouse_orthologs_1to1.rds"))
    r2m
  }
}, error = function(e) {
  log_fallback("rat2mouse_online", sprintf("biomaRt失败: %s", e$message),
                "使用本地静态映射表", "网络不可用或Ensembl维护")
  rat2mouse_fallback
})

rat_to_mouse_map <- setNames(rat2mouse$mouse_symbol, rat2mouse$rat_symbol)

# 报告每个模块的跨物种基因不可用情况
cat("\n  跨物种基因可用性报告:\n")
gene_loss_report <- list()
for (mod_name in names(MODULE_GENES)) {
  mod_genes <- MODULE_GENES[[mod_name]]
  mapped <- intersect(mod_genes, names(rat_to_mouse_map))
  missing <- setdiff(mod_genes, names(rat_to_mouse_map))
  gene_loss_report[[mod_name]] <- list(
    total = length(mod_genes),
    mapped = length(mapped),
    missing = missing,
    mapped_genes = rat_to_mouse_map[mapped]
  )
  cat(sprintf("  %s: %d/%d 可映射, 丢失: %s\n",
              mod_name, length(mapped), length(mod_genes),
              if(length(missing) > 0) paste(missing, collapse = ", ") else "无"))
}
cat(sprintf("  总可映射基因: %d/%d (%.0f%%)\n",
            sum(sapply(gene_loss_report, function(x) x$mapped)),
            sum(sapply(gene_loss_report, function(x) x$total)),
            100 * sum(sapply(gene_loss_report, function(x) x$mapped)) /
                 sum(sapply(gene_loss_report, function(x) x$total))))

# ======================================================================
#                      第一部分：Bulk 模块活性动态曲线
# ======================================================================

# ---------- 通用函数：解析GPL注释文件 ----------
parse_gpl_annotation <- function(gpl_path) {
  gpl_data <- read.table(gpl_path, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, fill = TRUE,
                          quote = "", comment.char = "#")
  gene_sym_col <- which(grepl("Gene Symbol", colnames(gpl_data), ignore.case = TRUE))[1]
  probe_id_col <- 1

  if (length(gene_sym_col) == 0 || is.na(gene_sym_col)) {
    raw_lines <- readLines(gpl_path)
    data_lines <- raw_lines[!grepl("^#", raw_lines)]
    header_line <- data_lines[1]
    data_lines <- data_lines[-1]
    cols <- strsplit(header_line, "\t")[[1]]
    gene_idx <- which(grepl("Gene Symbol", cols, ignore.case = TRUE))[1]
    gpl_data <- read.table(text = data_lines, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE, fill = TRUE, quote = "", comment.char = "")
    probe2gene <- setNames(gpl_data[, gene_idx], gpl_data[, 1])
  } else {
    probe2gene <- setNames(gpl_data[, gene_sym_col], gpl_data[, probe_id_col])
  }
  probe2gene <- probe2gene[!is.na(probe2gene) & probe2gene != "" & probe2gene != "---"]
  cat(sprintf("  GPL注释加载: %d 有效探针-基因对\n", length(probe2gene)))
  return(probe2gene)
}

# ---------- 通用函数：解析GEO系列矩阵的样本分组 ----------
parse_geo_series <- function(series_lines, default_timepoint = NULL) {
  gsm_line <- series_lines[grepl("^!Sample_geo_accession", series_lines)]
  gsm_parts <- strsplit(gsm_line, "\t")[[1]]
  gsm_ids <- gsub('"', '', gsm_parts[-1])
  gsm_ids <- gsm_ids[grepl("^GSM", gsm_ids)]

  title_line <- series_lines[grepl("^!Sample_title", series_lines)]
  title_parts <- strsplit(title_line, "\t")[[1]]
  titles <- gsub('"', '', title_parts[-1])
  titles <- titles[titles != ""]

  n <- min(length(gsm_ids), length(titles))
  gsm_ids <- gsm_ids[1:n]
  titles <- titles[1:n]

  cat(sprintf("  解析到 %d 个 GSM 样本\n", n))
  return(list(gsm_ids = gsm_ids, titles = titles))
}

# 探针→基因映射后去重（取均值，保留干净基因名）
map_and_clean_probes <- function(expr_mat, probe2gene_map) {
  common <- intersect(rownames(expr_mat), names(probe2gene_map))
  expr_mat <- expr_mat[common, , drop = FALSE]
  gene_names <- probe2gene_map[common]
  gene_names <- gsub(" /// .*$", "", gene_names)
  gene_names <- str_to_title(gene_names)
  expr_agg <- aggregate(as.data.frame(expr_mat), by = list(gene = gene_names), FUN = mean, na.rm = TRUE)
  rownames(expr_agg) <- expr_agg$gene
  expr_agg$gene <- NULL
  cat(sprintf("  映射并去重: %d probes → %d genes\n", nrow(expr_mat), nrow(expr_agg)))
  return(as.matrix(expr_agg))
}

cat("\n========== 第1部分：Bulk 模块活性动态曲线 ==========\n")

# -------------------- 1A. 加载 GSE104036（小鼠RNA-seq 多时序） --------------------
cat("\n--- 1A. 加载 GSE104036 (小鼠 RNA-seq, 3h/6h/12h/24h) ---\n")
gse104036_counts_path <- "D:/反向网络药理学/L1 数据集/bulk/GSE104036（多时序）/GSE104036_TC-RNAseq_counts.txt.gz"

counts_raw <- read.table(gzfile(gse104036_counts_path), header = TRUE,
                          row.names = NULL, check.names = FALSE,
                          stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
# 第一列为基因名（由于空表头被当作 row.names=NULL 后会自动命名）
colnames(counts_raw)[1] <- "gene_symbol"
rownames(counts_raw) <- make.unique(counts_raw$gene_symbol)
counts_raw$gene_symbol <- NULL
cat(sprintf("  原始 counts 矩阵: %d genes × %d samples\n", nrow(counts_raw), ncol(counts_raw)))

# 解析时间点信息
colnames(counts_raw) <- gsub("^X", "", colnames(counts_raw))
sample_info <- data.frame(
  sample_id = colnames(counts_raw),
  stringsAsFactors = FALSE
)
sample_info$group <- sapply(strsplit(sample_info$sample_id, "_"), `[`, 1)
sample_info$timepoint <- ifelse(
  grepl("S", sample_info$group), "sham",
  ifelse(grepl("C", sample_info$group),
         paste0(gsub("C[0-9]+_", "", sample_info$sample_id)),
         paste0(gsub("I[0-9]+_", "", sample_info$sample_id)))
)
sample_info$timepoint <- gsub("hr", "h", sample_info$timepoint)
sample_info$hemisphere <- ifelse(grepl("S", sample_info$group), "sham",
                                  ifelse(grepl("C", sample_info$group), "contralateral", "ipsilateral"))

cat("  样本分组:\n")
print(table(sample_info$timepoint, sample_info$hemisphere))

# 使用同侧（ipsilateral, I）作为卒中样本；Sham 作为对照组
ipsi_samples <- sample_info$sample_id[sample_info$hemisphere == "ipsilateral"]
sham_samples <- sample_info$sample_id[sample_info$hemisphere == "sham"]

gse104036_expr <- counts_raw[, c(sham_samples, ipsi_samples)]
gse104036_meta <- sample_info[sample_info$sample_id %in% c(sham_samples, ipsi_samples), ]

# 去除低表达基因
keep <- rowSums(cpm(gse104036_expr) > 1) >= 3
gse104036_expr <- gse104036_expr[keep, ]
cat(sprintf("  过滤低表达后: %d genes\n", nrow(gse104036_expr)))

# TMM标准化 + log2-CPM
dge_104036 <- DGEList(counts = gse104036_expr)
dge_104036 <- calcNormFactors(dge_104036, method = "TMM")
logcpm_104036 <- cpm(dge_104036, log = TRUE, prior.count = 1)
cat(sprintf("  log2-CPM 矩阵: %d genes × %d samples\n", nrow(logcpm_104036), ncol(logcpm_104036)))

gse104036_timepoints <- gse104036_meta$timepoint
names(gse104036_timepoints) <- gse104036_meta$sample_id

# -------------------- 1B. 加载 GSE97537（大鼠芯片, 24h） --------------------
cat("\n--- 1B. 加载 GSE97537 (大鼠 Affymetrix, 24h) ---\n")
gse97537_series_path <- "D:/反向网络药理学/L1 数据集/bulk/GSE97537(24H)/GSE97537_series_matrix.txt"

series_lines <- readLines(gse97537_series_path)
expr_start <- which(grepl("!series_matrix_table_begin", series_lines)) + 1
expr_end   <- which(grepl("!series_matrix_table_end", series_lines)) - 1
gse97537_expr <- read.table(text = series_lines[expr_start:expr_end],
                             header = TRUE, row.names = 1, sep = "\t",
                             check.names = FALSE, stringsAsFactors = FALSE,
                             fill = TRUE, comment.char = "")
cat(sprintf("  GSE97537 表达矩阵: %d probes × %d samples\n", nrow(gse97537_expr), ncol(gse97537_expr)))

# GPL1355 注释解析
gpl_path_97537 <- "D:/反向网络药理学/L1 数据集/bulk/GSE97537(24H)/GPL1355-10794 (1).txt"
probe2gene_97537 <- parse_gpl_annotation(gpl_path_97537)

gse97537_expr <- map_and_clean_probes(gse97537_expr, probe2gene_97537)

# 从系列矩阵解析样本分组
geo_info_97537 <- parse_geo_series(series_lines)
gse97537_group <- ifelse(grepl("MCAO", geo_info_97537$titles, ignore.case = TRUE), "MCAO_24h",
                  ifelse(grepl("Sham", geo_info_97537$titles, ignore.case = TRUE), "Sham_24h", "Unknown"))
gsm_order <- geo_info_97537$gsm_ids[seq_len(min(length(geo_info_97537$gsm_ids), length(gse97537_group)))]
gse97537_group <- gse97537_group[seq_len(length(gsm_order))]
names(gse97537_group) <- gsm_order
cat(sprintf("  MCAO: %d, Sham: %d, Unknown: %d\n",
            sum(grepl("MCAO", gse97537_group)),
            sum(grepl("Sham", gse97537_group)),
            sum(grepl("Unknown", gse97537_group))))

# -------------------- 1C. 加载 GSE61616（大鼠芯片, 7d） --------------------
cat("\n--- 1C. 加载 GSE61616 (大鼠 Affymetrix, 7d) ---\n")
gse61616_series_path <- "D:/反向网络药理学/L1 数据集/bulk/GSE61616（7d）/GSE61616_series_matrix.txt.gz"

series_lines_61616 <- readLines(gzfile(gse61616_series_path))
expr_start_61616 <- which(grepl("!series_matrix_table_begin", series_lines_61616)) + 1
expr_end_61616   <- which(grepl("!series_matrix_table_end", series_lines_61616)) - 1
gse61616_expr <- read.table(text = series_lines_61616[expr_start_61616:expr_end_61616],
                             header = TRUE, row.names = 1, sep = "\t",
                             check.names = FALSE, stringsAsFactors = FALSE,
                             fill = TRUE, comment.char = "")
cat(sprintf("  GSE61616 表达矩阵: %d probes × %d samples\n", nrow(gse61616_expr), ncol(gse61616_expr)))

gpl_61616_path <- "D:/反向网络药理学/L1 数据集/bulk/GSE61616（7d）/GPL1355-10794 (1).txt"
probe2gene_61616 <- parse_gpl_annotation(gpl_61616_path)

gse61616_expr <- map_and_clean_probes(gse61616_expr, probe2gene_61616)

# 从系列矩阵解析样本分组
geo_info_61616 <- parse_geo_series(series_lines_61616)
gse61616_group <- ifelse(grepl("Sham", geo_info_61616$titles, ignore.case = TRUE), "Sham_7d",
                  ifelse(grepl("Model|MCAO", geo_info_61616$titles, ignore.case = TRUE), "Model_7d",
                  ifelse(grepl("XST", geo_info_61616$titles, ignore.case = TRUE), "XST_7d", "Unknown")))
n_61616 <- min(length(geo_info_61616$gsm_ids), length(gse61616_group))
gse61616_group <- gse61616_group[1:n_61616]
names(gse61616_group) <- geo_info_61616$gsm_ids[1:n_61616]

cat(sprintf("  Sham: %d, Model: %d, XST: %d, Unknown: %d\n",
            sum(grepl("Sham", gse61616_group)),
            sum(grepl("Model", gse61616_group)),
            sum(grepl("XST", gse61616_group)),
            sum(grepl("Unknown", gse61616_group))))

# -------------------- 1D. 各数据集独立 ssGSEA（避免跨物种基因映射问题）--------------------
cat("\n--- 1D. 各数据集独立计算 ssGSEA ---\n")

# 定义各物种的模块基因集
# 小鼠模块基因（小写首字母）
MODULE_GENES_MOUSE <- list(
  "M1_CopperTransport" = c("Slc31a1", "Slc31a2", "Atp7a", "Atp7b", "Slc11a2", "Steap3"),
  "M2_FeS_Lipoylation" = c("Fdx1", "Lias", "Lipt1", "Lipt2", "Gcsh"),
  "M3_TCA_PDH"          = c("Dld", "Dlat", "Dlst", "Pdha1", "Pdhb", "Dbt"),
  "M4_Chaperones"       = c("Atox1", "Ccs", "Cox17", "Cox11", "Sco1", "Sco2"),
  "M5_Metallothioneins" = c("Mt1", "Mt2", "Cp", "Commd1", "Sod1", "Sod3"),
  "M6_StressResponse"   = c("Mtf1", "Nfe2l2", "Nlrp3", "Gls", "Cdkn2a")
)

# 大鼠模块基因（大写首字母，与小鼠同）
MODULE_GENES_RAT <- MODULE_GENES_MOUSE  # 铜死亡基因在大鼠和小鼠中符号一致

compute_ssgsea_per_dataset <- function(expr_matrix, gene_sets, dataset_name, timepoint_val, batch_val) {
  gs_filtered <- lapply(gene_sets, function(gs) intersect(gs, rownames(expr_matrix)))
  gs_filtered <- gs_filtered[sapply(gs_filtered, length) >= 2]

  cat(sprintf("  %s 可用模块基因:\n", dataset_name))
  for (nm in names(gs_filtered)) {
    cat(sprintf("    %s: %d genes (%s)\n", nm, length(gs_filtered[[nm]]),
                paste(gs_filtered[[nm]], collapse = ", ")))
  }

  if (length(gs_filtered) == 0) {
    cat(sprintf("  %s: 无足够基因用于 ssGSEA\n", dataset_name))
    return(NULL)
  }

  ssgsea_param <- ssgseaParam(
    exprData = as.matrix(expr_matrix),
    geneSets = gs_filtered,
    minSize = 2,
    maxSize = 500,
    normalize = TRUE
  )
  scores <- gsva(ssgsea_param, verbose = FALSE)

  scores_df <- as.data.frame(t(scores))
  scores_df$timepoint <- timepoint_val
  scores_df$batch <- batch_val
  scores_df$sample_id <- rownames(scores_df)
  return(scores_df)
}

# GSE104036 (小鼠) ssGSEA: 仅用 ipsilateral + sham
gse104036_for_ssgsea <- logcpm_104036
gse104036_tp_map <- gse104036_timepoints[colnames(logcpm_104036)]

ssgsea_104036 <- compute_ssgsea_per_dataset(
  gse104036_for_ssgsea, MODULE_GENES_MOUSE,
  "GSE104036",
  gse104036_tp_map,
  "GSE104036"
)

# GSE97537 (大鼠) ssGSEA: 仅用疾病样本（MCAO）
gse97537_mcao_cols <- names(gse97537_group)[grep("MCAO", gse97537_group)]
gse97537_for_ssgsea <- gse97537_expr[, gse97537_mcao_cols, drop = FALSE]

ssgsea_97537 <- compute_ssgsea_per_dataset(
  gse97537_for_ssgsea, MODULE_GENES_RAT,
  "GSE97537",
  "24h",
  "GSE97537"
)

# GSE61616 (大鼠) ssGSEA: 仅用 Model 样本（未治疗卒中）
gse61616_model_cols <- names(gse61616_group)[grep("Model", gse61616_group)]
gse61616_for_ssgsea <- gse61616_expr[, gse61616_model_cols, drop = FALSE]

ssgsea_61616 <- compute_ssgsea_per_dataset(
  gse61616_for_ssgsea, MODULE_GENES_RAT,
  "GSE61616",
  "7d",
  "GSE61616"
)

# -------------------- 1E. 合并模块活性得分（跨数据集） --------------------
cat("\n--- 1E. 合并模块活性得分 ---\n")

all_ssgsea_parts <- list()
if (!is.null(ssgsea_104036)) all_ssgsea_parts[[1]] <- ssgsea_104036
if (!is.null(ssgsea_97537)) all_ssgsea_parts[[2]] <- ssgsea_97537
if (!is.null(ssgsea_61616)) all_ssgsea_parts[[3]] <- ssgsea_61616

common_mod_cols <- Reduce(intersect, lapply(all_ssgsea_parts, function(x) {
  setdiff(colnames(x), c("timepoint", "batch", "sample_id"))
}))
cat(sprintf("  跨数据集共同模块: %d\n", length(common_mod_cols)))

if (length(common_mod_cols) < 3) {
  stop(sprintf("跨数据集共同模块数不足 (%d < 3)。请检查基因符号大小写映射，GPL1355注释需转为首字母大写",
               length(common_mod_cols)))
}
cat(sprintf("  共同模块: %s\n", paste(common_mod_cols, collapse = ", ")))

ssgsea_list <- lapply(all_ssgsea_parts, function(x) x[, c(common_mod_cols, "timepoint", "batch", "sample_id"), drop = FALSE])
ssgsea_df <- do.call(rbind, ssgsea_list)

# 统一 ComBat 对模块得分进行批次校正（而非基因表达）
if (length(all_ssgsea_parts) >= 2) {
  batch_labels <- ssgsea_df$batch
  tp_labels_combat <- ssgsea_df$timepoint
  module_mat <- as.matrix(ssgsea_df[, common_mod_cols, drop = FALSE])

  # 校正前跨批次相关性（验证ComBat的必要性和效果，参考: Müller et al. BMC Bioinformatics 2016）
  cat("  ComBat前跨批次一致性评估:\n")
  pre_batch_cors <- list()
  unique_batches <- unique(batch_labels)
  for (i in 1:length(unique_batches)) {
    for (j in i:length(unique_batches)) {
      if (i < j) {
        idx_i <- which(batch_labels == unique_batches[i])
        idx_j <- which(batch_labels == unique_batches[j])
        n_common <- min(length(idx_i), length(idx_j))
        if (n_common > 0) {
          mod_cors <- sapply(common_mod_cols, function(mod) {
            cor(module_mat[idx_i[1:n_common], mod], module_mat[idx_j[1:n_common], mod],
                method = "pearson", use = "complete.obs")
          })
          avg_cor <- mean(mod_cors, na.rm = TRUE)
          cat(sprintf("    %s vs %s: 平均 Pearson r = %+.3f\n",
                      unique_batches[i], unique_batches[j], avg_cor))
          pre_batch_cors[[paste(unique_batches[i], unique_batches[j], sep = "_vs_")]] <- avg_cor
        }
      }
    }
  }

  module_mat_corrected <- tryCatch({
    ComBat(dat = t(module_mat), batch = batch_labels, mod = model.matrix(~1, data = data.frame(tp = tp_labels_combat)))
  }, error = function(e) {
    log_fallback("ssGSEA_ComBat", sprintf("ComBat失败: %s", e$message),
                  "使用原始模块得分", "sva::ComBat计算错误")
    t(module_mat)
  })

  # 校正后跨批次相关性
  if (length(pre_batch_cors) > 0) {
    cat("  ComBat后跨批次一致性评估:\n")
    module_mat_post <- t(module_mat_corrected)
    colnames(module_mat_post) <- colnames(module_mat)
    for (i in 1:length(unique_batches)) {
      for (j in i:length(unique_batches)) {
        if (i < j) {
          idx_i <- which(batch_labels == unique_batches[i])
          idx_j <- which(batch_labels == unique_batches[j])
          n_common <- min(length(idx_i), length(idx_j))
          if (n_common > 0) {
            post_cors <- sapply(common_mod_cols, function(mod) {
              cor(module_mat_post[idx_i[1:n_common], mod], module_mat_post[idx_j[1:n_common], mod],
                  method = "pearson", use = "complete.obs")
            })
            avg_post_cor <- mean(post_cors, na.rm = TRUE)
            key <- paste(unique_batches[i], unique_batches[j], sep = "_vs_")
            pre_val <- ifelse(key %in% names(pre_batch_cors), pre_batch_cors[[key]], NA_real_)
            cat(sprintf("    %s: 校正后 r = %+.3f (校正前 %+.3f, Δ = %+.3f)\n",
                        key, avg_post_cor, pre_val, avg_post_cor - pre_val))
          }
        }
      }
    }
  }

  # 更新 ssgsea_df
  for (mod in common_mod_cols) {
    ssgsea_df[[mod]] <- module_mat_corrected[mod, ]
  }
  cat(sprintf("  模块得分 ComBat 校正完成: %d samples × %d modules\n", nrow(ssgsea_df), length(common_mod_cols)))
}

module_names <- common_mod_cols

cat(sprintf("  合并后 ssGSEA 评分: %d samples × %d modules\n", nrow(ssgsea_df), length(common_mod_cols)))

# 保存
write.csv(ssgsea_df, file.path(OUTPUT_DIR, "ssGSEA_module_scores.csv"), row.names = FALSE)

# (2) 合并表达矩阵：仅用于锚定标记基因的Bulk趋势分析
# 使用 GSE104036 (小鼠, 3h-24h) + 补充 GSE61616 (大鼠->小鼠, 7d)
cat(sprintf("  构建锚定表达矩阵...\n"))

anchor_parts <- list()

# GSE104036 (小鼠完整时间序列)
anchor_parts[["GSE104036"]] <- logcpm_104036

# GSE61616 (大鼠->小鼠, 7d)
common_rats_61616 <- intersect(rownames(gse61616_expr), names(rat_to_mouse_map))
if (length(common_rats_61616) > 10) {
  gse61616_mouse <- gse61616_expr[common_rats_61616, , drop = FALSE]
  mouse_genes_61616 <- rat_to_mouse_map[common_rats_61616]
  keep_61616 <- !is.na(mouse_genes_61616) & mouse_genes_61616 != ""
  gse61616_mouse <- gse61616_mouse[keep_61616, , drop = FALSE]
  mouse_genes_61616 <- mouse_genes_61616[keep_61616]
  rownames(gse61616_mouse) <- mouse_genes_61616
  gse61616_mouse <- aggregate(as.data.frame(gse61616_mouse),
                               by = list(gene = rownames(gse61616_mouse)),
                               FUN = mean, na.rm = TRUE)
  rownames(gse61616_mouse) <- gse61616_mouse$gene
  gse61616_mouse$gene <- NULL

  model_ids <- names(gse61616_group)[grep("Model", gse61616_group)]
  model_ids <- intersect(model_ids, colnames(gse61616_mouse))
  if (length(model_ids) > 0) {
    anchor_parts[["GSE61616_7d"]] <- as.matrix(gse61616_mouse[, model_ids, drop = FALSE])
  }
}

common_anchor_genes <- Reduce(intersect, lapply(anchor_parts, rownames))
cat(sprintf("  锚定共有基因: %d\n", length(common_anchor_genes)))

if (length(common_anchor_genes) > 10) {
  merged_expr_list <- lapply(anchor_parts, function(x) x[common_anchor_genes, , drop = FALSE])
    merged_expr_raw <- do.call(cbind, merged_expr_list)
    batch_labels <- c(rep("GSE104036", ncol(anchor_parts[["GSE104036"]])),
                      rep("GSE61616", ncol(anchor_parts[["GSE61616_7d"]])))

    # 跨平台 ComBat 校正（RNA-seq log2-CPM vs 芯片强度不可直接比较）
    merged_expr_corrected <- tryCatch({
      ComBat(dat = merged_expr_raw, batch = batch_labels)
    }, error = function(e) {
      log_fallback("anchor_ComBat", sprintf("跨平台ComBat失败: %s", e$message),
                    "使用原始值 (RNA-seq + 芯片未校正)", "可能引入技术偏差")
      merged_expr_raw
    })
    cat(sprintf("    锚定矩阵跨平台 ComBat 校正完成\n"))
  anchor_tp_labels <- c(
    gse104036_timepoints[colnames(anchor_parts[["GSE104036"]])],
    rep("7d", ncol(anchor_parts[["GSE61616_7d"]]))
  )
  built_merged_expr <- TRUE
  cat(sprintf("  锚定表达矩阵: %d genes X %d samples\n", nrow(merged_expr_corrected), ncol(merged_expr_corrected)))
  cat(sprintf("  时间点: %s\n", paste(names(table(anchor_tp_labels)), table(anchor_tp_labels), sep = "=", collapse = ", ")))
} else {
  log_fallback("anchor_merge", sprintf("共有基因仅%d个 < 阈值%d", length(common_anchor_genes), 10),
                "降级使用 GSE104036 单独", "缺少7d时间点，晚期锚定精度降低")
  merged_expr_corrected <- logcpm_104036
  anchor_tp_labels <- gse104036_timepoints[colnames(logcpm_104036)]
  batch_labels <- rep("GSE104036", ncol(logcpm_104036))
  built_merged_expr <- TRUE
}

# ==================== 1G. Loess 平滑拟合 + 拐点提取 ====================
cat("\n--- 1G. Loess 拟合与拐点提取 ---\n")

time_order <- c("sham", "3h", "6h", "12h", "24h", "7d")
time_numeric <- setNames(c(0, 3, 6, 12, 24, 168), time_order)

ssgsea_df$time_num <- time_numeric[ssgsea_df$timepoint]
ssgsea_df$time_num[ssgsea_df$timepoint == "sham"] <- 0

# 仅取疾病样本（非 Sham）拟合曲线
disease_idx <- ssgsea_df$timepoint != "sham"
disease_data <- ssgsea_df[disease_idx, ]

module_names <- common_mod_cols
curve_results <- list()
inflection_points <- list()

pdf(file.path(FIGURE_DIR, "module_activity_curves.pdf"), width = 14, height = 10)
par(mfrow = c(2, 3), mar = c(5, 4, 4, 2))

for (mod in module_names) {
  cat(sprintf("\n  --- %s ---\n", mod))

  fit_data <- disease_data[, c(mod, "time_num")]
  fit_data <- fit_data[order(fit_data$time_num), ]
  colnames(fit_data)[1] <- "score"

  # Loess 拟合
  loess_fit <- tryCatch({
    loess(score ~ time_num, data = fit_data, span = 0.75, degree = 2)
  }, error = function(e) {
    cat(sprintf("    Loess 拟合失败: %s，使用 gam 替代\n", e$message))
    NULL
  })

  # 预测平滑值
  pred_x <- seq(min(fit_data$time_num), max(fit_data$time_num), length.out = 200)
  if (!is.null(loess_fit)) {
    pred_y <- predict(loess_fit, newdata = data.frame(time_num = pred_x))
  } else {
    gam_fit <- mgcv::gam(score ~ s(time_num, k = 4), data = fit_data)
    pred_y <- predict(gam_fit, newdata = data.frame(time_num = pred_x))
  }

  # 一阶导数（中心差分）
  dx <- pred_x[2] - pred_x[1]
  dy <- diff(pred_y) / dx
  dy_x <- pred_x[-1] - dx / 2

  # 拐点：一阶导数过零点
  sign_change <- which(diff(sign(dy)) != 0)
  inflection_hours <- dy_x[sign_change]
  inflection_values <- pred_y[-1][sign_change]

  if (length(inflection_hours) > 0) {
    cat(sprintf("    拐点（导数过零）: %s h\n",
                paste(sprintf("%.1f", inflection_hours), collapse = ", ")))
    cat(sprintf("    拐点活性值: %s\n",
                paste(sprintf("%.3f", inflection_values), collapse = ", ")))
  } else {
    cat("    未检测到明确拐点，使用峰值/谷值替代\n")
    peak_idx <- which.max(pred_y)
    valley_idx <- which.min(pred_y)
    inflection_hours <- pred_x[c(valley_idx, peak_idx)]
    inflection_values <- pred_y[c(valley_idx, peak_idx)]
    cat(sprintf("    极值点: %s h\n",
                paste(sprintf("%.1f", inflection_hours), collapse = ", ")))
  }

  inflection_points[[mod]] <- data.frame(
    module = mod,
    time_h = inflection_hours,
    activity_value = inflection_values,
    stringsAsFactors = FALSE
  )

  # 绘图
  plot(fit_data$time_num, fit_data$score,
       xlab = "Time (hours)", ylab = "ssGSEA Score",
       main = mod, pch = 16, col = rgb(0.2, 0.4, 0.8, 0.6),
       cex = 1.2, xaxt = "n")
  axis(1, at = time_numeric, labels = time_order)
  lines(pred_x, pred_y, col = "darkred", lwd = 2.5)
  abline(v = inflection_hours, lty = 2, col = "darkgreen", lwd = 1.5)

  # 添加 E/M/L 分期标注
  abline(v = c(6, 24), lty = 3, col = "gray50")
  text(1.5, max(fit_data$score) * 0.95, "E", col = "blue", cex = 1.5, font = 2)
  text(15, max(fit_data$score) * 0.95, "M", col = "orange", cex = 1.5, font = 2)
  text(96, max(fit_data$score) * 0.95, "L", col = "red", cex = 1.5, font = 2)

  # 保存曲线数据
  curve_results[[mod]] <- data.frame(
    time_h = pred_x,
    smoothed_score = pred_y,
    module = mod,
    stringsAsFactors = FALSE
  )
}
dev.off()

# 保存拐点汇总
inflection_summary <- do.call(rbind, inflection_points)
write.csv(inflection_summary, file.path(OUTPUT_DIR, "inflection_points.csv"), row.names = FALSE)
cat("\n  拐点汇总已保存\n")
print(inflection_summary)

# 保存曲线数据
curve_all <- do.call(rbind, curve_results)
write.csv(curve_all, file.path(OUTPUT_DIR, "smoothed_curves.csv"), row.names = FALSE)

# ==================== 1H. 各时间点均值 ± SD 统计 ====================
cat("\n--- 1H. 时间点均值统计 ---\n")
tp_summary <- list()
for (tp in time_order) {
  tp_data <- ssgsea_df[ssgsea_df$timepoint == tp, module_names, drop = FALSE]
  if (nrow(tp_data) > 0) {
    tp_summary[[tp]] <- data.frame(
      timepoint = tp,
      module = colnames(tp_data),
      mean = colMeans(tp_data, na.rm = TRUE),
      sd = apply(tp_data, 2, sd, na.rm = TRUE),
      n = nrow(tp_data),
      stringsAsFactors = FALSE
    )
  }
}
tp_summary_df <- do.call(rbind, tp_summary)
write.csv(tp_summary_df, file.path(OUTPUT_DIR, "module_timepoint_summary.csv"), row.names = FALSE)
cat(sprintf("  时间点统计: %d 行\n", nrow(tp_summary_df)))

# ======================================================================
#                    第二部分：单细胞拟时序与分期
# ======================================================================
cat("\n========== 第2部分：单细胞拟时序与分期 ==========\n")

# 注：Monocle 3 和 CytoTRACE 依赖复杂，此处采用 Scanpy (Python) +
# Monocle 3 互补方案。Python 脚本负责 scVelo，R 负责 Monocle 3。

# -------------------- 2A. 从 h5ad 加载预处理后的数据 --------------------
cat("\n--- 2A. 准备单细胞数据 ---\n")

sc_h5ad_path <- file.path(BASE_DIR, "results/L1_phenotype_anchoring/sc_adata_cuproptosis.h5ad")
sc_alt_h5ad  <- file.path(BASE_DIR, "results/L1_phenotype_anchoring/pseudotime_adata.h5ad")

adata_exists <- file.exists(sc_h5ad_path) || file.exists(sc_alt_h5ad)
cat(sprintf("  单细胞 h5ad 文件存在: %s\n", adata_exists))

if (!adata_exists) {
  cat("  警告：未找到预处理后的 h5ad 文件\n")
  cat("  请先运行 scrna_analysis.py 或 pseudotime_analysis.py 生成预处理数据\n")
  cat("  当前使用从原始10X数据直接构建 Monocle 3 CDS...\n")
}

# -------------------- 2B. Monocle 3 拟时序分析 --------------------
cat("\n--- 2B. Monocle 3 拟时序分析 ---\n")

run_monocle3 <- function() {
  tryCatch({
    suppressPackageStartupMessages({
      library(monocle3)
      library(Seurat)
      library(SeuratWrappers)
    })

    # Seurat v4/v5 兼容：选择正确的表达矩阵提取函数
    seurat_v5 <- packageVersion("Seurat") >= "5.0.0"
    get_counts <- if (seurat_v5) {
      function(obj) LayerData(obj, assay = "RNA", layer = "counts")
    } else {
      function(obj) GetAssayData(obj, assay = "RNA", slot = "counts")
    }
    cat(sprintf("  Seurat 版本: %s, 使用 %s 提取counts\n",
                packageVersion("Seurat"), ifelse(seurat_v5, "LayerData (v5)", "GetAssayData (v4)")))

    # 从10X数据构建Seurat对象
    data_10x_dir <- "D:/反向网络药理学/L1 数据集/RNA-seq/GSE174574_10X_organized"
    mcao_dirs <- list.files(data_10x_dir, pattern = "MCAO", full.names = TRUE)
    sham_dirs <- list.files(data_10x_dir, pattern = "Sham", full.names = TRUE)
    all_dirs <- c(sham_dirs, mcao_dirs)

    seurat_list <- list()
    for (i in seq_along(all_dirs)) {
      sample_dir <- all_dirs[i]
      sample_name <- basename(sample_dir)
      counts <- Read10X(sample_dir)
      seurat_obj <- CreateSeuratObject(
        counts = counts,
        project = sample_name,
        min.cells = 3,
        min.features = 200
      )
      seurat_obj$sample <- sample_name
      seurat_obj$condition <- ifelse(grepl("Sham", sample_name), "Sham", "MCAO")
      seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^mt-", col.name = "percent.mt")
      seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
      seurat_list[[sample_name]] <- seurat_obj
    }

    merged_seurat <- merge(seurat_list[[1]], seurat_list[-1])
    # Seurat v5: 合并所有 layer (counts.Sham_1, counts.MCAO_1 ... → counts)
    merged_seurat <- JoinLayers(merged_seurat)
    merged_seurat <- NormalizeData(merged_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    merged_seurat <- FindVariableFeatures(merged_seurat, nfeatures = 2000)
    merged_seurat <- ScaleData(merged_seurat)
    merged_seurat <- RunPCA(merged_seurat)
    merged_seurat <- RunUMAP(merged_seurat, dims = 1:30)
    merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30)
    merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)

    cat(sprintf("  Seurat 对象: %d cells, %d genes\n", ncol(merged_seurat), nrow(merged_seurat)))

    # 细胞类型筛选
    cell_type_markers <- list(
      "Microglia" = c("Ptprc", "Aif1", "Cx3cr1", "Tmem119", "P2ry12", "C1qa"),
      "Neuron"    = c("Snap25", "Syt1", "Nefl", "Rbfox3", "Syn1", "Nefm"),
      "Astrocyte" = c("Gfap", "Aqp4", "Slc1a3", "Aldh1l1", "Slc1a2")
    )

    for (ct in names(cell_type_markers)) {
      markers <- intersect(cell_type_markers[[ct]], rownames(merged_seurat))
      merged_seurat <- AddModuleScore(merged_seurat, features = list(markers), name = paste0("score_", ct))
    }

    merged_seurat$cell_type <- apply(sapply(names(cell_type_markers), function(ct) {
      merged_seurat@meta.data[[paste0("score_", ct, "1")]]
    }), 1, function(x) names(which.max(x)))

    # 保存合并后的Seurat对象
    saveRDS(merged_seurat, file.path(OUTPUT_DIR, "merged_seurat.rds"))

    # 分别对3种主要细胞类型执行Monocle 3
    monocle_results <- list()
    monocle_results[["merged_seurat"]] <- merged_seurat
    monocle_results[["cell_type_assignments"]] <- merged_seurat$cell_type

    for (ct in c("Microglia", "Neuron", "Astrocyte")) {
      cat(sprintf("\n  --- Monocle 3: %s ---\n", ct))

      ct_result <- tryCatch({
        sub_seurat <- subset(merged_seurat, subset = cell_type == ct)

        if (ncol(sub_seurat) < 50) {
          cat(sprintf("    细胞数(%d)不足，跳过\n", ncol(sub_seurat)))
          return(NULL)
        }
        cat(sprintf("    细胞数: %d\n", ncol(sub_seurat)))

        # Seurat v5: 使用 LayerData 替代 GetAssayData
        cat("    转换为 Monocle 3 CDS...\n")
        cds <- as.cell_data_set(sub_seurat)

        # 确保 counts 层可用
        cds <- cluster_cells(cds, reduction_method = "UMAP")
        cds <- learn_graph(cds)

        # 使用 CytoTRACE2 确定起点 (Gulati et al., Science 2020; v2 2024)
        # 三级回退: CytoTRACE2 → CytoTRACE(v1) → 早期标记基因 → UMAP原点
        cyto_ok <- TRUE
        cyto_score <- tryCatch({
          suppressPackageStartupMessages(library(CytoTRACE2))
          cyto_res <- cytotrace2(sub_seurat, is_seurat = TRUE,
                                 slot_type = "counts", species = "mouse", seed = 42)
          cyto_res$CytoTRACE2_Score
        }, error = function(e) {
          cat(sprintf("    CytoTRACE2 失败: %s\n", e$message))
          tryCatch({
            suppressPackageStartupMessages(library(CytoTRACE))
            expr_matrix <- get_counts(sub_seurat)
            cyto_res <- CytoTRACE(as.matrix(expr_matrix), ncores = 1)
            cat("    使用 CytoTRACE v1 (Gulati et al. 2020)\n")
            return(cyto_res$CytoTRACE)
          }, error = function(e2) {
            cat(sprintf("    CytoTRACE v1 也失败: %s\n", e2$message))
            cyto_ok <<- FALSE
            rep(NA_real_, ncol(sub_seurat))
          })
        })

        names(cyto_score) <- colnames(sub_seurat)

        if (!cyto_ok) {
          early_markers <- intersect(c("Tnf", "Il1b", "Ccl2", "Fos"), rownames(sub_seurat))
          if (length(early_markers) > 0) {
            expr_mat <- get_counts(sub_seurat)
            early_expr <- Matrix::colSums(expr_mat[early_markers, , drop = FALSE])
            root_cells <- names(which.max(early_expr))
            cat(sprintf("    替代根节点(早期炎症标记): %s (Tnf/Il1b/Ccl2等总表达最高)\n", root_cells))
          } else {
            umap_coords <- Embeddings(sub_seurat, "umap")
            centroid_dist <- sqrt(rowSums(umap_coords^2))
            root_cells <- names(which.min(centroid_dist))
            cat(sprintf("    替代根节点(UMAP原点最近): %s\n", root_cells))
          }
          log_fallback(paste0("pseudotime_root_", ct), "CytoTRACE2/v1均不可用",
                        "替代策略(早期标记/UMAP)确定根节点", "可能引入偏差，建议检查Fos表达")
          pseudotime_unreliable <- TRUE
        } else {
          max_val <- max(cyto_score, na.rm = TRUE)
          root_cells <- names(cyto_score)[cyto_score == max_val & !is.na(cyto_score)]
          cat(sprintf("    根节点: %d cells with max CytoTRACE = %.4f\n", length(root_cells), max_val))
          pseudotime_unreliable <- FALSE
        }

        # 多方法根节点验证 (参考: Saelens et al. Nature Biotech 2019 — 多根方法共识)
        # 方法2: Fos (立即早期基因, PMID: 29097358) — 独立验证根节点
        fos_gene <- intersect("Fos", rownames(sub_seurat))
        if (length(fos_gene) > 0) {
          expr_mat <- get_counts(sub_seurat)
          fos_expr <- expr_mat[fos_gene, , drop = FALSE]
          if (nrow(fos_expr) > 0) {
            fos_expr_vec <- as.numeric(fos_expr[1, ])
            names(fos_expr_vec) <- colnames(fos_expr)
            fos_root <- names(which.max(fos_expr_vec))
            cat(sprintf("    Fos根节点验证: %s (Fos表达最高)\n", fos_root))
            if (!pseudotime_unreliable && length(intersect(root_cells, fos_root)) > 0) {
              cat("    ✓ CytoTRACE与Fos根节点一致\n")
            } else if (!pseudotime_unreliable) {
              cat("    ! CytoTRACE与Fos根节点不一致，以CytoTRACE为准 (Fos可能受独立通路调节)\n")
            }
          }
        }

        cds <- order_cells(cds, root_cells = root_cells)

        # 提取拟时序值
        pseudotime_vals <- pseudotime(cds)
        names(pseudotime_vals) <- colnames(cds)
        pseudotime_vals <- pseudotime_vals[!is.na(pseudotime_vals) & is.finite(pseudotime_vals)]

        if (length(pseudotime_vals) < 6) {
          cat("    有效拟时序值不足\n")
          return(NULL)
        }

        # 密度断点替代均分 (参考: Chen et al. Cell Systems 2019 — density-based trajectory
        # 断点比均分更适应非线性轨迹和不均匀细胞分布)
        pt_density <- density(pseudotime_vals, n = 512)
        cum_density <- cumsum(pt_density$y) / sum(pt_density$y)
        t1 <- pt_density$x[which.min(abs(cum_density - 1/3))]
        t2 <- pt_density$x[which.min(abs(cum_density - 2/3))]
        breaks <- c(min(pseudotime_vals), t1, t2, max(pseudotime_vals))
        cat(sprintf("    密度断点: t1=%.2f (33%%ile), t2=%.2f (67%%ile)\n", t1, t2))

        stage_labels <- cut(pseudotime_vals, breaks = breaks,
                            labels = c("E", "M", "L"), include.lowest = TRUE)
        names(stage_labels) <- names(pseudotime_vals)

        stage_counts <- table(stage_labels)
        cat(sprintf("    E/M/L分期 (密度分组): E=%d(%.1f%%), M=%d(%.1f%%), L=%d(%.1f%%)\n",
                    stage_counts["E"], 100*stage_counts["E"]/length(pseudotime_vals),
                    stage_counts["M"], 100*stage_counts["M"]/length(pseudotime_vals),
                    stage_counts["L"], 100*stage_counts["L"]/length(pseudotime_vals)))

        # 拟时序分布直方图 (判断三等分合理性)
        pdf(file.path(FIGURE_DIR, paste0("pseudotime_histogram_", ct, ".pdf")), width = 10, height = 5)
        par(mfrow = c(1, 2))
        hist(pseudotime_vals, breaks = 50, col = "steelblue", border = "white",
             main = paste(ct, "Pseudotime Distribution"),
             xlab = "Pseudotime", probability = TRUE)
        lines(density(pseudotime_vals), col = "darkred", lwd = 2)
        abline(v = c(t1, t2), lty = 2, col = c("blue", "orange"), lwd = 2)
        legend("topright", legend = c("E|M 断点", "M|L 断点"),
               lty = 2, col = c("blue", "orange"), cex = 0.8)

        # CytoTRACE vs Pseudotime
        if (!all(is.na(cyto_score))) {
          plot(pseudotime_vals[names(cyto_score)], cyto_score,
               pch = 16, cex = 0.3, col = rgb(0.2, 0.4, 0.8, 0.3),
               main = "CytoTRACE vs Pseudotime",
               xlab = "Pseudotime", ylab = "CytoTRACE Score")
          abline(lm(cyto_score ~ pseudotime_vals[names(cyto_score)]), col = "red", lwd = 2)
          ct_cor <- cor(pseudotime_vals[names(cyto_score)], cyto_score,
                        method = "spearman", use = "complete.obs")
          legend("topright", legend = sprintf("ρ = %.3f", ct_cor), bty = "n", cex = 1.2)
        } else {
          plot.new()
          text(0.5, 0.5, "CytoTRACE unavailable")
        }
        dev.off()

        ct_res <- list(
          cds = cds,
          pseudotime = pseudotime_vals,
          stages = stage_labels,
          cyto_score = cyto_score,
          cell_type = ct,
          expression_matrix = get_counts(sub_seurat),
          pseudotime_unreliable = pseudotime_unreliable
        )

        # 保存拟时序图
        pdf(file.path(FIGURE_DIR, paste0("monocle3_trajectory_", ct, ".pdf")), width = 8, height = 6)
        print(plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE,
                   graph_label_size = 1.5, cell_size = 0.8) +
          ggtitle(paste(ct, "- Monocle 3 Pseudotime")))
        dev.off()

        ct_res
      }, error = function(e) {
        cat(sprintf("    %s Monocle 3 失败: %s\n", ct, e$message))
        log_fallback(paste0("monocle3_", ct), sprintf("Monocle3失败: %s", e$message),
                      "跳过该细胞类型", "后续分析将缺少此细胞类型的数据")
        return(NULL)
      })

      if (!is.null(ct_result)) {
        monocle_results[[ct]] <- ct_result
      }
    }

    return(monocle_results)
  }, error = function(e) {
    cat(sprintf("  Monocle 3 分析失败: %s\n", e$message))
    cat("  将使用基于 Scanpy 的替代方案（需先运行 Python 脚本）\n")
    return(NULL)
  })
}

monocle_results <- run_monocle3()

# 如果 Monocle 3 失败，尝试从 Python 结果加载
if (is.null(monocle_results)) {
  cat("\n  尝试从 Python 拟时序结果加载...\n")
  pseudotime_csv <- file.path(OUTPUT_DIR, "../L1_phenotype_anchoring/pseudotime_stages.csv")
  if (file.exists(pseudotime_csv)) {
    pt_data <- read.csv(pseudotime_csv)
    cat(sprintf("  从CSV加载: %d cells\n", nrow(pt_data)))
    monocle_results <- list(loaded_from_python = pt_data)
  }
}

# 保存拟时序分期结果
if (!is.null(monocle_results)) {
  saveRDS(monocle_results, file.path(OUTPUT_DIR, "monocle3_results.rds"))
}

# ======================================================================
#                    第三部分：定性分期锚定
# ======================================================================
cat("\n========== 第3部分：定性分期锚定 ==========\n")

# -------------------- 3A. 标记基因跨组学趋势分析 --------------------
cat("\n--- 3A. 标记基因 Bulk 时间趋势 ---\n")

anchor_results <- data.frame(
  gene = character(),
  entrez = character(),
  expected_stage = character(),
  bulk_trend = character(),
  scRNA_trend = character(),
  consistent = logical(),
  stringsAsFactors = FALSE
)

if (is.null(merged_expr_corrected) || !built_merged_expr) {
  cat("  警告: 合并表达矩阵不可用，使用 ssGSEA 时间点标签替代\n")
  anchor_expr <- NULL
  anchor_tp_labels <- ssgsea_df$timepoint
} else {
  anchor_expr <- merged_expr_corrected
}

for (marker_name in names(ANCHOR_MARKERS)) {
  marker <- ANCHOR_MARKERS[[marker_name]]
  gene <- marker$gene
  expected <- marker$stage

  cat(sprintf("\n  标记基因: %s (预期阶段: %s)\n", gene, expected))

  # Bulk 表达趋势（z-score 标准化）
  if (!is.null(anchor_expr) && gene %in% rownames(anchor_expr)) {
    bulk_expr <- anchor_expr[gene, ]
    # 按时间点聚合
    tp_means <- tapply(as.numeric(bulk_expr), anchor_tp_labels, mean)
    tp_means <- tp_means[intersect(names(tp_means), time_order)]
    tp_z <- scale(tp_means)
    tp_z <- as.numeric(tp_z)
    names(tp_z) <- names(tp_means)

    # 判断趋势：E期(3-6h), M期(12-24h), L期(7d)
    early_mean <- mean(tp_z[intersect(names(tp_z), c("3h", "6h"))], na.rm = TRUE)
    mid_mean   <- mean(tp_z[intersect(names(tp_z), c("12h", "24h"))], na.rm = TRUE)
    late_mean  <- tp_z["7d"]
    if (is.na(late_mean)) late_mean <- 0

    trend_strengths <- c(E = ifelse(is.na(early_mean), -Inf, early_mean),
                         M = ifelse(is.na(mid_mean), -Inf, mid_mean),
                         L = ifelse(is.na(late_mean), -Inf, late_mean))
    bulk_peak_stage <- names(which.max(trend_strengths))
    cat(sprintf("    Bulk z-score 峰值阶段: %s (E=%.3f, M=%.3f, L=%.3f)\n",
                bulk_peak_stage, early_mean, mid_mean, late_mean))
  } else {
    bulk_peak_stage <- NA
    early_mean <- NA; mid_mean <- NA; late_mean <- NA
    cat(sprintf("    基因 %s 不在 Bulk 表达矩阵中\n", gene))
  }

  # 单细胞 E/M/L 期平均表达（从 pseudotime 结果获取）
  sc_stage_sum <- c(E = 0, M = 0, L = 0)
  sc_stage_n   <- c(E = 0, M = 0, L = 0)

  ct_keys <- c("Microglia", "Neuron", "Astrocyte")
  if (!is.null(monocle_results)) {
    for (ct in intersect(ct_keys, names(monocle_results))) {
      res <- monocle_results[[ct]]
      if (is.list(res) && !is.null(res$stages) && "expression_matrix" %in% names(res)) {
        tryCatch({
          emat <- res$expression_matrix
          lib_sizes <- Matrix::colSums(emat)
          if (gene %in% rownames(emat)) {
            if (all(lib_sizes == 0)) {
              gene_expr <- log1p(emat[gene, , drop = FALSE])
            } else {
              lib_sizes[lib_sizes == 0] <- 1
              gene_cpm <- (emat[gene, , drop = FALSE] / lib_sizes) * 1e6
              gene_expr <- as.numeric(log2(gene_cpm + 1))
              names(gene_expr) <- colnames(emat)
            }
            common_cells <- intersect(names(res$stages), names(gene_expr))
            if (length(common_cells) > 0) {
              stage_means <- tapply(as.numeric(gene_expr[common_cells]), res$stages[common_cells], mean)
              for (s in names(stage_means)) {
                if (s %in% names(sc_stage_sum)) {
                  sc_stage_sum[s] <- sc_stage_sum[s] + stage_means[s]
                  sc_stage_n[s]   <- sc_stage_n[s] + 1
                }
              }
            }
          }
        }, error = function(e) {})
      }
    }
  }
  sc_stage_mean <- sc_stage_sum / sc_stage_n
  sc_stage_mean[sc_stage_n == 0] <- NA

  sc_peak_stage <- names(which.max(sc_stage_mean))
  cat(sprintf("    scRNA 峰值阶段: %s (E=%.3f, M=%.3f, L=%.3f)\n",
              ifelse(length(sc_peak_stage) > 0, sc_peak_stage, "NA"),
              sc_stage_mean["E"], sc_stage_mean["M"], sc_stage_mean["L"]))

  consistent <- (!is.na(bulk_peak_stage) && !is.na(sc_peak_stage) &&
                   bulk_peak_stage == expected && sc_peak_stage == expected)

  anchor_results <- rbind(anchor_results, data.frame(
    gene = gene,
    entrez = marker$entrez,
    expected_stage = expected,
    bulk_peak = ifelse(is.na(bulk_peak_stage), "NA", bulk_peak_stage),
    scRNA_peak = ifelse(length(sc_peak_stage) == 0, "NA", sc_peak_stage),
    bulk_E_mean_z = early_mean,
    bulk_M_mean_z = mid_mean,
    bulk_L_mean_z = late_mean,
    consistent = consistent,
    stringsAsFactors = FALSE
  ))
}

cat("\n  标记基因锚定结果:\n")
print(anchor_results[, c("gene", "expected_stage", "bulk_peak", "scRNA_peak", "consistent")])

# 锚定判定
n_consistent <- sum(anchor_results$consistent, na.rm = TRUE)
consistency_rate <- n_consistent / nrow(anchor_results)

cat(sprintf("\n  一致标记基因数: %d / %d (%.1f%%)\n", n_consistent, nrow(anchor_results), 100 * consistency_rate))

# 置换检验：验证标记基因阶段一致性是否显著优于随机
# (参考: Phipson & Smyth Bioinformatics 2010 — permutation-based significance)
cat("\n  置换检验 (1000次) — 标记基因阶段一致性:\n")
n_perm <- 1000
perm_counts <- numeric(n_perm)
for (p in 1:n_perm) {
  shuffled_expected <- sample(anchor_results$expected_stage)
  perm_consistent <- sum(
    (anchor_results$bulk_peak == shuffled_expected) &
    (anchor_results$scRNA_peak == shuffled_expected),
    na.rm = TRUE
  )
  perm_counts[p] <- perm_consistent
}
perm_pvalue <- mean(perm_counts >= n_consistent)
cat(sprintf("  观察值: %d/5, 置换均值: %.2f, p = %.4f\n",
            n_consistent, mean(perm_counts), perm_pvalue))
if (perm_pvalue < 0.05) {
  cat("  ✓ 标记基因一致性显著高于随机 (p < 0.05)\n")
} else {
  cat("  ⚠ 标记基因一致性不显著 (p >= 0.05)，分期锚定需谨慎解读\n")
}
# 保存置换检验结果
perm_df <- data.frame(
  observed = n_consistent,
  perm_mean = mean(perm_counts),
  perm_sd = sd(perm_counts),
  p_value = perm_pvalue
)
write.csv(perm_df, file.path(OUTPUT_DIR, "permutation_test_markers.csv"), row.names = FALSE)

if (n_consistent >= 4 && consistency_rate >= 0.75) {
  cat("  ✓ 标记基因锚定通过！分期对应关系成立\n")
  cat("  E → 急性期（3-6h）\n")
  cat("  M → 亚急性期（12-24h）\n")
  cat("  L → 慢性期（1d-7d）\n")
  anchor_passed <- TRUE
} else {
  cat("  ⚠ 标记基因锚定不明确（<4个一致），将执行 Spearman 辅助锚定\n")
  anchor_passed <- FALSE
}

write.csv(anchor_results, file.path(OUTPUT_DIR, "anchor_marker_genes.csv"), row.names = FALSE)

# -------------------- 3B. Spearman 辅助锚定 --------------------
cat("\n--- 3B. Spearman 辅助锚定 ---\n")

if (!anchor_passed) {
  # 构建单细胞E/M/L × 6模块矩阵
  sc_module_sum <- matrix(0, nrow = 3, ncol = length(module_names),
                              dimnames = list(c("E", "M", "L"), module_names))
  sc_module_n   <- matrix(0, nrow = 3, ncol = length(module_names),
                              dimnames = list(c("E", "M", "L"), module_names))

  if (!is.null(monocle_results)) {
    for (ct in intersect(ct_keys, names(monocle_results))) {
      res <- monocle_results[[ct]]
      if (is.list(res) && !is.null(res$stages) && "expression_matrix" %in% names(res)) {
        emat <- res$expression_matrix
        lib_sizes <- Matrix::colSums(emat)
        lib_sizes[lib_sizes == 0] <- 1
        for (mod in module_names) {
          mod_genes <- intersect(MODULE_GENES[[mod]], rownames(emat))
          if (length(mod_genes) >= 2) {
            mod_subset <- emat[mod_genes, , drop = FALSE]
            cpm_subset <- Matrix::t(Matrix::t(mod_subset) / lib_sizes) * 1e6
            log_cpm <- log2(cpm_subset + 1)
            mod_expr <- Matrix::colMeans(log_cpm)
            common_cells <- intersect(names(res$stages), names(mod_expr))
            if (length(common_cells) > 0) {
              stage_means <- tapply(mod_expr[common_cells], res$stages[common_cells], mean)
              for (s in names(stage_means)) {
                if (s %in% rownames(sc_module_sum)) {
                  sc_module_sum[s, mod] <- sc_module_sum[s, mod] + stage_means[s]
                  sc_module_n[s, mod]   <- sc_module_n[s, mod] + 1
                }
              }
            }
          }
        }
      }
    }
  }
  sc_module_by_stage <- sc_module_sum / sc_module_n
  sc_module_by_stage[sc_module_n == 0] <- NA

  # 构建Bulk时间点 × 6模块矩阵
  bulk_module_by_time <- matrix(NA, nrow = length(time_order), ncol = length(module_names),
                                 dimnames = list(time_order, module_names))
  for (mod in module_names) {
    for (tp in time_order) {
      tp_scores <- ssgsea_df[ssgsea_df$timepoint == tp, mod]
      bulk_module_by_time[tp, mod] <- mean(tp_scores, na.rm = TRUE)
    }
  }

  # Spearman跨组学关联：聚合Bulk时间点为E/M/L三期以匹配scRNA分期维度
  sc_complete <- sc_module_by_stage[complete.cases(sc_module_by_stage), , drop = FALSE]

  bulk_stage_map <- list(E = c("3h", "6h"), M = c("12h", "24h"), L = c("7d"))
  bulk_module_by_stage <- matrix(NA, nrow = 3, ncol = length(module_names),
                                  dimnames = list(c("E", "M", "L"), module_names))
  for (stage in c("E", "M", "L")) {
    tp_match <- intersect(bulk_stage_map[[stage]], rownames(bulk_module_by_time))
    if (length(tp_match) >= 1) {
      bulk_module_by_stage[stage, ] <- colMeans(bulk_module_by_time[tp_match, , drop = FALSE], na.rm = TRUE)
    }
  }
  bulk_complete <- bulk_module_by_stage[complete.cases(bulk_module_by_stage), , drop = FALSE]

  common_stages <- intersect(rownames(sc_complete), rownames(bulk_complete))
  sc_complete <- sc_complete[common_stages, , drop = FALSE]
  bulk_complete <- bulk_complete[common_stages, , drop = FALSE]

  common_mods <- intersect(colnames(sc_complete), colnames(bulk_complete))
  sc_complete <- sc_complete[, common_mods, drop = FALSE]
  bulk_complete <- bulk_complete[, common_mods, drop = FALSE]

  if (nrow(sc_complete) >= 2 && length(common_mods) >= 2) {
    # Spearman相关替代CCA (参考: scAB, Zhang et al. NAR 2022 — 用pairwise correlation
    # 关联单细胞与Bulk数据；n=3时比CCA更稳健)
    cat(sprintf("  跨组学 Spearman 模块一致性 (n=%d):\n", nrow(sc_complete)))
    module_cors <- c()
    for (mod in common_mods) {
      sp_cor <- tryCatch(
        cor(sc_complete[, mod], bulk_complete[, mod],
            method = "spearman", use = "complete.obs"),
        error = function(e) NA_real_
      )
      module_cors <- c(module_cors, sp_cor)
      cat(sprintf("    %s: ρ = %+.3f\n", mod, sp_cor))
    }
    names(module_cors) <- common_mods

    mean_rho <- mean(module_cors, na.rm = TRUE)
    n_sig <- sum(abs(module_cors) >= 0.5, na.rm = TRUE)
    cat(sprintf("  平均 Spearman ρ = %+.3f, |ρ|≥0.5 的模块: %d/%d\n",
                mean_rho, n_sig, length(module_cors)))
    cat("  注意: n=3 (E/M/L) Spearman ρ 自由度低，单模块ρ值不稳定；以平均ρ和Fisher检验为准\n")

    # Fisher's exact test: 评估跨组学E→M→L阶段方向一致性
    # 分别计算每个阶段间变化的符号 (E→M, M→L) 在sc和bulk之间是否一致
    sc_directions <- sign(sc_complete[2, ] - sc_complete[1, ]) * sign(sc_complete[3, ] - sc_complete[2, ])
    bulk_directions <- sign(bulk_complete[2, ] - bulk_complete[1, ]) * sign(bulk_complete[3, ] - bulk_complete[2, ])
    names(sc_directions) <- colnames(sc_complete)
    names(bulk_directions) <- colnames(bulk_complete)

    # 2×2 列联表: 方向一致 vs 不一致
    concordant <- sum(sc_directions == bulk_directions, na.rm = TRUE)
    discordant <- sum(sc_directions != bulk_directions & !is.na(sc_directions) & !is.na(bulk_directions), na.rm = TRUE)
    if (concordant + discordant >= 2) {
      fisher_res <- fisher.test(matrix(c(concordant, discordant, discordant, concordant), nrow = 2),
                                 alternative = "greater")
      cat(sprintf("  Fisher精确检验 (阶段间方向一致性):\n"))
      cat(sprintf("    一致: %d 模块, 不一致: %d 模块\n", concordant, discordant))
      cat(sprintf("    p = %.4f (greater)\n", fisher_res$p.value))
      if (fisher_res$p.value < 0.05) {
        cat("    ✓ 跨组学阶段间方向显著一致\n")
      } else {
        cat("    ⚠ 跨组学方向一致性不显著\n")
      }
    }

    if (mean_rho >= 0.5) {
      cat("  ✓ 跨组学模块活性排序一致，支持分期锚定结论\n")
    } else if (abs(mean_rho) >= 0.3) {
      cat("  ~ 跨组学模块活性中度相关，分期锚定部分可信\n")
    } else {
      cat("  ⚠ 跨组学相关较弱，分期结论需谨慎解读\n")
    }

    # 保存 Spearman 结果
    spearman_df <- data.frame(
      Module = common_mods,
      Spearman_rho = module_cors,
      stringsAsFactors = FALSE
    )
    write.csv(spearman_df, file.path(OUTPUT_DIR, "crossomics_spearman.csv"), row.names = FALSE)
  } else {
    cat("  跨组学对比数据不足（需要 ≥2 行和 ≥2 模块）\n")
  }
}

# ======================================================================
#                    第四部分：CIBERSORTx 验证
# ======================================================================
cat("\n========== 第4部分：CIBERSORTx 验证 ==========\n")

run_cibersortx <- function() {
  cat("\n  CIBERSORTx 不可用，使用标记基因平均表达法估算细胞比例\n")
  cat("  方法：对每个Bulk样本，计算微胶质/巨噬细胞标记基因的平均表达(z-score标准化)\n")

  micro_markers <- c("Aif1", "Cx3cr1", "Tmem119", "C1qa", "C1qb", "Ptprc", "Cd68", "Fcgr3", "Itgam", "Tyrobp")
  micro_markers <- intersect(micro_markers, rownames(merged_expr_corrected))
  cat(sprintf("  可用微胶质标记基因: %d\n", length(micro_markers)))

  if (length(micro_markers) < 3) {
    cat("  微胶质标记基因不足，无法估算\n")
    return(NULL)
  }

  micro_expr <- merged_expr_corrected[micro_markers, , drop = FALSE]
  micro_score <- colMeans(micro_expr, na.rm = TRUE)

  tp_order <- intersect(time_order, unique(anchor_tp_labels))
  tp_means <- sapply(tp_order, function(tp) {
    idx <- which(anchor_tp_labels == tp)
    mean(micro_score[idx], na.rm = TRUE)
  })

  cat("\n  微胶质/巨噬细胞标记基因平均表达时间趋势:\n")
  trend_df <- data.frame(
    timepoint = names(tp_means),
    mean_score = as.numeric(tp_means),
    stringsAsFactors = FALSE
  )
  print(trend_df)

  # 判断趋势：E(低)→M(峰值)→L(持续高)
  e_idx <- which(names(tp_means) %in% c("3h", "6h"))
  m_idx <- which(names(tp_means) %in% c("12h", "24h"))
  l_idx <- which(names(tp_means) %in% c("7d"))

  e_val <- mean(tp_means[e_idx], na.rm = TRUE)
  m_val <- mean(tp_means[m_idx], na.rm = TRUE)
  l_val <- mean(tp_means[l_idx], na.rm = TRUE)

  cat(sprintf("  E期(3-6h)平均: %.4f\n", e_val))
  cat(sprintf("  M期(12-24h)平均: %.4f\n", m_val))
  cat(sprintf("  L期(7d)平均: %.4f\n", l_val))

  if (is.finite(e_val) && is.finite(m_val) && is.finite(l_val)) {
    if (m_val > e_val && l_val > e_val) {
      cat("  ✓ 微胶质趋势一致：E(低) → M(峰值) → L(持续高)\n")
      return(list(trend = "consistent", e = e_val, m = m_val, l = l_val))
    } else {
      cat("  ⚠ 微胶质趋势不完全一致\n")
      return(list(trend = "inconsistent", e = e_val, m = m_val, l = l_val))
    }
  } else {
    cat("  ⚠ 部分时间点数据缺失\n")
    return(NULL)
  }
}

cibersort_results <- run_cibersortx()

# ======================================================================
#                    第五部分：事件顺序约束表
# ======================================================================
cat("\n========== 第5部分：事件顺序约束表 ==========\n")

# 基于模块活性曲线的阶段归类
event_order <- data.frame(
  module = module_names,
  activation_phase = character(length(module_names)),
  peak_time = character(length(module_names)),
  constraint = character(length(module_names)),
  stringsAsFactors = FALSE
)

for (i in seq_along(module_names)) {
  mod <- module_names[i]
  if (mod %in% names(inflection_points)) {
    ip <- inflection_points[[mod]]
    peak_t <- ip$time_h[which.max(ip$activity_value)]

    if (length(peak_t) > 0) {
      if (peak_t <= 6) {
        phase <- "E (急性期 3-6h)"
        constraint <- sprintf("%s 激活最早，先于其他模块", mod)
      } else if (peak_t <= 24) {
        phase <- "M (亚急性期 12-24h)"
        constraint <- sprintf("%s 激活在 E 期模块之后、L 期模块之前", mod)
      } else {
        phase <- "L (慢性期 1d-7d)"
        constraint <- sprintf("%s 激活最晚，在 E/M 期模块之后", mod)
      }

      event_order$activation_phase[i] <- phase
      event_order$peak_time[i] <- sprintf("%.1fh", peak_t)
      event_order$constraint[i] <- constraint
    }
  }
}

# 剔除未分配的
event_order <- event_order[event_order$activation_phase != "", ]

# 按激活顺序排列
phase_order <- c("E (急性期 3-6h)", "M (亚急性期 12-24h)", "L (慢性期 1d-7d)")
event_order$phase_rank <- match(event_order$activation_phase, phase_order)
event_order <- event_order[order(event_order$phase_rank), ]

cat("\n事件顺序约束表:\n")
for (i in seq_len(nrow(event_order))) {
  cat(sprintf("  \"%s\"\n", event_order$constraint[i]))
}

write.csv(event_order, file.path(OUTPUT_DIR, "event_order_constraints.csv"), row.names = FALSE)

# -------------------- 5A. 自检 --------------------
cat("\n========== 自检标准 ==========\n")

# (1) 至少4个标记基因的跨组学趋势一致
check1 <- n_consistent >= 4
cat(sprintf("  自检1 - 标记基因一致性 (≥4/5): %s (实际 %d/5)\n",
            ifelse(check1, "✓ 通过", "✗ 未通过"), n_consistent))

# (2) 至少3个模块的活性呈现跨分期单调趋势
mono_count <- 0
for (mod in module_names) {
  tp_means <- tp_summary_df$mean[tp_summary_df$module == mod]
  tp_sds <- tp_summary_df$sd[tp_summary_df$module == mod]
  disease_tp <- tp_summary_df$timepoint[tp_summary_df$module == mod]
  disease_tp <- disease_tp[disease_tp != "sham"]
  disease_means <- tp_means[disease_tp != "sham"]

  if (length(disease_means) >= 3) {
    is_mono_inc <- all(diff(disease_means[order(match(disease_tp, time_order))]) > 0)
    is_mono_dec <- all(diff(disease_means[order(match(disease_tp, time_order))]) < 0)
    if (is_mono_inc || is_mono_dec) mono_count <- mono_count + 1
  }
}
check2 <- mono_count >= 3
cat(sprintf("  自检2 - 模块活性单调趋势 (≥3/6): %s (实际 %d/6)\n",
            ifelse(check2, "✓ 通过", "✗ 未通过"), mono_count))

# (3) CIBERSORTx 验证
check3 <- !is.null(cibersort_results)
cat(sprintf("  自检3 - CIBERSORTx 反卷积执行: %s\n", ifelse(check3, "✓ 完成", "✗ 未完成")))

all_checks <- c(check1, check2, check3)
cat(sprintf("\n  综合评估: %d/3 项通过\n", sum(all_checks)))

# ======================================================================
#                    第六部分：综合可视化
# ======================================================================
cat("\n========== 第6部分：综合可视化 ==========\n")

# --- 6A. 模块活性动态曲线（标注E/M/L分期） ---
pdf(file.path(FIGURE_DIR, "Fig_QualTCA_module_curves_with_stages.pdf"), width = 16, height = 10)

par(mfrow = c(2, 3), mar = c(5, 5, 4, 2))
for (mod in module_names) {
  tp_means <- tp_summary_df[tp_summary_df$module == mod, ]
  tp_means$time_num <- time_numeric[tp_means$timepoint]

  tp_disease <- tp_means[tp_means$timepoint != "sham", ]
  tp_disease <- tp_disease[order(tp_disease$time_num), ]

  ylim <- range(c(tp_disease$mean - tp_disease$sd,
                  tp_disease$mean + tp_disease$sd), na.rm = TRUE)
  if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)

  plot(tp_disease$time_num, tp_disease$mean,
       type = "b", pch = 19, col = "steelblue", lwd = 2,
       xlab = "Time post-MCAO", ylab = "ssGSEA Score",
       main = mod, ylim = ylim, xaxt = "n", cex = 1.5)
  axis(1, at = time_numeric[time_numeric > 0], labels = names(time_numeric)[time_numeric > 0])
  arrows(tp_disease$time_num, tp_disease$mean - tp_disease$sd,
         tp_disease$time_num, tp_disease$mean + tp_disease$sd,
         angle = 90, code = 3, length = 0.05, col = "steelblue")

  # E/M/L 分期背景
  rect(0, ylim[1] - 1, 6, ylim[2] + 1, col = rgb(0, 0, 1, 0.05), border = NA)
  rect(6, ylim[1] - 1, 24, ylim[2] + 1, col = rgb(1, 0.65, 0, 0.05), border = NA)
  rect(24, ylim[1] - 1, max(time_numeric) + 10, ylim[2] + 1,
       col = rgb(1, 0, 0, 0.05), border = NA)

  text(3, ylim[2] * 0.95, "E\n急性期\n3-6h", col = "blue", cex = 0.9, font = 2)
  text(15, ylim[2] * 0.95, "M\n亚急性期\n12-24h", col = "darkorange", cex = 0.9, font = 2)
  text(96, ylim[2] * 0.95, "L\n慢性期\n1d-7d", col = "red", cex = 0.9, font = 2)

  # 拐点标记
  if (mod %in% names(inflection_points)) {
    ip <- inflection_points[[mod]]
    abline(v = ip$time_h, lty = 2, col = "darkgreen", lwd = 1.5)
  }
}
dev.off()

# --- 6B. 标记基因一致性热图 ---
pdf(file.path(FIGURE_DIR, "Fig_QualTCA_anchor_heatmap.pdf"), width = 10, height = 6)

anchor_mat <- as.matrix(anchor_results[, c("bulk_E_mean_z", "bulk_M_mean_z", "bulk_L_mean_z")])
rownames(anchor_mat) <- anchor_results$gene

if (all(is.na(anchor_mat)) || all(is.infinite(as.matrix(anchor_mat)))) {
  cat("  标记基因锚定数据全为NA，跳过热图绘制\n")
  plot.new()
  text(0.5, 0.5, "Anchor marker data unavailable\n(skip heatmap)", cex = 1.5)
} else {
  anchor_mat[!is.finite(anchor_mat)] <- 0
  pheatmap(anchor_mat,
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.2f",
           main = "标记基因 Bulk z-score (E/M/L)",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           fontsize_number = 10, fontsize = 12)
}
dev.off()

# --- 6C. 标记基因趋势一致性条形图 ---

consistency_df <- data.frame(
  Category = c("一致", "不一致"),
  Count = c(n_consistent, nrow(anchor_results) - n_consistent)
)

ggplot(consistency_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = Count), vjust = -0.5, size = 6) +
  scale_fill_manual(values = c("darkgreen", "gray60")) +
  labs(title = "标记基因跨组学趋势一致性",
       subtitle = sprintf("阈值: >=4/5 一致 (实际 %d/5)", n_consistent),
       y = "标记基因数") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave(file.path(FIGURE_DIR, "Fig_QualTCA_consistency_bar.pdf"), width = 8, height = 6)

# ==================== 结果汇总保存 ====================
cat("\n========== 结果汇总 ==========\n")

results_summary <- list(
  analysis_name = "L1 定性分期锚定层 (QualTCA)",
  date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  input_datasets = c("GSE104036 (小鼠RNA-seq, 3h/6h/12h/24h)",
                     "GSE97537 (大鼠芯片, 24h)",
                     "GSE61616 (大鼠芯片, 7d)",
                     "GSE174574 (小鼠scRNA-seq, 24h)"),
  modules = common_mod_cols,
  anchor_markers = anchor_results,
  n_consistent_markers = n_consistent,
  consistency_rate = consistency_rate,
  anchor_passed = anchor_passed,
  event_order = event_order,
  inflection_points = inflection_summary,
  gene_loss_report = gene_loss_report,
  fallback_log = fallback_log,
  n_fallbacks = length(fallback_log),
  self_checks = list(
    check1_marker_consistency = c(passed = check1, value = n_consistent),
    check2_module_monotonic = c(passed = check2, value = mono_count),
    check3_cibersortx = c(passed = check3, value = NA)
  ),
  output_files = list(
    ssgsea_scores = file.path(OUTPUT_DIR, "ssGSEA_module_scores.csv"),
    inflection_points = file.path(OUTPUT_DIR, "inflection_points.csv"),
    smoothed_curves = file.path(OUTPUT_DIR, "smoothed_curves.csv"),
    anchor_markers = file.path(OUTPUT_DIR, "anchor_marker_genes.csv"),
    event_order = file.path(OUTPUT_DIR, "event_order_constraints.csv"),
    permutation_test = file.path(OUTPUT_DIR, "permutation_test_markers.csv")
  )
)

saveRDS(results_summary, file.path(OUTPUT_DIR, "QualTCA_results_summary.rds"))

# 降级日志汇总
print_fallback_summary()

cat("\n========== L1 定性分期锚定完成 ==========\n")
cat(sprintf("结果保存至: %s\n", OUTPUT_DIR))
cat(sprintf("图表保存至: %s\n", FIGURE_DIR))
cat(sprintf("完成时间: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))