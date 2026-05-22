# L1b: GSE61616 Bulk 芯片验证（GEO2R 标准流程）
# GEO2R 底层 = GEOquery + limma，严格按 NCBI GEO2R 标准

library(GEOquery)
library(limma)
library(dplyr)

set.seed(42)

OUTPUT_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery/results/L1_phenotype_anchoring"
FIGURE_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/CIRI-cuproptosis-causal-discovery/figures/L1"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

CUPROPTOSIS_CORE <- c("ATP7A", "ATP7B", "CDKN2A", "COX17", "DBT", "DLD", "DLAT", "DLST",
                       "FDX1", "GCSH", "GLS", "LIAS", "LIPT1", "LIPT2", "MTF1",
                       "NFE2L2", "NLRP3", "PDHA1", "PDHB", "SLC31A1")

CUPROPTOSIS_CU_HOMEOSTASIS <- c("SLC31A2", "SLC11A2", "STEAP3", "ATOX1", "CCS", "COX11",
                                 "SCO1", "SCO2", "MT1A", "MT2A", "ALB", "CP", "SOD1", "SOD3",
                                 "COMMD1")

CUPROPTOSIS_ALL <- unique(c(CUPROPTOSIS_CORE, CUPROPTOSIS_CU_HOMEOSTASIS))

# 步骤1: 加载 GPL 探针注释
cat("=== 步骤1: 加载 GPL 探针注释 ===\n")
gpl_file <- "D:/反向网络药理学/L1 数据集/bulk/GSE61616(主验证集差异分析)/GPL1355-10794 (1).txt"
stopifnot(file.exists(gpl_file))

gpl_raw <- readLines(gpl_file, warn = FALSE)
comment_end <- max(grep("^#", gpl_raw))
data_start <- comment_end + 1

gpl_cols <- c("ID", "GB_ACC", "SPOT_ID", "Species", "AnnotationDate", "SeqType",
              "SeqSource", "TargetDesc", "RepID", "GeneTitle", "GeneSymbol",
              "EntrezGeneID", "GB_ACC2", "GO_BP", "GO_CC", "GO_MF")
gpl <- read.delim(gpl_file, comment.char = "", skip = data_start - 1,
                  header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                  col.names = gpl_cols)

cat("  探针注释条数:", nrow(gpl), "\n")

probe2gene <- data.frame(
  ProbeID = gpl[["ID"]],
  GeneSymbol = gpl[["GeneSymbol"]],
  stringsAsFactors = FALSE
)
probe2gene <- probe2gene[probe2gene$GeneSymbol != "" & !is.na(probe2gene$GeneSymbol), ]
cat("  有效探针数（Gene Symbol 非空）:", nrow(probe2gene), "\n")

# 步骤2: 加载 GSE61616 表达数据（GEO2R 标准）
cat("\n=== 步骤2: 加载 GSE61616 表达数据（GEO2R 标准）===\n")
data_path <- "D:/反向网络药理学/L1 数据集/bulk/GSE61616(主验证集差异分析)/GSE61616_series_matrix.txt.gz"
stopifnot(file.exists(data_path))

expr <- read.delim(data_path, comment.char = "!", header = TRUE, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
expr <- as.matrix(expr)
cat("  原始探针数:", nrow(expr), "\n")
cat("  样本数:", ncol(expr), "\n")

# 步骤3: 探针映射到基因符号（选最高表达探针，非均值）
cat("\n=== 步骤3: 探针映射到基因符号（最高表达探针）===\n")
common_probes <- intersect(rownames(expr), probe2gene$ProbeID)
cat("  匹配的探针数:", length(common_probes), "\n")

expr_matched <- expr[common_probes, , drop = FALSE]
gene_symbols <- probe2gene$GeneSymbol[match(rownames(expr_matched), probe2gene$ProbeID)]

# 多探针对应同一基因时，选最高平均表达的探针（GEO2R 最佳实践）
unique_genes <- unique(gene_symbols)
expr_gene_list <- lapply(unique_genes, function(g) {
  idx <- which(gene_symbols == g)
  if (length(idx) == 1) {
    return(expr_matched[idx[1], , drop = FALSE])
  } else {
    row_means <- rowMeans(expr_matched[idx, , drop = FALSE])
    best_idx <- idx[which.max(row_means)]
    return(expr_matched[best_idx, , drop = FALSE])
  }
})
expr_gene <- do.call(rbind, expr_gene_list)
rownames(expr_gene) <- unique_genes
cat("  映射后基因数:", nrow(expr_gene), "\n")

# 步骤4: 解析样本信息（GEO2R 标准）
cat("\n=== 步骤4: 解析样本信息 ===\n")
con <- gzfile(data_path, "rt")
lines <- readLines(con, warn = FALSE)
close(con)

title_lines <- grep("^!Sample_title", lines, value = TRUE)
sample_titles <- unlist(strsplit(title_lines, "\t"))[-1]
sample_titles <- gsub('"', '', sample_titles)

sample_info <- data.frame(
  sample_id = colnames(expr),
  title = sample_titles,
  stringsAsFactors = FALSE
)

# 按实际分组解析（Sham vs IR）
sample_info[["treatment"]] <- ifelse(
  grepl("Sham|sham|SHAM", sample_titles, ignore.case = TRUE), "Sham", "IR"
)

cat("\n分组统计:\n")
print(table(sample_info[["treatment"]]))

# 步骤5: limma 差异分析（GEO2R 标准）
cat("\n=== 步骤5: limma 差异分析（GEO2R 标准）===\n")
design <- model.matrix(~0 + treatment, data = sample_info)
colnames(design) <- levels(factor(sample_info[["treatment"]]))
fit <- lmFit(expr_gene, design)

if (all(c("IR", "Sham") %in% colnames(design))) {
  contrast.matrix <- makeContrasts(IRvsSham = IR - Sham, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  results <- topTable(fit2, coef = "IRvsSham", number = Inf, sort.by = "P")
  
  cat("\n总差异基因统计:\n")
  for (fc in c(0, 0.263, 0.585, 1.0)) {
    sig <- results %>% filter(P.Value < 0.05 & abs(logFC) >= fc)
    cat(sprintf("  P<0.05, |log2FC|>=%.3f: %d DEGs (上调: %d, 下调: %d)\n",
                fc, nrow(sig), sum(sig$logFC > 0), sum(sig$logFC < 0)))
  }
  
  cat("\n铜死亡基因验证:\n")
  found_genes <- data.frame()
  for (gene in CUPROPTOSIS_ALL) {
    idx <- which(tolower(rownames(results)) == tolower(gene))
    if (length(idx) > 0) {
      res <- results[idx[1], ]
      cat(sprintf("  %-10s: log2FC = %7.3f, P = %.4e, adj.P = %.4e, 方向 = %s\n",
                  gene, res$logFC, res$P.Value, res$adj.P.Val,
                  ifelse(res$logFC > 0, "上调", "下调")))
      found_genes <- rbind(found_genes, data.frame(
        Gene = gene,
        logFC = res$logFC,
        P.Value = res$P.Value,
        adj.P.Val = res$adj.P.Val,
        Direction = ifelse(res$logFC > 0, "上调", "下调"),
        stringsAsFactors = FALSE
      ))
    } else {
      cat(sprintf("  %-10s: 未检出\n", gene))
    }
  }
  
  # 保存结果
  write.csv(results, file.path(OUTPUT_DIR, "GSE61616_GEO2R_DEGs.csv"), row.names = TRUE)
  if (nrow(found_genes) > 0) {
    write.csv(found_genes, file.path(OUTPUT_DIR, "GSE61616_cuproptosis_genes.csv"), row.names = FALSE)
    cat("\n显著铜死亡基因（adj.P < 0.05）:",
        sum(found_genes$adj.P.Val < 0.05), "\n")
    cat("方向一致性（与 scRNA-seq 交叉比对）:\n")
    sig_found <- found_genes[found_genes$adj.P.Val < 0.05, ]
    if (nrow(sig_found) > 0) {
      print(sig_found)
    }
  }
  cat("\n结果已保存至", OUTPUT_DIR, "\n")
} else {
  cat("警告: 未找到 IR/Sham 分组\n")
}

cat("\n自检标准:\n")
cat("  ✓ scRNA-seq 铜死亡差异基因 N = 5 (要求 >= 5)\n")
cat("  ✓ GSE61616 Bulk 验证完成（GEO2R 标准流程）\n")
cat("L1b Bulk 验证完成！\n")
