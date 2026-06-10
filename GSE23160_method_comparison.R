# 验证脚本：对比两种差异分析方法
# 方法1：所有时间点合并 vs Sham（Zhao原始文献方法）
# 方法2：分时间点pairwise比较（我们之前的方法）

setwd("D:/反向网络药理学/L1 数据集/bulk/GSE23160(主验证集时序差异分析，2h,8h,24h)")

library(limma)
library(ggplot2)
library(dplyr)

# ==================== 1. 读取数据 ====================
expr_data <- read.delim("GSE23160_series_matrix.txt.gz", 
                        comment.char = "!", 
                        header = TRUE, 
                        check.names = FALSE,
                        row.names = 1,
                        stringsAsFactors = FALSE)

expr_matrix <- as.matrix(expr_data)

# ==================== 2. 解析样本信息 ====================
con <- gzfile("GSE23160_series_matrix.txt.gz", "rt")
lines <- readLines(con, warn = FALSE)
close(con)

sample_titles_line <- grep("^!Sample_title", lines, value = TRUE)
sample_titles <- unlist(strsplit(sample_titles_line, "\t"))[-1]
sample_titles <- gsub('"', '', sample_titles)

sample_info <- data.frame(
  sample_id = colnames(expr_matrix),
  title = sample_titles,
  stringsAsFactors = FALSE
)

sample_info$brain_region <- ifelse(grepl("^Cortex-", sample_titles, ignore.case = TRUE), "Cortex",
                            ifelse(grepl("^Striatum-", sample_titles, ignore.case = TRUE), "Striatum", "Unknown"))

sample_info$treatment <- ifelse(grepl("-Sham-", sample_titles, ignore.case = TRUE), "Sham", "IR")

sample_info$time_group <- ifelse(sample_info$treatment == "Sham", "Ctrl",
                          ifelse(grepl("-2h-", sample_titles, ignore.case = TRUE), "2h",
                          ifelse(grepl("-8h-", sample_titles, ignore.case = TRUE), "8h",
                          ifelse(grepl("-24h-", sample_titles, ignore.case = TRUE), "24h", "Unknown"))))

# ==================== 3. 方法对比 ====================
compare_methods <- function(brain_region, expr_matrix, sample_info) {
  
  cat(paste0("\n", strrep("=", 60), "\n"))
  cat(paste0("脑区: ", brain_region, "\n"))
  
  region_idx <- sample_info$brain_region == brain_region
  expr_region <- expr_matrix[, region_idx]
  sample_region <- sample_info[region_idx, ]
  
  # ==================== 方法1：合并所有时间点 vs Sham（Zhao方法）====================
  cat("\n--- 方法1：合并所有I/R vs Sham（Zhao原始文献方法）---\n")
  
  # 选择Sham和所有IR样本
  ir_samples <- sample_region$treatment == "IR"
  sham_samples <- sample_region$treatment == "Sham"
  
  expr_compare <- expr_region[, ir_samples | sham_samples]
  sample_compare <- sample_region[ir_samples | sham_samples, ]
  
  sample_compare$group <- factor(sample_compare$treatment, levels = c("Sham", "IR"))
  
  design <- model.matrix(~0 + group, data = sample_compare)
  colnames(design) <- c("Sham", "IR")
  
  fit <- lmFit(expr_compare, design)
  contrast.matrix <- makeContrasts(IRvsSham = IR - Sham, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  results_zhao <- topTable(fit2, coef = "IRvsSham", number = Inf, sort.by = "P")
  
  sig_zhao <- results_zhao %>% filter(P.Value < 0.05 & abs(logFC) >= 0.585)
  
  cat(paste0("  总探针数: ", nrow(results_zhao), "\n"))
  cat(paste0("  DEGs (|log2FC|>=0.585, P<0.05): ", nrow(sig_zhao), "\n"))
  cat(paste0("    上调: ", sum(sig_zhao$logFC > 0), "\n"))
  cat(paste0("    下调: ", sum(sig_zhao$logFC < 0), "\n"))
  
  # 检查铜死亡基因
  copper_genes <- c("Fdx1", "Lias", "Lipt1", "Dlat", "Pdhb", "Pdhx", "Dld", "Dbt",
                    "Slc31a1", "Atp7a", "Atp7b", "Atox1", "Mtf1",
                    "Nfe2l2", "Hif1a", "Mtor", "Nfkb1", "Keap1",
                    "Cdkn2a", "Trp53", "Slc25a3", "Slc25a5",
                    "Sdha", "Sdhb", "Sdhc", "Sdhd",
                    "Cs", "Aco1", "Aco2", "Idh2", "Ogdh",
                    "Gcsh", "Gls", "Glud1")
  
  # 需要加载平台注释
  platform_annot <- read.delim("GPL6885-11608.txt", 
                               header = TRUE,
                               check.names = FALSE,
                               stringsAsFactors = FALSE,
                               quote = "",
                               fill = TRUE,
                               comment.char = "#")
  
  annot_table <- data.frame(
    ID = platform_annot[[1]],
    Gene_Symbol = if("Symbol" %in% names(platform_annot)) platform_annot[["Symbol"]] else
                  if("Gene Symbol" %in% names(platform_annot)) platform_annot[["Gene Symbol"]] else
                  if("GENE_SYMBOL" %in% names(platform_annot)) platform_annot[["GENE_SYMBOL"]] else NA,
    stringsAsFactors = FALSE
  )
  
  results_zhao$ProbeID <- rownames(results_zhao)
  results_zhao$Gene_Symbol <- annot_table$Gene_Symbol[match(results_zhao$ProbeID, annot_table$ID)]
  
  copper_found <- results_zhao %>%
    filter(!is.na(Gene_Symbol) & Gene_Symbol %in% copper_genes)
  
  cat(paste0("  铜死亡基因探针: ", nrow(copper_found), "\n"))
  copper_sig <- copper_found %>% filter(P.Value < 0.05 & abs(logFC) >= 0.585)
  if (nrow(copper_sig) > 0) {
    cat("  显著差异铜死亡基因:\n")
    print(copper_sig[, c("Gene_Symbol", "logFC", "P.Value")])
  }
  
  # ==================== 方法2：分时间点pairwise（我们的方法）====================
  cat("\n--- 方法2：分时间点pairwise比较 ---\n")
  
  time_points <- c("2h", "8h", "24h")
  
  for (tp in time_points) {
    tp_idx <- sample_region$time_group == tp | sample_region$time_group == "Ctrl"
    expr_tp <- expr_region[, tp_idx]
    sample_tp <- sample_region[tp_idx, ]
    
    sample_tp$group <- factor(sample_tp$treatment, levels = c("Sham", "IR"))
    
    design <- model.matrix(~0 + group, data = sample_tp)
    colnames(design) <- c("Sham", "IR")
    
    fit <- lmFit(expr_tp, design)
    contrast.matrix <- makeContrasts(IRvsSham = IR - Sham, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    results_tp <- topTable(fit2, coef = "IRvsSham", number = Inf, sort.by = "P")
    
    sig_tp <- results_tp %>% filter(P.Value < 0.05 & abs(logFC) >= 0.585)
    
    cat(paste0("  ", tp, ": ", nrow(sig_tp), " DEGs (上调: ", sum(sig_tp$logFC > 0), 
               ", 下调: ", sum(sig_tp$logFC < 0), ")\n"))
  }
  
  # 方法2总DEGs（去重后）
  all_sig <- list()
  for (tp in time_points) {
    tp_idx <- sample_region$time_group == tp | sample_region$time_group == "Ctrl"
    expr_tp <- expr_region[, tp_idx]
    sample_tp <- sample_region[tp_idx, ]
    
    sample_tp$group <- factor(sample_tp$treatment, levels = c("Sham", "IR"))
    design <- model.matrix(~0 + group, data = sample_tp)
    colnames(design) <- c("Sham", "IR")
    
    fit <- lmFit(expr_tp, design)
    contrast.matrix <- makeContrasts(IRvsSham = IR - Sham, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    results_tp <- topTable(fit2, coef = "IRvsSham", number = Inf, sort.by = "P")
    sig_tp <- results_tp %>% filter(P.Value < 0.05 & abs(logFC) >= 0.585)
    sig_tp$time_point <- tp
    all_sig[[tp]] <- sig_tp
  }
  
  all_combined <- do.call(rbind, all_sig)
  total_degs <- length(unique(all_combined$ProbeID))
  
  cat(paste0("  方法2总DEGs（去重后）: ", total_degs, "\n"))
  
  return(list(
    method1_degs = nrow(sig_zhao),
    method2_degs = total_degs,
    copper_sig_zhao = copper_sig
  ))
}

# ==================== 4. 运行对比 ====================
cortex_results <- compare_methods("Cortex", expr_matrix, sample_info)
striatum_results <- compare_methods("Striatum", expr_matrix, sample_info)

cat("\n", strrep("=", 60), "\n")
cat("方法对比总结\n")
cat(strrep("=", 60), "\n")

cat(paste0("Cortex:\n"))
cat(paste0("  方法1 (Zhao): ", cortex_results$method1_degs, " DEGs\n"))
cat(paste0("  方法2 (我们的): ", cortex_results$method2_degs, " DEGs\n"))
cat(paste0("  文献报道: ~2295 DEGs\n"))

cat(paste0("Striatum:\n"))
cat(paste0("  方法1 (Zhao): ", striatum_results$method1_degs, " DEGs\n"))
cat(paste0("  方法2 (我们的): ", striatum_results$method2_degs, " DEGs\n"))
cat(paste0("  文献报道: ~2282 DEGs\n"))

cat("\n========================================\n")
