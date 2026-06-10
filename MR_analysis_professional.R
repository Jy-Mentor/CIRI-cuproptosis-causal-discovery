# 专业级MR分析代码：多个eQTL数据集
# 分析9个hub基因与缺血性脑卒中的关系

library(TwoSampleMR)
library(readxl)
library(data.table)
library(ggplot2)

# 1. 设置参数
args <- commandArgs(trailingOnly = TRUE)

# 数据集配置
datasets <- list(
  list(
    path = "D:\\EQTL\\clump\\eQTLgen_allgene_p_5e-06_kb_1000_r2_0.01.xlsx",
    description = "eQTLgen (p=5e-06)",
    type = "excel"
  ),
  list(
    path = "D:\\EQTL\\clump\\eQTLgen_allgene_p_5e-08_kb_1000_r2_0.01.xlsx",
    description = "eQTLgen (p=5e-08)",
    type = "excel"
  ),
  list(
    path = "D:\\EQTL\\clump\\eQTLgen_allgene_p_1e-05_kb_1000_r2_0.01.xlsx",
    description = "eQTLgen (p=1e-05)",
    type = "excel"
  ),
  list(
    path = "C:\\Users\\Jy-Mentor-7\\Desktop\\生物信息学\\ETQL\\Whole_Blood.v11.eGenes.txt",
    description = "GTEx v11 Whole Blood",
    type = "gtex"
  )
)

# 分析的基因
genes <- c("NFKB1", "FDX1", "STAT3", "HIF1A", "HMOX1", "GPX4", "HSPA5", "AGER", "DLAT")

# 假设：哪些基因是风险基因，哪些是保护基因
hypotheses <- list(
  risk = c("NFKB1", "STAT3", "HIF1A"),  # 风险基因
  protective = c("FDX1", "HMOX1", "GPX4")  # 保护基因
)

# 2. 创建输出目录
output_dir <- paste0("D:/EQTL/MR_Multiple_Datasets_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3. 准备结局数据（缺血性脑卒中）
cat("======================================================================\n")
cat("Professional MR Analysis: Multiple eQTL Datasets\n")
cat("======================================================================\n\n")
cat(paste0("Analysis Date: ", Sys.time(), "\n"))
cat(paste0("Output Directory: ", output_dir, "\n\n"))

cat("----------------------------------------------------------------------\n")
cat("Step 1: Setting Up Outcome Data\n")
cat("----------------------------------------------------------------------\n")

# 设置CRAN镜像
options(repos = c(CRAN = "https://cran.r-project.org"))

# 结局数据文件路径
outcome_path <- "D:/EQTL/eqtlgen_ieu_outcome.csv"

# 检查文件是否存在
if (!file.exists(outcome_path)) {
  cat("  ERROR: Outcome data file not found\n")
  stop("Outcome data file not found")
}

cat(paste0("  Outcome: FinnGen R12 I9_STR \n"))
cat(paste0("  Outcome file: ", outcome_path, "\n\n"))

# 4. 读取eQTL数据的函数
read_eqtl_excel <- function(file_path, gene) {
  # 读取Excel文件
  eqtl_data <- read_excel(file_path)
  
  # 打印列名以调试
  cat(paste0("  Excel columns: ", paste(names(eqtl_data), collapse = ", "), "\n"))
  
  # 确定基因列名
  gene_col <- if ("gene" %in% names(eqtl_data)) "gene" else if ("Gene" %in% names(eqtl_data)) "Gene" else NULL
  
  if (is.null(gene_col)) {
    cat("  WARNING: No gene column found\n")
    return(NULL)
  }
  
  # 过滤特定基因的eQTL
  gene_eqtl <- eqtl_data[eqtl_data[[gene_col]] == gene, ]
  
  if (nrow(gene_eqtl) == 0) {
    return(NULL)
  }
  
  # 确定各列名（优先使用带.exposure后缀的列）
  snp_col <- if ("SNP" %in% names(gene_eqtl)) "SNP" else if ("rsid" %in% names(gene_eqtl)) "rsid" else NULL
  beta_col <- if ("beta.exposure" %in% names(gene_eqtl)) "beta.exposure" else if ("beta" %in% names(gene_eqtl)) "beta" else if ("Beta" %in% names(gene_eqtl)) "Beta" else NULL
  se_col <- if ("se.exposure" %in% names(gene_eqtl)) "se.exposure" else if ("se" %in% names(gene_eqtl)) "se" else if ("SE" %in% names(gene_eqtl)) "SE" else NULL
  pval_col <- if ("pval.exposure" %in% names(gene_eqtl)) "pval.exposure" else if ("pval" %in% names(gene_eqtl)) "pval" else if ("P" %in% names(gene_eqtl)) "P" else if ("p.value" %in% names(gene_eqtl)) "p.value" else NULL
  eaf_col <- if ("eaf.exposure" %in% names(gene_eqtl)) "eaf.exposure" else if ("eaf" %in% names(gene_eqtl)) "eaf" else if ("EAF" %in% names(gene_eqtl)) "EAF" else if ("maf" %in% names(gene_eqtl)) "maf" else NULL
  effect_allele_col <- if ("effect_allele.exposure" %in% names(gene_eqtl)) "effect_allele.exposure" else if ("effect_allele" %in% names(gene_eqtl)) "effect_allele" else if ("EA" %in% names(gene_eqtl)) "EA" else if ("A1" %in% names(gene_eqtl)) "A1" else NULL
  other_allele_col <- if ("other_allele.exposure" %in% names(gene_eqtl)) "other_allele.exposure" else if ("other_allele" %in% names(gene_eqtl)) "other_allele" else if ("OA" %in% names(gene_eqtl)) "OA" else if ("A2" %in% names(gene_eqtl)) "A2" else NULL
  
  # 检查必要列是否存在
  if (any(sapply(list(snp_col, beta_col, se_col, pval_col), is.null))) {
    cat("  WARNING: Missing required columns\n")
    return(NULL)
  }
  
  # 直接返回过滤后的eQTL数据（已经是TwoSampleMR格式）
  return(gene_eqtl)
}

read_eqtl_gtex <- function(file_path, gene) {
  # 读取GTEx eGenes文件
  gtex_data <- fread(file_path, header = TRUE)
  
  # 过滤特定基因
  gene_data <- gtex_data[gtex_data$gene_id == gene, ]
  
  if (nrow(gene_data) == 0) {
    return(NULL)
  }
  
  # 对于GTEx，我们需要从signif_pairs文件中获取SNP信息
  signif_pairs_path <- "C:\\Users\\Jy-Mentor-7\\Desktop\\生物信息学\\ETQL\\Whole_Blood.v11.eQTLs.signif_pairs.parquet"
  
  # 检查文件是否存在
  if (!file.exists(signif_pairs_path)) {
    cat("  WARNING: GTEx signif_pairs file not found\n")
    return(NULL)
  }
  
  # 读取parquet文件（需要arrow包）
  if (!requireNamespace("arrow", quietly = TRUE)) {
    install.packages("arrow", quiet = TRUE)
  }
  
  signif_pairs <- arrow::read_parquet(signif_pairs_path)
  
  # 过滤特定基因的SNP
  gene_snps <- signif_pairs[signif_pairs$gene_id == gene, ]
  
  if (nrow(gene_snps) == 0) {
    return(NULL)
  }
  
  # 转换为TwoSampleMR格式
  exposure_data <- data.frame(
    SNP = gene_snps$variant_id,
    beta.exposure = gene_snps$slope,
    se.exposure = gene_snps$slope_se,
    pval.exposure = gene_snps$pval_nominal,
    eaf.exposure = gene_snps$maf,
    effect_allele.exposure = gene_snps$alt,
    other_allele.exposure = gene_snps$ref,
    stringsAsFactors = FALSE
  )
  
  return(exposure_data)
}

# 5. 主分析函数
main_analysis <- function() {
  all_results <- list()
  
  # 分析每个数据集
  for (dataset in datasets) {
    cat("======================================================================\n")
    cat(paste0("Analyzing Dataset: ", dataset$description, "\n"))
    cat("======================================================================\n\n")
    
    dataset_results <- list()
    
    # 分析每个基因
    for (gene in genes) {
      cat("------------------------------------------------------------\n")
      cat(paste0("Gene: ", gene, "\n"))
      cat("------------------------------------------------------------\n")
      cat(paste0("  Processing ", gene, " ...\n"))
      
      # 读取eQTL数据
      if (dataset$type == "excel") {
        exposure_data <- read_eqtl_excel(dataset$path, gene)
      } else if (dataset$type == "gtex") {
        exposure_data <- read_eqtl_gtex(dataset$path, gene)
      }
      
      if (is.null(exposure_data) || nrow(exposure_data) == 0) {
        cat(paste0("  WARNING: No eQTLs found for ", gene, "\n"))
        next
      }
      
      cat(paste0("    Found ", nrow(exposure_data), " eQTLs for ", gene, "\n"))
      
      # 直接使用暴露数据（已经是TwoSampleMR格式）
      exp_fmt <- exposure_data
      
      # 提取结局数据
      snps <- exp_fmt$SNP
      cat(paste0("  Extracting outcome data for ", length(snps), " SNPs\n"))
      
      # 动态读取结局数据中匹配的SNP
      # 使用分块读取来处理大文件
      outcome_subset <- data.table()
      
      # 分块读取文件
      chunk_size <- 1000000  # 每块100万行
      con <- file(outcome_path, "r")
      
      # 读取表头
      header <- readLines(con, n = 1)
      header_cols <- strsplit(header, ",")[[1]]
      
      # 找到需要的列的索引
      cols_needed <- c("rsids", "beta", "sebeta", "pval", "af_alt", "alt", "ref")
      col_indices <- match(cols_needed, header_cols)
      
      if (any(is.na(col_indices))) {
        cat("  WARNING: Some required columns not found in outcome data\n")
        close(con)
        next
      }
      
      # 逐块读取并过滤
      while (TRUE) {
        chunk <- tryCatch({
          fread(con, nrows = chunk_size, select = cols_needed, header = FALSE)
        }, error = function(e) {
          NULL
        })
        
        if (is.null(chunk) || nrow(chunk) == 0) {
          break
        }
        
        # 过滤匹配的SNP
        chunk_match <- chunk[chunk$rsids %in% snps, ]
        if (nrow(chunk_match) > 0) {
          outcome_subset <- rbind(outcome_subset, chunk_match)
        }
      }
      
      close(con)
      
      if (nrow(outcome_subset) == 0) {
        cat("  WARNING: No matching SNPs in outcome data\n")
        next
      }
      
      cat(paste0("  Matched ", nrow(outcome_subset), " SNPs\n"))
      
      # 格式化结局数据
      outcome_fmt <- data.frame(
        SNP = outcome_subset$rsids,
        beta.outcome = outcome_subset$beta,
        se.outcome = outcome_subset$sebeta,
        pval.outcome = outcome_subset$pval,
        eaf.outcome = outcome_subset$af_alt,
        effect_allele.outcome = outcome_subset$alt,
        other_allele.outcome = outcome_subset$ref,
        outcome = "Ischemic Stroke",
        id.outcome = "Ischemic Stroke",
        stringsAsFactors = FALSE
      )
      
      #  Harmonize数据
      dat <- harmonise_data(exp_fmt, outcome_fmt)
      
      if (nrow(dat) == 0) {
        cat("  WARNING: No harmonized SNPs\n")
        next
      }
      
      cat(paste0("  Harmonized ", nrow(dat), " SNPs\n"))
      
      # 尝试不同的MR方法
      mr_results <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio"))
      
      # 优先选择IVW结果
      ivw_result <- mr_results[mr_results$method == "Inverse variance weighted", ]
      
      # 如果没有IVW结果，尝试Wald ratio
      if (nrow(ivw_result) == 0) {
        ivw_result <- mr_results[mr_results$method == "Wald ratio", ]
      }
      
      if (nrow(ivw_result) == 0) {
        cat("  WARNING: No MR results found\n")
        next
      }
      
      # 使用第一个结果
      ivw_result <- ivw_result[1, ]
      
      # 检查结果是否有效
      if (any(is.na(ivw_result$or), is.na(ivw_result$se), is.na(ivw_result$pval))) {
        cat("  WARNING: Invalid MR results\n")
        next
      }
      
      # 检查ivw_result的结构
      cat(paste0("  ivw_result structure: ", paste(names(ivw_result), collapse = ", "), "\n"))
      cat(paste0("  or length: ", length(ivw_result$or), "\n"))
      cat(paste0("  se length: ", length(ivw_result$se), "\n"))
      cat(paste0("  pval length: ", length(ivw_result$pval), "\n"))
      
      # 确保结果是标量
      or <- as.numeric(ivw_result$or[1])
      se <- as.numeric(ivw_result$se[1])
      pval <- as.numeric(ivw_result$pval[1])
      
      # 检查计算是否有效
      if (any(is.na(c(or, se, pval)))) {
        cat("  WARNING: Invalid MR results values\n")
        next
      }
      
      # 计算95% CI
      ci_lower <- exp(log(or) - 1.96 * se)
      ci_upper <- exp(log(or) + 1.96 * se)
      
      # 预期判断
      expectation <- "Unknown"
      if (gene %in% hypotheses$risk) {
        if (or > 1) {
          expectation <- "✓ Risk"
        } else {
          expectation <- "✗ Risk"
        }
      } else if (gene %in% hypotheses$protective) {
        if (or < 1) {
          expectation <- "✓ Protective"
        } else {
          expectation <- "✗ Protective"
        }
      }
      
      # 打印结果
      cat(paste0("  MR Result: OR = ", sprintf("%.4f", or), " [95% CI: ", sprintf("%.4f", ci_lower), " - ", sprintf("%.4f", ci_upper), "]    P = ", sprintf("%.6f", pval), "\n"))
      cat(paste0("  Expectation: ", expectation, "\n\n"))
      
      # 检查所有字段是否为标量
      if (all(c(length(gene) == 1, length(dataset$description) == 1, 
               length(or) == 1, length(se) == 1, length(ci_lower) == 1, 
               length(ci_upper) == 1, length(pval) == 1, nrow(dat) > 0, 
               length(expectation) == 1))) {
        # 构建结果
        result <- list(
          gene = as.character(gene),
          dataset = as.character(dataset$description),
          or = as.numeric(or),
          se = as.numeric(se),
          ci_lower = as.numeric(ci_lower),
          ci_upper = as.numeric(ci_upper),
          pval = as.numeric(pval),
          n_snps = as.integer(nrow(dat)),
          expectation = as.character(expectation),
          data = dat
        )
        
        dataset_results[[gene]] <- result
      } else {
        cat("  WARNING: Skipping result due to non-scalar values\n")
      }
      
      # 保存单基因结果
      gene_dir <- file.path(output_dir, dataset$description)
      dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)
      
      result_file <- file.path(gene_dir, paste0(gene, "_results.txt"))
      sink(result_file)
      cat(paste0("Gene: ", gene, "\n"))
      cat(paste0("Dataset: ", dataset$description, "\n"))
      cat(paste0("OR: ", sprintf("%.4f", or), "\n"))
      cat(paste0("SE: ", sprintf("%.4f", se), "\n"))
      cat(paste0("95% CI: ", sprintf("%.4f", ci_lower), " - ", sprintf("%.4f", ci_upper), "\n"))
      cat(paste0("P-value: ", sprintf("%.6f", as.numeric(ivw_result$pval)), "\n"))
      cat(paste0("Number of SNPs: ", nrow(dat), "\n"))
      cat(paste0("Expectation: ", expectation, "\n"))
      sink()
      
      # 保存SNP信息
      snp_file <- file.path(gene_dir, paste0(gene, "_snps.txt"))
      write.table(dat, snp_file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    
    all_results[[dataset$description]] <- dataset_results
  }
  
  # 生成汇总结果
  cat("======================================================================\n")
  cat("Generating Summary Results\n")
  cat("======================================================================\n")
  
  # 创建汇总数据框
  summary_data <- data.frame()
  
  for (dataset_name in names(all_results)) {
    dataset_results <- all_results[[dataset_name]]
    for (gene in names(dataset_results)) {
      result <- dataset_results[[gene]]
      # 检查结果是否有效
      if (!is.null(result) && all(c("gene", "dataset", "or", "se", "ci_lower", "ci_upper", "pval", "n_snps", "expectation") %in% names(result))) {
        row <- data.frame(
          Gene = as.character(result$gene),
          Dataset = as.character(result$dataset),
          OR = as.numeric(result$or),
          SE = as.numeric(result$se),
          CI_Lower = as.numeric(result$ci_lower),
          CI_Upper = as.numeric(result$ci_upper),
          P_Value = as.numeric(result$pval),
          N_SNPs = as.integer(result$n_snps),
          Expectation = as.character(result$expectation)
        )
        summary_data <- rbind(summary_data, row)
      }
    }
  }
  
  # 保存汇总结果
  summary_file <- file.path(output_dir, "summary_results.txt")
  write.table(summary_data, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # 保存汇总结果为Excel
  if (!requireNamespace("writexl", quietly = TRUE)) {
    install.packages("writexl", quiet = TRUE)
  }
  writexl::write_xlsx(summary_data, file.path(output_dir, "summary_results.xlsx"))
  
  cat("\nSummary results saved to:")
  cat(paste0("\n  ", summary_file))
  cat(paste0("\n  ", file.path(output_dir, "summary_results.xlsx")))
  
  # 生成森林图
  if (nrow(summary_data) > 0) {
    cat("\n\nGenerating Forest Plot...\n")
    
    # 准备森林图数据
    forest_data <- summary_data
    forest_data$label <- paste0(forest_data$Gene, " (", forest_data$Dataset, ")")
    
    # 按基因和数据集排序
    forest_data <- forest_data[order(forest_data$Gene, forest_data$Dataset), ]
    
    # 创建森林图
    p <- ggplot(forest_data, aes(x = label, y = OR, ymin = CI_Lower, ymax = CI_Upper)) +
      geom_pointrange() +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      coord_flip() +
      theme_bw() +
      labs(
        title = "MR Analysis Results: Genes vs Ischemic Stroke",
        x = "Gene (Dataset)",
        y = "Odds Ratio (95% CI)"
      ) +
      scale_y_continuous(trans = "log", breaks = c(0.1, 0.5, 1, 2, 5, 10))
    
    # 保存森林图
    forest_file <- file.path(output_dir, "forest_plot.png")
    ggsave(forest_file, p, width = 12, height = 8, dpi = 300)
    
    cat(paste0("Forest plot saved to: ", forest_file, "\n"))
  }
  
  cat("\n======================================================================\n")
  cat("Analysis Complete!\n")
  cat("======================================================================\n")
  
  return(all_results)
}

# 6. 运行分析
main_analysis()
