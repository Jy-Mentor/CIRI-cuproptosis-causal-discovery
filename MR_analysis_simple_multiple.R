#!/usr/bin/env Rscript
# 专业级MR分析代码 - 多eQTL数据集分析
# 作者: 生物信息学分析团队
# 日期: 2026-04-24
# 版本: 1.0.0

# 设置随机种子
set.seed(42)

# 加载必要的包
suppressPackageStartupMessages({
  if (!require("TwoSampleMR", quietly = TRUE)) {
    cat("Installing TwoSampleMR...\n")
    if (!require("devtools", quietly = TRUE)) {
      install.packages("devtools", repos = "https://cran.r-project.org")
      library(devtools, quietly = TRUE)
    }
    devtools::install_github("MRCIEU/TwoSampleMR")
    library(TwoSampleMR, quietly = TRUE)
  }
  
  if (!require("readxl", quietly = TRUE)) {
    install.packages("readxl", repos = "https://cran.r-project.org")
    library(readxl, quietly = TRUE)
  }
  
  if (!require("data.table", quietly = TRUE)) {
    install.packages("data.table", repos = "https://cran.r-project.org")
    library(data.table, quietly = TRUE)
  }
  
  if (!require("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2", repos = "https://cran.r-project.org")
    library(ggplot2, quietly = TRUE)
  }
})

# 打印包版本
cat("Package Versions:\n")
cat("  R:", R.version.string, "\n")
cat("  TwoSampleMR: Loaded\n")
cat("  readxl: Loaded\n")
cat("  data.table: Loaded\n")
cat("  ggplot2: Loaded\n")

# 基因列表
hub_genes <- c(
  "NFKB1", "FDX1", "STAT3",  # P0层核心基因
  "HIF1A", "HMOX1", "GPX4", "HSPA5", "AGER", "DLAT"  # P1层支持基因
)

# 数据集配置
datasets <- list(
  list(
    name = "eQTLgen_p5e-06",
    type = "excel",
    path = "D:/EQTL/clump/eQTLgen_allgene_p_5e-06_kb_1000_r2_0.01.xlsx",
    description = "eQTLgen (p=5e-06)"
  ),
  list(
    name = "eQTLgen_p5e-08",
    type = "excel",
    path = "D:/EQTL/clump/eQTLgen_allgene_p_5e-08_kb_1000_r2_0.01.xlsx",
    description = "eQTLgen (p=5e-08)"
  ),
  list(
    name = "eQTLgen_p1e-05",
    type = "excel",
    path = "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.01.xlsx",
    description = "eQTLgen (p=1e-05)"
  ),
  list(
    name = "GTEx_v11_whole_blood",
    type = "text",
    path = "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eGenes.txt",
    description = "GTEx v11 Whole Blood"
  )
)

# 结局数据配置
outcome_config <- list(
  name = "Ischemic Stroke",
  file = "D:/EQTL/eqtlgen_ieu_outcome.csv",
  description = "FinnGen R12 I9_STR"
)

# 科学假设
hypotheses <- list(
  risk = c("NFKB1", "FDX1", "STAT3", "AGER"),  # 高表达增加风险
  protective = c("GPX4", "HMOX1"),  # 高表达降低风险
  unknown = c("HSPA5", "DLAT", "HIF1A")  # 未知
)

# 函数：从Excel文件读取eQTL数据
read_eqtl_excel <- function(file_path, gene_symbol) {
  tryCatch({
    eqtl_data <- read_excel(file_path)
    
    # 检查必要列
    required_cols <- c("gene", "SNP", "beta.exposure", "se.exposure", "pval.exposure", 
                      "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")
    
    missing_cols <- setdiff(required_cols, names(eqtl_data))
    if (length(missing_cols) > 0) {
      cat("    ERROR: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
      return(NULL)
    }
    
    # 筛选基因
    gene_data <- eqtl_data[eqtl_data$gene == gene_symbol, ]
    
    if (nrow(gene_data) == 0) {
      cat("    WARNING: No eQTLs found for", gene_symbol, "\n")
      return(NULL)
    }
    
    cat("    Found", nrow(gene_data), "eQTLs for", gene_symbol, "\n")
    
    # 转换格式
    exp_dat <- data.frame(
      SNP = gene_data$SNP,
      beta.exposure = gene_data$beta.exposure,
      se.exposure = gene_data$se.exposure,
      pval.exposure = gene_data$pval.exposure,
      eaf.exposure = gene_data$eaf.exposure,
      effect_allele.exposure = gene_data$effect_allele.exposure,
      other_allele.exposure = gene_data$other_allele.exposure,
      exposure = paste0(gene_symbol, " (eQTLgen)"),
      stringsAsFactors = FALSE
    )
    
    # 过滤无rsID的SNP
    exp_dat <- exp_dat[exp_dat$SNP != "", ]
    
    if (nrow(exp_dat) == 0) {
      cat("    WARNING: No SNPs with rsID found\n")
      return(NULL)
    }
    
    return(exp_dat)
  }, error = function(e) {
    cat("    ERROR reading Excel:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# 函数：从GTEx v11文本文件读取eQTL数据
read_eqtl_gtex <- function(file_path, gene_symbol) {
  tryCatch({
    eqtl_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # 筛选基因
    gene_data <- eqtl_data[eqtl_data$gene_name == gene_symbol, ]
    
    if (nrow(gene_data) == 0) {
      cat("    WARNING: No eQTLs found for", gene_symbol, "\n")
      return(NULL)
    }
    
    cat("    Found", nrow(gene_data), "eQTLs for", gene_symbol, "\n")
    
    # 转换格式
    exp_dat <- data.frame(
      SNP = gene_data$rs_id_dbSNP157_GRCh38p14,
      beta.exposure = gene_data$slope,
      se.exposure = gene_data$slope_se,
      pval.exposure = gene_data$pval_nominal,
      eaf.exposure = gene_data$af,
      effect_allele.exposure = gene_data$alt,
      other_allele.exposure = gene_data$ref,
      exposure = paste0(gene_symbol, " (GTEx v11 Whole Blood)"),
      stringsAsFactors = FALSE
    )
    
    # 过滤无rsID的SNP
    exp_dat <- exp_dat[exp_dat$SNP != "", ]
    
    if (nrow(exp_dat) == 0) {
      cat("    WARNING: No SNPs with rsID found\n")
      return(NULL)
    }
    
    return(exp_dat)
  }, error = function(e) {
    cat("    ERROR reading GTEx file:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# 函数：获取eQTL数据
get_eqtl_data <- function(dataset, gene_symbol) {
  cat("  Processing", gene_symbol, "...\n")
  
  if (dataset$type == "excel") {
    return(read_eqtl_excel(dataset$path, gene_symbol))
  } else if (dataset$type == "text") {
    return(read_eqtl_gtex(dataset$path, gene_symbol))
  } else {
    cat("    ERROR: Unknown dataset type\n")
    return(NULL)
  }
}

# 主分析函数
main_analysis <- function() {
  cat("======================================================================\n")
  cat("Professional MR Analysis: Multiple eQTL Datasets\n")
  cat("======================================================================\n\n")
  
  # 创建输出目录
  output_dir <- file.path("D:/EQTL", paste0("MR_Multiple_Datasets_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(output_dir, recursive = TRUE)
  
  cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Output Directory:", output_dir, "\n\n")
  
  # 读取结局数据
  cat("----------------------------------------------------------------------\n")
  cat("Step 1: Setting Up Outcome Data\n")
  cat("----------------------------------------------------------------------\n")
  
  cat("  Outcome:", outcome_config$description, "\n")
  cat("  File:", outcome_config$file, "\n")
  
  # 检查文件是否存在
  if (!file.exists(outcome_config$file)) {
    cat("  ERROR: Outcome data file not found\n")
    return()
  }
  
  cat("  Outcome data file ready for chunked processing\n")
  
  # 运行分析
  all_results <- list()
  
  for (dataset in datasets) {
    cat("\n======================================================================\n")
    cat("Analyzing Dataset:", dataset$description, "\n")
    cat("======================================================================\n")
    
    dataset_results <- list()
    
    for (gene in hub_genes) {
      cat("\n------------------------------------------------------------\n")
      cat("Gene:", gene, "\n")
      cat("------------------------------------------------------------\n")
      
      # 获取eQTL数据
      exp_dat <- get_eqtl_data(dataset, gene)
      
      if (is.null(exp_dat)) {
        cat("  Skipping MR analysis due to missing eQTL data\n")
        next
      }
      
      # 格式转换
      exp_fmt <- format_data(
        exp_dat,
        type = "exposure",
        snp_col = "SNP",
        beta_col = "beta.exposure",
        se_col = "se.exposure",
        pval_col = "pval.exposure",
        eaf_col = "eaf.exposure",
        effect_allele_col = "effect_allele.exposure",
        other_allele_col = "other_allele.exposure"
      )
      
      # 提取结局数据
      snps <- exp_fmt$SNP
      cat("  Extracting outcome data for", length(snps), "SNPs\n")
      
      # 动态读取结局数据中匹配的SNP（分块读取）
      outcome_subset <- data.table()
      
      # 分块读取文件
      chunk_size <- 1000000  # 每块100万行
      con <- file(outcome_config$file, "r")
      
      # 读取表头
      header <- readLines(con, n = 1)
      header_cols <- strsplit(header, ",")[[1]]
      
      # 找到需要的列
      cols_needed <- c("rsids", "beta", "sebeta", "pval", "af_alt", "alt", "ref")
      
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
      
      cat("  Matched", nrow(outcome_subset), "SNPs\n")
      
      # 转换结局数据格式
      outcome_fmt <- data.frame(
        SNP = outcome_subset$rsids,
        beta.outcome = outcome_subset$beta,
        se.outcome = outcome_subset$sebeta,
        pval.outcome = outcome_subset$pval,
        eaf.outcome = outcome_subset$af_alt,
        effect_allele.outcome = outcome_subset$alt,
        other_allele.outcome = outcome_subset$ref,
        outcome = outcome_config$name,
        id.outcome = outcome_config$name,
        stringsAsFactors = FALSE
      )
      
      # 数据协调
      dat <- harmonise_data(
        exposure_dat = exp_fmt,
        outcome_dat = outcome_fmt
      )
      
      if (nrow(dat) == 0) {
        cat("  WARNING: No harmonized SNPs found\n")
        next
      }
      
      cat("  Harmonized", nrow(dat), "SNPs\n")
      
      # 初始化结果变量
      ivw_result <- NULL
      
      # 尝试使用IVW方法
      tryCatch({
        ivw_result <- mr_ivw(dat)
      }, error = function(e) {
        ivw_result <- NULL
      })
      
      # 如果IVW失败，尝试Wald ratio（单SNP）
      if (is.null(ivw_result) || (nrow(ivw_result) == 0)) {
        if (nrow(dat) == 1) {
          cat("  Trying Wald ratio for single SNP\n")
          tryCatch({
            ivw_result <- mr_wald_ratio(dat)
            if (!is.null(ivw_result) && nrow(ivw_result) > 0) {
              cat("  Wald ratio analysis successful\n")
            }
          }, error = function(e) {
            cat("  Error in Wald ratio: ", e$message, "\n")
            ivw_result <- NULL
          })
        } else {
          cat("  Skipping Wald ratio: ", nrow(dat), " SNPs available\n")
        }
      }
      
      # 检查是否有有效结果
      if (is.null(ivw_result) || nrow(ivw_result) == 0) {
        cat("  WARNING: No MR results found\n")
        next
      }
      
      # 检查结果是否有效
      if (!all(c("or", "se", "pval") %in% names(ivw_result)) || 
          any(is.na(ivw_result$or), is.na(ivw_result$se), is.na(ivw_result$pval))) {
        cat("  WARNING: Invalid MR results\n")
        next
      }
      
      # 预期判断
      expectation <- "Unknown"
      or_value <- as.numeric(ivw_result$or)
      
      if (gene %in% hypotheses$risk) {
        if (or_value > 1) {
          expectation <- "✓ Risk"
        } else {
          expectation <- "✗ Risk"
        }
      } else if (gene %in% hypotheses$protective) {
        if (or_value < 1) {
          expectation <- "✓ Protective"
        } else {
          expectation <- "✗ Protective"
        }
      }
      
      # 构建结果
      result <- list(
        gene = gene,
        dataset = dataset$description,
        or = as.numeric(ivw_result$or),
        se = as.numeric(ivw_result$se),
        pval = as.numeric(ivw_result$pval),
        n_snps = nrow(dat),
        expectation = expectation,
        data = dat
      )
      
      dataset_results[[gene]] <- result
      
      # 保存单基因结果
      write.csv(result$data, 
                file.path(output_dir, paste0("MR_", gene, "_", dataset$name, ".csv")), 
                row.names = FALSE)
      
      # 打印结果
      cat("  MR Result: OR =", sprintf("%.4f", result$or), 
          " [95%% CI: %.4f - %.4f]", result$or - 1.96 * result$se, result$or + 1.96 * result$se, 
          " P =", sprintf("%.6f", result$pval), "\n")
      cat("  Expectation:", result$expectation, "\n")
    }
    
    all_results[[dataset$name]] <- dataset_results
  }
  
  # 生成汇总结果
  cat("\n======================================================================\n")
  cat("Generating Summary Results\n")
  cat("======================================================================\n")
  
  summary_data <- data.frame()
  
  for (dataset_name in names(all_results)) {
    dataset_results <- all_results[[dataset_name]]
    
    for (gene in names(dataset_results)) {
      result <- dataset_results[[gene]]
      
      summary_data <- rbind(summary_data, data.frame(
        Gene = result$gene,
        Dataset = result$dataset,
        OR = result$or,
        SE = result$se,
        P_Value = result$pval,
        N_SNPs = result$n_snps,
        Expectation = result$expectation,
        Significant = ifelse(result$pval < 0.05, "Yes", "No")
      ))
    }
  }
  
  # 保存汇总结果
  write.csv(summary_data, 
            file.path(output_dir, "MR_Summary_All_Datasets.csv"), 
            row.names = FALSE)
  
  # 生成森林图
  if (nrow(summary_data) > 0) {
    cat("Generating forest plot...\n")
    
    # 计算95% CI
    summary_data$ci_lower <- summary_data$OR - 1.96 * summary_data$SE
    summary_data$ci_upper <- summary_data$OR + 1.96 * summary_data$SE
    
    # 排序
    summary_data <- summary_data[order(summary_data$Gene, summary_data$Dataset), ]
    
    # 创建森林图
    p <- ggplot(summary_data, aes(x = OR, y = interaction(Gene, Dataset)))
    p <- p + geom_point(aes(color = Dataset), size = 3)
    p <- p + geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = Dataset), height = 0.2)
    p <- p + geom_vline(xintercept = 1, linetype = "dashed", color = "gray50")
    p <- p + scale_x_log10(breaks = c(0.8, 1, 1.2, 1.4))
    p <- p + labs(
      title = "MR Analysis Results: Hub Genes vs Ischemic Stroke",
      x = "Odds Ratio (log scale)",
      y = "Gene + Dataset",
      color = "Dataset"
    )
    p <- p + theme_minimal()
    p <- p + theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      legend.position = "top"
    )
    
    # 保存
    ggsave(file.path(output_dir, "MR_Forest_All_Datasets.png"), p, width = 12, height = 10, dpi = 300)
    ggsave(file.path(output_dir, "MR_Forest_All_Datasets.pdf"), p, width = 12, height = 10)
    
    cat("Forest plot saved\n")
  } else {
    cat("WARNING: No data for forest plot\n")
  }
  
  # 打印显著结果
  cat("\n======================================================================\n")
  cat("Significant Results (P < 0.05)\n")
  cat("======================================================================\n")
  
  significant_results <- summary_data[summary_data$Significant == "Yes", ]
  
  if (nrow(significant_results) > 0) {
    print(significant_results[, c("Gene", "Dataset", "OR", "P_Value", "Expectation")])
  } else {
    cat("No significant results found\n")
  }
  
  cat("\n======================================================================\n")
  cat("Analysis Complete!\n")
  cat("======================================================================\n")
  cat("Output Directory:", output_dir, "\n")
  cat("Files Generated:\n")
  cat("  - MR_Summary_All_Datasets.csv\n")
  cat("  - MR_*_*.csv (per gene per dataset)\n")
  cat("  - MR_Forest_All_Datasets.png/pdf\n")
  cat("======================================================================\n")
}

# 运行分析
main_analysis()
