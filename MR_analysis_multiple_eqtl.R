#!/usr/bin/env Rscript
# 专业级MR分析代码 - 多eQTL数据集分析
# 作者: 生物信息学分析团队
# 日期: 2026-04-24
# 版本: 1.0.0

# 设置随机种子
set.seed(42)

# 加载必要的包
# 首先安装TwoSampleMR
if (!require("TwoSampleMR", quietly = TRUE)) {
  cat("Installing TwoSampleMR...\n")
  # 从GitHub安装最新版本
  if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cran.r-project.org")
    library(devtools, quietly = TRUE)
  }
  devtools::install_github("MRCIEU/TwoSampleMR")
  library(TwoSampleMR, quietly = TRUE)
}

# 安装其他包
if (!require("readxl", quietly = TRUE)) {
  cat("Installing readxl...\n")
  install.packages("readxl", repos = "https://cran.r-project.org")
  library(readxl, quietly = TRUE)
}

if (!require("dplyr", quietly = TRUE)) {
  cat("Installing dplyr...\n")
  install.packages("dplyr", repos = "https://cran.r-project.org")
  library(dplyr, quietly = TRUE)
}

if (!require("ggplot2", quietly = TRUE)) {
  cat("Installing ggplot2...\n")
  install.packages("ggplot2", repos = "https://cran.r-project.org")
  library(ggplot2, quietly = TRUE)
}

if (!require("data.table", quietly = TRUE)) {
  cat("Installing data.table...\n")
  install.packages("data.table", repos = "https://cran.r-project.org")
  library(data.table, quietly = TRUE)
}

# 打印包版本
cat("Package Versions:\n")
cat("  R:", R.version.string, "\n")
cat("  TwoSampleMR: Loaded\n")
cat("  readxl: Loaded\n")
cat("  dplyr: Loaded\n")
cat("  ggplot2: Loaded\n")
cat("  data.table: Loaded\n")

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
  file = "D:/EQTL/finngen_R12_I9_STR",
  description = "FinnGen R12 I9_STR"
)

# 科学假设
hypotheses <- list(
  risk = c("NFKB1", "FDX1", "STAT3", "AGER"),  # 高表达增加风险
  protective = c("GPX4", "HMOX1"),  # 高表达降低风险
  unknown = c("HSPA5", "DLAT", "HIF1A")  # 未知
)

# 核心函数：从Excel文件读取eQTL数据
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

# 核心函数：从GTEx v11文本文件读取eQTL数据
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

# 核心函数：获取eQTL数据
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

# 核心函数：运行MR分析
run_mr_analysis <- function(exp_dat, outcome_dat, gene_symbol, dataset_name) {
  tryCatch({
    # 格式转换 - 转换为data.frame
    exp_dat_df <- as.data.frame(exp_dat)
    
    # 检查exposure数据的列名
    exp_cols <- names(exp_dat_df)
    cat("    Exposure columns:", paste(exp_cols, collapse = ", "), "\n")
    
    # 格式转换 - exposure
    exp_fmt <- format_data(
      exp_dat_df,
      type = "exposure"
    )
    
    # 提取结局数据 - 转换为data.frame
    outcome_dat_df <- as.data.frame(outcome_dat)
    
    # 检查outcome数据的列名
    out_cols <- names(outcome_dat_df)
    cat("    Outcome columns:", paste(out_cols, collapse = ", "), "\n")
    
    # 格式转换 - outcome
    # 尝试不同的列名格式
    tryCatch({
      outcome_fmt <- format_data(
        outcome_dat_df,
        type = "outcome"
      )
    }, error = function(e) {
      cat("    Trying alternative outcome column names...\n")
      # 尝试直接使用列名
    outcome_fmt <- format_data(
      outcome_dat_df,
      type = "outcome",
      snp_col = ifelse("SNP" %in% out_cols, "SNP", out_cols[1]),
      beta_col = ifelse("beta" %in% out_cols, "beta", 
                      ifelse("BETA" %in% out_cols, "BETA", 
                             ifelse("OR" %in% out_cols, "OR", out_cols[2]))),
      se_col = ifelse("se" %in% out_cols, "se", 
                    ifelse("SE" %in% out_cols, "SE", 
                           ifelse("std.error" %in% out_cols, "std.error", out_cols[3]))),
      pval_col = ifelse("pval" %in% out_cols, "pval", 
                      ifelse("P" %in% out_cols, "P", 
                             ifelse("p.value" %in% out_cols, "p.value", out_cols[4])))
    )
    })
    
    # 数据协调
    dat <- harmonise_data(
      exposure_dat = exp_fmt,
      outcome_dat = outcome_fmt
    )
    
    if (nrow(dat) == 0) {
      cat("    WARNING: No harmonized SNPs found\n")
      return(NULL)
    }
    
    cat("    Harmonized", nrow(dat), "SNPs\n")
    
    # Steiger过滤
    tryCatch({
      dat_steiger <- steiger_filtering(dat)
      if (!is.null(dat_steiger) && nrow(dat_steiger) > 0 && "steiger_dir" %in% names(dat_steiger)) {
        n_pass <- sum(dat_steiger$steiger_dir == TRUE, na.rm = TRUE)
        cat("    Steiger pass:", n_pass, "/", nrow(dat_steiger), "\n")
        
        if (n_pass > 0) {
          dat <- dat_steiger[dat_steiger$steiger_dir == TRUE | is.na(dat_steiger$steiger_dir), ]
          cat("    Using", nrow(dat), "SNPs after Steiger filtering\n")
        }
      }
    }, error = function(e) {
      cat("    WARNING: Steiger error:", conditionMessage(e), "\n")
    })
    
    # MR分析
    mr_results <- mr(dat)
    
    # 获取IVW结果
    ivw_result <- mr_results[mr_results$method == "Inverse variance weighted", ]
    
    if (nrow(ivw_result) == 0) {
      cat("    WARNING: No IVW results found\n")
      return(NULL)
    }
    
    # 预期判断
    expectation <- "Unknown"
    if (gene_symbol %in% hypotheses$risk) {
      expectation <- ifelse(ivw_result$or > 1, "✓ Risk", "✗ Risk")
    } else if (gene_symbol %in% hypotheses$protective) {
      expectation <- ifelse(ivw_result$or < 1, "✓ Protective", "✗ Protective")
    }
    
    # 构建结果
    result <- list(
      gene = gene_symbol,
      dataset = dataset_name,
      or = ivw_result$or,
      se = ivw_result$se,
      pval = ivw_result$pval,
      n_snps = nrow(dat),
      expectation = expectation,
      data = dat
    )
    
    return(result)
  }, error = function(e) {
    cat("    ERROR in MR analysis:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# 函数：生成森林图
generate_forest_plot <- function(results_list, output_dir) {
  cat("Generating forest plot...\n")
  
  # 整理数据
  plot_data <- data.frame()
  
  for (dataset_results in results_list) {
    for (result in dataset_results) {
      if (!is.null(result)) {
        plot_data <- rbind(plot_data, data.frame(
          gene = result$gene,
          dataset = result$dataset,
          or = result$or,
          se = result$se,
          pval = result$pval,
          n_snps = result$n_snps,
          expectation = result$expectation
        ))
      }
    }
  }
  
  if (nrow(plot_data) == 0) {
    cat("WARNING: No data for forest plot\n")
    return()
  }
  
  # 计算95% CI
  plot_data$ci_lower <- plot_data$or - 1.96 * plot_data$se
  plot_data$ci_upper <- plot_data$or + 1.96 * plot_data$se
  
  # 排序
  plot_data <- plot_data[order(plot_data$gene, plot_data$dataset), ]
  
  # 创建森林图
  p <- ggplot(plot_data, aes(x = or, y = interaction(gene, dataset)))
  p <- p + geom_point(aes(color = dataset), size = 3)
  p <- p + geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = dataset), height = 0.2)
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
  cat("Step 1: Reading Outcome Data\n")
  cat("----------------------------------------------------------------------\n")
  
  tryCatch({
    outcome_dat <- fread(outcome_config$file)
    cat("  Outcome:", outcome_config$description, "\n")
    cat("  SNPs:", nrow(outcome_dat), "\n")
    cat("  Columns:", paste(names(outcome_dat), collapse = ", "), "\n")
    
    # 检查并调整列名
    if ("SNP" %in% names(outcome_dat)) {
      cat("  SNP column found: SNP\n")
    } else if ("rsid" %in% names(outcome_dat)) {
      cat("  SNP column found: rsid, renaming to SNP\n")
      setnames(outcome_dat, "rsid", "SNP")
    } else if ("variant" %in% names(outcome_dat)) {
      cat("  SNP column found: variant, renaming to SNP\n")
      setnames(outcome_dat, "variant", "SNP")
    } else {
      cat("  WARNING: No SNP column found, using first column as SNP\n")
      setnames(outcome_dat, 1, "SNP")
    }
    
    # 确保必要的列存在
    required_cols <- c("SNP", "beta", "se", "pval", "eaf", "effect_allele", "other_allele")
    for (col in required_cols) {
      if (!col %in% names(outcome_dat)) {
        cat("  WARNING: Column", col, "not found in outcome data\n")
      }
    }
    
  }, error = function(e) {
    cat("  ERROR reading outcome data:", conditionMessage(e), "\n")
    return()
  })
  
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
      
      # 运行MR分析
      result <- run_mr_analysis(exp_dat, outcome_dat, gene, dataset$description)
      
      if (!is.null(result)) {
        cat("  MR Result: OR =", sprintf("%.4f", result$or), 
            "[95%% CI: %.4f - %.4f]", result$or - 1.96 * result$se, result$or + 1.96 * result$se, 
            "P =", sprintf("%.6f", result$pval), "\n")
        cat("  Expectation:", result$expectation, "\n")
        
        dataset_results[[gene]] <- result
        
        # 保存单基因结果
        write.csv(result$data, 
                  file.path(output_dir, paste0("MR_", gene, "_", dataset$name, ".csv")), 
                  row.names = FALSE)
      }
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
  generate_forest_plot(all_results, output_dir)
  
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
