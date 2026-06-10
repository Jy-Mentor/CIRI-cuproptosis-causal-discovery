#!/usr/bin/env Rscript

# 双样本孟德尔随机化分析：铜死亡核心基因与缺血性脑卒中
# 兼容版：ieugwasr 1.1.0 with old API endpoint

options(stringsAsFactors = FALSE, encoding = "UTF-8")
set.seed(2026)

# ========== 关键修复：使用旧版API端点 ==========
# 在加载ieugwasr之前设置环境变量
Sys.setenv(IEUGWASR_API = "https://api.opengwas.io")
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3NTMxNTM2NywiZXhwIjoxNzc2NTI0OTY3fQ.j7i0kagMHKPkoJBh8w3Xzyz6BEkp7QCd0mjdT1tC6qGIfFhs1hO1Z6awmzGB17lpgiMTyUwZ2M6sJSo0PNlvNb2UIGKAAqqDSrKEHqmlA-zLvrws0JQ9mcec56cKxdSF3Qn-tP_deKnmMaKOq-eYA7bgMgkEcQRo_7LRAOQDB1-MLCKP8ffo4oxAPUzLnGEyZPtQIpmhfxOM-nvSAhZCpcZqgrZfS1QrpLUjgvvDlAuVjpuKcgDSmhnkmXImPJeCvxgwMLZCte0MBIl7ATvwaEve_2yKkqH8xYlmmO0HFVtua_2HWCoT4JujHpbP7jKamglM-Mc7oy-8Douat7AOGg")
cat("Set IEUGWASR_API to: https://api.opengwas.io\n")
cat("Set OPENGWAS_JWT via environment variable\n")

# 设置CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# ========== 包加载 ==========
tryCatch({
  if (!require("remotes", character.only = TRUE)) {
    install.packages("remotes")
    library(remotes)
  }

  if (!require("TwoSampleMR", character.only = TRUE)) {
    remotes::install_github("MRCIEU/TwoSampleMR")
    library(TwoSampleMR)
  }

  if (!require("ieugwasr", character.only = TRUE)) {
    install.packages("ieugwasr")
    library(ieugwasr)
  }

  if (!require("dplyr", character.only = TRUE)) {
    install.packages("dplyr")
    library(dplyr)
  }

  if (!require("ggplot2", character.only = TRUE)) {
    install.packages("ggplot2")
    library(ggplot2)
  }

  if (!require("MRPRESSO", character.only = TRUE)) {
    remotes::install_github("rondolab/MRPRESSO")
    library(MRPRESSO)
  }

  # 尝试设置token
  cat("Setting API token...\n")
  token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3NTMxNTM2NywiZXhwIjoxNzc2NTI0OTY3fQ.j7i0kagMHKPkoJBh8w3Xzyz6BEkp7QCd0mjdT1tC6qGIfFhs1hO1Z6awmzGB17lpgiMTyUwZ2M6sJSo0PNlvNb2UIGKAAqqDSrKEHqmlA-zLvrws0JQ9mcec56cKxdSF3Qn-tP_deKnmMaKOq-eYA7bgMgkEcQRo_7LRAOQDB1-MLCKP8ffo4oxAPUzLnGEyZPtQIpmhfxOM-nvSAhZCpcZqgrZfS1QrpLUjgvvDlAuVjpuKcgDSmhnkmXImPJeCvxgwMLZCte0MBIl7ATvwaEve_2yKkqH8xYlmmO0HFVtua_2HWCoT4JujHpbP7jKamglM-Mc7oy-8Douat7AOGg"
  
  # 尝试不同的方式设置token
  tryCatch({
    # 尝试1: 查看get_query_content函数的源代码
    cat("Checking get_query_content function...\n")
    cat("Function body:", paste(deparse(body(get_query_content)), collapse = "\n"), "\n")
    
    # 尝试2: 查看get_opengwas_jwt函数的源代码
    cat("\nChecking get_opengwas_jwt function...\n")
    cat("Function body:", paste(deparse(body(get_opengwas_jwt)), collapse = "\n"), "\n")
    
    # 尝试3: 直接设置环境变量
    cat("\nSetting token via environment variables...\n")
    Sys.setenv(IEUGWASR_TOKEN = token)
    Sys.setenv(OPENGWAS_TOKEN = token)
    cat("Set IEUGWASR_TOKEN and OPENGWAS_TOKEN environment variables\n")
    
    # 尝试4: 查看当前环境变量
    cat("\nChecking current environment variables...\n")
    cat("IEUGWASR_API:", Sys.getenv("IEUGWASR_API"), "\n")
    cat("IEUGWASR_TOKEN:", substr(Sys.getenv("IEUGWASR_TOKEN"), 1, 20), "...\n")
    cat("OPENGWAS_TOKEN:", substr(Sys.getenv("OPENGWAS_TOKEN"), 1, 20), "...\n")
    
  }, error = function(e) {
    cat(sprintf("Error setting token: %s\n", e$message))
  })

  # 测试连接
  cat("Testing API connection...\n")
  tryCatch({
    # 使用gwasinfo函数测试连接
    test_info <- gwasinfo(id = "ieu-a-1239")
    if (!is.null(test_info)) {
      cat("✓ API connection successful (MEGASTROKE accessible)\n")
    } else {
      cat("✗ API connection failed\n")
    }
  }, error = function(e) {
    cat(sprintf("API connection test failed: %s\n", e$message))
  })

  cat("Packages loaded successfully\n")

}, error = function(e) {
  cat(sprintf("✗ Initialization failed: %s\n", e$message))
  stop("Cannot proceed")
})

# ========== 分析代码 ==========

# 设置输出目录
output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/MR_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# 基因列表（人类Gene Symbol）
# 临时只测试前3个基因，验证脚本是否正常运行
genes <- c("FDX1", "SLC31A1", "DLAT")

# 结局GWAS ID（MEGASTROKE all ischemic stroke）
outcome_id <- "ieu-a-1239"
outcome_name <- "Ischemic Stroke"

# 存储所有结果
results_list <- list()
snp_details_list <- list()

# 预先提取卒中GWAS的instrument（用于反向MR）
cat("Pre-extracting stroke GWAS instruments for reverse MR...\n")
stroke_instruments <- extract_instruments(
  outcomes = outcome_id,
  p1 = 5e-8,  # GWAS显著性阈值
  clump = TRUE,
  r2 = 0.001,  # 严格LD clumping
  kb = 10000  # 10Mb
)
if(nrow(stroke_instruments) > 0) {
  cat(sprintf("Extracted %d stroke GWAS instruments\n", nrow(stroke_instruments)))
} else {
  cat("Warning: No stroke GWAS instruments found\n")
  stroke_instruments <- data.frame()
}

# 循环分析每个基因
for(gene in genes) {
  cat(sprintf("\n===== Analyzing gene: %s =====\n", gene))
  
  tryCatch({
    # 1. 获取eQTL工具变量（GTEx Brain Cortex）
    cat("1. Extracting eQTL instruments...\n")

    # 使用variants_gene函数获取特定基因的eQTL数据
    cat(sprintf("Getting eQTL data for gene: %s\n", gene))
    eQTL_data <- variants_gene(gene = gene, pval = 5e-6)
    
    if(!is.null(eQTL_data) && nrow(eQTL_data) > 0) {
      # 筛选Brain Cortex组织的数据
      brain_cortex_data <- eQTL_data[grep("cortex|brain", eQTL_data$tissue, ignore.case = TRUE), ]
      
      if(nrow(brain_cortex_data) > 0) {
        # 提取SNP信息
        exposure_dat <- data.frame(
          SNP = brain_cortex_data$rsid,
          effect_allele.exposure = brain_cortex_data$ea,
          other_allele.exposure = brain_cortex_data$nea,
          eaf.exposure = brain_cortex_data$eaf,
          beta.exposure = brain_cortex_data$beta,
          se.exposure = brain_cortex_data$se,
          pval.exposure = brain_cortex_data$pval
        )
        cat(sprintf("Extracted %d eQTLs for %s in Brain Cortex\n", nrow(exposure_dat), gene))
      } else {
        cat("No Brain Cortex eQTL data found, skipping...\n")
        next
      }
    } else {
      cat("No eQTL data found for the gene, skipping...\n")
      next
    }
    
    # 方法3：如果仍然失败，跳过该基因
    if(nrow(exposure_dat) == 0) {
      cat("eQTL dataset extraction failed, skipping...\n")
      next
    }
    
    if(nrow(exposure_dat) == 0) {
      cat("No eQTL instruments found, skipping...\n")
      next
    }
    
    # 计算F统计量
    exposure_dat$F_statistic <- (exposure_dat$beta.exposure)^2 / (exposure_dat$se.exposure)^2
    cat(sprintf("Extracted %d eQTLs (F > 10: %d)\n", 
                nrow(exposure_dat), sum(exposure_dat$F_statistic > 10)))
    
    # 移除弱工具变量（F < 10）
    exposure_dat <- exposure_dat[exposure_dat$F_statistic >= 10, ]
    if(nrow(exposure_dat) == 0) {
      cat("No strong instruments after filtering, skipping...\n")
      next
    }
    
    # 保存SNP详情
    snp_details <- exposure_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
                                   "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "F_statistic")]
    colnames(snp_details) <- c("SNP", "effect_allele", "other_allele", 
                              "eaf", "beta.exposure", "se.exposure", "pval.exposure", "F_statistic")
    snp_details_list[[gene]] <- snp_details
    
    # 2. 获取结局数据（MEGASTROKE）
    cat("2. Extracting outcome data...\n")
    outcome_dat <- extract_outcome_data(
      snps = exposure_dat$SNP,
      outcomes = outcome_id
    )
    
    if(is.null(outcome_dat) || nrow(outcome_dat) == 0) {
      cat("No outcome data found, skipping...\n")
      next
    }
    
    # 3. 协调数据
    cat("3. Harmonizing data...\n")
    dat <- harmonise_data(exposure_dat, outcome_dat)
    
    # 4. MR分析
    cat("4. Performing MR analysis...\n")
    
    # 选择合适的方法
    if(nrow(dat) == 1) {
      # 只有1个SNP，使用Wald ratio
      method_list <- c("mr_wald_ratio")
    } else {
      # 多个SNP，使用多种方法
      method_list <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", 
                       "mr_weighted_mode", "mr_simple_mode")
    }
    
    res <- mr(dat, method_list = method_list)
    
    # 5. 敏感性分析
    cat("5. Performing sensitivity analysis...\n")
    
    # 异质性检验
    hetero <- mr_heterogeneity(dat)
    
    # 多效性检验
    pleio <- mr_pleiotropy_test(dat)
    
    # 留一法分析
    loo <- mr_leaveoneout(dat)
    
    # MR-PRESSO（如果有多个SNP）
    presso_res <- NULL
    if(nrow(dat) > 1) {
      tryCatch({
        presso_res <- mr_presso(BetaOutcome = "beta.outcome", 
                               BetaExposure = "beta.exposure", 
                               SdOutcome = "se.outcome", 
                               SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE, 
                               DISTORTIONtest = TRUE, 
                               data = dat, 
                               nboot = 1000)
      }, error = function(e) {
        cat(sprintf("MR-PRESSO failed: %s\n", e$message))
      })
    }
    
    # 6. 可视化
    cat("6. Generating visualizations...\n")
    
    # 森林图
    p_forest <- mr_forest_plot(res, dat)
    
    # 散点图
    p_scatter <- mr_scatter_plot(res, dat)
    
    # 漏斗图
    p_funnel <- mr_funnel_plot(res, dat)
    
    # 留一法图
    p_loo <- mr_leaveoneout_plot(loo)
    
    # 7. 整理结果
    cat("7. Compiling results...\n")
    
    # 合并结果与敏感性分析
    result_full <- res %>%
      left_join(hetero, by = "method") %>%
      left_join(pleio, by = "method")
    
    # 添加基因信息
    result_full$gene <- gene
    result_full$outcome <- outcome_name
    
    # 计算OR值（如果beta可用）
    if("b" %in% colnames(result_full)) {
      result_full$or <- exp(result_full$b)
      result_full$or_lci95 <- exp(result_full$b - 1.96 * result_full$se)
      result_full$or_uci95 <- exp(result_full$b + 1.96 * result_full$se)
    }
    
    # 保存结果
    results_list[[gene]] <- list(
      result = result_full,
      heterogeneity = hetero,
      pleiotropy = pleio,
      leaveoneout = loo,
      presso = presso_res,
      plots = list(forest = p_forest, scatter = p_scatter, 
                   funnel = p_funnel, leaveoneout = p_loo)
    )
    
    # 6. 反向MR分析（卒中风险 → 基因表达）
    cat("6. Performing reverse MR analysis...\n")
    tryCatch({
      if(nrow(stroke_instruments) > 0 && !is.null(eqtl_id) && eqtl_id %in% all_datasets$id) {
        # 获取基因表达作为结局
        # 使用之前确定的eQTL ID
        reverse_outcome_dat <- extract_outcome_data(
          snps = stroke_instruments$SNP,
          outcomes = eqtl_id  # 使用之前确定的eQTL ID
        )
        
        if(!is.null(reverse_outcome_dat) && nrow(reverse_outcome_dat) > 0) {
          # 协调数据
          reverse_dat <- harmonise_data(
            exposure_dat = stroke_instruments,
            outcome_dat = reverse_outcome_dat
          )
          
          # 执行反向MR分析
          reverse_res <- mr(reverse_dat, method_list = c("mr_ivw"))
          
          # 保存反向MR结果
          results_list[[gene]]$reverse_result <- reverse_res
          
          cat(sprintf("Reverse MR analysis completed for %s\n", gene))
        } else {
          cat("No reverse outcome data found, skipping...\n")
        }
      } else {
        cat("Stroke instruments not available or eqtl_id invalid, skipping reverse MR...\n")
      }
    }, error = function(e) {
      cat(sprintf("Reverse MR analysis failed: %s\n", e$message))
    })
    
    cat(sprintf("Analysis completed for %s\n", gene))
    
  }, error = function(e) {
    cat(sprintf("Error analyzing %s: %s\n", gene, e$message))
  })
}

# 导出结果
cat("\n===== Exporting results =====\n")

# 1. 保存每个基因的详细结果
for(gene in names(results_list)) {
  res <- results_list[[gene]]$result
  if(!is.null(res)) {
    output_file <- file.path(output_dir, sprintf("MR_results_%s_%s.csv", gene, gsub(" ", "_", outcome_name)))
    write.csv(res, output_file, row.names = FALSE)
    cat(sprintf("Saved detailed results for %s\n", gene))
  }
}

# 2. 保存SNP详情
for(gene in names(snp_details_list)) {
  snp_dat <- snp_details_list[[gene]]
  if(!is.null(snp_dat)) {
    output_file <- file.path(output_dir, sprintf("SNP_details_%s.csv", gene))
    write.csv(snp_dat, output_file, row.names = FALSE)
    cat(sprintf("Saved SNP details for %s\n", gene))
  }
}

# 3. 生成汇总表格
cat("Generating summary table...\n")
summary_list <- list()
for(gene in names(results_list)) {
  res <- results_list[[gene]]$result
  if(!is.null(res)) {
    # 提取IVW结果
    ivw_res <- res %>% filter(method == "Inverse variance weighted")
    if(nrow(ivw_res) > 0) {
      summary_list[[gene]] <- ivw_res
    }
  }
}

if(length(summary_list) > 0) {
  summary_df <- bind_rows(summary_list)
  
  # 标记显著性（Bonferroni校正）
  alpha_bonferroni <- 0.05 / length(genes)
  summary_df$significant <- ifelse(summary_df$pval < alpha_bonferroni, "Yes", "No")
  
  # 标记异质性和多效性
  summary_df$heterogeneity <- ifelse(summary_df$Q_pval < 0.05, "High", "Low")
  summary_df$pleiotropy <- ifelse(summary_df$egger_intercept_pval < 0.05, "Significant", "Not significant")
  
  # 保存汇总表格
  output_file <- file.path(output_dir, "MR_summary_all_genes.csv")
  write.csv(summary_df, output_file, row.names = FALSE)
  cat("Saved summary table\n")
}

# 4. 保存可视化
cat("Saving visualizations...\n")

# 合并森林图
pdf(file.path(output_dir, "forest_plots_combined.pdf"), width = 12, height = 8)
tryCatch({
  for(gene in names(results_list)) {
    p <- results_list[[gene]]$plots$forest
    if(!is.null(p)) {
      print(p + ggtitle(sprintf("Forest plot: %s vs %s", gene, outcome_name)) + theme_bw(base_size = 12))
    }
  }
}, error = function(e) {}
) 
dev.off()

# 合并散点图
pdf(file.path(output_dir, "scatter_plots_combined.pdf"), width = 10, height = 8)
tryCatch({
  for(gene in names(results_list)) {
    p <- results_list[[gene]]$plots$scatter
    if(!is.null(p)) {
      print(p + ggtitle(sprintf("Scatter plot: %s vs %s", gene, outcome_name)) + theme_bw(base_size = 12))
    }
  }
}, error = function(e) {}
) 
dev.off()

# 合并漏斗图
pdf(file.path(output_dir, "funnel_plots_combined.pdf"), width = 10, height = 8)
tryCatch({
  for(gene in names(results_list)) {
    p <- results_list[[gene]]$plots$funnel
    if(!is.null(p)) {
      print(p + ggtitle(sprintf("Funnel plot: %s vs %s", gene, outcome_name)) + theme_bw(base_size = 12))
    }
  }
}, error = function(e) {}
) 
dev.off()

# 合并留一法图
pdf(file.path(output_dir, "loo_plots_combined.pdf"), width = 12, height = 8)
tryCatch({
  for(gene in names(results_list)) {
    p <- results_list[[gene]]$plots$leaveoneout
    if(!is.null(p)) {
      print(p + ggtitle(sprintf("Leave-one-out plot: %s vs %s", gene, outcome_name)) + theme_bw(base_size = 12))
    }
  }
}, error = function(e) {}
) 
dev.off()

# 单个基因单独成图
for(gene in names(results_list)) {
  plots <- results_list[[gene]]$plots
  if(!is.null(plots)) {
    # 森林图
    pdf(file.path(output_dir, sprintf("forest_plot_%s.pdf", gene)), width = 12, height = 8)
    tryCatch(print(plots$forest + ggtitle(sprintf("Forest plot: %s vs %s", gene, outcome_name)) + theme_bw(base_size = 12)), error = function(e) {})
    dev.off()
    
    # 散点图
    pdf(file.path(output_dir, sprintf("scatter_plot_%s.pdf", gene)), width = 10, height = 8)
    tryCatch(print(plots$scatter + ggtitle(sprintf("Scatter plot: %s vs %s", gene, outcome_name)) + theme_bw(base_size = 12)), error = function(e) {})
    dev.off()
    
    # 漏斗图
    pdf(file.path(output_dir, sprintf("funnel_plot_%s.pdf", gene)), width = 10, height = 8)
    tryCatch(print(plots$funnel + ggtitle(sprintf("Funnel plot: %s vs %s", gene, outcome_name)) + theme_bw(base_size = 12)), error = function(e) {})
    dev.off()
    
    # 留一法图
    pdf(file.path(output_dir, sprintf("loo_plot_%s.pdf", gene)), width = 12, height = 8)
    tryCatch(print(plots$leaveoneout + ggtitle(sprintf("Leave-one-out plot: %s vs %s", gene, outcome_name)) + theme_bw(base_size = 12)), error = function(e) {})
    dev.off()
  }
}

cat("\nAnalysis completed! All results saved to:", output_dir, "\n")
