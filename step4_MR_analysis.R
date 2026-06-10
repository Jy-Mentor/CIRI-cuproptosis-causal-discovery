# ============================================================================
# 第4步：孟德尔随机化（MR）人群验证
# 验证BCP核心基因与缺血性脑卒中的因果关系
# ============================================================================

options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

# 安装TwoSampleMR包
if(!"TwoSampleMR" %in% installed.packages()){
  install.packages("remotes")
  remotes::install_github("MRCIEU/TwoSampleMR")
}
if(!"ieugwasr" %in% installed.packages()){
  remotes::install_github("MRCIEU/ieugwasr")
}

library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(ggplot2)

cat("=== 第4步：孟德尔随机化（MR）人群验证 ===\n\n")

result_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/causal_analysis_results"
mr_dir <- file.path(result_dir, "MR_results")
dir.create(mr_dir, showWarnings=FALSE, recursive=TRUE)

# ============================================================================
# 1. 定义分析基因
# ============================================================================
cat("【4.1】定义BCP核心基因...\n")

# 8个BCP核心基因
bcp_genes <- c("IL6", "STAT3", "NFKB1", "TGFB1", "AGER", "PTGS2", "TLR4", "FDX1")

# 检查哪些基因已经在D:/EQTL中分析过
existing_results <- c("STAT3", "TLR4", "FDX1")  # 根据之前的csv文件
new_genes <- setdiff(bcp_genes, existing_results)

cat(sprintf("  总共 %d 个BCP基因\n", length(bcp_genes)))
cat(sprintf("  已分析: %s\n", paste(existing_results, collapse=", ")))
cat(sprintf("  新分析: %s\n", paste(new_genes, collapse=", ")))

# ============================================================================
# 2. 读取eQTL数据（从D:/EQTL）
# ============================================================================
cat("\n【4.2】读取eQTL数据...\n")

# 读取clump后的eQTL数据
eqtl_file <- "D:/EQTL/clump/eQTLgen_allgene_p_5e-08_kb_10000_r2_0.01.xlsx"
if(!file.exists(eqtl_file)) {
  # 尝试其他文件
  eQTL_files <- list.files("D:/EQTL/clump", pattern="eQTLgen.*\.xlsx", full.names=TRUE)
  if(length(eQTL_files) > 0) {
    eqtl_file <- eQTL_files[1]
  }
}

if(!file.exists(eqtl_file)) {
  cat("  [警告] eQTL文件不存在，尝试使用在线IEU GWAS数据...\n")
  use_online <- TRUE
} else {
  cat(sprintf("  读取eQTL文件: %s\n", basename(eqtl_file)))
  library(readxl)
  eqtl_data <- read_excel(eqtl_file)
  cat(sprintf("  eQTL数据: %d 行\n", nrow(eqtl_data)))
  use_online <- FALSE
}

# ============================================================================
# 3. 读取FinnGen GWAS数据（缺血性脑卒中）
# ============================================================================
cat("\n【4.3】读取FinnGen GWAS数据（缺血性脑卒中）...\n")

finngen_file <- "D:/EQTL/finngen_R12_I9_STR"
if(file.exists(paste0(finngen_file, ".rds"))) {
  gwas_data <- readRDS(paste0(finngen_file, ".rds"))
  cat(sprintf("  GWAS数据: %d 行\n", nrow(gwas_data)))
} else if(file.exists(paste0(finngen_file, ".csv"))) {
  gwas_data <- fread(paste0(finngen_file, ".csv"))
  cat(sprintf("  GWAS数据: %d 行\n", nrow(gwas_data)))
} else {
  cat("  [警告] FinnGen数据不存在，尝试使用在线数据...\n")
  # 使用IEU GWAS的卒中数据
  gwas_data <- NULL
}

# ============================================================================
# 4. MR分析函数
# ============================================================================
cat("\n【4.4】执行MR分析...\n")

mr_analysis_gene <- function(gene_name, eqtl_data, gwas_data) {
  cat(sprintf("\n  分析基因: %s...\n", gene_name))
  
  tryCatch({
    # 方法1: 使用IEU GWAS在线数据
    if(use_online || is.null(eqtl_data)) {
      cat("    使用IEU GWAS在线数据...\n")
      
      # 获取eQTL数据（血液）
      exposure_dat <- extract_instruments(
        outcomes = paste0("eqtl-a-", gene_name),
        p1 = 5e-08,
        clump = TRUE,
        r2 = 0.001,
        kb = 10000
      )
      
      if(is.null(exposure_dat) || nrow(exposure_dat) == 0) {
        cat(sprintf("    [跳过] %s 没有找到eQTL instruments\n", gene_name))
        return(NULL)
      }
      
      cat(sprintf("    找到 %d 个eQTL SNPs\n", nrow(exposure_dat)))
      
      # 获取结局数据（缺血性脑卒中）
      # FinnGen R12 I9_STR
      outcome_dat <- extract_outcome_data(
        snps = exposure_dat$SNP,
        outcomes = "finn-b-I9_STR"  # FinnGen缺血性脑卒中
      )
      
    } else {
      # 方法2: 使用本地数据
      cat("    使用本地eQTL数据...\n")
      
      # 提取该基因的eQTL数据
      gene_eqtl <- eqtl_data[eqtl_data$gene == gene_name, ]
      
      if(nrow(gene_eqtl) == 0) {
        cat(sprintf("    [跳过] %s 在eQTL数据中未找到\n", gene_name))
        return(NULL)
      }
      
      cat(sprintf("    找到 %d 个eQTL SNPs\n", nrow(gene_eqtl)))
      
      # 格式化exposure数据
      exposure_dat <- data.frame(
        SNP = gene_eqtl$SNP,
        beta = gene_eqtl$beta.exposure,
        se = gene_eqtl$se.exposure,
        effect_allele = gene_eqtl$effect_allele.exposure,
        other_allele = gene_eqtl$other_allele.exposure,
        pval = gene_eqtl$pval.exposure,
        eaf = gene_eqtl$eaf.exposure,
        exposure = gene_name,
        id.exposure = gene_name,
        stringsAsFactors = FALSE
      )
      
      # 获取结局数据
      if(!is.null(gwas_data)) {
        # 从本地GWAS数据中提取
        outcome_snps <- gwas_data[gwas_data$rsids %in% exposure_dat$SNP, ]
        
        if(nrow(outcome_snps) == 0) {
          cat(sprintf("    [跳过] 没有找到共同SNPs\n"))
          return(NULL)
        }
        
        outcome_dat <- data.frame(
          SNP = outcome_snps$rsids,
          beta = outcome_snps$beta,
          se = outcome_snps$sebeta,
          effect_allele = outcome_snps$alt,
          other_allele = outcome_snps$ref,
          pval = outcome_snps$pval,
          eaf = outcome_snps$af_alt,
          outcome = "Ischemic Stroke",
          id.outcome = "finngen_R12_I9_STR",
          stringsAsFactors = FALSE
        )
      } else {
        outcome_dat <- NULL
      }
    }
    
    if(is.null(outcome_dat) || nrow(outcome_dat) == 0) {
      cat(sprintf("    [跳过] 没有结局数据\n"))
      return(NULL)
    }
    
    # 协调数据
    cat("    协调exposure和outcome数据...\n")
    dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat,
      action = 2
    )
    
    if(nrow(dat) < 3) {
      cat(sprintf("    [跳过] 协调后SNP数量不足 (%d < 3)\n", nrow(dat)))
      return(NULL)
    }
    
    cat(sprintf("    协调后: %d 个SNPs\n", nrow(dat)))
    
    # 执行MR分析
    cat("    执行MR分析...\n")
    res <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw", "mr_weighted_median", "mr_egger_regression"))
    
    # 敏感性分析
    cat("    敏感性分析...\n")
    
    # Cochran's Q检验（异质性）
    het <- mr_heterogeneity(dat)
    
    # Egger截距（水平多效性）
    pleio <- mr_pleiotropy_test(dat)
    
    # 留一法
    single <- mr_leaveoneout(dat)
    
    # 保存结果
    results <- list(
      gene = gene_name,
      mr_results = res,
      heterogeneity = het,
      pleiotropy = pleio,
      leaveoneout = single,
      n_snps = nrow(dat),
      dat = dat
    )
    
    # 输出主要结果
    cat("    MR结果:\n")
    for(i in 1:nrow(res)) {
      method <- res$method[i]
      b <- res$b[i]
      p <- res$pval[i]
      sig <- ifelse(p < 0.05, "*", "")
      cat(sprintf("      %s: b=%.3f, p=%.3f %s\n", method, b, p, sig))
    }
    
    # 保存到文件
    write.csv(res, file.path(mr_dir, paste0("MR_results_", gene_name, ".csv")), row.names=FALSE)
    
    cat(sprintf("    ✅ %s 分析完成\n", gene_name))
    return(results)
    
  }, error = function(e) {
    cat(sprintf("    [错误] %s: %s\n", gene_name, conditionMessage(e)))
    return(NULL)
  })
}

# ============================================================================
# 5. 对所有BCP基因进行MR分析
# ============================================================================

mr_results_all <- list()

for(gene in bcp_genes) {
  res <- mr_analysis_gene(gene, eqtl_data, gwas_data)
  if(!is.null(res)) {
    mr_results_all[[gene]] <- res
  }
}

# ============================================================================
# 6. 汇总结果
# ============================================================================
cat("\n【4.5】MR分析汇总...\n")

if(length(mr_results_all) > 0) {
  cat(sprintf("\n  成功分析 %d 个基因:\n", length(mr_results_all)))
  
  # 汇总表
  summary_df <- data.frame(
    gene = names(mr_results_all),
    n_snps = sapply(mr_results_all, function(x) x$n_snps),
    ivw_beta = NA,
    ivw_p = NA,
    wm_beta = NA,
    wm_p = NA,
    stringsAsFactors = FALSE
  )
  
  for(i in 1:nrow(summary_df)) {
    gene <- summary_df$gene[i]
    res <- mr_results_all[[gene]]$mr_results
    
    # IVW
    ivw_row <- res[res$method == "Inverse variance weighted", ]
    if(nrow(ivw_row) > 0) {
      summary_df$ivw_beta[i] <- ivw_row$b
      summary_df$ivw_p[i] <- ivw_row$pval
    }
    
    # Weighted median
    wm_row <- res[res$method == "Weighted median", ]
    if(nrow(wm_row) > 0) {
      summary_df$wm_beta[i] <- wm_row$b
      summary_df$wm_p[i] <- wm_row$pval
    }
  }
  
  cat("\n  MR分析结果汇总:\n")
  print(summary_df)
  
  write.csv(summary_df, file.path(mr_dir, "MR_summary_results.csv"), row.names=FALSE)
  
  # 识别显著基因
  sig_genes <- summary_df$gene[summary_df$ivw_p < 0.05 | summary_df$wm_p < 0.05]
  if(length(sig_genes) > 0) {
    cat(sprintf("\n  ✅ 显著因果关联基因: %s\n", paste(sig_genes, collapse=", ")))
  } else {
    cat("\n  ⚠️ 没有发现显著因果关联 (p < 0.05)\n")
  }
  
  # 保存所有结果
  saveRDS(mr_results_all, file.path(mr_dir, "all_MR_results.rds"))
  
} else {
  cat("  [警告] 没有成功完成任何MR分析\n")
}

# ============================================================================
# 7. 可视化
# ============================================================================
cat("\n【4.6】生成可视化...\n")

if(length(mr_results_all) > 0) {
  pdf(file.path(mr_dir, "MR_forest_plot.pdf"), width=12, height=10)
  
  # 森林图
  res_all <- bind_rows(lapply(mr_results_all, function(x) x$mr_results))
  
  if(nrow(res_all) > 0) {
    p <- ggplot(res_all, aes(x=exposure, y=b, ymin=b-1.96*se, ymax=b+1.96*se, color=method)) +
      geom_pointrange(position=position_dodge(width=0.5)) +
      geom_hline(yintercept=0, linetype="dashed", color="gray") +
      coord_flip() +
      labs(title="MR Analysis - BCP Genes vs Ischemic Stroke",
           x="Gene", y="Effect Size (Beta)") +
      theme_bw() +
      theme(legend.position="bottom")
    print(p)
  }
  
  dev.off()
  cat("  ✅ 森林图已保存\n")
}

cat("\n"); cat(rep("=", 60), sep=""); cat("\n")
cat("第4步完成！\n")
cat("结果目录:", mr_dir, "\n")
cat("\n关键发现:\n")
if(exists("sig_genes") && length(sig_genes) > 0) {
  cat(sprintf("  - 显著因果关联基因: %d 个\n", length(sig_genes)))
  cat(sprintf("  - %s\n", paste(sig_genes, collapse=", ")))
} else {
  cat("  - 没有发现显著因果关联\n")
}
cat("  - 在人群水平验证了BCP基因与缺血性脑卒中的关系\n")
cat(rep("=", 60), sep=""); cat("\n")
