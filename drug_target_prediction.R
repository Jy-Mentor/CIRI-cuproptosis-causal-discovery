#!/usr/bin/env Rscript
# ================================================================================
# 药物靶点预测脚本 - MR 分析优化路线 B
# 数据源：DrugBank, DGIdb, Open Targets
# ================================================================================

# 包安装与加载
install_and_load_packages <- function() {
  packages <- c(
    "dplyr",
    "readr",
    "stringr",
    "httr",
    "jsonlite",
    "tibble",
    "writexl",
    "openxlsx"
  )
  
  message("正在检查并安装所需包...")
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("正在安装:", pkg))
      tryCatch({
        install.packages(pkg, repos = "https://cloud.r-project.org/")
      }, error = function(e) {
        message(paste("安装失败:", pkg, "-", e$message))
      })
    }
  }
  
  invisible(lapply(packages, library, character.only = TRUE))
  message("所有包加载完成")
}

# 从 MR 结果中提取显著基因
extract_significant_genes <- function(mr_results_file) {
  message("\n=== 提取显著基因 ===")
  
  results <- read.csv(mr_results_file, stringsAsFactors = FALSE)
  
  # 筛选成功的基因
  successful_genes <- results %>%
    filter(status == "SUCCESS" | status == "HETEROGENEITY") %>%
    filter(!is.na(discovery_pval))
  
  # FDR 显著基因
  fdr_sig_genes <- successful_genes %>%
    filter(!is.na(fdr_qval)) %>%
    filter(fdr_qval < 0.05) %>%
    pull(gene) %>%
    unique()
  
  # P 值显著基因
  pval_sig_genes <- successful_genes %>%
    filter(discovery_pval < 0.05) %>%
    pull(gene) %>%
    unique()
  
  message(paste("FDR 显著基因:", length(fdr_sig_genes)))
  message(paste("P 值显著基因:", length(pval_sig_genes)))
  
  return(list(
    fdr_significant = fdr_sig_genes,
    pval_significant = pval_sig_genes,
    all_successful = unique(successful_genes$gene)
  ))
}

# 查询 DrugBank (通过 web scraping)
query_drugbank <- function(gene_list) {
  message("\n=== 查询 DrugBank ===")
  
  results <- list()
  
  for (i in seq_along(gene_list)) {
    gene <- gene_list[i]
    message(paste("  [", i, "/", length(gene_list), "] 查询:", gene, sep = ""))
    
    tryCatch({
      # DrugBank 搜索 URL
      search_url <- paste0("https://go.drugbank.com/uniprot/search?query=", gene)
      
      # 使用 httr 获取页面
      response <- httr::GET(search_url, 
                           httr::set_cookies(.cookies = c()),
                           httr::timeout(10))
      
      if (response$status_code == 200) {
        # 简单解析（实际应该用 rvest）
        content <- httr::content(response, "text", encoding = "UTF-8")
        
        # 检查是否有药物关联
        if (grepl("drug", content, ignore.case = TRUE)) {
          results[[gene]] <- list(
            found = TRUE,
            url = search_url,
            note = "可能存在药物关联"
          )
        } else {
          results[[gene]] <- list(
            found = FALSE,
            url = search_url,
            note = "未找到药物关联"
          )
        }
      } else {
        results[[gene]] <- list(
          found = FALSE,
          url = search_url,
          note = paste("HTTP 错误:", response$status_code)
        )
      }
      
    }, error = function(e) {
      results[[gene]] <<- list(
        found = FALSE,
        error = e$message
      )
    })
    
    # 避免请求过快
    Sys.sleep(1)
  }
  
  return(results)
}

# 查询 DGIdb (Drug-Gene Interaction Database)
query_dgidb <- function(gene_list) {
  message("\n=== 查询 DGIdb (Drug-Gene Interaction Database) ===")
  
  all_results <- list()
  
  # DGIdb API URL
  base_url <- "https://dgidb.org/api/interactions.json"
  
  for (i in seq_along(gene_list)) {
    gene <- gene_list[i]
    message(paste("  [", i, "/", length(gene_list), "] 查询:", gene, sep = ""))
    
    tryCatch({
      # 构建请求
      request_url <- paste0(base_url, "?genes=", gene, "&drug_types=small+molecule,biotech")
      
      response <- httr::GET(request_url, httr::timeout(10))
      
      if (response$status_code == 200) {
        data <- jsonlite::fromJSON(httr::content(response, "text"))
        
        if (!is.null(data$interactions) && nrow(data$interactions) > 0) {
          interactions <- data$interactions %>%
            select(gene_name, drug_name, interaction_type, 
                   drug_chembl_id, drug_type, source) %>%
            distinct()
          
          all_results[[gene]] <- list(
            found = TRUE,
            interactions = interactions,
            count = nrow(interactions)
          )
          
          message(paste("    找到", nrow(interactions), "个药物相互作用"))
        } else {
          all_results[[gene]] <- list(
            found = FALSE,
            interactions = NULL,
            count = 0
          )
          message("    未找到药物相互作用")
        }
      } else {
        all_results[[gene]] <- list(
          found = FALSE,
          error = paste("HTTP", response$status_code)
        )
      }
      
    }, error = function(e) {
      all_results[[gene]] <<- list(
        found = FALSE,
        error = e$message
      )
      message(paste("    查询失败:", e$message))
    })
    
    Sys.sleep(0.5)
  }
  
  return(all_results)
}

# 查询 Open Targets Platform
query_opentargets <- function(gene_list, disease = "stroke") {
  message("\n=== 查询 Open Targets Platform ===")
  message(paste("疾病:", disease))
  
  all_results <- list()
  
  # Open Targets GraphQL API
  api_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  for (i in seq_along(gene_list)) {
    gene <- gene_list[i]
    message(paste("  [", i, "/", length(gene_list), "] 查询:", gene, sep = ""))
    
    tryCatch({
      # GraphQL 查询
      query <- list(
        query = paste0('
          {
            geneInfo(geneId: "', gene, '") {
              id
              symbol
              approvedSymbol
            }
            targetsForGene(geneId: "', gene, '") {
              id
              approvedSymbol
              tractability {
                level
                source
              }
              knownDrugs {
                id
                name
              }
            }
          }
        ')
      )
      
      response <- httr::POST(
        api_url,
        body = jsonlite::toJSON(query, auto_unbox = TRUE),
        encode = "json",
        httr::timeout(10)
      )
      
      if (response$status_code == 200) {
        data <- jsonlite::fromJSON(httr::content(response, "text"))
        
        if (!is.null(data$data) && !is.null(data$data$targetsForGene)) {
          targets <- data$data$targetsForGene
          
          all_results[[gene]] <- list(
            found = TRUE,
            tractability = if(!is.null(targets$tractability)) targets$tractability$level else NA,
            known_drugs = if(!is.null(targets$knownDrugs)) nrow(targets$knownDrugs) else 0
          )
          
          message(paste("    可成药性:", all_results[[gene]]$tractability))
          message(paste("    已知药物:", all_results[[gene]]$known_drugs))
        } else {
          all_results[[gene]] <- list(
            found = FALSE,
            tractability = NA,
            known_drugs = 0
          )
          message("    未找到靶点信息")
        }
      } else {
        all_results[[gene]] <- list(
          found = FALSE,
          error = paste("HTTP", response$status_code)
        )
      }
      
    }, error = function(e) {
      all_results[[gene]] <<- list(
        found = FALSE,
        error = e$message
      )
      message(paste("    查询失败:", e$message))
    })
    
    Sys.sleep(0.5)
  }
  
  return(all_results)
}

# 整合所有药物靶点信息
integrate_drug_targets <- function(dgidb_results, opentargets_results, gene_list) {
  message("\n=== 整合药物靶点信息 ===")
  
  target_table <- data.frame(
    gene = gene_list,
    dgidb_found = FALSE,
    dgidb_count = 0,
    dgidb_drugs = "",
    opentargets_found = FALSE,
    tractability = NA_character_,
    known_drugs_count = 0,
    drug_priority = "Unknown",
    stringsAsFactors = FALSE
  )
  
  for (gene in gene_list) {
    # DGIdb 结果
    if (!is.null(dgidb_results[[gene]]) && dgidb_results[[gene]]$found) {
      target_table[target_table$gene == gene, "dgidb_found"] <- TRUE
      target_table[target_table$gene == gene, "dgidb_count"] <- dgidb_results[[gene]]$count
      
      if (!is.null(dgidb_results[[gene]]$interactions)) {
        drugs <- unique(dgidb_results[[gene]]$interactions$drug_name)
        target_table[target_table$gene == gene, "dgidb_drugs"] <- paste(drugs, collapse = "; ")
      }
    }
    
    # Open Targets 结果
    if (!is.null(opentargets_results[[gene]]) && opentargets_results[[gene]]$found) {
      target_table[target_table$gene == gene, "opentargets_found"] <- TRUE
      target_table[target_table$gene == gene, "tractability"] <- 
        opentargets_results[[gene]]$tractability
      target_table[target_table$gene == gene, "known_drugs_count"] <- 
        opentargets_results[[gene]]$known_drugs
    }
    
    # 药物优先级
    priority <- "Unknown"
    if (target_table[target_table$gene == gene, "dgidb_count"] > 0) {
      priority <- "High"
    } else if (!is.na(target_table[target_table$gene == gene, "tractability"]) && 
               target_table[target_table$gene == gene, "tractability"] %in% c("Tchem", "Tbio")) {
      priority <- "Medium"
    } else if (target_table[target_table$gene == gene, "known_drugs_count"] > 0) {
      priority <- "Low"
    }
    
    target_table[target_table$gene == gene, "drug_priority"] <- priority
  }
  
  return(target_table)
}

# 保存药物靶点结果
save_drug_targets <- function(target_table, output_dir) {
  message("\n=== 保存药物靶点结果 ===")
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 保存为 CSV
  csv_file <- file.path(output_dir, "drug_targets_summary.csv")
  write.csv(target_table, csv_file, row.names = FALSE)
  message(paste("靶点摘要表:", csv_file))
  
  # 保存为 Excel
  xlsx_file <- file.path(output_dir, "drug_targets_summary.xlsx")
  openxlsx::write.xlsx(target_table, xlsx_file)
  message(paste("Excel 文件:", xlsx_file))
  
  # 高优先级靶点
  high_priority <- target_table %>%
    filter(drug_priority == "High") %>%
    arrange(desc(dgidb_count))
  
  if (nrow(high_priority) > 0) {
    high_file <- file.path(output_dir, "high_priority_targets.csv")
    write.csv(high_priority, high_file, row.names = FALSE)
    message(paste("高优先级靶点:", high_file, "(", nrow(high_priority), "个基因)"))
  }
  
  # 生成可视化
  message("生成可视化图表...")
  
  # 1. 药物靶点优先级分布
  priority_counts <- table(target_table$drug_priority)
  png(file.path(output_dir, "drug_priority_distribution.png"), 
      width = 800, height = 600, res = 100)
  barplot(priority_counts,
          main = "药物靶点优先级分布",
          xlab = "优先级",
          ylab = "基因数量",
          col = c("red", "orange", "yellow", "gray"))
  dev.off()
  message("  优先级分布图：drug_priority_distribution.png")
  
  # 2. DGIdb 药物数量分布
  if (sum(target_table$dgidb_found) > 0) {
    png(file.path(output_dir, "dgidb_drug_count.png"), 
        width = 800, height = 600, res = 100)
    hist(target_table$dgidb_count[target_table$dgidb_found],
         main = "DGIdb 药物相互作用数量分布",
         xlab = "药物数量",
         ylab = "基因频率",
         col = "steelblue",
         breaks = 10)
    dev.off()
    message("  DGIdb 药物数量分布图：dgidb_drug_count.png")
  }
  
  message(paste("所有结果已保存至:", output_dir))
}

# 生成药物靶点报告
generate_drug_target_report <- function(target_table, output_file) {
  message("\n=== 生成药物靶点报告 ===")
  
  report <- c(
    "# 药物靶点预测报告",
    "",
    "## 分析概述",
    paste("分析日期:", Sys.time()),
    paste("分析基因数:", nrow(target_table)),
    "",
    "## 高优先级靶点",
    ""
  )
  
  high_priority <- target_table %>%
    filter(drug_priority == "High") %>%
    arrange(desc(dgidb_count))
  
  if (nrow(high_priority) > 0) {
    for (i in 1:nrow(high_priority)) {
      gene <- high_priority$gene[i]
      report <- c(report, paste0(
        "### ", gene, "  \n",
        "- **DGIdb 药物数**: ", high_priority$dgidb_count[i], "  \n",
        "- **已知药物**: ", high_priority$dgidb_drugs[i], "  \n",
        "- **Open Targets 可成药性**: ", high_priority$tractability[i], "  \n",
        "- **Open Targets 已知药物数**: ", high_priority$known_drugs_count[i], "  \n"
      ))
    }
  } else {
    report <- c(report, "无高优先级靶点")
  }
  
  report <- c(report, "", "## 所有靶点详情", "")
  report <- c(report, "| Gene | DGIdb Found | DGIdb Count | Drugs | OpenTargets | Tractability | Priority |")
  report <- c(report, "|------|-------------|-------------|-------|-------------|--------------|----------|")
  
  for (i in 1:nrow(target_table)) {
    report <- c(report, paste0(
      "| ", target_table$gene[i], 
      " | ", ifelse(target_table$dgidb_found[i], "Yes", "No"),
      " | ", target_table$dgidb_count[i],
      " | ", ifelse(nchar(target_table$dgidb_drugs[i]) > 30, 
                   paste0(substr(target_table$dgidb_drugs[i], 1, 27), "..."),
                   target_table$dgidb_drugs[i]),
      " | ", ifelse(target_table$opentargets_found[i], "Yes", "No"),
      " | ", target_table$tractability[i],
      " | ", target_table$drug_priority[i],
      " |"
    ))
  }
  
  writeLines(report, output_file)
  message(paste("报告已保存:", output_file))
}

# 主函数
main <- function() {
  message(rep("=", 60))
  message("药物靶点预测 - MR 分析优化路线 B")
  message(rep("=", 60))
  
  # 1. 加载包
  install_and_load_packages()
  
  # 2. 提取显著基因
  mr_results_file <- "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/MR_results_main_optimized.csv"
  
  if (!file.exists(mr_results_file)) {
    stop("找不到 MR 结果文件")
  }
  
  genes <- extract_significant_genes(mr_results_file)
  
  # 使用 FDR 显著 + P 值显著基因
  target_gene_list <- unique(c(genes$fdr_significant, genes$pval_significant))
  
  if (length(target_gene_list) == 0) {
    message("无显著基因，使用所有成功基因")
    target_gene_list <- genes$all_successful
  }
  
  message(paste("\n用于药物靶点预测的基因数:", length(target_gene_list)))
  message("基因列表:", paste(target_gene_list, collapse = ", "))
  
  # 3. 查询 DGIdb
  dgidb_results <- query_dgidb(target_gene_list)
  
  # 4. 查询 Open Targets
  opentargets_results <- query_opentargets(target_gene_list, disease = "stroke")
  
  # 5. 整合结果
  target_table <- integrate_drug_targets(dgidb_results, opentargets_results, target_gene_list)
  
  # 6. 保存结果
  output_dir <- "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/drug_targets"
  save_drug_targets(target_table, output_dir)
  
  # 7. 生成报告
  report_file <- file.path(output_dir, "drug_target_report.md")
  generate_drug_target_report(target_table, report_file)
  
  message("\n", rep("=", 60))
  message("药物靶点预测完成!")
  message(paste("结果目录:", output_dir))
  message(rep("=", 60))
}

# 运行主函数
if (!interactive()) {
  main()
} else {
  message("交互模式下运行，请手动调用 main()")
}
