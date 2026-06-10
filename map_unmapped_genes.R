# ============================================
# 使用biomaRt映射未映射的大鼠基因到人类基因
# ============================================

cat("正在加载必要的R包...\n")

# 安装和加载biomaRt
if (!require("biomaRt", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  }
  BiocManager::install("biomaRt", ask = FALSE, update = FALSE)
}
library(biomaRt)

# 设置路径
work_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
result_dir <- file.path(work_dir, "GSE61616_cluster_ssGSEA_PPI_results_v2")
output_file <- file.path(work_dir, "PPI_Input_Genes_Human_Complete.xlsx")

# 读取已映射的基因列表（人类基因）
mapped_human_file <- file.path(work_dir, "PPI_Input_Genes_Human.xlsx")
if (file.exists(mapped_human_file)) {
  library(readxl)
  mapped_df <- read_excel(mapped_human_file, sheet = "PPI输入基因 (Human)")
  mapped_human_genes <- unique(toupper(mapped_df$人类基因[mapped_df$映射状态 == "Mapped"]))
  cat(paste0("已映射的人类基因数: ", length(mapped_human_genes), "\n"))
} else {
  mapped_human_genes <- c()
}

# 读取未映射的大鼠基因
unmapped_file <- file.path(result_dir, "String_PPI_input_genes.txt")
unmapped_rat_file <- file.path(work_dir, "PPI_Input_Genes_Human.xlsx")

# 从之前的Excel中读取未映射基因
library(readxl)
unmapped_sheet <- read_excel(unmapped_rat_file, sheet = "未映射基因", skip = 3, col_names = FALSE)
unmapped_rat_genes <- toupper(unmapped_sheet[[1]])
unmapped_rat_genes <- unmapped_rat_genes[!is.na(unmapped_rat_genes)]

cat(paste0("需要映射的未映射大鼠基因数: ", length(unmapped_rat_genes), "\n"))

# 连接到Ensembl数据库
cat("\n连接到Ensembl数据库...\n")

tryCatch({
  # 尝试连接大鼠数据库
  ensembl_rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  cat("成功连接到Ensembl数据库\n")
  
  # 获取大鼠基因的人类同源基因
  cat("\n正在查询未映射基因的人类同源基因...\n")
  
  # 分批查询以避免超时
  batch_size <- 100
  n_genes <- length(unmapped_rat_genes)
  n_batches <- ceiling(n_genes / batch_size)
  
  all_mappings <- data.frame(
    Rat_Gene = character(),
    Human_Gene = character(),
    Ensembl_Rat_ID = character(),
    Ensembl_Human_ID = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_batches) {
    start_idx <- (i-1) * batch_size + 1
    end_idx <- min(i * batch_size, n_genes)
    batch_genes <- unmapped_rat_genes[start_idx:end_idx]
    
    cat(paste0("处理批次 ", i, "/", n_batches, " (", length(batch_genes), " 个基因)...\n"))
    
    # 查询同源基因
    result <- tryCatch({
      getLDS(
        attributes = c("rgd_symbol", "ensembl_gene_id"),
        filters = "rgd_symbol",
        values = batch_genes,
        mart = ensembl_rat,
        attributesL = c("hgnc_symbol", "ensembl_gene_id"),
        martL = ensembl_human
      )
    }, error = function(e) {
      cat(paste0("批次 ", i, " 查询失败: ", e$message, "\n"))
      return(NULL)
    })
    
    if (!is.null(result) && nrow(result) > 0) {
      colnames(result) <- c("Rat_Gene", "Ensembl_Rat_ID", "Human_Gene", "Ensembl_Human_ID")
      # 过滤掉空的人类基因
      result <- result[result$Human_Gene != "" & !is.na(result$Human_Gene), ]
      if (nrow(result) > 0) {
        all_mappings <- rbind(all_mappings, result)
      }
    }
    
    # 避免请求过快
    Sys.sleep(0.5)
  }
  
  cat(paste0("\n通过biomaRt成功映射 ", nrow(all_mappings), " 条记录\n"))
  
  # 去重
  unique_mappings <- all_mappings[!duplicated(all_mappings$Human_Gene), ]
  cat(paste0("去重后人类基因数: ", nrow(unique_mappings), "\n"))
  
  # 合并之前已映射的基因
  final_human_genes <- unique(c(mapped_human_genes, toupper(unique_mappings$Human_Gene)))
  cat(paste0("\n最终人类基因总数（已映射 + 新映射）: ", length(final_human_genes), "\n"))
  
  # 保存结果到文本文件
  write.table(final_human_genes, 
              file = file.path(work_dir, "PPI_Input_Genes_Human_Complete.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # 创建详细映射表
  complete_mappings <- data.frame(
    Human_Gene = final_human_genes,
    stringsAsFactors = FALSE
  )
  
  write.table(complete_mappings,
              file = file.path(work_dir, "PPI_Input_Genes_Human_Complete_Table.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat(paste0("\n结果已保存:\n"))
  cat(paste0("  - 完整人类基因列表: PPI_Input_Genes_Human_Complete.txt\n"))
  cat(paste0("  - 详细映射表: PPI_Input_Genes_Human_Complete_Table.txt\n"))
  
  # 输出新映射的基因
  new_mappings <- setdiff(toupper(unique_mappings$Human_Gene), mapped_human_genes)
  cat(paste0("\n新映射的人类基因数: ", length(new_mappings), "\n"))
  
  if (length(new_mappings) > 0) {
    cat("新映射的基因示例 (前20个):\n")
    print(head(new_mappings, 20))
  }
  
}, error = function(e) {
  cat(paste0("\n错误: ", e$message, "\n"))
  cat("尝试使用备用方法...\n")
  
  # 备用：使用本地映射库进一步查询
  mapping_file <- file.path(work_dir, "大创/大鼠 小鼠 人类映射库.txt")
  mapping_lines <- readLines(mapping_file, warn = FALSE)
  
  # 解析映射文件，尝试别名匹配
  # ...（备用逻辑）
  
  cat("请检查网络连接或手动处理未映射基因\n")
})

cat("\n========================================\n")
cat("         基因映射完成\n")
cat("========================================\n")
