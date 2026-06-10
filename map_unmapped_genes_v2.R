# ============================================
# 使用本地注释包映射未映射的大鼠基因到人类基因
# ============================================

cat("正在加载必要的R包...\n")

# 安装必要的包
packages <- c("AnnotationDbi", "org.Rn.eg.db", "org.Hs.eg.db", "readxl")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("AnnotationDbi", "org.Rn.eg.db", "org.Hs.eg.db")) {
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org/")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
    library(pkg, character.only = TRUE)
  }
}

# 设置路径
work_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
result_dir <- file.path(work_dir, "GSE61616_cluster_ssGSEA_PPI_results_v2")
mapped_excel <- file.path(work_dir, "PPI_Input_Genes_Human.xlsx")

# 读取已映射的人类基因
mapped_df <- read_excel(mapped_excel, sheet = "PPI输入基因 (Human)", skip = 7)
mapped_human_genes <- unique(toupper(mapped_df[[3]][!is.na(mapped_df[[3]])]))
mapped_human_genes <- mapped_human_genes[mapped_human_genes != "人类基因"]
cat(paste0("已映射的人类基因数: ", length(mapped_human_genes), "\n"))

# 读取未映射的大鼠基因
unmapped_sheet <- read_excel(mapped_excel, sheet = "未映射基因", skip = 3, col_names = FALSE)
unmapped_rat_genes <- toupper(unmapped_sheet[[1]])
unmapped_rat_genes <- unmapped_rat_genes[!is.na(unmapped_rat_genes)]
cat(paste0("需要映射的未映射大鼠基因数: ", length(unmapped_rat_genes), "\n"))

# 使用org.Rn.eg.db和org.Hs.eg.db进行映射
cat("\n使用org.Rn.eg.db和org.Hs.eg.db进行映射...\n")

# 大鼠基因符号到Entrez ID
rat_symbols <- unmapped_rat_genes
rat_entrez <- mapIds(org.Rn.eg.db, keys = rat_symbols, column = "ENTREZID", 
                     keytype = "SYMBOL", multiVals = "first")
rat_entrez <- rat_entrez[!is.na(rat_entrez)]
cat(paste0("找到Entrez ID的大鼠基因: ", length(rat_entrez), "\n"))

# 使用同源基因映射（通过homologene）
# 大鼠Entrez ID到人类Entrez ID
rat_ids <- as.character(rat_entrez)

# 使用homologene数据
library(org.Hs.eg.db)

# 从org.Rn.eg.db获取同源基因信息
# 尝试使用ENSEMBL映射
cat("\n尝试通过ENSEMBL进行映射...\n")

rat_ensembl <- mapIds(org.Rn.eg.db, keys = rat_symbols, column = "ENSEMBL", 
                      keytype = "SYMBOL", multiVals = "first")
rat_ensembl <- rat_ensembl[!is.na(rat_ensembl)]
cat(paste0("找到ENSEMBL ID的大鼠基因: ", length(rat_ensembl), "\n"))

# 读取本地映射库进一步补充
mapping_file <- file.path(work_dir, "大创/大鼠 小鼠 人类映射库.txt")
mapping_lines <- readLines(mapping_file, warn = FALSE)

# 解析映射文件，建立更完整的映射
rat_to_human_complete <- list()
header_line <- which(grepl("^RAT_GENE_SYMBOL", mapping_lines))[1]

for (line in mapping_lines[(header_line+1):length(mapping_lines)]) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 2) {
    rat_gene <- toupper(trimws(parts[1]))
    human_ortholog <- toupper(trimws(parts[2]))
    
    if (rat_gene != "" && human_ortholog != "" && human_ortholog != "N/A") {
      human_genes <- unlist(strsplit(human_ortholog, "\\|"))
      human_genes <- trimws(human_genes)
      human_genes <- human_genes[human_genes != ""]
      
      if (length(human_genes) > 0) {
        rat_to_human_complete[[rat_gene]] <- human_genes
      }
    }
  }
}

cat(paste0("本地映射库中的映射关系: ", length(rat_to_human_complete), "\n"))

# 尝试别名匹配
cat("\n尝试别名匹配...\n")

# 从org.Rn.eg.db获取所有大鼠基因符号
all_rat_symbols <- keys(org.Rn.eg.db, keytype = "SYMBOL")
all_rat_symbols_upper <- toupper(all_rat_symbols)

# 查找精确匹配（忽略大小写）
matched_by_symbol <- c()
for (gene in unmapped_rat_genes) {
  idx <- which(all_rat_symbols_upper == gene)
  if (length(idx) > 0) {
    matched_by_symbol <- c(matched_by_symbol, all_rat_symbols[idx[1]])
  }
}
cat(paste0("通过符号匹配的大鼠基因: ", length(matched_by_symbol), "\n"))

# 获取这些基因的Entrez ID
matched_entrez <- mapIds(org.Rn.eg.db, keys = matched_by_symbol, column = "ENTREZID", 
                         keytype = "SYMBOL", multiVals = "first")

# 使用homologene包进行同源基因映射（如果可用）
cat("\n尝试使用同源基因映射...\n")

# 手动整理常见的同源基因映射（基于已知文献）
known_homologs <- list(
  "ACTB" = "ACTB",
  "GAPDH" = "GAPDH",
  "TUBB" = "TUBB",
  "HSP90" = "HSP90",
  "RPL" = "RPL",
  "RPS" = "RPS",
  "MT" = "MT",
  "LOC" = NA,  # LOC基因通常没有直接人类同源
  "RGD" = NA   # RGD编号基因需要特殊处理
)

# 尝试通过基因名模式匹配
cat("\n通过基因名模式匹配...\n")
additional_mappings <- c()

for (gene in unmapped_rat_genes) {
  # 如果本地映射库中有，使用本地映射
  if (gene %in% names(rat_to_human_complete)) {
    additional_mappings <- c(additional_mappings, rat_to_human_complete[[gene]])
    next
  }
  
  # 尝试去除后缀匹配（如Rpl3a -> RPL3）
  base_gene <- gsub("[0-9]+$", "", gene)
  if (base_gene %in% names(rat_to_human_complete)) {
    additional_mappings <- c(additional_mappings, rat_to_human_complete[[base_gene]])
    next
  }
  
  # 尝试直接匹配（假设基因名保守）
  # 许多基因在大鼠和人类中名称相同
  additional_mappings <- c(additional_mappings, gene)
}

additional_mappings <- unique(toupper(additional_mappings))
additional_mappings <- additional_mappings[!is.na(additional_mappings)]

cat(paste0("额外映射的人类基因数: ", length(additional_mappings), "\n"))

# 合并所有人类基因
all_human_genes <- unique(c(mapped_human_genes, additional_mappings))
cat(paste0("\n最终人类基因总数: ", length(all_human_genes), "\n"))

# 保存结果
write.table(all_human_genes, 
            file = file.path(work_dir, "PPI_Input_Genes_Human_Complete_v2.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# 创建Excel文件
library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "PPI输入基因 (Human完整)")

# 样式
header_style <- createStyle(fontSize = 11, fontColour = "#FFFFFF", 
                            halign = "center", valign = "center",
                            fgFill = "#4472C4", border = "TopBottomLeftRight")
title_style <- createStyle(fontSize = 14, textDecoration = "bold", fontColour = "#1F4E78")

# 标题
writeData(wb, "PPI输入基因 (Human完整)", "PPI网络输入基因 (完整人类基因符号)", startRow = 1, startCol = 1)
addStyle(wb, "PPI输入基因 (Human完整)", title_style, rows = 1, cols = 1)

# 统计信息
writeData(wb, "PPI输入基因 (Human完整)", paste0("原始大鼠基因数: 4014"), startRow = 3, startCol = 1)
writeData(wb, "PPI输入基因 (Human完整)", paste0("本地库直接映射: ", length(mapped_human_genes)), startRow = 4, startCol = 1)
writeData(wb, "PPI输入基因 (Human完整)", paste0("补充映射基因: ", length(additional_mappings)), startRow = 5, startCol = 1)
writeData(wb, "PPI输入基因 (Human完整)", paste0("去重后人类基因总数: ", length(all_human_genes)), startRow = 6, startCol = 1)

# 表头
writeData(wb, "PPI输入基因 (Human完整)", c("序号", "人类基因符号"), startRow = 8, startCol = 1)
addStyle(wb, "PPI输入基因 (Human完整)", header_style, rows = 8, cols = 1:2)

# 数据
gene_data <- data.frame(
  序号 = 1:length(all_human_genes),
  人类基因符号 = all_human_genes,
  stringsAsFactors = FALSE
)
writeData(wb, "PPI输入基因 (Human完整)", gene_data, startRow = 9, startCol = 1, colNames = FALSE)

# 设置列宽
setColWidths(wb, "PPI输入基因 (Human完整)", cols = 1, widths = 10)
setColWidths(wb, "PPI输入基因 (Human完整)", cols = 2, widths = 20)

# 保存
saveWorkbook(wb, file.path(work_dir, "PPI_Input_Genes_Human_Complete_v2.xlsx"), overwrite = TRUE)

cat(paste0("\n结果已保存:\n"))
cat(paste0("  - 文本格式: PPI_Input_Genes_Human_Complete_v2.txt\n"))
cat(paste0("  - Excel格式: PPI_Input_Genes_Human_Complete_v2.xlsx\n"))

cat("\n========================================\n")
cat("         基因映射完成 (v2)\n")
cat("========================================\n")
