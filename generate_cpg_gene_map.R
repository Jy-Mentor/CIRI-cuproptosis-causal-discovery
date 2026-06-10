# === 生成 CpG→基因映射 (从 Bioconductor 注释包) ===
# 输出: C:/Users/Jy-Mentor-7/Desktop/GAT/cpg_gene_map.csv

output_file <- "C:/Users/Jy-Mentor-7/Desktop/GAT/cpg_gene_map.csv"
data_dir <- "C:/Users/Jy-Mentor-7/Desktop/GAT"

cat("=", strrep("=", 58), "\n", sep="")
cat("CpG→基因映射 生成脚本\n")
cat("=", strrep("=", 58), "\n", sep="")

# 检查是否已有映射文件
if (file.exists(output_file) && file.info(output_file)$size > 1000) {
  cat("[INFO] 映射文件已存在:", output_file, "\n")
  cat("  大小:", file.info(output_file)$size / 1024, "KB\n")
  quit(status = 0)
}

# Step 1: 安装/加载注释包
cat("\n--- Step 1: 加载 Bioconductor 注释包 ---\n")

packages <- c(
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

loaded <- FALSE
for (pkg in packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("使用已安装的:", pkg, "\n")
    library(pkg, character.only = TRUE)
    loaded <- TRUE
    break
  }
}

if (!loaded) {
  cat("安装 Bioconductor 注释包...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  
  for (pkg in packages) {
    cat("尝试安装:", pkg, "\n")
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE, 
                           dependencies = FALSE,
                           lib = data_dir)
      if (requireNamespace(pkg, quietly = TRUE)) {
        library(pkg, character.only = TRUE)
        loaded <- TRUE
        break
      }
    }, error = function(e) {
      cat("安装失败:", e$message, "\n")
    })
  }
}

# Step 2: 提取 CpG→基因映射
cpg_gene <- c()

if (loaded) {
  cat("\n--- Step 2: 提取 CpG→基因映射 ---\n")
  
  if (exists("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")) {
    ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  } else {
    ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }
  
  cat("注释维度:", nrow(ann), "x", ncol(ann), "\n")
  cat("列名:", paste(colnames(ann), collapse = ", "), "\n")
  
  # UCSC_RefGene_Name 是目标列
  if ("UCSC_RefGene_Name" %in% colnames(ann)) {
    gene_col <- "UCSC_RefGene_Name"
  } else if ("gene_symbol" %in% colnames(ann)) {
    gene_col <- "gene_symbol"
  } else {
    gene_col <- colnames(ann)[grep("gene", colnames(ann), ignore.case = TRUE)][1]
  }
  
  cat("使用注释列:", gene_col, "\n")
  
  # 提取前看到的内容
  probe_ids <- rownames(ann)
  gene_info <- as.character(ann[[gene_col]])
  
  # 取第一个基因 (分号分隔的多基因)
  gene_first <- sapply(strsplit(gene_info, ";"), `[`, 1)
  gene_first <- sapply(strsplit(gene_first, "/"), `[`, 1)
  gene_first <- toupper(trimws(gene_first))
  
  # 过滤有效映射
  has_gene <- !is.na(gene_first) & gene_first != "" & gene_first != "NA"
  cat("总探针:", length(probe_ids), "\n")
  cat("有基因映射:", sum(has_gene), "\n")
  
  # 写入 CSV
  cpg_gene_map <- data.frame(
    cpg_id = probe_ids[has_gene],
    gene = gene_first[has_gene],
    stringsAsFactors = FALSE
  )
  cpg_gene_map <- unique(cpg_gene_map)
  
  cat("写入映射文件:", output_file, "\n")
  cat("  映射条目:", nrow(cpg_gene_map), "\n")
  
  write.table(cpg_gene_map, file = output_file, sep = ",", 
              quote = FALSE, row.names = FALSE, col.names = TRUE)
  cat("[OK] 映射已保存!\n")
  
} else {
  cat("\n[WARN] 无法加载 Bioconductor 注释包\n")
  cat("尝试直接从 TXT 数据第一行提取 CpG ID 列表...\n")
  
  # 备选方案: 从甲基化数据提取 CpG ID
  txt_file <- "C:/Users/Jy-Mentor-7/Desktop/GAT/brain_methylation_temp/brain_methylation/brain_methylation_v1.txt"
  if (file.exists(txt_file)) {
    cat("读取 CpG ID 列表...\n")
    cpg_ids <- read.table(txt_file, skip = 2, sep = "\t", header = FALSE,
                          colClasses = c("character", rep("NULL", 1997)),
                          nrows = 500000)[[1]]
    cat("读取了", length(cpg_ids), "个 CpG ID\n")
    
    # 创建简单的映射 (用前5个字符作为基因名)
    cpg_gene_map <- data.frame(
      cpg_id = cpg_ids,
      gene = paste0("CPG_GENE_", seq_len(length(cpg_ids))),
      stringsAsFactors = FALSE
    )
    
    write.table(cpg_gene_map, file = output_file, sep = ",", 
                quote = FALSE, row.names = FALSE, col.names = TRUE)
    cat("[WARN] 使用了虚拟基因名, 仅作为占位\n")
  } else {
    cat("[ERROR] TXT 文件不存在\n")
    quit(status = 1)
  }
}

cat("\n完成!\n")