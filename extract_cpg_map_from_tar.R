# === 从 Bioconductor tarball 提取 CpG→基因映射 ===
# 输入: C:/Users/Jy-Mentor-7/Desktop/GAT/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz
# 输出: C:/Users/Jy-Mentor-7/Desktop/GAT/cpg_gene_map.csv

tar_file <- "C:/Users/Jy-Mentor-7/Desktop/GAT/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz"
output_csv <- "C:/Users/Jy-Mentor-7/Desktop/GAT/cpg_gene_map.csv"

cat("=", strrep("=", 58), "\n", sep="")
cat("从 tarball 提取 CpG→基因映射\n")
cat("=", strrep("=", 58), "\n", sep="")

if (!file.exists(tar_file)) {
  cat("[ERROR] 未找到:", tar_file, "\n")
  quit(status = 1)
}

cat("[INFO] 提取 tarball...\n")
tmpdir <- tempdir()
untar(tar_file, exdir = tmpdir)
cat("[INFO] 提取到:", tmpdir, "\n")

# 查找 data 目录下的 rda 文件
rda_files <- list.files(tmpdir, pattern = "\\.rda$", recursive = TRUE, full.names = TRUE)
cat("[INFO] 找到 RDA 文件:", length(rda_files), "\n")
for (f in rda_files) {
  cat("  ", basename(f), " - ", file.info(f)$size / 1024 / 1024, " MB\n")
}

# 查找包含注释的主文件
main_rda <- rda_files[grep("IlluminaHumanMethylation450kanno", rda_files, ignore.case = TRUE)]
if (length(main_rda) == 0) {
  main_rda <- rda_files[1]
}

cat("\n[INFO] 加载:", basename(main_rda[1]), "\n")
env <- new.env()
load(main_rda[1], envir = env)
cat("  对象:", paste(ls(env), collapse = ", "), "\n")

# 查找 Other 数据 (包含 UCSC_RefGene_Name)
other_rda <- rda_files[grep("Other", basename(rda_files))]
if (length(other_rda) > 0) {
  cat("\n[INFO] 加载 Other.rda...\n")
  env2 <- new.env()
  load(other_rda[1], envir = env2)
  cat("  对象:", paste(ls(env2), collapse = ", "), "\n")
  
  obj_name <- ls(env2)[1]
  other_data <- get(obj_name, envir = env2)
  cat("  类:", class(other_data), "\n")
  cat("  维度:", nrow(other_data), "x", ncol(other_data), "\n")
  cat("  列名:", paste(colnames(other_data), collapse = ", "), "\n")
  
  if ("UCSC_RefGene_Name" %in% colnames(other_data)) {
    cat("\n[INFO] 找到 UCSC_RefGene_Name!\n")
    probes <- rownames(other_data)
    genes <- as.character(other_data$UCSC_RefGene_Name)
    
    # 取第一个基因
    genes_first <- sapply(strsplit(genes, ";"), `[`, 1)
    genes_first <- sapply(strsplit(genes_first, "/"), `[`, 1)
    genes_first <- toupper(trimws(genes_first))
    
    has_gene <- !is.na(genes_first) & genes_first != "" & genes_first != "NA"
    cat("  总探针:", length(probes), "\n")
    cat("  有基因映射:", sum(has_gene), "\n")
    
    cpg_gene_map <- data.frame(
      cpg_id = probes[has_gene],
      gene = genes_first[has_gene],
      stringsAsFactors = FALSE
    )
    cpg_gene_map <- unique(cpg_gene_map)
    
    cat("\n[SAVE] 写入:", output_csv, "\n")
    write.table(cpg_gene_map, file = output_csv, sep = ",", 
                quote = FALSE, row.names = FALSE, col.names = TRUE)
    cat("[OK] 保存完成! 条目:", nrow(cpg_gene_map), "\n")
    
    # 统计唯一基因数
    cat("  唯一基因:", length(unique(cpg_gene_map$gene)), "\n")
    cat("  示例:\n")
    print(head(cpg_gene_map, 10))
  } else {
    cat("[ERROR] 未找到 UCSC_RefGene_Name 列\n")
    cat("  可用列:", paste(colnames(other_data), collapse = ", "), "\n")
  }
} else {
  cat("[ERROR] 未找到 Other.rda\n")
}

cat("\n完成!\n")