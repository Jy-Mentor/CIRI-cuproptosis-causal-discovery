#!/usr/bin/env Rscript

# 设置编码为UTF-8
options(encoding = "UTF-8")

# 测试读取本地文件（使用英文路径）
test_file_eng <- "C:/Users/Jy-Mentor-7/Desktop/test/test.txt"

# 创建测试目录和文件
if (!dir.exists(dirname(test_file_eng))) {
  dir.create(dirname(test_file_eng), recursive = TRUE)
}
writeLines("Hello, World!", test_file_eng)

cat("尝试读取英文路径文件: ", test_file_eng, "\n")

tryCatch({
  content <- readLines(test_file_eng)
  cat("成功读取内容:\n")
  print(content)
}, error = function(e) {
  cat("读取失败: ", e$message, "\n")
})

# 测试读取10x数据（使用英文路径）
cat("\n尝试读取10x数据...\n")
tryCatch({
  library(Seurat)
  # 复制10x文件到英文路径
  source_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE/GSM5319987_sham1"
  target_dir <- "C:/Users/Jy-Mentor-7/Desktop/test/GSM5319987_sham1"
  
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  
  # 复制文件
  files <- c("barcodes.tsv", "genes.tsv", "matrix.mtx")
  for (file in files) {
    source_file <- file.path(source_dir, file)
    target_file <- file.path(target_dir, file)
    if (file.exists(source_file)) {
      file.copy(source_file, target_file, overwrite = TRUE)
      cat("复制文件: ", file, "\n")
    }
  }
  
  # 尝试读取
  counts <- Read10X(target_dir)
  cat("成功读取，基因数: ", nrow(counts), " 细胞数: ", ncol(counts), "\n")
}, error = function(e) {
  cat("读取失败: ", e$message, "\n")
})

# 清理测试文件
unlink("C:/Users/Jy-Mentor-7/Desktop/test", recursive = TRUE)
