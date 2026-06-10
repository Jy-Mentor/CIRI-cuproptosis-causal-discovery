#!/usr/bin/env Rscript

# 测试文件读取权限
setwd("C:/Users/Jy-Mentor-7/Desktop")

# 测试文件路径
test_file <- "GSE/GSM5319987_sham1/barcodes.tsv"

cat("测试文件读取权限...\n")
cat("文件路径: ", test_file, "\n")

# 测试1: 直接读取
tryCatch({
  cat("\n1. 尝试直接读取...\n")
  content <- readLines(test_file)
  cat("成功读取! 前5行:\n")
  print(head(content, 5))
}, error = function(e) {
  cat("失败: ", e$message, "\n")
})

# 测试2: 使用file()函数
tryCatch({
  cat("\n2. 尝试使用file()函数...\n")
  con <- file(test_file, "r")
  content <- readLines(con, n = 5)
  close(con)
  cat("成功读取! 前5行:\n")
  print(content)
}, error = function(e) {
  cat("失败: ", e$message, "\n")
})

# 测试3: 使用scan()函数
tryCatch({
  cat("\n3. 尝试使用scan()函数...\n")
  content <- scan(test_file, what = "character", n = 5, sep = "\n")
  cat("成功读取! 前5行:\n")
  print(content)
}, error = function(e) {
  cat("失败: ", e$message, "\n")
})

# 测试4: 检查文件权限
cat("\n4. 检查文件权限...\n")
cat("文件存在: ", file.exists(test_file), "\n")
cat("文件可读: ", file.access(test_file, mode = 4) == 0, "\n")
cat("文件可写: ", file.access(test_file, mode = 2) == 0, "\n")
cat("文件可执行: ", file.access(test_file, mode = 1) == 0, "\n")

# 测试5: 尝试使用绝对路径
tryCatch({
  cat("\n5. 尝试使用绝对路径...\n")
  absolute_path <- "C:/Users/Jy-Mentor-7/Desktop/GSE/GSM5319987_sham1/barcodes.tsv"
  content <- readLines(absolute_path)
  cat("成功读取! 前5行:\n")
  print(head(content, 5))
}, error = function(e) {
  cat("失败: ", e$message, "\n")
})
