# 解压GSE97537 series_matrix.gz文件
gz_file <- "D:/反向网络药理学/L1 数据集/bulk/GSE97537(24H)/GSE97537_series_matrix (1).txt.gz"
out_file <- "D:/反向网络药理学/L1 数据集/bulk/GSE97537(24H)/GSE97537_series_matrix.txt"

cat("Reading gz file...\n")
data <- gzfile(gz_file, "r")
lines <- readLines(data)
close(data)

cat("Writing uncompressed file...\n")
writeLines(lines, out_file)

cat("Done! Output:", out_file, "\n")
cat("Lines:", length(lines), "\n")
