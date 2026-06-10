#!/usr/bin/env Rscript
# 将 GTEx parquet 转换为 CSV
library(arrow)
library(data.table)

cat("转换 GTEx parquet 为 CSV\n")

# 加载脑组织数据
cat("加载脑皮层数据...\n")
brain_eqtl <- read_parquet("C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Brain_Cortex.v11.eQTLs.signif_pairs.parquet")
cat(sprintf("✓ 脑皮层：%d 行\n", nrow(brain_eqtl)))

# 保存为 CSV
cat("保存为 CSV...\n")
fwrite(brain_eqtl, "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/brain_eqtl.csv")
cat("✓ 已保存：brain_eqtl.csv\n\n")

# 加载全血数据
cat("加载全血数据...\n")
blood_eqtl <- read_parquet("C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eQTLs.signif_pairs.parquet")
cat(sprintf("✓ 全血：%d 行\n", nrow(blood_eqtl)))

# 保存为 CSV
cat("保存为 CSV...\n")
fwrite(blood_eqtl, "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/blood_eqtl.csv")
cat("✓ 已保存：blood_eqtl.csv\n\n")

cat("完成！现在可以用 Python 处理 CSV 文件了\n")
