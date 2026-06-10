# 检查Excel文件结构
library(readxl)

bulk_file <- "C:/Users/Jy-Mentor-7/Downloads/GSE163614_mRNA_Expression_Profiling.xlsx"

# 读取前10行查看结构
bulk_df <- read_excel(bulk_file, n_max=10)
print("列名:")
print(colnames(bulk_df))
print("\n前5行:")
print(head(bulk_df, 5))
print("\n数据维度:")
print(dim(bulk_df))
