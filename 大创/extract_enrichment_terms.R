# 安装并加载必要的包
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
library(readxl)

# 定义文件路径
dir_path <- "c:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创\\nb\\"

# 读取GO BP结果
go_bp_file <- paste0(dir_path, "GO_BP_enrichment_results.xlsx")
go_bp <- read_excel(go_bp_file)

# 筛选显著富集的条目 (p < 0.05 且 FDR < 0.05)
go_bp_significant <- go_bp[go_bp$pvalue < 0.05 & go_bp$`p.adjust` < 0.05, ]

# 按Adjusted P-value排序
go_bp_significant <- go_bp_significant[order(go_bp_significant$`p.adjust`), ]

# 提取Top3条目
go_bp_top3 <- head(go_bp_significant, 3)

# 读取GO CC结果
go_cc_file <- paste0(dir_path, "GO_CC_enrichment_results.xlsx")
go_cc <- read_excel(go_cc_file)

# 筛选显著富集的条目
go_cc_significant <- go_cc[go_cc$pvalue < 0.05 & go_cc$`p.adjust` < 0.05, ]

# 按Adjusted P-value排序
go_cc_significant <- go_cc_significant[order(go_cc_significant$`p.adjust`), ]

# 提取Top4条目
go_cc_top4 <- head(go_cc_significant, 4)

# 读取GO MF结果
go_mf_file <- paste0(dir_path, "GO_MF_enrichment_results.xlsx")
go_mf <- read_excel(go_mf_file)

# 筛选显著富集的条目
go_mf_significant <- go_mf[go_mf$pvalue < 0.05 & go_mf$`p.adjust` < 0.05, ]

# 按Adjusted P-value排序
go_mf_significant <- go_mf_significant[order(go_mf_significant$`p.adjust`), ]

# 提取Top4条目
go_mf_top4 <- head(go_mf_significant, 4)

# 读取KEGG结果
kegg_file <- paste0(dir_path, "KEGG_enrichment_results.xlsx")
kegg <- read_excel(kegg_file)

# 筛选显著富集的通路
kegg_significant <- kegg[kegg$pvalue < 0.05 & kegg$`p.adjust` < 0.05, ]

# 按Adjusted P-value排序
kegg_significant <- kegg_significant[order(kegg_significant$`p.adjust`), ]

# 提取Top4通路
kegg_top4 <- head(kegg_significant, 4)

# 输出结果
cat("\n===== GO Biological Process Top 3 =====\n")
print(go_bp_top3[, c("ID", "Description", "pvalue", "p.adjust")])

cat("\n===== GO Cellular Component Top 4 =====\n")
print(go_cc_top4[, c("ID", "Description", "pvalue", "p.adjust")])

cat("\n===== GO Molecular Function Top 4 =====\n")
print(go_mf_top4[, c("ID", "Description", "pvalue", "p.adjust")])

cat("\n===== KEGG Pathway Top 4 =====\n")
print(kegg_top4[, c("ID", "Description", "pvalue", "p.adjust")])

# 保存结果为CSV文件
write.csv(go_bp_top3, paste0(dir_path, "GO_BP_top3.csv"), row.names = FALSE)
write.csv(go_cc_top4, paste0(dir_path, "GO_CC_top4.csv"), row.names = FALSE)
write.csv(go_mf_top4, paste0(dir_path, "GO_MF_top4.csv"), row.names = FALSE)
write.csv(kegg_top4, paste0(dir_path, "KEGG_top4.csv"), row.names = FALSE)

cat("\n结果已保存到CSV文件中。\n")
