# 设置工作目录
setwd("c:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\大创")

# 读取文件
bp_file <- "GO_BP_enrichment_results.tsv"
cc_file <- "GO_CC_enrichment_results.tsv"
mf_file <- "GO_MF_enrichment_results.tsv"
kegg_file <- "KEGG_enrichment_results.tsv"

df_bp <- read.delim(bp_file, sep="\t", stringsAsFactors=FALSE)
df_cc <- read.delim(cc_file, sep="\t", stringsAsFactors=FALSE)
df_mf <- read.delim(mf_file, sep="\t", stringsAsFactors=FALSE)
df_kegg <- read.delim(kegg_file, sep="\t", stringsAsFactors=FALSE)

# 过滤显著富集的条目 (p<0.05 且 FDR<0.05)
df_bp_sig <- df_bp[df_bp$pvalue < 0.05 & df_bp$p.adjust < 0.05, ]
df_cc_sig <- df_cc[df_cc$pvalue < 0.05 & df_cc$p.adjust < 0.05, ]
df_mf_sig <- df_mf[df_mf$pvalue < 0.05 & df_mf$p.adjust < 0.05, ]
df_kegg_sig <- df_kegg[df_kegg$pvalue < 0.05 & df_kegg$p.adjust < 0.05, ]

# 按 Adjusted P-value 排序
df_bp_sig_sorted <- df_bp_sig[order(df_bp_sig$p.adjust), ]
df_cc_sig_sorted <- df_cc_sig[order(df_cc_sig$p.adjust), ]
df_mf_sig_sorted <- df_mf_sig[order(df_mf_sig$p.adjust), ]
df_kegg_sig_sorted <- df_kegg_sig[order(df_kegg_sig$p.adjust), ]

# 提取 Top 条目
top3_bp <- head(df_bp_sig_sorted$Description, 3)
top4_cc <- head(df_cc_sig_sorted$Description, 4)
top4_mf <- head(df_mf_sig_sorted$Description, 4)
top4_kegg <- head(df_kegg_sig_sorted$Description, 4)

# 输出结果
cat("Top 3 GO BP 显著富集条目:\n")
for (i in 1:length(top3_bp)) {
  cat(paste(i, ". ", top3_bp[i], "\n", sep=""))
}

cat("\nTop 4 GO CC 显著富集条目:\n")
for (i in 1:length(top4_cc)) {
  cat(paste(i, ". ", top4_cc[i], "\n", sep=""))
}

cat("\nTop 4 GO MF 显著富集条目:\n")
for (i in 1:length(top4_mf)) {
  cat(paste(i, ". ", top4_mf[i], "\n", sep=""))
}

cat("\nTop 4 KEGG 显著富集通路:\n")
for (i in 1:length(top4_kegg)) {
  cat(paste(i, ". ", top4_kegg[i], "\n", sep=""))
}

# 保存结果到文件
output_file <- "../enrichment_terms_summary.txt"

cat("Top 3 GO BP 显著富集条目:\n", file=output_file)
cat(top3_bp, file=output_file, append=TRUE, sep="\n")

cat("\nTop 4 GO CC 显著富集条目:\n", file=output_file, append=TRUE)
cat(top4_cc, file=output_file, append=TRUE, sep="\n")

cat("\nTop 4 GO MF 显著富集条目:\n", file=output_file, append=TRUE)
cat(top4_mf, file=output_file, append=TRUE, sep="\n")

cat("\nTop 4 KEGG 显著富集通路:\n", file=output_file, append=TRUE)
cat(top4_kegg, file=output_file, append=TRUE, sep="\n")

cat("\n结果已保存到 enrichment_terms_summary.txt 文件\n")