# 国内镜像设置----
options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

# 安装scTenifoldKnk----
if(!"BiocManager" %in% installed.packages()){install.packages('BiocManager')}
if(!"remotes" %in% installed.packages()){install.packages('remotes')}

# 加载scTenifoldKnk----
scTenifoldKnk_path <- "./scTenifoldKnk/scTenifoldKnk-master"
if(file.exists(scTenifoldKnk_path)) {
  if(!"scTenifoldKnk" %in% installed.packages()) {
    cat("安装scTenifoldKnk...\n")
    install.packages(scTenifoldKnk_path, repos = NULL, type = "source")
  }
}

library(scTenifoldKnk)
library(openxlsx)

# 读取表达矩阵----
sc_Matrix <- readRDS('./result/knk_input/sc_Matrix.rds')
cat(sprintf("表达矩阵: %d 基因 x %d 细胞\n", nrow(sc_Matrix), ncol(sc_Matrix)))

# BCP轴基因----
knk_genes <- c("Ager", "Nfkb1", "Fdx1")

# 结果保存目录----
dir_save <- './result/knk_results_new'
if(!dir.exists(dir_save)){dir.create(dir_save, recursive = T)}

# 批量虚拟敲除----
for (gKO in knk_genes) {
  cat(sprintf("\n=== 虚拟敲除: %s ===\n", gKO))

  set.seed(666)
  knk_res <- scTenifoldKnk(
    countMatrix = sc_Matrix,
    gKO = gKO,
    qc = TRUE,
    qc_mtThreshold = 0.1,
    qc_minLSize = 1000,
    nc_lambda = 0,
    nc_nNet = 10,
    nc_nCells = min(500, ncol(sc_Matrix)),
    nc_nComp = 3,
    nc_scaleScores = TRUE,
    nc_symmetric = FALSE,
    nc_q = 0.9,
    td_K = 3,
    td_maxIter = 1000,
    td_maxError = 1e-05,
    td_nDecimal = 3,
    ma_nDim = 2,
    nCores = 3
  )

  # 差异基因结果----
  diff_df <- knk_res$diffRegulation
  colnames(diff_df) <- c("gene", "FC", "log2FC", "p.value", "adj.p.value", "AUC", "delta")

  cat(sprintf("差异基因数: %d\n", sum(diff_df$p.value < 0.05 & abs(diff_df$log2FC) > 0.25)))

  # 保存完整结果----
  save_path <- paste0(dir_save, '/', gKO, '_diffRegulation.xlsx')
  write.xlsx(diff_df, save_path)
  cat(sprintf("结果保存至: %s\n", save_path))

  # 显著差异基因----
  sig_df <- diff_df[diff_df$p.value < 0.05 & abs(diff_df$log2FC) > 0.25, ]
  if (nrow(sig_df) > 0) {
    sig_path <- paste0(dir_save, '/', gKO, '_significant.xlsx')
    write.xlsx(sig_df, sig_path)
  }

  # BCP轴相关基因检查----
  bcp_genes <- c("Ager", "Nfkb1", "Fdx1", "Tlr4", "Stat1", "Stat3")
  axis_check <- diff_df[diff_df$gene %in% bcp_genes, ]
  if (nrow(axis_check) > 0) {
    cat("BCP轴基因变化:\n")
    print(axis_check)
    axis_path <- paste0(dir_save, '/', gKO, '_BCPaxis.xlsx')
    write.xlsx(axis_check, axis_path)
  }
}

cat("\n=== 虚拟敲除完成 ===\n")
