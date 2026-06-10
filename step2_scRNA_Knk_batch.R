options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

options(future.globals.maxSize = 5000 * 1024^2)

if(!"BiocManager" %in% installed.packages()){install.packages('BiocManager')}
if(!"remotes" %in% installed.packages()){install.packages('remotes')}

scTenifoldKnk_path <- "scTenifoldKnk/scTenifoldKnk-master"
if(!file.exists(scTenifoldKnk_path)) {
  stop("scTenifoldKnk包路径不存在: ", scTenifoldKnk_path)
}

if(!"scTenifoldNet" %in% installed.packages()) {
  cat("安装scTenifoldNet...\n")
  remotes::install_github('cailab-tamu/scTenifoldNet')
}

if(!"scTenifoldKnk" %in% installed.packages()) {
  cat("安装scTenifoldKnk从本地...\n")
  install.packages(scTenifoldKnk_path, repos = NULL, type = "source")
}

if(!"dplyr" %in% installed.packages()){install.packages('dplyr')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
if(!"ggrepel" %in% installed.packages()){install.packages('ggrepel')}
if(!"openxlsx" %in% installed.packages()){install.packages('openxlsx')}
if(!"parallel" %in% installed.packages()){install.packages('parallel')}

library(scTenifoldKnk)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(parallel)

dir_data <- "./result"
dir_save <- "./result/knk_results"
if(!dir.exists(dir_save)){dir.create(dir_save, recursive = T)}

knk_genes <- c("Ager", "Nfkb1", "Fdx1")

sc_obj <- readRDS(file.path(dir_data, "sc_annotated.rds"))
sc_obj <- JoinLayers(sc_obj)

cat("=== 准备虚拟敲除数据 ===\n")
cat("细胞类型分布:\n")
print(table(sc_obj$cell_type))

for (knk_gene in knk_genes) {
  cat(sprintf("\n=== 虚拟敲除: %s ===\n", knk_gene))

  for (cell_type in levels(sc_obj$cell_type)) {
    cat(sprintf("\n--- 处理 %s 中的 %s ---\n", cell_type, knk_gene))

    cells_ct <- subset(sc_obj, cell_type == !!cell_type)
    n_cells <- ncol(cells_ct)
    cat(sprintf("该类型细胞数: %d\n", n_cells))

    if (n_cells < 100) {
      cat("细胞数不足100，跳过\n")
      next
    }

    if (!knk_gene %in% rownames(cells_ct)) {
      cat(sprintf("%s 不在数据中，跳过\n", knk_gene))
      next
    }

    set.seed(666)
    if (n_cells > 2000) {
      cells_ct <- subset(cells_ct, downsample = 2000)
      cat(sprintf("下采样至 %d 细胞\n", ncol(cells_ct)))
    }

    cat("提取表达矩阵...\n")
    sc_Matrix <- SeuratObject::LayerData(cells_ct, assay = "RNA", layer = "counts")

    if (!knk_gene %in% rownames(sc_Matrix)) {
      cat(sprintf("%s 不在count矩阵中，跳过\n", knk_gene))
      next
    }

    cat(sprintf("执行虚拟敲除: %s (矩阵维度: %d x %d)...\n",
                knk_gene, nrow(sc_Matrix), ncol(sc_Matrix)))

    tryCatch({
      set.seed(666)
      knk_res <- scTenifoldKnk(
        countMatrix = sc_Matrix,
        gKO = knk_gene,
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
        nCores = parallel::detectCores() - 1
      )

      knk_df <- knk_res$diffRegulation %>%
        dplyr::mutate(logFC = log2(knk_res$diffRegulation$FC)) %>%
        dplyr::filter(gene != knk_gene) %>%
        dplyr::filter(p.value < 0.05)

      if (nrow(knk_df) > 0) {
        knk_df[, 2:7] <- sapply(knk_df[, 2:7], as.numeric)
      }

      save_path <- file.path(dir_save, paste0(cell_type, "_", knk_gene, "_knk.xlsx"))
      write.xlsx(knk_df, save_path)
      cat(sprintf("结果已保存: %s (%d 个差异基因)\n", basename(save_path), nrow(knk_df)))

      p1 <- ggplot(knk_df, aes(x = reorder(gene, logFC), y = logFC)) +
        geom_bar(stat = 'identity', fill = 'orange') +
        coord_flip() +
        labs(title = paste0("Differentially Regulated Genes after ", knk_gene, " KO (", cell_type, ")"),
             x = "", y = "logFC") +
        theme_bw() +
        theme(axis.text = element_text(size = 8))

      ggsave(plot = p1,
             device = 'png',
             width = 8,
             height = max(6, nrow(knk_df) * 0.15),
             dpi = 300,
             filename = file.path(dir_save, paste0(cell_type, "_", knk_gene, "_bar.png")))

      df <- knk_res$diffRegulation %>%
        dplyr::mutate(label_genes = ifelse(gene %in% knk_df$gene, gene, NA))

      p2 <- ggplot(df, aes(x = log2(FC), y = -log10(p.value))) +
        geom_point(alpha = ifelse(df$gene %in% df$label_genes, 0.9, 0.4),
                   fill = ifelse(df$gene %in% df$label_genes, '#ffd6a5', 'gray'),
                   shape = 21, color = 'white',
                   aes(size = -log10(p.value))) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        geom_text_repel(data = df, aes(label = label_genes),
                        size = 2, parse = F, fontface = 'bold.italic',
                        color = '#c7522a', max.overlaps = 30) +
        labs(title = paste0(knk_gene, " KO in ", cell_type),
             x = "logFC", y = "-log10(Pvalue)") +
        theme_bw() +
        guides(size = 'none')

      ggsave(plot = p2,
             device = 'png',
             width = 6,
             height = 6,
             dpi = 300,
             filename = file.path(dir_save, paste0(cell_type, "_", knk_gene, "_volcano.png")))

      cat(sprintf("可视化已保存\n"))

    }, error = function(e) {
      cat(sprintf("错误: %s\n", conditionMessage(e)))
    })
  }
}

cat("\n=== 虚拟敲除分析完成 ===\n")
cat(sprintf("结果保存在: %s\n", dir_save))
