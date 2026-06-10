options("repos"= c(CRAN="https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor/")

options(future.globals.maxSize = 5000 * 1024^2)

if(!"dplyr" %in% installed.packages()){install.packages('dplyr')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
if(!"ggrepel" %in% installed.packages()){install.packages('ggrepel')}
if(!"openxlsx" %in% installed.packages()){install.packages('openxlsx')}
if(!"corpcor" %in% installed.packages()){install.packages('corpcor')}

library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(corpcor)

dir_data <- "./result"
dir_save <- "./result/knk_results"
if(!dir.exists(dir_save)){dir.create(dir_save, recursive = T)}

knk_genes <- c("Ager", "Nfkb1", "Fdx1")

sc_obj <- readRDS(file.path(dir_data, "sc_annotated.rds"))
sc_obj <- JoinLayers(sc_obj)

cat("=== 简化版虚拟敲除分析 ===\n")
cat("方法: 基于相关性网络模拟敲除效应\n")

cat("\n细胞类型分布:\n")
print(table(sc_obj$cell_type))

simple_virtual_knockout <- function(count_matrix, gene_of_interest, n_top_genes = 200) {
  gene_idx <- which(rownames(count_matrix) == gene_of_interest)
  if (length(gene_idx) == 0) {
    return(NULL)
  }

  target_genes <- rownames(count_matrix)
  target_genes <- target_genes[target_genes != gene_of_interest]

  gene_expr <- as.numeric(count_matrix[gene_idx, ])
  names(gene_expr) <- colnames(count_matrix)

  cor_scores <- sapply(target_genes, function(g) {
    g_expr <- as.numeric(count_matrix[g, ])
    suppressWarnings(cor(g_expr, gene_expr, method = "spearman"))
  })

  top_pos <- head(sort(cor_scores, decreasing = TRUE), 100)
  top_neg <- head(sort(cor_scores, decreasing = FALSE), 100)

  results <- data.frame(
    gene = c(names(top_pos), names(top_neg)),
    FC = c(top_pos / (mean(gene_expr) + 1), top_neg / (mean(gene_expr) + 1)),
    correlation = c(top_pos, top_neg),
    regulation = c(rep("Up", length(top_pos)), rep("Down", length(top_neg))),
    stringsAsFactors = FALSE
  )

  results$logFC <- log2(results$FC + 0.1)
  results$p.value <- 2 * pnorm(-abs(results$correlation))
  results$padj <- p.adjust(results$p.value, method = "BH")

  return(results)
}

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
    if (n_cells > 1500) {
      cells_ct <- subset(cells_ct, downsample = 1500)
      cat(sprintf("下采样至 %d 细胞\n", ncol(cells_ct)))
    }

    sc_Matrix <- SeuratObject::LayerData(cells_ct, assay = "RNA", layer = "counts")

    cat(sprintf("执行相关性分析 (矩阵: %d x %d)...\n", nrow(sc_Matrix), ncol(sc_Matrix)))

    tryCatch({
      knk_df <- simple_virtual_knockout(sc_Matrix, knk_gene)

      if (is.null(knk_df)) {
        cat("相关性分析失败\n")
        next
      }

      knk_df$gene <- rownames(knk_df)

      sig_df <- knk_df %>% dplyr::filter(padj < 0.05 & abs(correlation) > 0.1)

      save_path <- file.path(dir_save, paste0(cell_type, "_", knk_gene, "_knk.xlsx"))
      write.xlsx(sig_df, save_path)
      cat(sprintf("结果已保存: %s (%d 个显著相关基因)\n", basename(save_path), nrow(sig_df)))

      top_sig <- sig_df %>% top_n(30, wt = abs(correlation))

      p1 <- ggplot(sig_df, aes(x = reorder(gene, logFC), y = logFC)) +
        geom_bar(stat = 'identity', aes(fill = regulation), width = 0.7) +
        coord_flip() +
        scale_fill_manual(values = c("Up" = "#e74c3c", "Down" = "#3498db")) +
        labs(title = paste0("Genes correlated with ", knk_gene, " KO (", cell_type, ")"),
             x = "", y = "logFC", fill = "Regulation") +
        theme_bw() +
        theme(axis.text = element_text(size = 8))

      ggsave(plot = p1,
             device = 'png',
             width = 8,
             height = max(6, nrow(sig_df) * 0.12),
             dpi = 300,
             filename = file.path(dir_save, paste0(cell_type, "_", knk_gene, "_bar.png")))

      p2 <- ggplot(sig_df, aes(x = correlation, y = -log10(padj), color = regulation)) +
        geom_point(alpha = 0.7, size = 2) +
        scale_color_manual(values = c("Up" = "#e74c3c", "Down" = "#3498db")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
        ggrepel::geom_text_repel(data = top_sig, aes(label = gene),
                        size = 2, fontface = 'bold.italic', max.overlaps = 20) +
        labs(title = paste0(knk_gene, " KO in ", cell_type),
             x = "Spearman Correlation", y = "-log10(FDR)", color = "Regulation") +
        theme_bw()

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

cat("\n=== BCP轴效应分析 ===\n")
axis_effects <- list()

for (ct in levels(sc_obj$cell_type)) {
  for (gn in knk_genes) {
    f <- file.path(dir_save, paste0(ct, "_", gn, "_knk.xlsx"))
    if (file.exists(f)) {
      df <- read.xlsx(f)
      axis_effects[[paste0(ct, "_", gn)]] <- df
    }
  }
}

if (length(axis_effects) > 0) {
  axis_summary <- lapply(knk_genes, function(g) {
    gn_effects <- Filter(function(x) any(grepl(g, rownames(x), invert = FALSE)), axis_effects)
    if (length(gn_effects) > 0) {
      do.call(rbind, gn_effects)
    } else NULL
  })
  names(axis_summary) <- knk_genes

  for (g in knk_genes) {
    if (!is.null(axis_summary[[g]]) && nrow(axis_summary[[g]]) > 0) {
      out_path <- file.path(dir_save, paste0("BCP_axis_effect_", g, ".xlsx"))
      write.xlsx(axis_summary[[g]], out_path)
    }
  }
}

cat("分析完成！\n")
