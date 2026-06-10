# PC因果推断分析 - 针对30个炎症/信号通路基因

library(Seurat)
library(pcalg)
library(data.table)

target_genes <- c(
  "IL6", "STAT3", "NFKB1", "PPARG", "TGFB1", "CCL2", "TLR4", "PTGS2",
  "CCND1", "STAT1", "ICAM1", "PTPRC", "RELA", "CASP8", "CXCR4",
  "NOTCH1", "MAPK1", "MDM2", "HSPA5", "PARP1", "JAK1", "CREBBP",
  "MMP2", "SREBF1", "CDC42", "STAT5A", "NFE2L2", "IRF1", "HMOX1"
)

cat("=== PC Causal Inference Analysis ===\n")

cache_file <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/analysis_cache.rds"
cache_data <- readRDS(cache_file)
neurons <- cache_data$neurons

cat(sprintf("Seurat object: %d cells x %d genes\n", ncol(neurons), nrow(neurons)))

all_genes <- rownames(neurons)
target_upper <- toupper(target_genes)
genes_upper <- toupper(all_genes)
matched_idx <- which(genes_upper %in% target_upper)
matched_genes <- all_genes[matched_idx]
matched_genes_upper <- target_genes[target_upper %in% genes_upper]

cat(sprintf("\nMatched %d/%d target genes\n", length(matched_genes), length(target_genes)))

# 提取表达矩阵
cat("\nExtracting expression matrix (Seurat v5)...\n")
expr_data <- neurons[["RNA"]]$data
cat(sprintf("Raw assay data: %d genes x %d cells\n", nrow(expr_data), ncol(expr_data)))

# 提取目标基因
target_exp <- expr_data[matched_genes, , drop = FALSE]
cat(sprintf("Target genes extracted: %d x %d\n", nrow(target_exp), ncol(target_exp)))

# 转置
target_exp <- t(as.matrix(target_exp))
cat(sprintf("Transposed: %d cells x %d genes\n", nrow(target_exp), ncol(target_exp)))

# 过滤低表达基因
gene_means <- colMeans(target_exp)
gene_counts <- colSums(target_exp > 0)
min_cells <- as.integer(max(5, 0.005 * nrow(target_exp)))

cat(sprintf("\nGene filtering parameters:\n"))
cat(sprintf("  Min cells threshold: %d\n", min_cells))
cat(sprintf("  Gene counts range: %d - %d\n", min(gene_counts), max(gene_counts)))

high_genes_idx <- gene_counts >= min_cells & gene_means > 0.01
high_genes <- colnames(target_exp)[high_genes_idx]
target_exp <- target_exp[, high_genes, drop = FALSE]

cat(sprintf("After filtering: %d genes retained\n", ncol(target_exp)))

if (ncol(target_exp) < 3) {
  stop("Insufficient genes after filtering")
}

# Z-score标准化
target_exp_scaled <- scale(target_exp)
target_exp_scaled[is.na(target_exp_scaled)] <- 0

# PC算法
cat("\n=== Running PC Algorithm ===\n")
corr_matrix <- cor(target_exp_scaled)
n_samples <- nrow(target_exp_scaled)

cat(sprintf("Correlation matrix: %d x %d\n", nrow(corr_matrix), ncol(corr_matrix)))
cat(sprintf("Sample size: %d\n", n_samples))

# 尝试不同的alpha值
found_edges <- FALSE
for (alpha_val in c(0.01, 0.05, 0.1)) {
  cat(sprintf("\nTrying alpha = %g...\n", alpha_val))
  
  suffStat <- list(C = corr_matrix, n = n_samples)
  indepTest <- gaussCItest
  
  pc_result <- tryCatch({
    pc(suffStat, indepTest = indepTest, alpha = alpha_val,
       labels = colnames(target_exp_scaled), verbose = FALSE)
  }, error = function(e) {
    cat(sprintf("  PC failed: %s\n", conditionMessage(e)))
    NULL
  })
  
  if (!is.null(pc_result)) {
    adj_matrix <- as(pc_result@graph, "matrix")
    cat(sprintf("  Adjacency matrix sum: %d\n", sum(adj_matrix)))
    
    n <- nrow(adj_matrix)
    edge_list <- data.frame()
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        if (adj_matrix[i,j] == 1 || adj_matrix[j,i] == 1) {
          from <- colnames(adj_matrix)[i]
          to <- colnames(adj_matrix)[j]
          
          if (adj_matrix[i,j] == 1 && adj_matrix[j,i] == 1) {
            direction <- "<->"
          } else if (adj_matrix[i,j] == 1) {
            direction <- "->"
          } else if (adj_matrix[j,i] == 1) {
            direction <- "->"
          } else {
            direction <- "?->"
          }
          
          edge_list <- rbind(edge_list, data.frame(
            from = from, to = to, direction = direction, stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    cat(sprintf("  Found %d edges\n", nrow(edge_list)))
    
    if (nrow(edge_list) > 0) {
      cat("\n=== PC Causal Network Results (alpha =", alpha_val, ") ===\n")
      cat(sprintf("Total edges: %d\n", nrow(edge_list)))
      print(edge_list)
      
      # 保存结果
      output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创"
      write.csv(edge_list, file.path(output_dir, "PC_30genes_pc_edges.csv"), row.names = FALSE)
      write.csv(data.frame(gene = matched_genes_upper), file.path(output_dir, "PC_30genes_matched.csv"), row.names = FALSE)
      write.csv(corr_matrix, file.path(output_dir, "PC_30genes_correlation.csv"), row.names = TRUE)
      
      cat("\n=== Summary ===\n")
      cat(sprintf("Target genes: %d\n", length(target_genes)))
      cat(sprintf("Matched: %d\n", length(matched_genes)))
      cat(sprintf("After filtering: %d\n", ncol(target_exp)))
      cat(sprintf("PC edges: %d\n", nrow(edge_list)))
      
      cat("\nResults saved:\n")
      cat("  - PC_30genes_pc_edges.csv\n")
      cat("  - PC_30genes_matched.csv\n")
      cat("  - PC_30genes_correlation.csv\n")
      
      found_edges <- TRUE
      break
    }
  }
}

if (!found_edges) {
  cat("\nNo edges found with PC algorithm. Saving correlation matrix instead.\n")
  
  output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创"
  write.csv(data.frame(gene = matched_genes_upper), file.path(output_dir, "PC_30genes_matched.csv"), row.names = FALSE)
  write.csv(corr_matrix, file.path(output_dir, "PC_30genes_correlation.csv"), row.names = TRUE)
  
  # 保存强相关性对
  strong_corr <- which(corr_matrix > 0.3 & upper.tri(corr_matrix), arr.ind = TRUE)
  if (nrow(strong_corr) > 0) {
    edge_list <- data.frame()
    for (i in 1:nrow(strong_corr)) {
      from <- rownames(corr_matrix)[strong_corr[i,1]]
      to <- colnames(corr_matrix)[strong_corr[i,2]]
      cor_val <- corr_matrix[strong_corr[i,1], strong_corr[i,2]]
      edge_list <- rbind(edge_list, data.frame(
        from = from, to = to, direction = "corr", correlation = cor_val, stringsAsFactors = FALSE
      ))
    }
    cat(sprintf("\nFound %d strong correlations (>0.3):\n", nrow(edge_list)))
    print(edge_list)
    write.csv(edge_list, file.path(output_dir, "PC_30genes_correlations.csv"), row.names = FALSE)
  }
}
