#!/usr/bin/env Rscript
# Exposure Data Optimization Script v2.0
# Goal: Improve exposure rate from 82% to 95%+

library(readxl)
library(dplyr)

# Configuration
input_file <- "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_10000_r2_0.01.xlsx"
output_dir <- "D:/下载/MR_batch_results/exposure_optimized"
gene_list_file <- "D:/下载/MR_batch_results/gene_list_optimized.txt"

cat("\n============================================================\n")
cat("Exposure Data Optimization Script v2.0\n")
cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")

# Read data
cat("\nReading file:", input_file, "\n")
data <- read_excel(input_file)
cat("  Completed:", nrow(data), "rows,", ncol(data), "columns\n")

# Convert column names to uppercase
colnames(data) <- toupper(colnames(data))

# Filter candidate genes
gene_list_file_original <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/MR_batch_summary_20260506/10_gene_list.txt"
gene_list <- read.table(gene_list_file_original, header=TRUE, sep="\t", stringsAsFactors=FALSE)
candidate_genes <- gene_list$GENE

cat("\nFiltering candidate genes...\n")
data <- data %>% filter(GENE %in% candidate_genes)
cat("  Genes after filtering:", n_distinct(data$GENE), "\n")
cat("  SNPs after filtering:", nrow(data), "\n")

# Create output directory
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

# Export by gene
cat("\nExporting by gene...\n")
gene_stats <- data.frame(GENE=character(), N_SNPS=integer(), STATUS=character(), stringsAsFactors=FALSE)

for (gene in candidate_genes) {
  gene_df <- data %>% filter(GENE == gene)
  
  if (nrow(gene_df) == 0) {
    cat(sprintf("  %-10s: NO DATA\n", gene))
    gene_stats <- rbind(gene_stats, data.frame(GENE=gene, N_SNPS=0, STATUS="NO_DATA"))
    next
  }
  
  # Add samplesize column
  gene_df$samplesize <- 31684
  
  # Save
  output_file <- file.path(output_dir, paste0(gene, ".exposure.csv"))
  write.csv(gene_df, output_file, row.names=FALSE)
  
  cat(sprintf("  %-10s: %d SNPs\n", gene, nrow(gene_df)))
  gene_stats <- rbind(gene_stats, data.frame(GENE=gene, N_SNPS=nrow(gene_df), STATUS="OK"))
}

# Save gene list
cat("\nSaving gene list:", gene_list_file, "\n")
write.table(gene_stats, gene_list_file, sep="\t", row.names=FALSE, quote=FALSE)

# Summary statistics
n_with_data <- sum(gene_stats$N_SNPS > 0)
n_without_data <- nrow(gene_stats) - n_with_data
total_snps <- sum(gene_stats$N_SNPS)

cat("\n============================================================\n")
cat("Optimization Complete - Summary\n")
cat("============================================================\n")
cat("Total candidate genes:", length(candidate_genes), "\n")
cat("Genes with exposure data:", n_with_data, sprintf("(%.1f%%)", n_with_data/length(candidate_genes)*100), "\n")
cat("Genes without exposure data:", n_without_data, sprintf("(%.1f%%)", n_without_data/length(candidate_genes)*100), "\n")
cat("Total SNPs:", total_snps, "\n")
cat("Average SNPs per gene:", sprintf("%.2f", total_snps/n_with_data), "\n")
cat("\nOutput directory:", output_dir, "\n")
cat("Gene list file:", gene_list_file, "\n")

# Compare with original
original_rate <- 107/130*100
optimized_rate <- n_with_data/length(candidate_genes)*100
improvement <- optimized_rate - original_rate

cat("\n============================================================\n")
cat("Optimization Effect Comparison\n")
cat("============================================================\n")
cat(sprintf("Original exposure rate:  %.1f%% (107/130)\n", original_rate))
cat(sprintf("Optimized exposure rate: %.1f%% (%d/%d)\n", optimized_rate, n_with_data, length(candidate_genes)))
cat(sprintf("Improvement:             +%.1f%%\n", improvement))
cat(sprintf("New genes added:         %d\n", n_with_data - 107))

if (n_with_data > 107) {
  cat("\nSUCCESS! Exposure rate improved from", sprintf("%.1f%%", original_rate), "to", sprintf("%.1f%%", optimized_rate), "\n")
} else {
  cat("\nWarning: Limited improvement, please check data source\n")
}

cat("\nNext steps:\n")
cat("1. Check optimized exposure data in:", output_dir, "\n")
cat("2. Run MR analysis: Rscript mr_analysis_batch.R", output_dir, "./outcome", gene_list_file, "\n")
cat("3. Compare results before and after optimization\n")
cat("\nDone!\n\n")
