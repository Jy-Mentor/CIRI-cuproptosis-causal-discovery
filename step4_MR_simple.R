# Step 4: MR Analysis Summary
library(dplyr)

cat("=== Step 4: Mendelian Randomization (MR) Analysis ===\n\n")

# Read MR results
mr_files <- c(
  "D:/EQTL/MR_10000kb_Results/mr_main_results.csv",
  "D:/EQTL/MR_100kbcis_Results/mr_main_results.csv",
  "D:/EQTL/MR_Fixed_Results/mr_main_results.csv"
)

existing_files <- mr_files[file.exists(mr_files)]
all_results <- lapply(existing_files, read.csv, stringsAsFactors=FALSE)
combined_results <- bind_rows(all_results)

# BCP genes
bcp_genes <- c("IL6", "STAT3", "NFKB1", "TGFB1", "AGER", "PTGS2", "TLR4", "FDX1")
bcp_results <- combined_results[combined_results$gene %in% bcp_genes, ]

cat("BCP Genes MR Results:\n\n")

for(gene in unique(bcp_results$gene)) {
  gene_res <- bcp_results[bcp_results$gene == gene, ]
  cat(sprintf("[%s]\n", gene))
  for(i in 1:nrow(gene_res)) {
    sig <- ifelse(gene_res$pval[i] < 0.05, "***", ifelse(gene_res$pval[i] < 0.1, "*", ""))
    cat(sprintf("  %s: b=%.3f, p=%.4f %s\n", 
                gene_res$method[i], gene_res$b[i], gene_res$pval[i], sig))
  }
  cat("\n")
}

# Summary
cat("="); cat(rep("=", 59), sep=""); cat("\n")
cat("Summary:\n\n")

gene_best <- bcp_results %>% group_by(gene) %>% slice_min(pval, n=1)
sig_genes <- gene_best$gene[gene_best$pval < 0.05]
sugg_genes <- gene_best$gene[gene_best$pval >= 0.05 & gene_best$pval < 0.1]

if(length(sig_genes) > 0) {
  cat("Significant (p<0.05):", paste(sig_genes, collapse=", "), "\n")
}
if(length(sugg_genes) > 0) {
  cat("Suggestive (p<0.1):", paste(sugg_genes, collapse=", "), "\n")
}

missing <- setdiff(bcp_genes, unique(bcp_results$gene))
if(length(missing) > 0) {
  cat("Missing:", paste(missing, collapse=", "), "\n")
}

cat("="); cat(rep("=", 59), sep=""); cat("\n")
cat("Step 4 Complete!\n")
