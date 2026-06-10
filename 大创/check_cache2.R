cache <- readRDS('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/analysis_cache.rds')
cat("Cache keys:\n")
print(names(cache))

cat("\n\nData dimensions:\n")
for (nm in names(cache)) {
  if (is.matrix(cache[[nm]])) {
    cat(nm, ":", dim(cache[[nm]])[1], "x", dim(cache[[nm]])[2], "\n")
  } else if (is.vector(cache[[nm]])) {
    cat(nm, ": length =", length(cache[[nm]]), "\n")
  } else {
    cat(nm, ":", class(cache[[nm]]), "\n")
  }
}

# Check if there's a gene list
if ("target_genes" %in% names(cache)) {
  cat("\n\nTarget genes:\n")
  print(cache$target_genes)
}

# Try to find expression matrix
possible_expr_names <- c("target_dense", "expr_matrix", "expression", "data", "target_data")
for (nm in possible_expr_names) {
  if (nm %in% names(cache) && is.matrix(cache[[nm]])) {
    cat("\n\nExpression matrix found:", nm, "\n")
    cat("Genes (first 50):\n")
    print(colnames(cache[[nm]])[1:min(50, ncol(cache[[nm]]))])
    break
  }
}
