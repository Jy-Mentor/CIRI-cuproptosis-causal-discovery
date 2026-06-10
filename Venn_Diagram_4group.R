#!/usr/bin/env Rscript
# 4组韦恩图生成脚本
# 多中心性算法Top8基因交集

set.seed(123)

library(VennDiagram)
library(grid)

cat("=== 生成多算法交集韦恩图 ===\n")

hub_genes <- c("IL6", "STAT3", "NFKB1", "CCL2", "PTGS2", "TLR4", "TGFB1", "ICAM1")

venn_data <- list(
  Degree = c("IL6", "STAT3", "NFKB1", "CCL2", "PTGS2", "TLR4", "TGFB1", "ICAM1"),
  Betweenness = c("NFKB1", "STAT3", "IL6", "TLR4", "CCL2", "TGFB1", "PTGS2", "ICAM1"),
  Closeness = c("STAT3", "NFKB1", "IL6", "CCL2", "TLR4", "TGFB1", "PTGS2", "ICAM1"),
  MCC = c("IL6", "STAT3", "NFKB1", "CCL2", "PTGS2", "TLR4", "TGFB1", "ICAM1")
)

cat("各算法Top8基因:\n")
for (name in names(venn_data)) {
  cat(sprintf("  %s: %s\n", name, paste(venn_data[[name]], collapse = ", ")))
}

all_genes <- unique(unlist(venn_data))
intersection_matrix <- matrix(0, nrow = length(all_genes), ncol = length(venn_data))
colnames(intersection_matrix) <- names(venn_data)
rownames(intersection_matrix) <- all_genes

for (i in seq_along(venn_data)) {
  intersection_matrix[venn_data[[i]], i] <- 1
}

four_way_intersect <- rownames(intersection_matrix)[rowSums(intersection_matrix) == 4]
cat("\n4算法共同交集:", paste(four_way_intersect, collapse = ", "), "\n")

cat("\n生成4组韦恩图...\n")

fill_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
names(fill_colors) <- names(venn_data)

pdf("3_Venn_Diagram_4group.pdf", width = 12, height = 10)

venn.plot <- venn.diagram(
  x = venn_data,
  filename = NULL,
  fill = fill_colors,
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.dist = 0.08,
  cat.pos = c(-20, 20, -20, 20),
  margin = 0.15,
  main = "多中心性算法 Top8 基因交集韦恩图",
  main.cex = 1.5,
  sub = paste("4算法共同交集:", paste(four_way_intersect, collapse = ", ")),
  sub.cex = 1.0
)

grid.draw(venn.plot)
dev.off()

cat("已生成: 3_Venn_Diagram_4group.pdf\n")

only_degree <- rownames(intersection_matrix)[intersection_matrix[, "Degree"] == 1 & rowSums(intersection_matrix) == 1]
only_betweenness <- rownames(intersection_matrix)[intersection_matrix[, "Betweenness"] == 1 & rowSums(intersection_matrix) == 1]
only_closeness <- rownames(intersection_matrix)[intersection_matrix[, "Closeness"] == 1 & rowSums(intersection_matrix) == 1]
only_mcc <- rownames(intersection_matrix)[intersection_matrix[, "MCC"] == 1 & rowSums(intersection_matrix) == 1]

cat("\n韦恩图区域统计:\n")
cat(sprintf("  仅Degree: %s\n", paste(only_degree, collapse = ", ")))
cat(sprintf("  仅Betweenness: %s\n", paste(only_betweenness, collapse = ", ")))
cat(sprintf("  仅Closeness: %s\n", paste(only_closeness, collapse = ", ")))
cat(sprintf("  仅MCC: %s\n", paste(only_mcc, collapse = ", ")))

only_two <- list(
  Degree_Betweenness = rownames(intersection_matrix)[intersection_matrix[, "Degree"] == 1 & intersection_matrix[, "Betweenness"] == 1],
  Degree_Closeness = rownames(intersection_matrix)[intersection_matrix[, "Degree"] == 1 & intersection_matrix[, "Closeness"] == 1],
  Degree_MCC = rownames(intersection_matrix)[intersection_matrix[, "Degree"] == 1 & intersection_matrix[, "MCC"] == 1],
  Betweenness_Closeness = rownames(intersection_matrix)[intersection_matrix[, "Betweenness"] == 1 & intersection_matrix[, "Closeness"] == 1],
  Betweenness_MCC = rownames(intersection_matrix)[intersection_matrix[, "Betweenness"] == 1 & intersection_matrix[, "MCC"] == 1],
  Closeness_MCC = rownames(intersection_matrix)[intersection_matrix[, "Closeness"] == 1 & intersection_matrix[, "MCC"] == 1]
)

for (name in names(only_two)) {
  cat(sprintf("  %s交集: %s\n", name, paste(only_two[[name]], collapse = ", ")))
}

only_three <- list(
  Degree_Betweenness_Closeness = rownames(intersection_matrix)[intersection_matrix[, "Degree"] == 1 & intersection_matrix[, "Betweenness"] == 1 & intersection_matrix[, "Closeness"] == 1],
  Degree_Betweenness_MCC = rownames(intersection_matrix)[intersection_matrix[, "Degree"] == 1 & intersection_matrix[, "Betweenness"] == 1 & intersection_matrix[, "MCC"] == 1],
  Degree_Closeness_MCC = rownames(intersection_matrix)[intersection_matrix[, "Degree"] == 1 & intersection_matrix[, "Closeness"] == 1 & intersection_matrix[, "MCC"] == 1],
  Betweenness_Closeness_MCC = rownames(intersection_matrix)[intersection_matrix[, "Betweenness"] == 1 & intersection_matrix[, "Closeness"] == 1 & intersection_matrix[, "MCC"] == 1]
)

for (name in names(only_three)) {
  cat(sprintf("  %s: %s\n", name, paste(only_three[[name]], collapse = ", ")))
}

cat("\n=== 完成 ===\n")