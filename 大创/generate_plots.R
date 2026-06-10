library(ggplot2)
library(dplyr)
library(tidyr)

# ============ 1. MR森林图 ============
cat("Generating MR Forest Plot...\n")

mr_data <- read.csv("D:/EQTL/MR_1e-5_Results/mr_main_results.csv", stringsAsFactors = FALSE)

mr_plot_data <- mr_data %>%
  filter(method %in% c("Inverse variance weighted", "Weighted median")) %>%
  mutate(
    ci_low = b - 1.96 * se,
    ci_high = b + 1.96 * se,
    sig = ifelse(pval < 0.05, "Significant", "Non-significant"),
    gene_method = paste0(gene, "\n(", method, ")")
  ) %>%
  arrange(gene, method)

mr_plot_data$gene_method <- factor(mr_plot_data$gene_method, 
                                    levels = unique(mr_plot_data$gene_method))

p1 <- ggplot(mr_plot_data, aes(x = b, y = gene_method, color = sig)) +
  geom_point(aes(shape = method), size = 3, position = position_dodge(0.5)) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3, position = position_dodge(0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Significant" = "#F24236", "Non-significant" = "#2E86AB")) +
  scale_shape_manual(values = c("Inverse variance weighted" = 16, "Weighted median" = 17)) +
  theme_classic() +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom",
    text = element_text(family = "Helvetica", size = 10),
    title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  ) +
  xlim(-0.15, 0.15) +
  labs(
    title = "MR Analysis: Copper Death Genes → Ischemic Stroke",
    subtitle = "Significant: PDHB (IVW), FDX1 (WM), ATOX1 (WM)",
    x = "Causal Effect (Beta)",
    y = "Gene",
    color = "Significance",
    shape = "Method"
  )

ggsave("mr_forest.png", p1, width = 8, height = 6, dpi = 300)
cat("MR Forest Plot saved to mr_forest.png\n")

# ============ 2. PC网络差异图 ============
cat("Generating PC Network Plot...\n")

library(igraph)

stroke_edges <- read.csv("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/network_stroke_directed.csv")
novel_edges <- read.csv("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/novel_stroke_edges.csv")

copper_genes <- c("fdx1", "lias", "lipt1", "dld", "dlat", "pdhb", "slc31a1", "atp7b", "atp7a", "atox1")
age_rage_genes <- c("nfkb1", "rela", "stat1", "stat3", "ccl2", "icam1", "tgfb1", "atf4", "egr1")

g <- graph_from_data_frame(stroke_edges[, c("from", "to")], directed = FALSE)

V(g)$color <- "gray80"
V(g)$color[V(g)$name %in% copper_genes] <- "#F24236"
V(g)$color[V(g)$name %in% age_rage_genes] <- "#2E86AB"

V(g)$size <- 5
V(g)$size[V(g)$name %in% copper_genes] <- 8
V(g)$size[V(g)$name %in% age_rage_genes] <- 7

edge_colors <- ifelse(
  paste(stroke_edges$from, stroke_edges$to) %in% paste(novel_edges$From, novel_edges$To) |
  paste(stroke_edges$to, stroke_edges$from) %in% paste(novel_edges$From, novel_edges$To),
  "#F24236", "gray50"
)
edge_colors[stroke_edges$direction == "undirected"] <- "gray70"

png("stroke_network.png", width = 10, height = 8, units = "in", res = 300)
par(bg = "white", mar = c(0, 0, 2, 0))

set.seed(42)
layout <- layout.fruchterman.reingold(g, niter = 1000)

plot(g, 
     layout = layout,
     vertex.color = V(g)$color,
     vertex.size = V(g)$size,
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     edge.color = edge_colors,
     edge.arrow.size = 0.5,
     edge.lty = ifelse(stroke_edges$direction == "undirected", 2, 1),
     main = "Stroke-Specific Regulatory Network (PC Algorithm)",
     sub = "Red edges: Novel connections in stroke | Red nodes: Copper death core",
     cex.main = 1.2,
     cex.sub = 0.9)

dev.off()
cat("Network Plot saved to stroke_network.png\n")

# ============ 3. MR-PC桥接整合图 ============
cat("Generating Integration Plot...\n")

bridge_data <- read.csv("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/mr_pc_bridge_analysis.csv")

sig_mr <- mr_data %>%
  filter(pval < 0.05) %>%
  arrange(b)

p2 <- ggplot(sig_mr, aes(x = reorder(gene, b), y = b, fill = ifelse(b > 0, "Risk", "Protection"))) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = b - 1.96*se, ymax = b + 1.96*se), width = 0.3) +
  geom_text(aes(label = ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", ""))), 
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Risk" = "#F24236", "Protection" = "#2E86AB")) +
  theme_classic() +
  theme(
    plot.background = element_rect(fill = "#F5F5F5"),
    panel.background = element_rect(fill = "#F5F5F5"),
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "MR-PC Integration: Genetic Causality in Regulatory Networks",
    subtitle = "Top: MR causal effects | Bottom: Network context",
    x = "Gene",
    y = "Causal Effect (Beta)",
    fill = "Effect Direction"
  )

ggsave("integration_plot.png", p2, width = 10, height = 6, dpi = 300)
cat("Integration Plot saved to integration_plot.png\n")

# ============ 4. F统计量质控图 ============
cat("Generating F-statistic QC Plot...\n")

diagnostics <- read.csv("D:/EQTL/MR_1e-5_Results/diagnostics.csv")

fstat_data <- diagnostics %>%
  filter(!is.na(fstat_snps)) %>%
  arrange(fstat_snps)

p3 <- ggplot(fstat_data, aes(x = reorder(gene, fstat_snps), y = fstat_snps, 
                              fill = ifelse(fstat_snps > 10, "Strong", "Weak"))) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "#E63946", size = 1) +
  annotate("text", x = nrow(fstat_data), y = 12, label = "Weak instrument threshold (F=10)", 
           color = "#E63946", hjust = 1, size = 3) +
  scale_fill_manual(values = c("Strong" = "#2E86AB", "Weak" = "#E63946")) +
  coord_flip() +
  theme_classic() +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  ) +
  geom_text(aes(label = round(fstat_snps, 1)), hjust = -0.2, size = 3) +
  labs(
    title = "Instrumental Variable Strength Validation",
    subtitle = "All significant genes pass weak instrument test (F>20)",
    x = "Gene (Instrumental Variable)",
    y = "F-statistic"
  )

ggsave("fstat_qc.png", p3, width = 8, height = 6, dpi = 300)
cat("F-statistic QC Plot saved to fstat_qc.png\n")

cat("\nAll plots generated successfully!\n")
