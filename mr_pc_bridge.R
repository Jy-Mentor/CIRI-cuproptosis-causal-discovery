library(data.table)
library(igraph)
library(stringr)

tryCatch({
  cat("=== MR-PC Bridge Analysis Starting ===\n")
  
  output_dir <- "D:/EQTL/MR_PC_Bridge"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Output directory created:", output_dir, "\n")
  
  mr_main_path <- "D:/EQTL/MR_1e-5_Results/mr_main_results.csv"
  mr_sens_path <- "D:/EQTL/MR_1e-5_Results/mr_sensitivity.csv"
  mr_sig_path <- "D:/EQTL/MR_1e-5_Results/mr_significant_results.csv"
  mr_steiger_path <- "D:/EQTL/MR_1e-5_Results/mr_steiger.csv"
  
  pc_stroke_edges <- "c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/network_stroke_edges.csv"
  pc_stroke_directed <- "c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/network_stroke_directed.csv"
  pc_control_edges <- "c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/network_control_edges.csv"
  pc_novel <- "c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/novel_stroke_edges.csv"
  
  cat("=== Step 1: Reading MR Results ===\n")
  mr_main <- fread(mr_main_path)
  mr_sens <- fread(mr_sens_path)
  mr_sig <- fread(mr_sig_path)
  mr_steiger <- fread(mr_steiger_path)
  
  cat("MR main results rows:", nrow(mr_main), "\n")
  cat("MR sensitivity rows:", nrow(mr_sens), "\n")
  cat("MR significant rows:", nrow(mr_sig), "\n")
  
  ivw_data <- mr_main[method == "Inverse variance weighted", .(gene, beta = b, pval)]
  setnames(ivw_data, c("beta", "pval"), c("ivw_beta", "ivw_pval"))
  
  cat("IVW data extracted:", nrow(ivw_data), "genes\n")
  
  mr_combined <- merge(ivw_data, mr_sens[, .(gene, egger_p, presso_p)], by = "gene", all.x = TRUE)
  
  mr_combined[, pph4 := 0.5]
  mr_combined[, steiger_p := NA_real_]
  
  steiger_summary <- mr_steiger[, .(steiger_p = mean(as.numeric(pval.outcome), na.rm = TRUE)), by = gene]
  mr_combined <- merge(mr_combined, steiger_summary, by = "gene", all.x = TRUE)
  
  mr_combined[, gene_lower := tolower(gene)]
  
  mr_combined[, evidence_level := "None"]
  mr_combined[ivw_pval < 0.05 & pph4 > 0.8 & (is.na(egger_p) | egger_p > 0.05), evidence_level := "Strong"]
  mr_combined[ivw_pval < 0.05 & pph4 > 0.5 & evidence_level == "None", evidence_level := "Moderate"]
  mr_combined[ivw_pval < 0.05 & evidence_level == "None", evidence_level := "Weak"]
  
  cat("Evidence levels:\n")
  print(mr_combined[, .N, by = evidence_level])
  
  cat("\n=== Step 2: Reading PC Network ===\n")
  
  read_pc_edges <- function(path) {
    if (!file.exists(path)) {
      cat("Warning: File not found:", path, "\n")
      return(NULL)
    }
    dt <- fread(path)
    if ("from" %in% names(dt) && "to" %in% names(dt)) {
      return(dt[, .(source = from, target = to, weight)])
    }
    if ("regulator" %in% names(dt) && "target" %in% names(dt)) {
      return(dt[, .(source = regulator, target = target, weight)])
    }
    return(NULL)
  }
  
  edges_stroke <- read_pc_edges(pc_stroke_edges)
  edges_directed <- read_pc_edges(pc_stroke_directed)
  edges_control <- read_pc_edges(pc_control_edges)
  edges_novel <- read_pc_edges(pc_novel)
  
  all_edges <- rbindlist(list(edges_stroke, edges_directed, edges_control, edges_novel), fill = TRUE, use.names = TRUE)
  all_edges <- unique(all_edges, by = c("source", "target"))
  cat("Total unique edges:", nrow(all_edges), "\n")
  
  g <- graph_from_data_frame(all_edges[, .(source, target, weight)], directed = FALSE)
  
  degree_cent <- degree(g, normalized = FALSE)
  betweenness_cent <- betweenness(g, normalized = TRUE)
  
  degree_quantile <- quantile(degree_cent, 0.75)
  cat("Degree 75% quantile:", degree_quantile, "\n")
  
  node_metrics <- data.table(
    gene = names(degree_cent),
    degree = as.integer(degree_cent),
    degree_normalized = degree_cent / max(degree_cent),
    betweenness = as.numeric(betweenness_cent),
    is_hub = degree_cent > degree_quantile
  )
  
  cat("Hub nodes (degree > 75% quantile):", sum(node_metrics$is_hub), "\n")
  
  cat("\n=== Step 3: Bridge Logic ===\n")
  
  source_mr <- mr_combined[, .(gene_lower, gene, ivw_beta, ivw_pval, egger_p, pph4, evidence_level)]
  setnames(source_mr, c("gene_lower", "gene", "ivw_beta", "ivw_pval", "egger_p", "pph4", "evidence_level"),
           c("source", "source_mr_gene", "source_beta", "source_pval", "source_egger_p", "source_pph4", "source_evidence"))
  
  target_mr <- mr_combined[, .(gene_lower, gene, ivw_beta, ivw_pval, egger_p, pph4, evidence_level)]
  setnames(target_mr, c("gene_lower", "gene", "ivw_beta", "ivw_pval", "egger_p", "pph4", "evidence_level"),
           c("target", "target_mr_gene", "target_beta", "target_pval", "target_egger_p", "target_pph4", "target_evidence"))
  
  bridge_annotated <- merge(all_edges, source_mr, by = "source", all.x = TRUE)
  bridge_annotated <- merge(bridge_annotated, target_mr, by = "target", all.x = TRUE)
  
  bridge_annotated[, mr_supported := FALSE]
  bridge_annotated[(!is.na(source_pval) & source_pval < 0.05) | (!is.na(target_pval) & target_pval < 0.05), mr_supported := TRUE]
  
  bridge_annotated[, effect_concordance := "Unknown"]
  bridge_annotated[!is.na(source_beta) & !is.na(target_beta), 
                  effect_concordance := ifelse(sign(source_beta) == sign(target_beta), "Synergistic", "Antagonistic")]
  
  bridge_annotated[, bridge_score := 0]
  bridge_annotated[source_evidence == "Strong", bridge_score := bridge_score + 2]
  bridge_annotated[source_evidence == "Moderate", bridge_score := bridge_score + 1]
  bridge_annotated[target_evidence == "Strong", bridge_score := bridge_score + 2]
  bridge_annotated[target_evidence == "Moderate", bridge_score := bridge_score + 1]
  
  cat("MR-supported edges:", sum(bridge_annotated$mr_supported), "\n")
  cat("Synergistic edges:", sum(bridge_annotated$effect_concordance == "Synergistic"), "\n")
  
  cat("\n=== Step 4: Node Priority ===\n")
  
  mr_for_node <- mr_combined[, .(gene_lower, ivw_beta, ivw_pval, pph4, evidence_level)]
  setnames(mr_for_node, "gene_lower", "gene")
  
  node_metrics <- merge(node_metrics, mr_for_node, by = "gene", all.x = TRUE)
  
  node_metrics[, evidence_score := 0]
  node_metrics[evidence_level == "Strong", evidence_score := 3]
  node_metrics[evidence_level == "Moderate", evidence_score := 2]
  node_metrics[evidence_level == "Weak", evidence_score := 1]
  
  node_metrics[, hub_bonus := as.integer(is_hub)]
  node_metrics[, priority_rank := evidence_score + hub_bonus]
  
  node_metrics[, recommendation := "Exploratory"]
  node_metrics[priority_rank >= 3 | (evidence_level == "Strong" & is_hub), recommendation := "Priority_Target"]
  node_metrics[priority_rank >= 2 & recommendation == "Exploratory", recommendation := "Secondary"]
  
  cat("Recommendations:\n")
  print(node_metrics[, .N, by = recommendation])
  
  cat("\n=== Step 5: Mechanism Classification ===\n")
  
  copper_genes <- c("slc31a1", "atox1", "atp7a", "atp7b")
  lipoic_genes <- c("dlat", "dld", "lias", "fdx1")
  inflammation_genes <- c("nfkb1", "stat1", "stat3", "ccl2", "icam1")
  
  node_metrics[, mechanism := "Other"]
  node_metrics[tolower(gene) %in% copper_genes, mechanism := "Copper_Transport"]
  node_metrics[tolower(gene) %in% lipoic_genes, mechanism := "Lipoic_Acid"]
  node_metrics[tolower(gene) %in% inflammation_genes, mechanism := "Inflammation"]
  
  cat("Mechanism distribution:\n")
  print(node_metrics[, .N, by = mechanism])
  
  cat("\n=== Step 6: Generating Outputs ===\n")
  
  fwrite(bridge_annotated, file.path(output_dir, "bridge_annotated_edges.csv"))
  cat("Saved: bridge_annotated_edges.csv\n")
  
  fwrite(node_metrics, file.path(output_dir, "bridge_node_priority.csv"))
  cat("Saved: bridge_node_priority.csv\n")
  
  high_conf <- node_metrics[recommendation == "Priority_Target" | pph4 > 0.8]
  fwrite(high_conf, file.path(output_dir, "high_confidence_targets.csv"))
  cat("Saved: high_confidence_targets.csv\n")
  
  mechanism_summary <- node_metrics[mechanism != "Other", .(
    strong_evidence = sum(evidence_level == "Strong", na.rm = TRUE),
    hub_nodes = sum(is_hub, na.rm = TRUE),
    total_genes = .N
  ), by = mechanism]
  fwrite(mechanism_summary, file.path(output_dir, "mechanism_summary.csv"))
  cat("Saved: mechanism_summary.csv\n")
  
  cat("\n=== Console Summary ===\n")
  
  total_edges <- nrow(bridge_annotated)
  mr_supported_count <- sum(bridge_annotated$mr_supported)
  mr_ratio <- mr_supported_count / total_edges * 100
  cat("MR-supported edge ratio:", round(mr_ratio, 2), "%\n")
  
  top5 <- node_metrics[!is.na(ivw_pval)][order(ivw_pval)]
  if (nrow(top5) > 0) {
    top5 <- head(top5, 5)
    cat("\nTOP5 Targets (by p-value):\n")
    print(top5[, .(gene, beta = ivw_beta, pval = ivw_pval, pph4)])
  }
  
  cat("\n=== Analysis Complete ===\n")
  
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  cat("Stack trace:\n")
  print(sys.calls())
})
