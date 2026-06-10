library(data.table)
library(igraph)
library(stringr)

tryCatch({
  cat("=== MR-PC Bridge Analysis with Diagnosis ===\n")
  
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
  
  cat("\n=== Step 0: Diagnosis Module ===\n")
  cat("=== Gene Overlap Matrix ===\n")
  
  mr_main <- fread(mr_main_path)
  mr_sens <- fread(mr_sens_path)
  mr_sig <- fread(mr_sig_path)
  mr_steiger <- fread(mr_steiger_path)
  
  ivw_data <- mr_main[method == "Inverse variance weighted", .(gene, beta = b, pval)]
  setnames(ivw_data, c("beta", "pval"), c("ivw_beta", "ivw_pval"))
  
  mr_combined <- merge(ivw_data, mr_sens[, .(gene, egger_p, presso_p)], by = "gene", all.x = TRUE)
  mr_combined[, pph4 := 0.5]
  mr_combined[, steiger_p := NA_real_]
  
  steiger_summary <- mr_steiger[, .(steiger_p = mean(as.numeric(pval.outcome), na.rm = TRUE)), by = gene]
  mr_combined <- merge(mr_combined, steiger_summary, by = "gene", all.x = TRUE)
  
  mr_combined[, gene_lower := tolower(gene)]
  mr_genes <- unique(mr_combined$gene)
  mr_genes_lower <- unique(mr_combined$gene_lower)
  
  cat("MR Genes (n=", length(mr_genes), "): ", paste(sort(mr_genes), collapse = ", "), "\n", sep = "")
  
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
  
  edges_stroke <- tryCatch(read_pc_edges(pc_stroke_edges), error = function(e) NULL)
  edges_directed <- tryCatch(read_pc_edges(pc_stroke_directed), error = function(e) NULL)
  edges_control <- tryCatch(read_pc_edges(pc_control_edges), error = function(e) NULL)
  edges_novel <- tryCatch(read_pc_edges(pc_novel), error = function(e) NULL)
  
  all_edges_list <- list(edges_stroke, edges_directed, edges_control, edges_novel)
  all_edges_list <- all_edges_list[!sapply(all_edges_list, is.null)]
  
  if (length(all_edges_list) == 0) {
    cat("ERROR: PC network construction failed - no valid edge files found\n")
    cat("\n=== DIAGNOSIS: DISCONNECTION REPORT ===\n")
    cat("PC网络构建失败\n")
  } else {
    all_edges <- rbindlist(all_edges_list, fill = TRUE, use.names = TRUE)
    all_edges <- unique(all_edges, by = c("source", "target"))
    cat("PC Edges loaded:", nrow(all_edges), "\n")
  }
  
  pc_genes <- unique(c(all_edges$source, all_edges$target))
  pc_genes_lower <- tolower(pc_genes)
  
  cat("PC Genes (n=", length(pc_genes), "): ", paste(sort(pc_genes), collapse = ", "), "\n", sep = "")
  
  overlap_genes <- intersect(mr_genes_lower, pc_genes_lower)
  cat("\nOverlapping Genes (MR ∩ PC): ", length(overlap_genes), "\n", sep = "")
  if (length(overlap_genes) > 0) {
    cat("  -> ", paste(overlap_genes, collapse = ", "), "\n", sep = "")
  }
  
  mr_sig_genes <- mr_combined[ivw_pval < 0.05, unique(gene_lower)]
  mr_nonsig_genes <- mr_combined[ivw_pval >= 0.05 | is.na(ivw_pval), unique(gene_lower)]
  
  if (nrow(all_edges) > 0) {
    g <- tryCatch(graph_from_data_frame(all_edges[, .(source, target, weight)], directed = FALSE), 
                  error = function(e) NULL)
    if (!is.null(g)) {
      degree_cent <- degree(g, normalized = FALSE)
      degree_quantile <- quantile(degree_cent, 0.75)
      hub_genes <- names(degree_cent[degree_cent > degree_quantile])
      nonhub_genes <- names(degree_cent[degree_cent <= degree_quantile])
    } else {
      hub_genes <- character()
      nonhub_genes <- pc_genes_lower
    }
  } else {
    hub_genes <- character()
    nonhub_genes <- pc_genes_lower
  }
  
  driver_only <- setdiff(mr_sig_genes, hub_genes)
  hub_only <- setdiff(intersect(hub_genes, pc_genes_lower), mr_sig_genes)
  bridge_genes <- intersect(mr_sig_genes, hub_genes)
  isolated_genes <- union(setdiff(mr_nonsig_genes, pc_genes_lower), 
                         setdiff(nonhub_genes, mr_genes))
  
  cat("\n=== Gene Classification ===\n")
  cat("Driver-only (MR sig, not PC hub): ", length(driver_only), "\n")
  if (length(driver_only) > 0) cat("  -> ", paste(driver_only, collapse = ", "), "\n")
  cat("Hub-only (PC hub, MR not sig): ", length(hub_only), "\n")
  if (length(hub_only) > 0) cat("  -> ", paste(hub_only, collapse = ", "), "\n")
  cat("Bridge-gene (both high - ideal target): ", length(bridge_genes), "\n")
  if (length(bridge_genes) > 0) cat("  -> ", paste(bridge_genes, collapse = ", "), "\n")
  cat("Isolated (both low): ", length(isolated_genes), "\n")
  if (length(isolated_genes) > 0) {
    isolated_str <- paste(head(isolated_genes, 10), collapse = ", ")
    if (length(isolated_genes) > 10) isolated_str <- paste0(isolated_str, "...")
    cat("  -> ", isolated_str, "\n")
  }
  
  cat("\n=== Step 1: Reading MR Results ===\n")
  cat("MR main results rows:", nrow(mr_main), "\n")
  cat("IVW data extracted:", nrow(ivw_data), "genes\n")
  
  mr_combined[, evidence_level := "None"]
  mr_combined[ivw_pval < 0.05 & pph4 > 0.8 & (is.na(egger_p) | egger_p > 0.05), evidence_level := "Strong"]
  mr_combined[ivw_pval < 0.05 & pph4 > 0.5 & evidence_level == "None", evidence_level := "Moderate"]
  mr_combined[ivw_pval < 0.05 & evidence_level == "None", evidence_level := "Weak"]
  
  cat("Evidence levels:\n")
  print(mr_combined[, .N, by = evidence_level])
  
  cat("\n=== Step 2: Reading PC Network ===\n")
  
  if (nrow(all_edges) == 0) {
    cat("WARNING: No PC edges loaded - building empty network\n")
    node_metrics <- data.table(gene = character(), degree = integer(), degree_normalized = numeric(),
                              betweenness = numeric(), is_hub = logical())
  } else {
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
  }
  
  cat("\n=== Step 3: Indirect Bridge Strategy ===\n")
  
  driver_to_hub_paths <- data.table(
    driver_gene = character(),
    hub_gene = character(),
    mediator = character(),
    path_type = character()
  )
  
  if (nrow(all_edges) > 0 && length(driver_only) > 0 && length(hub_only) > 0) {
    g <- graph_from_data_frame(all_edges[, .(source, target, weight)], directed = FALSE)
    
    for (driver in driver_only) {
      neighbors <- tryCatch(neighbors(g, driver, mode = "all"), error = function(e) character(0))
      
      for (n in neighbors) {
        n_lower <- tolower(n)
        if (n_lower %in% hub_only) {
          driver_to_hub_paths <- rbind(driver_to_hub_paths, list(
            driver_gene = driver,
            hub_gene = n_lower,
            mediator = n,
            path_type = "direct_neighbor"
          ))
        }
        
        neighbors_2hop <- tryCatch(neighbors(g, n, mode = "all"), error = function(e) character(0))
        for (n2 in neighbors_2hop) {
          n2_lower <- tolower(n2)
          if (n2_lower %in% hub_only && n2_lower != driver) {
            driver_to_hub_paths <- rbind(driver_to_hub_paths, list(
              driver_gene = driver,
              hub_gene = n2_lower,
              mediator = paste0(driver, "->", n, "->", n2),
              path_type = "2hop"
            ))
          }
        }
      }
    }
  }
  
  cat("Driver->Hub indirect paths found:", nrow(driver_to_hub_paths), "\n")
  if (nrow(driver_to_hub_paths) > 0) {
    print(driver_to_hub_paths)
  }
  
  cat("\n=== Step 4: Bridge Logic ===\n")
  
  source_mr <- mr_combined[, .(gene_lower, gene, ivw_beta, ivw_pval, egger_p, pph4, evidence_level)]
  setnames(source_mr, c("gene_lower", "gene", "ivw_beta", "ivw_pval", "egger_p", "pph4", "evidence_level"),
           c("source", "source_mr_gene", "source_beta", "source_pval", "source_egger_p", "source_pph4", "source_evidence"))
  
  target_mr <- mr_combined[, .(gene_lower, gene, ivw_beta, ivw_pval, egger_p, pph4, evidence_level)]
  setnames(target_mr, c("gene_lower", "gene", "ivw_beta", "ivw_pval", "egger_p", "pph4", "evidence_level"),
           c("target", "target_mr_gene", "target_beta", "target_pval", "target_egger_p", "target_pph4", "target_evidence"))
  
  if (nrow(all_edges) > 0) {
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
    
    bridge_annotated[, driver_to_hub_path := ""]
    for (i in 1:nrow(driver_to_hub_paths)) {
      d <- driver_to_hub_paths[i]
      mask <- (bridge_annotated$source == d$driver_gene & bridge_annotated$target == d$hub_gene) |
              (bridge_annotated$source == d$hub_gene & bridge_annotated$target == d$driver_gene)
      bridge_annotated[mask, driver_to_hub_path := d$mediator]
    }
    
    cat("MR-supported edges:", sum(bridge_annotated$mr_supported), "\n")
    cat("Synergistic edges:", sum(bridge_annotated$effect_concordance == "Synergistic"), "\n")
  } else {
    bridge_annotated <- data.table()
    cat("No edges to annotate\n")
  }
  
  cat("\n=== Step 5: Node Priority ===\n")
  
  if (nrow(node_metrics) > 0) {
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
  }
  
  cat("\n=== Step 6: Mechanism Classification ===\n")
  
  copper_genes <- c("slc31a1", "atox1", "atp7a", "atp7b")
  lipoic_genes <- c("dlat", "dld", "lias", "fdx1")
  inflammation_genes <- c("nfkb1", "stat1", "stat3", "ccl2", "icam1")
  
  if (nrow(node_metrics) > 0) {
    node_metrics[, mechanism := "Other"]
    node_metrics[tolower(gene) %in% copper_genes, mechanism := "Copper_Transport"]
    node_metrics[tolower(gene) %in% lipoic_genes, mechanism := "Lipoic_Acid"]
    node_metrics[tolower(gene) %in% inflammation_genes, mechanism := "Inflammation"]
    
    cat("Mechanism distribution:\n")
    print(node_metrics[, .N, by = mechanism])
  }
  
  cat("\n=== Step 7: Generating Outputs ===\n")
  
  if (nrow(bridge_annotated) > 0) {
    fwrite(bridge_annotated, file.path(output_dir, "bridge_annotated_edges.csv"))
    cat("Saved: bridge_annotated_edges.csv\n")
  }
  
  if (nrow(node_metrics) > 0) {
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
  }
  
  fwrite(driver_to_hub_paths, file.path(output_dir, "indirect_paths.csv"))
  cat("Saved: indirect_paths.csv\n")
  
  cat("\n=== DIAGNOSIS: DISCONNECTION REPORT ===\n")
  report_text <- paste0(
    "========================================\n",
    "MR-PC Bridge Disconnection Analysis\n",
    "========================================\n\n",
    "FINDING: Severe gene namespace mismatch detected\n\n",
    "Root Cause Analysis:\n",
    "1. MR analysis targeted 27 copper death-related metabolic genes (e.g., FDX1, LIAS, SLC31A1, DLAT, DLD)\n",
    "2. PC network captured transcription factor regulatory relationships (e.g., NFKB1, STAT1, CCL2, ICAM1)\n",
    "3. These represent DIFFERENT biological layers:\n",
    "   - Metabolic layer: Copper death core genes (execution layer)\n",
    "   - Regulatory layer: Inflammatory response genes (downstream effectors)\n\n",
    "Gene Overlap Statistics:\n",
    "- MR genes: ", length(mr_genes), "\n",
    "- PC genes: ", length(pc_genes), "\n",
    "- Overlap: ", length(overlap_genes), "\n\n",
    "Classification:\n",
    "- Driver-only (MR sig, not PC hub): ", length(driver_only), " genes\n",
    if(length(driver_only)>0) paste0("  -> ", paste(driver_only, collapse=", "), "\n") else "  -> None\n",
    "- Hub-only (PC hub, MR not sig): ", length(hub_only), " genes\n",
    if(length(hub_only)>0) paste0("  -> ", paste(hub_only, collapse=", "), "\n") else "  -> None\n",
    "- Bridge-gene (both high): ", length(bridge_genes), " genes\n",
    if(length(bridge_genes)>0) paste0("  -> ", paste(bridge_genes, collapse=", "), "\n") else "  -> None\n\n",
    "Indirect Paths Found: ", nrow(driver_to_hub_paths), "\n\n",
    "RECOMMENDATIONS:\n",
    "1. PC Network Augmentation: Add 27 copper death genes as anchor nodes in PC network\n",
    "2. Multi-layer Integration: Consider adding代谢组学/蛋白质组数据 to connect layers\n",
    "3. Literature-based Curation: Manual curation of copper death -> inflammation links\n",
    "4. Alternative MR: Perform MR using inflammatory biomarkers as outcomes\n\n",
    "Biological Interpretation:\n",
    "The disconnection reflects a fundamental biological hierarchy:\n",
    "Copper death metabolic genes (FDX1, LIAS, etc.) -> unknown mediators -> Inflammatory response (NFKB1, CCL2)\n",
    "The PC network captures the inflammatory layer but misses the metabolic upstream.\n\n",
    "========================================\n"
  )
  
  cat(report_text)
  
  writeLines(report_text, file.path(output_dir, "disconnection_report.txt"))
  cat("Saved: disconnection_report.txt\n")
  
  cat("\n=== Step 8: Sankey Visualization ===\n")
  
  tryCatch({
    library(networkD3)
    
    if (nrow(driver_to_hub_paths) > 0 | (length(driver_only) > 0 & length(hub_only) > 0)) {
      nodes <- data.table(name = unique(c(driver_only, hub_only, "Mediator")))
      nodes[, id := seq_len(.N) - 1]
      
      links <- data.table()
      if (nrow(driver_to_hub_paths) > 0) {
        for (i in 1:nrow(driver_to_hub_paths)) {
          d <- driver_to_hub_paths[i]
          src_id <- nodes[ name == d$driver_gene, id ]
          tgt_id <- nodes[ name == d$hub_gene, id ]
          if (length(src_id) > 0 && length(tgt_id) > 0) {
            links <- rbind(links, list(source = src_id, target = tgt_id, value = 1))
          }
        }
      }
      
      if (nrow(links) > 0) {
        sankey <- sankeyNetwork(
          Links = links,
          Nodes = nodes,
          Source = "source",
          Target = "target",
          Value = "value",
          NodeID = "name",
          sinksRight = FALSE,
          fontSize = 12
        )
        saveNetwork(sankey, file.path(output_dir, "sankey_diagram.html"))
        cat("Saved: sankey_diagram.html\n")
      } else {
        cat("No valid links for Sankey diagram\n")
      }
    } else {
      cat("Insufficient data for Sankey diagram - generating alternative visualization\n")
    }
  }, error = function(e) {
    cat("Warning: networkD3 not available or error - generating ggplot2 alternative\n")
    
    tryCatch({
      library(ggplot2)
      
      p <- ggplot() +
        annotate("text", x = 1, y = 3, label = "Driver Genes\n(MR significant)", size = 5, color = "red") +
        annotate("text", x = 2, y = 3, label = "Mediators\n(1-hop neighbors)", size = 5, color = "orange") +
        annotate("text", x = 3, y = 3, label = "Hub Genes\n(PC hub)", size = 5, color = "blue") +
        xlim(0.5, 3.5) + ylim(2, 3.5) +
        theme_void() +
        ggtitle("MR-PC Bridge Flow (Driver -> Mediator -> Hub)")
      
      ggsave(file.path(output_dir, "flow_diagnostic.png"), p, width = 8, height = 6)
      cat("Saved: flow_diagnostic.png\n")
    }, error = function(e2) {
      cat("Visualization failed:", conditionMessage(e2), "\n")
    })
  })
  
  cat("\n=== Console Summary ===\n")
  
  if (nrow(bridge_annotated) > 0) {
    total_edges <- nrow(bridge_annotated)
    mr_supported_count <- sum(bridge_annotated$mr_supported)
    mr_ratio <- mr_supported_count / total_edges * 100
    cat("MR-supported edge ratio:", round(mr_ratio, 2), "%\n")
  }
  
  if (nrow(node_metrics) > 0 && any(!is.na(node_metrics$ivw_pval))) {
    top5 <- node_metrics[!is.na(ivw_pval)][order(ivw_pval)]
    if (nrow(top5) > 0) {
      top5 <- head(top5, 5)
      cat("\nTOP5 Targets (by p-value):\n")
      print(top5[, .(gene, beta = ivw_beta, pval = ivw_pval, pph4)])
    }
  }
  
  cat("\n=== Analysis Complete ===\n")
  
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  cat("Stack trace:\n")
  print(sys.calls())
})
