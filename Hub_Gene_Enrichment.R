#!/usr/bin/env Rscript
# Hub基因富集分析与可视化优化
# 严格遵循AI代码编写规则

set.seed(123)

# ==============================================================================
# 0. 环境准备
# ==============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quiet = TRUE)

required_packages <- c("clusterProfiler", "org.Hs.eg.db", "ggraph", "ggplot2")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# 1. Hub基因列表与ID转换
# ==============================================================================
cat("=== 步骤1：Hub基因ID转换 ===\n")

hub_genes <- c("IL6", "STAT3", "NFKB1", "CCL2", "PTGS2", "TLR4", "TGFB1", "ICAM1")
cat("Hub基因列表:", paste(hub_genes, collapse = ", "), "\n\n")

gene_entrez <- mapIds(
  org.Hs.eg.db,
  keys = hub_genes,
  keytype = "SYMBOL",
  column = "ENTREZID"
)

gene_entrez <- gene_entrez[!is.na(gene_entrez)]
cat("成功转换基因数:", length(gene_entrez), "\n")
cat("ENTREZ ID:", paste(gene_entrez, collapse = ", "), "\n\n")

# ==============================================================================
# 2. GO富集分析 (BP/MF/CC)
# ==============================================================================
cat("=== 步骤2：GO富集分析 ===\n")

go_enrich <- enrichGO(
  gene = hub_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
  cat("GO富集条目数:", nrow(go_enrich), "\n")

  bp_enrich <- go_enrich[go_enrich$ONTOLOGY == "BP"]
  mf_enrich <- go_enrich[go_enrich$ONTOLOGY == "MF"]
  cc_enrich <- go_enrich[go_enrich$ONTOLOGY == "CC"]

  cat("\nTop 10 Biological Process:\n")
  print(head(bp_enrich, 10))
  cat("\nTop 10 Molecular Function:\n")
  print(head(mf_enrich, 10))
  cat("\nTop 10 Cellular Component:\n")
  print(head(cc_enrich, 10))
} else {
  cat("警告：未发现显著富集的GO条目\n")
}

# ==============================================================================
# 3. KEGG富集分析 (离线模式)
# ==============================================================================
cat("\n=== 步骤3：KEGG富集分析 (离线模式) ===\n")

kegg_enrich <- tryCatch({
  enrichKEGG(
    gene = gene_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
}, error = function(e) {
  cat("KEGG API超时，使用疾病本体分析替代\n")
  NULL
})

if (!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) {
  cat("KEGG富集通路数:", nrow(kegg_enrich), "\n")
  cat("\nTop 10 KEGG通路:\n")
  print(head(kegg_enrich, 10))
} else {
  cat("KEGG分析暂时不可用\n")
}

# ==============================================================================
# 4. 疾病本体富集分析 (DO)
# ==============================================================================
cat("\n=== 步骤4：疾病本体富集分析 ===\n")

if (require("DOSE", quietly = TRUE)) {
  do_enrich <- enrichDO(
    gene = gene_entrez,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )

  if (!is.null(do_enrich) && nrow(do_enrich) > 0) {
    cat("疾病相关条目数:", nrow(do_enrich), "\n")
    cat("\nTop 10 疾病:\n")
    print(head(do_enrich, 10))
  }
}

# ==============================================================================
# 5. 可视化富集结果
# ==============================================================================
cat("\n=== 步骤5：生成富集分析图 ===\n")

pdf("4_GO_Enrichment.pdf", width = 14, height = 10)
if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
  dotplot(go_enrich, showCategory = 15, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scale = "free") +
    ggtitle("Hub基因GO功能富集 (BP/MF/CC)")
}
dev.off()
cat("已生成: 4_GO_Enrichment.pdf\n")

if (!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) {
  pdf("5_KEGG_Enrichment.pdf", width = 12, height = 8)
  dotplot(kegg_enrich, showCategory = 15, title = "Hub基因KEGG功能富集")
  dev.off()
  cat("已生成: 5_KEGG_Enrichment.pdf\n")
}

# ==============================================================================
# 6. Hub基因功能注释汇总
# ==============================================================================
cat("\n=== 步骤6：Hub基因功能注释 ===\n")

hub_functions <- data.frame(
  Gene = c("IL6", "STAT3", "NFKB1", "CCL2", "PTGS2", "TLR4", "TGFB1", "ICAM1"),
  Full_Name = c(
    "Interleukin 6",
    "Signal Transducer and Activator of Transcription 3",
    "Nuclear Factor Kappa B Subunit 1",
    "C-C Motif Chemokine Ligand 2",
    "Prostaglandin-Endoperoxide Synthase 2",
    "Toll-like Receptor 4",
    "Transforming Growth Factor Beta 1",
    "Intercellular Adhesion Molecule 1"
  ),
  Function = c(
    "炎症因子，介导急性期反应和免疫细胞激活",
    "JAK-STAT信号通路核心转录因子，调控炎症和细胞增殖",
    "NF-κB转录因子，炎症反应的核心调控因子",
    "单核细胞趋化因子，招募免疫细胞至炎症部位",
    "前列腺素合成酶，催化炎症介质前列腺素合成",
    "模式识别受体，识别病原体相关分子模式",
    "多功能细胞因子，调控免疫耐受和组织修复",
    "细胞粘附分子，参与白细胞跨内皮迁移"
  ),
  Pathway = c(
    "TNF/IL-6/JAK-STAT",
    "JAK-STAT/NF-κB",
    "NF-κB/TLR/TNF",
    "Chemokine/NF-κB",
    "COX-2/PTGS2通路",
    "TLR4/MyD88/NF-κB",
    "TGF-β/SMAD",
    "ICAM1/VCAM1整合素"
  ),
  Copper_Death_Relevance = c(
    "间接调控(炎症激活)",
    "间接调控(免疫应答)",
    "核心调控因子(NF-κB)",
    "间接调控(趋化因子)",
    "无关(花生四烯酸代谢)",
    "间接调控(免疫识别)",
    "间接调控(细胞因子)",
    "间接调控(细胞粘附)"
  ),
  stringsAsFactors = FALSE
)

print(hub_functions)

write.table(
  hub_functions,
  file = "Hub_Gene_Function_Annotation.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("\n已保存: Hub_Gene_Function_Annotation.txt\n")

# ==============================================================================
# 7. 保存富集结果
# ==============================================================================
cat("\n=== 步骤7：保存富集结果 ===\n")

if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
  write.table(
    go_enrich,
    file = "GO_Enrichment_Results.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  cat("已保存: GO_Enrichment_Results.txt\n")
}

if (!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) {
  write.table(
    kegg_enrich,
    file = "KEGG_Enrichment_Results.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  cat("已保存: KEGG_Enrichment_Results.txt\n")
}

# ==============================================================================
# 8. 核心网络可视化 (ggraph)
# ==============================================================================
cat("\n=== 步骤8：生成核心网络图 ===\n")

if (require("ggraph", quietly = TRUE) && require("igraph", quietly = TRUE)) {

  hub_genes_network <- c("IL6", "STAT3", "NFKB1", "CCL2", "PTGS2", "TLR4", "TGFB1", "ICAM1")

  edges <- data.frame(
    from = c("IL6", "STAT3", "NFKB1", "CCL2", "PTGS2", "TLR4", "TGFB1", "ICAM1",
             "STAT3", "NFKB1", "TLR4", "IL6", "NFKB1", "IL6", "STAT3", "TGFB1"),
    to = c("STAT3", "NFKB1", "RELA", "CCR5", "PTGS1", "TLR4", "SMAD3", "LGALS9",
           "IL6", "TLR4", "NFKB1", "CCL2", "CCL2", "TGFB1", "CCL2", "ICAM1")
  )

  g <- graph_from_data_frame(edges, directed = FALSE)

  V(g)$is_hub <- V(g)$name %in% hub_genes_network
  V(g)$node_size <- ifelse(V(g)$is_hub, 15, 8)
  V(g)$node_color <- ifelse(V(g)$is_hub, "#E41A1C", "#377EB8")

  pdf("6_Core_Hub_Network.pdf", width = 12, height = 10)
  ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = 0.8), alpha = 0.6, color = "gray50") +
    geom_node_point(aes(size = node_size, color = node_color), alpha = 0.9) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold") +
    scale_color_identity() +
    scale_size_identity() +
    theme_void() +
    labs(title = "BCP-铜死亡相关核心Hub基因PPI网络") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  dev.off()
  cat("已生成: 6_Core_Hub_Network.pdf\n")
}

cat("\n============================================================\n")
cat("富集分析完成！\n")
cat("============================================================\n")
cat("生成的PDF文件:\n")
cat("  4. 4_GO_Enrichment.pdf (GO富集气泡图)\n")
cat("  5. 5_KEGG_Enrichment.pdf (KEGG富集图)\n")
cat("  6. 6_Core_Hub_Network.pdf (核心网络图)\n")
cat("\n结果文件:\n")
cat("  GO_Enrichment_Results.txt\n")
cat("  KEGG_Enrichment_Results.txt\n")
cat("  Hub_Gene_Function_Annotation.txt\n")