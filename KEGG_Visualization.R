#!/usr/bin/env Rscript
# KEGG富集结果可视化
# 基于用户提供的KEGG分析结果

library(ggplot2)
library(dplyr)

set.seed(123)

cat("=== 生成KEGG富集图 ===\n")

kegg_results <- data.frame(
  ID = c("hsa04933", "hsa05417", "hsa05171", "hsa05163", "hsa05164",
         "hsa04668", "hsa05321", "hsa05167", "hsa05144", "hsa05200",
         "hsa05323", "hsa05142", "hsa05161", "hsa04066", "hsa05145",
         "hsa04659", "hsa05166", "hsa04621", "hsa05162", "hsa05152",
         "hsa05135", "hsa05140", "hsa04657", "hsa05169", "hsa04064",
         "hsa05146", "hsa04932", "hsa05418", "hsa05235", "hsa04151",
         "hsa05133", "hsa04068", "hsa04060", "hsa04625", "hsa04620",
         "hsa04931", "hsa05168", "hsa05022", "hsa05130", "hsa05206",
         "hsa04936", "hsa05134", "hsa05132", "hsa04218", "hsa05212",
         "hsa05010", "hsa04062", "hsa05205", "hsa01521", "hsa05207",
         "hsa04920", "hsa01523", "hsa05220", "hsa05221", "hsa04061",
         "hsa04672", "hsa05410", "hsa04010", "hsa04917", "hsa05202",
         "hsa04630", "hsa05222", "hsa05131", "hsa04613", "hsa05170",
         "hsa04518", "hsa04926", "hsa04380", "hsa04217", "hsa05415",
         "hsa05165", "hsa04148", "hsa05143", "hsa04623", "hsa05160",
         "hsa05203", "hsa04390", "hsa05211", "hsa05223", "hsa05208",
         "hsa05210", "hsa04923", "hsa04722", "hsa04622", "hsa04071",
         "hsa04913", "hsa04517", "hsa04662", "hsa05416", "hsa04519",
         "hsa04081", "hsa04514", "hsa04024", "hsa04921", "hsa04935",
         "hsa05226", "hsa05225", "hsa05120", "hsa04370", "hsa04145",
         "hsa04670", "hsa00590", "hsa04014", "hsa04726", "hsa05332",
         "hsa04110", "hsa05414", "hsa04550", "hsa05030", "hsa04640",
         "hsa04350", "hsa05215", "hsa04650", "hsa01100", "hsa05150",
         "hsa04210", "hsa05204", "hsa04658", "hsa05020", "hsa04660",
         "hsa04723", "hsa04211"),
  Description = c(
    "AGE-RAGE signaling pathway in diabetic complications",
    "Lipid and atherosclerosis",
    "Coronavirus disease - COVID-19",
    "Human cytomegalovirus infection",
    "Influenza A",
    "TNF signaling pathway",
    "Inflammatory bowel disease",
    "Kaposi sarcoma-associated herpesvirus infection",
    "Malaria",
    "Pathways in cancer",
    "Rheumatoid arthritis",
    "Chagas disease",
    "Hepatitis B",
    "HIF-1 signaling pathway",
    "Toxoplasmosis",
    "Th17 cell differentiation",
    "Human T-cell leukemia virus 1 infection",
    "NOD-like receptor signaling pathway",
    "Measles",
    "Tuberculosis",
    "Yersinia infection",
    "Leishmaniasis",
    "IL-17 signaling pathway",
    "Epstein-Barr virus infection",
    "NF-kappa B signaling pathway",
    "Amoebiasis",
    "Non-alcoholic fatty liver disease",
    "Fluid shear stress and atherosclerosis",
    "PD-L1 expression and PD-1 checkpoint pathway in cancer",
    "PI3K-Akt signaling pathway",
    "Pertussis",
    "FoxO signaling pathway",
    "Cytokine-cytokine receptor interaction",
    "C-type lectin receptor signaling pathway",
    "Toll-like receptor signaling pathway",
    "Insulin resistance",
    "Herpes simplex virus 1 infection",
    "Pathways of neurodegeneration - multiple diseases",
    "Pathogenic Escherichia coli infection",
    "MicroRNAs in cancer",
    "Alcoholic liver disease",
    "Legionellosis",
    "Salmonella infection",
    "Cellular senescence",
    "Pancreatic cancer",
    "Alzheimer disease",
    "Chemokine signaling pathway",
    "Proteoglycans in cancer",
    "EGFR tyrosine kinase inhibitor resistance",
    "Chemical carcinogenesis - receptor activation",
    "Adipocytokine signaling pathway",
    "Antifolate resistance",
    "Chronic myeloid leukemia",
    "Acute myeloid leukemia",
    "Viral protein interaction with cytokine and cytokine receptor",
    "Intestinal immune network for IgA production",
    "Hypertrophic cardiomyopathy",
    "MAPK signaling pathway",
    "Prolactin signaling pathway",
    "Transcriptional misregulation in cancer",
    "JAK-STAT signaling pathway",
    "Small cell lung cancer",
    "Shigellosis",
    "Neutrophil extracellular trap formation",
    "Human immunodeficiency virus 1 infection",
    "Integrin signaling",
    "Relaxin signaling pathway",
    "Osteoclast differentiation",
    "Necroptosis",
    "Diabetic cardiomyopathy",
    "Human papillomavirus infection",
    "Efferocytosis",
    "African trypanosomiasis",
    "Cytosolic DNA-sensing pathway",
    "Hepatitis C",
    "Viral carcinogenesis",
    "Hippo signaling pathway",
    "Renal cell carcinoma",
    "Non-small cell lung cancer",
    "Chemical carcinogenesis - reactive oxygen species",
    "Colorectal cancer",
    "Regulation of lipolysis in adipocytes",
    "Neurotrophin signaling pathway",
    "RIG-I-like receptor signaling pathway",
    "Sphingolipid signaling pathway",
    "Ovarian steroidogenesis",
    "IgSF CAM signaling",
    "B cell receptor signaling pathway",
    "Viral myocarditis",
    "Cadherin signaling",
    "Hormone signaling",
    "Cell adhesion molecule (CAM) interaction",
    "cAMP signaling pathway",
    "Oxytocin signaling pathway",
    "Growth hormone synthesis, secretion and action",
    "Gastric cancer",
    "Hepatocellular carcinoma",
    "Epithelial cell signaling in Helicobacter pylori infection",
    "VEGF signaling pathway",
    "Phagosome",
    "Leukocyte transendothelial migration",
    "Arachidonic acid metabolism",
    "Ras signaling pathway",
    "Serotonergic synapse",
    "Graft-versus-host disease",
    "Cell cycle",
    "Dilated cardiomyopathy",
    "Signaling pathways regulating pluripotency of stem cells",
    "Cocaine addiction",
    "Hematopoietic cell lineage",
    "TGF-beta signaling pathway",
    "Prostate cancer",
    "Natural killer cell mediated cytotoxicity",
    "Metabolic pathways",
    "Staphylococcus aureus infection",
    "Apoptosis",
    "Chemical carcinogenesis - DNA adducts",
    "Th1 and Th2 cell differentiation",
    "Prion disease",
    "T cell receptor signaling pathway",
    "Retrograde endocannabinoid signaling",
    "Longevity regulating pathway"
  ),
  Count = c(6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1),
  stringsAsFactors = FALSE
)

top_kegg <- kegg_results[order(-kegg_results$Count), ]
top_kegg <- head(top_kegg, 30)

top_kegg$Description <- factor(top_kegg$Description,
                               levels = rev(top_kegg$Description))

pdf("5_KEGG_Enrichment_Barplot.pdf", width = 14, height = 10)
ggplot(top_kegg, aes(x = Count, y = Description, fill = Count)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "#377EB8", high = "#E41A1C") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    title = "Hub Gene KEGG Enrichment Analysis (Top 30)",
    x = "Gene Count",
    y = "KEGG Pathway",
    fill = "Count"
  )
dev.off()
cat("Generated: 5_KEGG_Enrichment_Barplot.pdf\n")

pdf("5_KEGG_Enrichment_Bubble.pdf", width = 12, height = 12)
ggplot(top_kegg, aes(x = Count, y = Description, size = Count, color = Count)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 12)) +
  scale_color_gradient(low = "#377EB8", high = "#E41A1C") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    title = "Hub Gene KEGG Enrichment Bubble Chart (Top 30)",
    x = "Gene Count",
    y = "KEGG Pathway",
    size = "Count",
    color = "Count"
  )
dev.off()
cat("Generated: 5_KEGG_Enrichment_Bubble.pdf\n")

inflammation_keywords <- "inflammatory|TNF|IL-|NF-kappa|TLR|JAK|chemokine|cytokine|bowel|arthritis|virus|COVID|influenza|cytomegalovirus|measles|herpes|infection|disease"
inflammation_kegg <- kegg_results[grep(inflammation_keywords, kegg_results$Description, ignore.case = TRUE), ]
inflammation_kegg <- inflammation_kegg[order(-inflammation_kegg$Count), ]
inflammation_kegg <- head(inflammation_kegg, 20)
inflammation_kegg$Description <- factor(inflammation_kegg$Description,
                                       levels = rev(inflammation_kegg$Description))

pdf("5_KEGG_Inflammation_Pathways.pdf", width = 12, height = 8)
ggplot(inflammation_kegg, aes(x = Count, y = Description, fill = Count)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "#4DAF4A", high = "#E41A1C") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    title = "Hub Gene Inflammation/Immune Related Pathways",
    x = "Gene Count",
    y = "KEGG Pathway",
    fill = "Count"
  )
dev.off()
cat("Generated: 5_KEGG_Inflammation_Pathways.pdf\n")

write.table(kegg_results, file = "KEGG_Enrichment_Full_Results.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: KEGG_Enrichment_Full_Results.txt\n")

cat("\n=== Complete ===\n")
cat("Generated files:\n")
cat("  5_KEGG_Enrichment_Barplot.pdf (Bar plot)\n")
cat("  5_KEGG_Enrichment_Bubble.pdf (Bubble chart)\n")
cat("  5_KEGG_Inflammation_Pathways.pdf (Inflammation pathways)\n")
cat("  KEGG_Enrichment_Full_Results.txt (Full results)\n")