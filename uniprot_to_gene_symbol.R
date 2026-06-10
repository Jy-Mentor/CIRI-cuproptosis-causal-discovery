#!/usr/bin/env Rscript
library(org.Hs.eg.db)
library(annotate)

uniprot_ids <- c(
  "ITAL_HUMAN", "BMP2_HUMAN", "P49137", "AK1C2_HUMAN", "APOA2_HUMAN",
  "STS_HUMAN", "ALBU_HUMAN", "KIF11_HUMAN", "P10828", "MK01_HUMAN",
  "CASP7_HUMAN", "TTHY_HUMAN", "VTDB_HUMAN", "PRGR_HUMAN", "P49638",
  "PIM1_HUMAN", "MCR_HUMAN", "ADA17_HUMAN", "ANDR_HUMAN", "P08069",
  "ALDR_HUMAN", "P14061", "Q06520", "AMPM2_HUMAN", "ESR1_HUMAN",
  "EST1_HUMAN", "P62508", "DUS6_HUMAN", "P53041", "CHLE_HUMAN",
  "Q16539", "CAH2_HUMAN", "BACE1_HUMAN", "EPHB4_HUMAN", "SRC_HUMAN",
  "AOFB_HUMAN", "RORA_HUMAN", "PDPK1_HUMAN", "PDE4B_HUMAN", "SHBG_HUMAN",
  "RARG_HUMAN", "PYRD_HUMAN", "PPARG_HUMAN", "P24941", "DHI1_HUMAN",
  "Q9BY41", "DPP4_HUMAN", "THRB_HUMAN", "O00204", "P55263",
  "MDM2_HUMAN", "P49841", "NONE", "Q92731", "PTN11_HUMAN",
  "GSTP1_HUMAN", "Q04828", "EGFR_HUMAN", "AK1C3_HUMAN", "PNMT_HUMAN",
  "NR1I3_HUMAN", "CALM_HUMAN", "PLGF_HUMAN", "RXRA_HUMAN", "P07900",
  "PDE4D_HUMAN", "P08581", "P18031", "P00517", "FABP6_HUMAN",
  "CTNA1_HUMAN", "ADHX_HUMAN", "FABPH_HUMAN", "O76054", "NQO1_HUMAN",
  "PTGD2_HUMAN", "FABP5_HUMAN", "GCR_HUMAN", "RARB_HUMAN", "FNTA_HUMAN",
  "MMP13_HUMAN", "FABP7_HUMAN", "ERBB4_HUMAN", "P08263", "NR1H2_HUMAN",
  "LCK_HUMAN", "PDK2_HUMAN", "O43617", "P16442", "DPEP1_HUMAN",
  "RET4_HUMAN", "A1AT_HUMAN", "JAK3_HUMAN", "WASP_HUMAN", "NR1H3_HUMAN",
  "KSYK_HUMAN", "Q08881", "CCNA2_HUMAN", "P14555", "RXRB_HUMAN",
  "PPARA_HUMAN", "TGFR1_HUMAN", "CASP3_HUMAN", "FA10_HUMAN", "TIE2_HUMAN",
  "Q14541", "PPARD_HUMAN", "RENI_HUMAN", "NR1I2_HUMAN", "P53779",
  "P10827", "ERG7_HUMAN", "MMP2_HUMAN", "P07602", "DCK_HUMAN",
  "KPCT_HUMAN", "S10A9_HUMAN", "ZAP70_HUMAN", "P00517", "NR1H4_HUMAN",
  "CP2C9_HUMAN", "ST1E1_HUMAN", "P22830", "P35968", "P29373",
  "MMP3_HUMAN", "Q13126", "NGAL_HUMAN", "MP2K1_HUMAN", "VDR_HUMAN",
  "P50135", "HMDH_HUMAN", "O14965", "FKB1A_HUMAN", "TGM3_HUMAN",
  "P00519", "BRAF1_HUMAN", "P08238", "CHK1_HUMAN", "RARA_HUMAN",
  "IL2_HUMAN", "FGFR1_HUMAN", "P85A_HUMAN", "CATB_HUMAN", "KIT_HUMAN",
  "P36873", "Q07817", "GRB2_HUMAN", "Q9UKL6", "PUR2_HUMAN",
  "GSTT2_HUMAN", "EPCR_HUMAN", "AOFA_HUMAN"
)

result_df <- data.frame(
  Original = uniprot_ids,
  GeneSymbol = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(uniprot_ids)) {
  uid <- uniprot_ids[i]

  if (uid == "NONE" || uid == "") {
    result_df$GeneSymbol[i] <- "NONE"
    next
  }

  is_entry_name <- grepl("_HUMAN$", uid)
  is_pure_uniprot <- nchar(uid) == 6 && substr(uid, 1, 1) %in% c("P", "Q", "O")

  if (is_entry_name) {
    clean <- gsub("_HUMAN$", "", uid)
    mapped <- tryCatch({
      mapIds(org.Hs.eg.db, keys=clean, keytype="UNIPROT", column="SYMBOL")
    }, error = function(e) NA)

    if (!is.na(mapped) && mapped != clean) {
      result_df$GeneSymbol[i] <- as.character(mapped)
    } else {
      result_df$GeneSymbol[i] <- clean
    }
  } else if (is_pure_uniprot) {
    mapped <- tryCatch({
      mapIds(org.Hs.eg.db, keys=uid, keytype="UNIPROT", column="SYMBOL")
    }, error = function(e) NA)

    if (!is.na(mapped)) {
      result_df$GeneSymbol[i] <- as.character(mapped)
    } else {
      result_df$GeneSymbol[i] <- uid
    }
  } else {
    result_df$GeneSymbol[i] <- uid
  }
}

write.table(result_df, file = "c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/uniprot_to_genesymbol.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

unmapped_ids <- result_df$Original[result_df$GeneSymbol == result_df$Original &
                                    result_df$Original != "NONE" &
                                    !grepl("_HUMAN$", result_df$Original)]
if (length(unmapped_ids) > 0) {
  writeLines(unmapped_ids, "c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/unmapped_uniprot.txt")
}

cat("=== 转换结果 ===\n\n")
for (i in 1:nrow(result_df)) {
  cat(sprintf("%s\t%s\n", result_df$Original[i], result_df$GeneSymbol[i]))
}

cat("\n=== 统计信息 ===\n")
cat(sprintf("总ID数: %d\n", nrow(result_df)))
mapped_count <- sum(result_df$GeneSymbol != "NONE" &
                     (result_df$GeneSymbol != result_df$Original |
                      grepl("_HUMAN$", result_df$Original)))
cat(sprintf("成功映射: %d\n", mapped_count))
cat(sprintf("保持原样: %d\n", nrow(result_df) - mapped_count))

cat("\n结果已保存到: c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/uniprot_to_genesymbol.txt\n")
if (length(unmapped_ids) > 0) {
  cat("未映射ID已保存到: c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/unmapped_uniprot.txt\n")
}
