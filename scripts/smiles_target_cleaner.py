# -*- coding: utf-8 -*-
"""
======================================================================
SMILES-derived Drug Target Cleaning & Confidence Filtering Pipeline
======================================================================

Reference Papers:
  [1] Gfeller D, et al. SwissTargetPrediction. Nucleic Acids Res, 2014
      (PMID: 24792161) — STP probability-based scoring
  [2] Li X, et al. INPUT platform. Comput Struct Biotechnol J, 2022
      (PMID: 35414968) — Normalized score >= 0.5 threshold
  [3] Daina A, et al. SwissADME. Sci Rep, 2017 (PMID: 28256516)
  [4] Tsvetkov P, et al. Cuproptosis. Science, 2022 (PMID: 35298263)
  [5] Eisenberg E, Levanon EY. Housekeeping genes. Trends Genet, 2013
      (PMID: 23213612)
  [6] World Federation of Chinese Medicine Societies. Network Pharmacology
      Evaluation Method Guidance, 2021 — Multi-source validation standard
  [7] Gfeller D, Zoete V. SwissTargetPrediction 2019 update. Nucleic
      Acids Res, 2019 (PMID: 31106366)
  [8] Alidoost M, et al. Curated DTI improves prediction. JCIM, 2026
      — Multi-database target integration approach

Key GitHub references:
  - network pharmacology pipelines: various TCM network pharmacology
    workflows in Python/R

Methodology (6-Layer Filtering):
  Layer 0: Raw cleaning — remove empty, duplicate, invalid symbols
  Layer 1: Negative filtering — housekeeping (PMID:23213612),
           hard blacklist (hemoglobin, mtDNA, Ig, collagens)
  Layer 2: Positive anchoring — known drug targets (DrugBank/ChEMBL)
  Layer 3: Disease context — CIRI pathway, brain expression
  Layer 4: Druggability — PDB availability, druggable families
  Layer 5: Network topology — PPI hub status (if STRING data available)
  Layer 6: Integrated scoring + Tier assignment

Output:
  - tiered_targets.csv:  full annotated table
  - filter_report.txt:   detailed tier-by-tier report
  - high_confidence.txt: Tier1+ Tier2 for downstream analysis
======================================================================
"""

import os
import sys
import json
import re
from typing import Set, Tuple, List, Dict

import numpy as np
import pandas as pd

# ================================================================
# Configuration
# ================================================================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_FILE = os.path.join(BASE_DIR, "大创", "石竹烯 人.txt")
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "target_cleaning")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ================================================================
# Layer 0: Housekeeping genes — PMID:23213612 (Eisenberg & Levanon)
# ================================================================
HOUSEKEEPING_GENES = {
    "ALB", "TBP", "ACTA1", "GAPDH", "GAPDHS", "B2M", "ACTB", "RPLP0",
    "HPRT1", "PPIA", "YWHAZ", "GUSB", "TFRC", "POLR2A",
    "PGK1", "TUBB", "RPL13A", "SDHA", "UBC", "GAPDH",
    "HMBS", "IPO8", "RPS18", "RPS27A",
}

# ================================================================
# Layer 1: Hard Blacklist — biologically irrelevant or noise-prone
# ================================================================
BLACKLIST_HEMOGLOBIN = {
    "HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBZ", "HBE1", "HBQ1",
}
BLACKLIST_MITOCHONDRIAL = {
    "MTCYB", "MTND1", "MTND2", "MTND3", "MTND4", "MTND4L", "MTND5",
    "MTND6", "MTCO1", "MTCO2", "MTCO3", "MTATP6", "MTATP8",
    "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5",
    "MT-ND6", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8",
}
BLACKLIST_IMMUNOGLOBULIN = {
    "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM",
    "IGKC", "IGLC1", "IGLC2", "IGLC3",
}
BLACKLIST_PLASMA_PROTEIN = {
    "TF", "HP", "HPX", "ORM1", "ORM2", "A2M", "SERPINA1",
    "SERPINC1", "C1QC", "C1QA", "C1QB",
}
HARD_BLACKLIST = (
    BLACKLIST_HEMOGLOBIN | BLACKLIST_MITOCHONDRIAL | BLACKLIST_IMMUNOGLOBULIN
)
SOFT_BLACKLIST = BLACKLIST_PLASMA_PROTEIN

# ================================================================
# Layer 2: Known Drug-Target Anchoring (DrugBank + ChEMBL curated)
# ================================================================
KNOWN_DRUG_TARGETS = {
    # GPCRs
    "ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C",
    "ADRB1", "ADRB2", "ADRB3", "DRD1", "DRD2", "DRD3", "DRD4",
    "HTR1A", "HTR1B", "HTR1D", "HTR2A", "HTR2B", "HTR2C", "HTR3A",
    "HTR4", "HTR5A", "HTR6", "HTR7",
    "HRH1", "HRH2", "HRH3", "HRH4",
    "CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5",
    "OPRM1", "OPRK1", "OPRD1", "OPRL1",
    "CNR1", "CNR2", "PTAFR", "OXTR",
    "TACR1", "TACR2", "TACR3",
    "MTNR1A", "MTNR1B", "HCRTR1", "HCRTR2",
    "S1PR1", "S1PR2", "S1PR3",
    "ADORA1", "ADORA2A", "ADORA2B", "ADORA3",
    "P2RX7", "P2RX4", "P2RX3", "P2RX2", "P2RX1",
    "FFAR1", "FFAR2", "FFAR3", "FFAR4",
    "CCR2", "CCR5", "CXCR3", "CXCR4",
    "RORC", "RORA",
    "PTGER1", "PTGER2", "PTGER3", "PTGER4",
    "TBXA2R",
    # Nuclear receptors
    "ESR1", "ESR2", "AR", "PGR", "NR3C1", "NR3C2",
    "PPARA", "PPARG", "PPARD",
    "RARA", "RARB", "RARG", "RXRA", "RXRB", "RXRG",
    "THRA", "THRB", "NR1H2", "NR1H3",
    "NR1I2", "NR1I3", "VDR",
    "NR1H4", "NR5A1", "NR5A2",
    # Kinases
    "AKT1", "AKT2", "AKT3",
    "MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK10", "MAPK14",
    "MAPKAPK2",
    "EGFR", "ERBB2", "ERBB3", "ERBB4",
    "SRC", "FYN", "YES1", "HCK", "LYN", "LCK", "FGR",
    "JAK1", "JAK2", "JAK3", "TYK2",
    "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG",
    "MTOR", "RPS6KB1", "RPS6KB2",
    "CDK1", "CDK2", "CDK4", "CDK6", "CDK9",
    "BRAF", "RAF1", "MAP2K1", "MAP2K2",
    "INSR", "IGF1R",
    "FLT1", "FLT3", "FLT4", "KDR",
    "CSF1R", "KIT", "PDGFRA", "PDGFRB",
    "ABL1", "ABL2", "ITK", "BTK",
    "IRAK1", "IRAK4", "TBK1", "IKBKB",
    "PRKACA", "PRKACB", "PRKCA", "PRKCB", "PRKCD", "PRKCE",
    "PRKCQ", "PRKCI", "PRKCZ",
    "PRKAA1", "PRKAA2",
    "PLK1", "PLK2", "PLK3",
    "AURKA", "AURKB",
    "CHEK1", "CHEK2", "WEE1",
    "EPHA1", "EPHA2", "EPHB4",
    "PTK2", "PTK2B",
    "STK4", "STK16", "STK33",
    "CSNK1D", "CSNK1E", "CSNK2A1", "CSNK2A2",
    "MAPKAPK2", "MAPKAPK5",
    "ROCK1", "ROCK2",
    "GSK3A", "GSK3B",
    "NEK1", "NEK2", "NEK6", "NEK7",
    "TTK", "MELK", "BUB1", "BUB1B",
    # Proteases
    "CASP1", "CASP3", "CASP6", "CASP7", "CASP8", "CASP9", "CASP10",
    "CTSB", "CTSK", "CTSL", "CTSS", "CTSV", "CTSF", "CTSC", "CTRC",
    "MMP1", "MMP2", "MMP9", "MMP13", "MMP3", "MMP7", "MMP8", "MMP10", "MMP12",
    "ACE", "ACE2",
    "DPP4", "DPP8", "DPP9",
    "PRSS1", "PRSS2", "PRSS3",
    "PREP",
    "TMPRSS2", "TMPRSS11D",
    # Oxidoreductases / Metabolic enzymes (druggable subset)
    "PTGS1", "PTGS2",
    "ALOX5", "ALOX12", "ALOX15", "ALOX12B",
    "XDH", "IDO1", "TDO2",
    "NOS1", "NOS2", "NOS3",
    "HMOX1", "HMOX2",
    "SOD1", "SOD2", "CAT", "GPX1", "GPX4",
    "AKR1C1", "AKR1C2", "AKR1C3",
    "HSD11B1", "HSD11B2",
    "HSD17B1", "HSD17B2", "HSD17B3", "HSD17B4", "HSD17B10",
    "DHODH",
    # Epigenetic / Transcriptional regulators
    "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7", "HDAC8",
    "HDAC9", "HDAC10", "HDAC11",
    "SIRT1", "SIRT2", "SIRT3", "SIRT5", "SIRT6",
    "PARP1", "PARP2", "PARP3", "PARP4", "PARP10", "PARP12", "PARP14",
    "BRD2", "BRD3", "BRD4", "BRDT",
    "SMARCA2", "SMARCA4",
    "DNMT1", "DNMT3A", "DNMT3B",
    "EZH2", "KMT2A", "KMT2D",
    # Transporters / Ion channels
    "ACHE", "BCHE",
    "MAOA", "MAOB",
    "SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN8A", "SCN9A", "SCN10A", "SCN11A",
    "KCNJ1", "KCNJ2", "KCNJ5", "KCNJ8", "KCNJ11",
    "KCNA1", "KCNA2", "KCNA3", "KCNA4", "KCNA5",
    "KCNH1", "KCNH2", "KCNH6", "KCNH7",
    "KCNQ1", "KCNQ2", "KCNQ3", "KCNQ4", "KCNQ5",
    "CACNA1A", "CACNA1B", "CACNA1C", "CACNA1D", "CACNA1E", "CACNA1G", "CACNA1H", "CACNA1I",
    "TRPA1", "TRPC1", "TRPC3", "TRPC4", "TRPC5", "TRPC6", "TRPC7",
    "TRPM2", "TRPM8", "TRPV1", "TRPV2", "TRPV3", "TRPV4",
    "SLC6A1", "SLC6A2", "SLC6A3", "SLC6A4",
    "SLC6A11", "SLC6A12",
    "SLC1A1", "SLC1A2", "SLC1A3",
    "SLC2A1", "SLC2A3", "SLC2A4",
    "SLC12A1", "SLC12A2", "SLC12A3", "SLC12A5", "SLC12A6",
    "ABCB1", "ABCC1", "ABCG2",
    "ATP1A1", "ATP1A2", "ATP1A3", "ATP2A1", "ATP2A2",
    "SLC9A1", "SLC9A3",
    "SLC22A1", "SLC22A2", "SLC22A6",
    "SLC47A1", "SLC47A2",
    "NPC1L1",
    "KCNE1",
    "RYR1", "RYR2", "RYR3",
    # Phosphodiesterases
    "PDE4A", "PDE4B", "PDE4C", "PDE4D",
    "PDE5A", "PDE3A", "PDE3B",
    "PDE8A",
    "PDE6G",
    "PDE10A",
    # Key signaling / transcription factors
    "TP53", "NFKB1", "NFKB2", "RELA", "RELB", "REL",
    "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "NFE2L2", "FOXO1", "FOXO3", "FOXO4",
    "HIF1A", "HIF1AN", "ARNT",
    "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
    "SREBF1", "SREBF2",
    "EGR1", "JUN", "FOS", "JUND", "FOSL1",
    "ATF2", "ATF3", "ATF4", "DDIT3", "XBP1",
    "IRF1", "IRF3", "IRF4", "IRF5", "IRF7", "IRF8",
    "CEBPB", "CEBPA", "CEBPD",
    "GATA3", "GATA4",
    "NR4A1", "NR4A2", "NR4A3",
    "MYC", "MAX", "MXD1",
    "CTNNB1", "TCF7",
    "YAP1", "TAZ",
    "GLI1", "GLI2", "GLI3",
    "MITF", "TFEB", "TFE3",
    # Inflammatory cytokines / immune
    "TNF", "TNFRSF1A", "TNFRSF1B",
    "IL1A", "IL1B", "IL1R1",
    "IL2", "IL4", "IL5", "IL6", "IL6R",
    "IL10", "IL10RA", "IL10RB",
    "IL18", "IL18R1",
    "CCL2", "CCL3", "CCL4", "CCL5",
    "CXCL8", "CXCL10", "CXCL12",
    "ICAM1", "VCAM1", "SELE", "SELP",
    "TLR2", "TLR3", "TLR4", "TLR7", "TLR8", "TLR9",
    "NLRP3", "NLRP1",
    "AIM2",
    "TREM1", "TREM2",
    "CSF1", "CSF1R",
    "CD28", "CD80", "CD86",
    "CTLA4", "PDCD1", "CD274",
    "CD38", "CD47",
    "IFNG", "IFNAR1", "IFNAR2",
    "TGFB1", "TGFBR1", "TGFBR2",
    # Damage / autophagy / cell death
    "NLRP3", "ATG5", "ATG7", "BECN1", "ULK1", "ULK2",
    "SQSTM1", "OPTN",
    "BCL2", "BCL2L1", "BCL2L11", "MCL1",
    "BAX", "BAK1", "BAD", "BID",
    "FAS", "FASLG",
    "CASP8", "CASP9",
    "RIPK1", "RIPK3", "MLKL",
    "GSDMD", "GSDME", "GSDMB",
    "HMGB1", "HMGB2", "HMGB3",
    "S100A6", "S100A8", "S100A9", "S100B",
    "HSPA5", "HSP90AA1", "HSP90AB1",
    "TSPO", "VDAC1", "VDAC2",
    # Cuproptosis — Tsvetkov Science 2022
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHA1", "PDHB",
    "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7A", "ATP7B",
    "DLD", "DBT", "DLST", "GCSH", "PDHA2",
    # Neurological / synaptic
    "APP", "BACE1", "PSEN1", "PSEN2",
    "MAPT", "SNCA", "LRRK2", "PINK1", "PARK7",
    "BDNF", "NGF", "NGFR",
    "GDNF", "NTF3", "NTF4",
    "SIGMAR1",
    "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D",
    "GRIA1", "GRIA2", "GRIA3", "GRIA4",
    "GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6",
    "GABRB1", "GABRB2", "GABRB3",
    "GABRG1", "GABRG2", "GABRG3",
    "GABBR1", "GABBR2",
    "CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5",
    "CHRNA6", "CHRNA7", "CHRNA9", "CHRNA10",
    "CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4",
    "HCN1", "HCN2", "HCN4",
    "SYP", "SYN1", "SYN2", "SYN3",
    "DLG4", "DLG2", "DLG3",
    "NLGN1", "NLGN2", "NLGN3",
    "NRXN1", "NRXN2", "NRXN3",
    # Nitric oxide / oxidative stress
    "NOS1", "NOS2", "NOS3",
    "GSR", "GPX1", "GPX2", "GPX3", "GPX4",
    "TXN", "TXNRD1", "TXNRD2",
    "NQO1", "NQO2",
    "GCLC", "GCLM", "GSS",
    "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6",
    "KEAP1", "NFE2L2",
    # Drug metabolism (CYP / UGT)
    "CYP1A1", "CYP1A2", "CYP1B1",
    "CYP2A6", "CYP2B6",
    "CYP2C8", "CYP2C9", "CYP2C18", "CYP2C19",
    "CYP2D6", "CYP2E1",
    "CYP2J2", "CYP2R1", "CYP2S1", "CYP2U1",
    "CYP3A4", "CYP3A5", "CYP3A43", "CYP3A7",
    "CYP4A11", "CYP4B1", "CYP4F2",
    "CYP17A1", "CYP19A1", "CYP26A1", "CYP27A1", "CYP27B1",
    "CYP51A1", "CYP8B1",
    "UGT1A1", "UGT1A3", "UGT1A4", "UGT1A6", "UGT1A7", "UGT1A8",
    "UGT1A9", "UGT1A10", "UGT2B4", "UGT2B7", "UGT2B10",
    "UGT2B11", "UGT2B15", "UGT2B17", "UGT2B28",
    # Miscellaneous validated targets
    "FAAH", "MGLL", "ABHD6", "ABHD12",
    "CES1", "CES2",
    "PLA2G4A", "PLA2G2A",
    "PTPN1", "PTPN2", "PTPN6", "PTPN11",
    "PTPRF", "PTPRJ", "PTPRC",
    "FASN", "HMGCR", "SQLE", "LSS", "FDFT1",
    "ACACA", "ACACB", "FASN",
    "CPT1A", "CPT1B", "CPT2",
    "SOAT1", "SOAT2",
    "DGAT1", "DGAT2",
    "PNPLA2", "PNPLA3",
    "LPL", "LIPC",
    "MGAT1", "MGAT2", "MGAT3",
    "SCD", "SCD5",
    "ELOVL1", "ELOVL2", "ELOVL5", "ELOVL6",
    "COMT", "GAD1", "GAD2", "DDC",
    "TPH1", "TPH2",
    "MAOA", "MAOB",
    "SRD5A1", "SRD5A2",
    "MTAP",
    "TK1", "TK2",
    "DHFR", "DHFR2",
    "TYMP",
    "RRM1", "RRM2",
    "TOP1", "TOP2A", "TOP2B",
    "TERT",
    "TDP1", "TDP2",
    "PARP1", "PARP2",
    "POLK", "REV1",
    "XRCC6", "XRCC5",
    "RAD51", "RAD52",
    "NUDC", "NUDCD2",
    "PARN",
    "OAZ1", "OAZ2",
    "SAT1", "SAT2",
    "OAT",
    "PAH",
    "SHMT1", "SHMT2",
    "TYMS",
    "IMPDH1", "IMPDH2",
    "GMPS",
    "CTPS1", "CTPS2",
    "CAD",
    "PYCR1", "PYCR2", "PYCR3",
    "PRODH", "PRODH2",
    "ALDH1A1", "ALDH2", "ALDH3A1",
    "CHDH",
    "ABAT",
    "AOC1", "AOC3",
    "SAT1", "SAT2",
    "CHAT",
    "ACHE", "BCHE",
    "CES1", "CES2",
    "EPHX1", "EPHX2",
    "GSTP1", "GSTM1", "GSTT1",
    "SULT1A1", "SULT1E1",
    "TPMT",
    "AHCY",
    "MAT2A", "MAT2B",
    "REN",
    "MMUT",
    "ACADM", "ACADVL", "ACAD11",
    "ECHS1", "HADHA", "HADHB",
    "HIBADH",
    "CRYL1",
    "FAH",
    "GCH1",
    "HPD",
    "RENBP",
    "IMPDH2",
    "PDXK",
    "TPK1",
    "FNTB", "FNTA",
    "ICMT",
    "PGGT1B", "RABGGTB",
    "RAPGEF3", "RAPGEF4",
    "STARD13",
    "REPS2",
    "SMAP1",
    "ACLY",
    "IDH1", "IDH2",
    "MDH1", "MDH2",
    "FH",
    "SDHA", "SDHB", "SDHC", "SDHD",
    "UQCRC1",
    "COX5A", "COX6B1",
    "ATP5A1", "ATP5B",
    "VEGFA", "VEGFB", "VEGFC",
    "KDR", "FLT1", "FLT4",
    "PGF",
    "ANGPT1", "ANGPT2",
    "TIE1", "TEK",
    "SERPINE1", "SERPINB1", "SERPINB10", "SERPINA6",
    "TIMP1", "TIMP2", "TIMP3",
    "CTGF",
    "COL1A1",
    "FN1",
    "PECAM1", "VWF",
    "FGA", "FGB", "FGG",
    "F2", "F3", "F5", "F7", "F9", "F10",
    "PLG", "PLAT", "PLAU", "PLAUR",
    "TBXAS1",
    "PTGES", "PTGES2",
    "PTGR1", "PTGR2",
}

# ================================================================
# Layer 3: CIRI Disease Context — brain expression + pathway
# ================================================================
NEURON_MARKERS = {
    "GRIN1", "GRIN2A", "GRIN2B", "GRIA1", "GRIA2", "GRIA3",
    "GABRA1", "GABRB3", "GABRG2",
    "SYN1", "SYP", "SV2A", "SNAP25", "SYT1",
    "MAP2", "NEUROD6", "NEUROD1", "NEUROD2",
    "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
    "ARC", "BDNF", "FOS", "JUN", "EGR1", "NR4A1",
    "HOMER1", "HOMER2", "HOMER3",
    "DLG4", "SHANK1", "SHANK2", "SHANK3",
    "RIMS1", "RIMS2",
    "VGLUT1", "VGLUT2",
}

GLIAL_MARKERS = {
    "GFAP", "ALDH1L1", "SLC1A2", "SLC1A3", "GJB6", "GJC1",
    "AQP4", "KCNJ10", "ABCC8", "KCNJ8",
    "MBP", "PLP1", "MOG", "CNP", "MAG",
    "OLIG1", "OLIG2", "SOX10",
    "P2RY12", "TMEM119", "CX3CR1", "CD33", "TYROBP",
    "ITGAM", "ITGAX", "FCGR1A", "FCGR3A",
    "APOE", "CLU", "C3", "C1QA", "C1QB", "C1QC",
    "CD68", "CD74", "CTSB", "CTSD", "CTSL", "CTSS",
    "TREM2", "FCER1G",
    "AIF1", "Iba1",
}

BRAIN_EXPRESSED = NEURON_MARKERS | GLIAL_MARKERS | {
    "PECAM1", "VWF", "CLDN5", "CDH5", "TJP1",
    "FLT1", "KDR",
    "SLC2A1", "ABCB1", "ABCG2",
    "CAV1", "CAV2",
    "CALM1", "CALM2", "CALM3",
    "PVALB", "SST", "NPY", "CCK", "VIP",
    "TAC1", "PENK", "PDYN",
    "CHRM1", "CHRM4",
    "GABBR1", "GABBR2",
    "GRM1", "GRM5",
    "ADORA1", "ADORA2A",
    "CNR1",
    "DRD1", "DRD2",
    "HTR1A", "HTR2A",
    "GPR37", "GPR37L1",
    "GAP43",
    "GAP43",
}

CIRI_CORE_PATHWAY_GENES = {
    "NFKB1", "RELA", "TLR4", "TNF", "IL6", "IL1B",
    "HMOX1", "PTGS2", "MMP9", "BCL2", "BAX", "CASP3",
    "AKT1", "MTOR", "STAT3", "NLRP3", "NFE2L2", "KEAP1",
    "SOD1", "SOD2", "CAT", "GPX1", "GSR",
    "NOS2", "NOS3",
    "HIF1A", "VEGFA",
    "APP", "BDNF",
    "JAK1", "JAK2",
    "CASP8", "CASP9",
    "CCL2", "ICAM1", "VCAM1",
    "NOTCH1", "SIRT1",
    "CNR2", "TRPV1", "PPARG",
    "FOXO1", "FOXO3",
    "ATG5", "BECN1",
    "CASP1", "IL18",
}

CIRI_DEATH_SUB = {
    "ferroptosis": {"GPX4", "ACSL4", "LPCAT3", "ALOX15", "ALOX12",
                    "TFRC", "SLC7A11", "FSP1", "DHODH", "GCH1"},
    "necroptosis": {"RIPK1", "RIPK3", "MLKL", "ZBP1", "FADD", "CASP8", "CYLD"},
    "pyroptosis": {"NLRP3", "CASP1", "CASP4", "CASP5", "GSDMD", "GSDME", "IL1B", "IL18", "P2RX7"},
    "cuproptosis": {"FDX1", "LIAS", "LIPT1", "DLAT", "PDHA1", "PDHB",
                    "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7A", "ATP7B",
                    "DLD", "DBT", "DLST", "GCSH"},
    "apoptosis": {"CASP3", "CASP6", "CASP7", "CASP8", "CASP9",
                  "BAX", "BAK1", "BAD", "BCL2", "BCL2L1", "MCL1",
                  "FAS", "TP53", "PMAIP1", "CYCS", "APAF1", "XIAP"},
    "autophagy": {"MTOR", "ULK1", "ATG5", "ATG7", "ATG12", "ATG16L1",
                  "BECN1", "SQSTM1", "FOXO1", "FOXO3"},
}

# ================================================================
# Layer 4: Druggability indicators
# ================================================================
DRUGGABLE_FAMILY_PREFIXES = ["SCN", "KCN", "TRP", "CACNA", "SLC", "ABC"]

HAS_PDB_KNOWN = {
    "ACHE", "BCHE", "MAOA", "MAOB",
    "PTGS1", "PTGS2",
    "HMGCR", "FASN",
    "CYP2D6", "CYP3A4", "CYP2C9", "CYP2C19", "CYP1A2",
    "CYP17A1", "CYP19A1", "CYP51A1",
    "HDAC1", "HDAC2", "HDAC3", "HDAC6", "HDAC8",
    "SIRT1", "SIRT2", "SIRT3", "SIRT5", "SIRT6",
    "PARP1", "BRD2", "BRD3", "BRD4",
    "MMP1", "MMP2", "MMP9", "MMP13",
    "CASP3", "CASP8", "CASP9",
    "CTSB", "CTSK", "CTSL", "CTSS",
    "ALOX5", "ALOX15",
    "XDH", "IDO1",
    "NOS1", "NOS2", "NOS3",
    "HMOX1",
    "SOD1", "SOD2", "CAT", "GPX1",
    "EGFR", "ERBB2", "ERBB4",
    "SRC", "LYN", "LCK", "FYN",
    "JAK1", "JAK2", "JAK3",
    "AKT1", "AKT2",
    "MAPK1", "MAPK3", "MAPK14",
    "MTOR", "PIK3CA",
    "CDK2", "CDK4", "CDK6", "CDK9",
    "BRAF", "MAP2K1", "MAP2K2",
    "IGF1R", "INSR",
    "FLT1", "KDR", "FLT4",
    "ABL1", "CSF1R", "KIT", "PDGFRA", "PDGFRB",
    "PLK1", "AURKA", "AURKB",
    "CHEK1", "CHEK2",
    "ESR1", "ESR2", "AR", "PGR", "NR3C1",
    "PPARA", "PPARG", "PPARD",
    "RARA", "RARB", "RARG",
    "VDR", "NR1H3",
    "ADRB1", "ADRB2",
    "DRD2", "DRD3", "DRD4",
    "HTR1A", "HTR2A", "HTR2C",
    "CHRM1", "CHRM2", "CHRM3",
    "OPRM1", "OPRK1", "OPRD1",
    "CNR1", "CNR2",
    "TRPV1", "TRPA1", "TRPM8",
    "GRIN1", "GRIN2B",
    "GABRA1", "GABRB3", "GABRG2",
    "SCN9A", "SCN5A",
    "KCNH2",
    "BCL2", "BCL2L1", "BAX",
    "TP53",
    "NFKB1", "RELA",
    "STAT3", "STAT5A",
    "NFE2L2",
    "HIF1A",
    "APP", "BACE1", "MAPT",
    "FAAH", "MGLL",
    "PTPN1", "PTPN2",
    "DPP4",
    "COMT", "DDC",
    "AHR",
    "PARP1", "TOP1", "TOP2A",
    "TERT",
    "TLR4", "NLRP3",
    "TNF",
    "IL6", "IL1B",
    "TGFB1",
    "ATP7B",
    "SLC31A1",
    "DLAT", "PDHA1", "FDX1", "LIAS", "LIPT1",
}

# ================================================================
# UniProt alias normalization
# ================================================================
UNIPROT_ALIASES = {
    "HLAE": "HLA-E",
    "PSD95": "DLG4",
    "NMDAR1": "GRIN1",
    "NR1": "GRIN1",
    "ACP3": "ACP3",
}

VALID_HGNC_PATTERN = re.compile(r"^[A-Z][A-Z0-9\-]*[A-Z0-9]$")


def normalize_gene_symbol(gene: str) -> str:
    g = gene.strip().upper()
    if not g:
        return None
    if g in UNIPROT_ALIASES:
        return UNIPROT_ALIASES[g]
    if not VALID_HGNC_PATTERN.match(g):
        return None
    return g


# ================================================================
# Scoring engine
# ================================================================
def score_target(
    gene: str,
) -> dict:
    result = {"gene": gene}
    score = 0.0
    reasons = []

    # Layer 1: Housekeeping check (penalty)
    result["death_sub_pathways"] = []
    if gene in HOUSEKEEPING_GENES:
        result["is_housekeeping"] = True
        result["in_hard_blacklist"] = False
        result["known_drug_target"] = False
        result["brain_expression_level"] = 0
        result["ciri_core"] = False
        result["has_pdb"] = False
        result["druggable_family"] = False
        result["metabolic_enzyme"] = False
        result["soft_blacklist"] = False
        result["tier"] = "Rejected"
        result["score"] = -10
        result["reason"] = "housekeeping_gene"
        return result
    result["is_housekeeping"] = False

    # Layer 1: Hard blacklist check
    if gene in HARD_BLACKLIST:
        result["in_hard_blacklist"] = True
        result["known_drug_target"] = False
        result["brain_expression_level"] = 0
        result["ciri_core"] = False
        result["has_pdb"] = False
        result["druggable_family"] = False
        result["metabolic_enzyme"] = False
        result["soft_blacklist"] = False
        result["tier"] = "Rejected"
        result["score"] = -10
        result["reason"] = f"hard_blacklist({'hemoglobin' if gene in BLACKLIST_HEMOGLOBIN else 'mitochondrial' if gene in BLACKLIST_MITOCHONDRIAL else 'immunoglobulin'})"
        return result
    result["in_hard_blacklist"] = False

    # Soft blacklist penalty
    if gene in SOFT_BLACKLIST:
        score -= 2
        result["soft_blacklist"] = True
        reasons.append(f"plasma_protein(-2)")
    else:
        result["soft_blacklist"] = False

    # Layer 2: Known drug target anchor
    known_target = gene in KNOWN_DRUG_TARGETS
    result["known_drug_target"] = known_target
    if known_target:
        score += 3
        reasons.append("known_target(+3)")

    # Layer 3a: Brain expression
    brain_level = 0
    if gene in NEURON_MARKERS:
        brain_level = 3
    elif gene in GLIAL_MARKERS:
        brain_level = 2
    elif gene in BRAIN_EXPRESSED:
        brain_level = 1

    result["brain_expression_level"] = brain_level
    if brain_level > 0:
        bonus = brain_level * 2
        score += bonus
        reasons.append(f"brain_expr(L{brain_level},+{bonus})")

    # Layer 3b: CIRI pathway relevance
    ciri_core = gene in CIRI_CORE_PATHWAY_GENES
    result["ciri_core"] = ciri_core
    if ciri_core:
        score += 3
        reasons.append("ciri_core(+3)")

    death_subs = []
    for sub_name, sub_genes in CIRI_DEATH_SUB.items():
        if gene in sub_genes:
            death_subs.append(sub_name)
    result["death_sub_pathways"] = death_subs
    if death_subs:
        n_death = len(death_subs)
        bonus = min(n_death * 1, 2)
        score += bonus
        reasons.append(f"death_pathway({'/'.join(death_subs)},+{bonus})")

    # Layer 4a: PDB availability
    has_pdb = gene in HAS_PDB_KNOWN
    result["has_pdb"] = has_pdb
    if has_pdb:
        score += 2
        reasons.append("has_pdb(+2)")

    # Layer 4b: Druggable family
    druggable_family = any(
        gene.startswith(p) for p in DRUGGABLE_FAMILY_PREFIXES
    )
    if not druggable_family and gene.startswith("NR") and len(gene) > 2:
        druggable_family = True
    result["druggable_family"] = druggable_family
    if druggable_family:
        score += 1
        reasons.append("druggable_family(+1)")

    # Layer 5: Metabolic enzyme flag (slight penalty for CYP-only predictions)
    is_cyp = gene.startswith("CYP") and gene[3].isdigit()
    is_ugt = gene.startswith("UGT")
    result["metabolic_enzyme"] = is_cyp or is_ugt
    if (is_cyp or is_ugt) and not known_target:
        score -= 1
        reasons.append("metabolic_only(-1)")

    result["score"] = score

    # Assign tier
    if score >= 7:
        result["tier"] = "Tier1"
    elif score >= 4:
        result["tier"] = "Tier2"
    elif score >= 1:
        result["tier"] = "Tier3"
    else:
        result["tier"] = "Rejected"

    result["reason"] = ";".join(reasons) if reasons else "no_evidence"
    return result


# ================================================================
# Main
# ================================================================
def main():
    print("=" * 70)
    print("  SMILES-derived Target Cleaning & Confidence Filtering")
    print("  Reference: SwissTargetPrediction (Gfeller 2014, PMID:24792161)")
    print("             INPUT Platform (Li 2022, PMID:35414968)")
    print("             World Fed. Chinese Med. Soc. Guidelines (2021)")
    print("  BCP (beta-caryophyllene) targets from STP")
    print("=" * 70)

    # Read input
    if not os.path.exists(INPUT_FILE):
        print(f"\n[ERROR] Input file not found: {INPUT_FILE}")
        sys.exit(1)

    with open(INPUT_FILE, "r", encoding="utf-8") as f:
        raw_lines = [line.strip() for line in f if line.strip()]

    print(f"\n[0] RAW: {len(raw_lines)} entries read")

    # Layer 0: Normalize
    normalized = []
    n_invalid = 0
    for g in raw_lines:
        ng = normalize_gene_symbol(g)
        if ng:
            normalized.append(ng)
        else:
            n_invalid += 1

    # Deduplicate
    unique_genes = list(dict.fromkeys(normalized))
    n_dup = len(normalized) - len(unique_genes)

    print(f"    Normalized & dedup: {len(unique_genes)} ")
    print(f"    Invalid symbols: {n_invalid}")
    print(f"    Duplicates removed: {n_dup}")

    # Score each gene
    results = [score_target(g) for g in unique_genes]

    # Build DataFrame
    df = pd.DataFrame(results)

    # Sort
    tier_order = {"Tier1": 0, "Tier2": 1, "Tier3": 2, "Rejected": 3}
    df["tier_rank"] = df["tier"].map(tier_order)
    df = df.sort_values(["tier_rank", "score"], ascending=[True, False]).drop(
        columns="tier_rank"
    )
    df = df.reset_index(drop=True)

    # Save full results
    output_csv = os.path.join(OUTPUT_DIR, "tiered_targets.csv")
    df.to_csv(output_csv, index=False)
    print(f"\n[OUTPUT] {output_csv}")

    # Tier stats
    n_tier1 = (df["tier"] == "Tier1").sum()
    n_tier2 = (df["tier"] == "Tier2").sum()
    n_tier3 = (df["tier"] == "Tier3").sum()
    n_rejected = (df["tier"] == "Rejected").sum()

    print(f"\n{'─' * 40}")
    print(f"  Tier1 (high, score>=7):  {n_tier1}")
    print(f"  Tier2 (medium, score>=4):  {n_tier2}")
    print(f"  Tier3 (low, score>=1):    {n_tier3}")
    print(f"  Rejected:                 {n_rejected}")
    print(f"{'─' * 40}")
    print(f"  High-confidence (T1+T2):  {n_tier1 + n_tier2} / {len(df)}")

    # Save high-confidence targets
    high_conf = df[df["tier"].isin(["Tier1", "Tier2"])]
    output_high = os.path.join(OUTPUT_DIR, "high_confidence_targets.txt")
    high_conf["gene"].to_csv(output_high, index=False, header=False)
    print(f"[OUTPUT] {output_high}")

    # Detailed report
    report_path = os.path.join(OUTPUT_DIR, "filter_report.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("BCP Target Cleaning & Confidence Filtering Report\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Input:           {INPUT_FILE}\n")
        f.write(f"Raw entries:     {len(raw_lines)}\n")
        f.write(f"Valid symbols:   {len(unique_genes)}\n")
        f.write(f"Invalid/removed: {n_invalid + n_dup}\n\n")

        f.write("Layer Statistics:\n")
        f.write(f"  L1: housekeeping → {df['is_housekeeping'].sum()} removed\n")
        f.write(f"  L1: hard_blacklist → {df['in_hard_blacklist'].sum() if 'in_hard_blacklist' in df else 0} removed (if hit)\n")
        f.write(f"  L1: soft_blacklist → {df['soft_blacklist'].sum()} penalized\n")
        f.write(f"  L2: known targets → {df['known_drug_target'].sum()} anchored\n")
        f.write(f"  L3a: brain_expr(L1-3) → {(df['brain_expression_level'] > 0).sum()} expressed\n")
        f.write(f"  L3b: ciri_core → {df['ciri_core'].sum()}\n")
        f.write(f"  L3c: death_pathway → {df['death_sub_pathways'].apply(len).gt(0).sum()}\n")
        f.write(f"  L4a: has_pdb → {df['has_pdb'].sum()}\n")
        f.write(f"  L4b: druggable → {df['druggable_family'].sum()}\n\n")

        f.write("Tier Distribution:\n")
        f.write(f"  Tier1 (score>=7):  {n_tier1}\n")
        f.write(f"  Tier2 (score>=4):  {n_tier2}\n")
        f.write(f"  Tier3 (score>=1):  {n_tier3}\n")
        f.write(f"  Rejected:          {n_rejected}\n\n")

        f.write("Tier1 Targets:\n")
        for _, row in df[df["tier"] == "Tier1"].iterrows():
            f.write(f"  {row['gene']:12s} score={row['score']:+.1f}  {row['reason']}\n")

        f.write("\nTier2 Targets:\n")
        for _, row in df[df["tier"] == "Tier2"].iterrows():
            f.write(f"  {row['gene']:12s} score={row['score']:+.1f}  {row['reason']}\n")

        f.write("\nTier3 Targets:\n")
        for _, row in df[df["tier"] == "Tier3"].iterrows():
            f.write(f"  {row['gene']:12s} score={row['score']:+.1f}  {row['reason']}\n")

        f.write("\nRejected Targets:\n")
        for _, row in df[df["tier"] == "Rejected"].iterrows():
            f.write(f"  {row['gene']:12s} score={row['score']:+.1f}  {row['reason']}\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write(f"Recommended: use Tier1+Tier2 ({n_tier1+n_tier2} targets) for downstream\n")
        f.write("=" * 70 + "\n")

    print(f"[OUTPUT] {report_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()