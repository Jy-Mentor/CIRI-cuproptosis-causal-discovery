# -*- coding: utf-8 -*-
"""
BCP靶点分级清洗: 6层交叉验证 -> Tier1/2/3/Rejected (v5)
=================================================
v5修复 (对照v4审查报告):
  P1: 修复PROBE_BLACKLIST大小写不一致 (1368518_at → 1368518_AT)
  P1: pathway_synergy改为铜死亡+死亡相关子通路协同
  P2: bc_p90预计算移到循环外，避免循环内重复计算
"""

import os
import sys
import json
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import BASE_DIR, RESULTS_DIR, CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED

OUTPUT_DIR = os.path.join(RESULTS_DIR, "target_filter")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# 内置数据库
# ============================================================

KNOWN_TARGET_SET = {
    "ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C",
    "ADRB1", "ADRB2", "ADRB3", "DRD1", "DRD2", "DRD3", "DRD4", "DRD5",
    "HTR1A", "HTR1B", "HTR1D", "HTR2A", "HTR2B", "HTR2C", "HTR3A",
    "HTR4", "HTR5A", "HTR6", "HTR7",
    "HRH1", "HRH2", "HRH3", "HRH4",
    "CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5",
    "OPRM1", "OPRK1", "OPRD1", "OPRL1",
    "CNR1", "CNR2", "PTAFR", "OXTR",
    "GABRA1", "GABRA2", "GABRA3", "GABRA5",
    "GABRB1", "GABRB2", "GABRB3",
    "GABRG1", "GABRG2", "GABRG3",
    "GRIN1", "GRIN2A", "GRIN2B",
    "P2RX1", "P2RX2", "P2RX3", "P2RX4", "P2RX7",
    "ESR1", "ESR2", "AR", "PGR", "NR3C1", "NR3C2",
    "PPARA", "PPARG", "PPARD",
    "RARA", "RARB", "RARG", "RXRA", "RXRB", "RXRG",
    "THRA", "THRB", "NR1H2", "NR1H3",
    "NR1I2", "NR1I3", "VDR",
    "AKT1", "AKT2", "AKT3",
    "MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK10", "MAPK14",
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
    "PRKACA", "PRKACB",
    "PRKAA1", "PRKAA2",
    "PRKCQ", "PRKCI", "PRKCB",
    "PLK1", "AURKA", "AURKB",
    "CHEK1", "CHEK2", "WEE1",
    "EPHA1", "EPHA2",
    "PTK2", "PTK2B",
    "ACHE", "BCHE", "MAOA", "MAOB",
    "PDE4A", "PDE4B", "PDE5A",
    "DPP4", "FASN", "HMGCR", "CYP1A2", "CYP2C9",
    "CYP2C19", "CYP2D6", "CYP3A4", "CYP17A1", "CYP51A1",
    "CYP19A1", "HSD11B1", "HSD17B1", "HSD17B3",
    "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6",
    "SIRT1", "SIRT2", "SIRT5", "SIRT6",
    "PARP1", "BRD2", "BRD3", "BRD4",
    "MMP1", "MMP2", "MMP9", "MMP13",
    "CASP1", "CASP3", "CASP8", "CASP9",
    "CTSB", "CTSK", "CTSL", "CTSS", "CASP6",
    "PTGS1", "PTGS2", "ALOX5", "ALOX12", "ALOX15",
    "XDH", "IDO1", "TDO2",
    "HSPA5", "HSP90AA1",
    "NOS1", "NOS2", "NOS3",
    "HMOX1", "SOD1", "SOD2", "CAT", "GPX1", "GSR",
    "TSPO", "COMT", "GAD1", "GAD2", "DDC",
    "FAAH", "MGLL", "CES1", "CES2", "PLA2G4A",
    "PTPN1", "PTPN2", "PTPN6", "PTPN11",
    "BCL2", "BAX", "BCL2L1", "MCL1",
    "SCN1A", "SCN2A", "SCN5A", "SCN9A",
    "KCNJ1", "KCNJ2", "KCNJ8", "KCNJ11",
    "KCNA1", "KCNA2", "KCNA5",
    "KCNH1", "KCNH2",
    "KCNQ1", "KCNQ2", "KCNQ3", "KCNQ4", "KCNQ5",
    "CACNA1A", "CACNA1C", "CACNA1D", "CACNA1G",
    "TRPA1", "TRPC1", "TRPC4", "TRPC5", "TRPC6",
    "TRPM8", "TRPV1", "TRPV2", "TRPV3", "TRPV4",
    "SLC6A1", "SLC6A2", "SLC6A3", "SLC6A4",
    "SLC2A1", "SLC2A3", "SLC2A4",
    "ABCB1", "ABCC1", "ABCG2",
    "ATP1A1", "ATP1A2", "ATP1A3", "ATP2A1", "ATP2A2",
    "TP53", "NFKB1", "RELA", "RELB", "REL",
    "STAT1", "STAT2", "STAT3", "STAT5A", "STAT6",
    "NFE2L2", "FOXO1", "FOXO3", "HIF1A",
    "NOTCH1", "NOTCH2", "NOTCH3",
    "SREBF1", "SREBF2", "EGR1", "JUN", "FOS",
    "ATF4", "DDIT3", "XBP1", "IRF1", "IRF3",
    "CEBPB", "CEBPD", "NR4A1", "NR4A2",
    "TNF", "IL1A", "IL1B", "IL2", "IL4", "IL5", "IL6", "IL8", "IL10",
    "CCL2", "CCL3", "CCL4", "CCL5",
    "CXCL8", "CXCL10", "CXCL12",
    "ICAM1", "VCAM1", "SELE", "SELP",
    "TLR2", "TLR3", "TLR4", "TLR7", "TLR9",
    "NLRP3", "ATG5", "ATG7", "BECN1",
    "APP", "MAPT", "PSEN1", "PSEN2",
    "BDNF", "NGF", "GDNF", "NT3",
    "SIGMAR1",
    "GFAP", "AIF1", "AHR", "ALDH1A1", "ARG1",
    "CCND1", "CXCR3", "KEAP1",
    "PDGFRA", "TGFB1", "TIMP1", "VEGFA",
    "RORC", "RORA",
}

BRAIN_EXPRESSED_SET = {
    "兴奋性神经元": {
        "GRIN1", "GRIN2A", "GRIN2B", "GRIA1", "GRIA2", "GRIA3",
        "GABRA1", "GABRB3", "GABRG2",
        "SYN1", "SYP", "SV2A", "SNAP25", "SYT1",
        "MAP2", "NEUROD6", "NEUROD1", "NEUROD2",
        "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
        "ARC", "BDNF", "FOS", "JUN", "EGR1", "NR4A1",
        "HOMER1", "HOMER2", "HOMER3",
        "DLG4", "SHANK1", "SHANK2", "SHANK3",
    },
    "胶质细胞": {
        "GFAP", "ALDH1L1", "SLC1A2", "SLC1A3", "GJB6", "GJC1",
        "AQP4", "KCNJ10", "ABCC8", "KCNJ8",
        "MBP", "PLP1", "MOG", "CNP", "MAG",
        "OLIG1", "OLIG2", "SOX10",
        "P2RY12", "TMEM119", "CX3CR1", "CD33", "TYROBP",
        "ITGAM", "ITGAX", "FCGR1A", "FCGR3A",
        "APOE", "CLU", "C3", "C1QA", "C1QB", "C1QC",
        "CD68", "CD74", "CTSB", "CTSD", "CTSL", "CTSS",
        "TREM2", "FCER1G",
    },
    "血管内皮": {
        "PECAM1", "VWF", "CLDN5", "CDH5", "TJP1", "TJP2",
        "FLT1", "KDR", "PDGFRB",
        "SLC2A1", "ABCB1", "ABCG2",
        "CAV1", "CAV2",
    },
    "脑区富集": {
        "CALM1", "CALM2", "CALM3", "SYN2", "SNAP23",
        "GAD1", "GAD2", "SLC32A1", "SLC6A11",
        "PVALB", "SST", "NPY", "CCK", "VIP",
        "TAC1", "TAC3", "PENK", "PDYN", "POMC",
        "NPY1R", "NPY2R",
        "DRD1", "DRD2", "DRD3", "HTR1A", "HTR2A",
        "ADORA2A", "CNR1", "CHRNA4", "CHRNB2",
        "GRM1", "GRM2", "GRM3", "GRM4", "GRM5",
        "GRM7", "GRM8",
        "CHRM1", "CHRM4", "GABBR1", "GABBR2",
        "GPR37", "GPR37L1", "GPR88", "GPR139",
    },
    "突触相关": {
        "DLG4", "DLG2", "DLG3",
        "SHANK1", "SHANK2", "SHANK3",
        "HOMER1", "HOMER2", "HOMER3",
        "CAMK2A", "CAMK2B", "SYNGAP1",
        "NRXN1", "NRXN2", "NRXN3",
        "NLGN1", "NLGN2", "NLGN3", "NLGN4X",
        "CNTNAP2", "LGI1",
    },
}

NEURON_MARKERS = BRAIN_EXPRESSED_SET.get("兴奋性神经元", set())
GLIAL_MARKERS = BRAIN_EXPRESSED_SET.get("胶质细胞", set())
VASCULAR_MARKERS = BRAIN_EXPRESSED_SET.get("血管内皮", set())
BRAIN_EXPRESSED_FLAT = {g for v in BRAIN_EXPRESSED_SET.values() for g in v}

BLACKLIST_PLASMA = {"ALB", "TF", "TFRC", "HP", "HPX", "ORM1", "ORM2", "A2M", "SERPINA1",
                    "SERPINC1", "C1QC", "C1QA", "C1QB"}
BLACKLIST_STRUCTURAL = {"ACTA1", "ACTA2", "ACTG1", "ACTG2",
                        "MYH1", "MYH2", "MYH3", "MYH4", "MYH6", "MYH7", "MYH9", "MYH10", "MYH11",
                        "MYBPC1", "MYBPC2", "MYBPC3",
                        "TNNC1", "TNNC2", "TNNT1", "TNNT2", "TNNT3",
                        "TPM1", "TPM2", "TPM3", "TPM4",
                        "FLNA", "FLNB", "FLNC", "VIM", "DES", "NEB", "TTN",
                        "SNTA1", "SNTB1", "SNTB2"}
BLACKLIST_HEMOGLOBIN = {"HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBZ", "HBE1", "HBQ1"}
BLACKLIST_MITODNA = {"MTCYB", "MTND1", "MTND2", "MTND3", "MTND4", "MTND4L", "MTND5",
                     "MTND6", "MTCO1", "MTCO2", "MTCO3", "MTATP6", "MTATP8",
                     "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5",
                     "MT-ND6", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8",
                     "MT-CYB"}
BLACKLIST_COLLAGEN = {"COL1A1", "COL1A2", "COL2A1", "COL3A1", "COL4A1", "COL4A2",
                      "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
                      "COL7A1", "COL8A1", "COL9A1", "COL10A1", "COL11A1",
                      "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1",
                      "COL17A1", "COL18A1", "COL19A1", "COL20A1",
                      "FN1", "LAMA1", "LAMA2", "LAMA3", "LAMB1", "LAMC1",
                      "LTBP1", "FBLN1", "FBLN2", "FBLN5",
                      "CTHRC1", "C1QTNF1",
                      "BGN", "DCN", "FMOD", "LUM",
                      "ELN", "EFEMP1", "EFEMP2"}
BLACKLIST_IG = {"IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM", "IGKC", "IGLC1",
                "IGLC2", "IGLC3",
                "B2M", "CD1B", "CD1C", "CD1D"}

HARD_BLACKLIST = BLACKLIST_HEMOGLOBIN | BLACKLIST_MITODNA | BLACKLIST_IG

METABOLIC_ENZYMES = {
    "CYP1A1", "CYP1A2", "CYP1B1", "CYP2A6", "CYP2A13", "CYP2B6",
    "CYP2C8", "CYP2C9", "CYP2C18", "CYP2C19", "CYP2D6", "CYP2E1",
    "CYP2J2", "CYP2R1", "CYP2S1", "CYP2U1",
    "CYP3A4", "CYP3A5", "CYP3A43", "CYP3A7",
    "CYP4A11", "CYP4B1", "CYP4F2", "CYP4F3",
    "CYP17A1", "CYP19A1", "CYP27A1", "CYP27B1",
    "CYP51A1", "CYP8B1", "CYP26A1", "CYP26B1",
    "CYP24A1",
    "ALDH1A1", "ALDH1A2", "ALDH1A3", "ALDH1B1",
    "ALDH2", "ALDH3A1", "ALDH3A2", "ALDH3B1", "ALDH3B2",
    "ALDH4A1", "ALDH5A1", "ALDH6A1", "ALDH7A1", "ALDH8A1",
    "ALDH9A1", "ALDH16A1", "ALDH18A1", "ALDH1L1", "ALDH1L2",
}

DRUGGABLE_PREFIXES = [
    "SCN", "KCN", "TRP", "CACNA", "CACNB", "CACNG",
    "SLC", "ABC",
]

DRUGGABLE_EXACT_PREFIXES = [
    "NR",
]

INTERMEDIATE_FILAMENTS = {"GFAP", "VIM", "KRT1", "KRT2", "KRT3", "KRT4", "KRT5", "KRT6A",
                          "KRT6B", "KRT6C", "KRT7", "KRT8", "KRT9", "KRT10", "KRT12",
                          "KRT13", "KRT14", "KRT15", "KRT16", "KRT17", "KRT18", "KRT19",
                          "KRT20", "NEFL", "NEFM", "NEFH", "SYNM",
                          "DES", "BFSP1", "BFSP2"}

HAS_PDB_SET = {
    "ACHE", "BCHE", "MAOA", "MAOB", "PDE4A", "PDE4B", "PDE5A",
    "DPP4", "FASN", "HMGCR", "CYP1A2", "CYP2C9", "CYP2C19",
    "CYP2D6", "CYP3A4", "CYP17A1", "CYP51A1", "CYP19A1",
    "HSD11B1", "HSD17B1", "HDAC1", "HDAC2", "HDAC3", "HDAC4",
    "HDAC5", "HDAC6", "HDAC7", "HDAC8", "HDAC9", "HDAC10",
    "SIRT1", "SIRT2", "SIRT5", "SIRT6",
    "PARP1", "BRD2", "BRD3", "BRD4",
    "MMP1", "MMP2", "MMP9", "MMP13",
    "CASP1", "CASP3", "CASP8", "CASP9",
    "CTSB", "CTSK", "CTSL", "CTSS",
    "PTGS1", "PTGS2", "ALOX5", "ALOX12", "ALOX15",
    "XDH", "IDO1", "TDO2",
    "HSPA5", "HSP90AA1",
    "NOS1", "NOS2", "NOS3",
    "HMOX1", "SOD1", "SOD2", "CAT", "GPX1",
    "TSPO", "COMT", "GAD1", "GAD2", "DDC",
    "FAAH", "MGLL", "CES1", "CES2",
    "PLA2G4A", "PTPN1", "PTPN2", "PTPN6", "PTPN11",
    "BCL2", "BAX", "BCL2L1", "MCL1",
    "EGFR", "ERBB2", "ERBB3", "ERBB4",
    "SRC", "FYN", "LYN", "LCK", "FGR",
    "JAK1", "JAK2", "JAK3", "TYK2",
    "PIK3CA", "MTOR", "RPS6KB1",
    "CDK1", "CDK2", "CDK4", "CDK6", "CDK9",
    "BRAF", "RAF1", "MAP2K1", "MAP2K2",
    "INSR", "IGF1R",
    "FLT1", "FLT3", "FLT4", "KDR",
    "CSF1R", "KIT", "PDGFRA", "PDGFRB",
    "ABL1", "ABL2", "ITK", "BTK",
    "IRAK1", "IRAK4", "TBK1", "IKBKB",
    "AKT1", "AKT2", "AKT3",
    "MAPK1", "MAPK3", "MAPK14",
    "PRKACA", "PRKACB",
    "PRKAA1", "PRKAA2",
    "PLK1", "AURKA", "AURKB",
    "CHEK1", "CHEK2", "WEE1",
    "ESR1", "ESR2", "AR", "PGR", "NR3C1", "NR3C2",
    "PPARA", "PPARG", "PPARD",
    "RARA", "RARB", "RARG", "RXRA", "RXRB",
    "THRA", "THRB", "NR1H2", "NR1H3", "VDR",
    "ADRA1A", "ADRA1B", "ADRA2A", "ADRA2B", "ADRA2C",
    "ADRB1", "ADRB2", "ADRB3",
    "DRD1", "DRD2", "DRD3", "DRD4",
    "HTR1A", "HTR2A", "HTR2C", "HTR3A", "HTR6", "HTR7",
    "HRH1", "HRH2", "HRH3",
    "CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5",
    "OPRM1", "OPRK1", "OPRD1",
    "CNR1", "CNR2",
    "P2RX7", "GABRA1", "GABRA2", "GABRB3", "GABRG2",
    "GRIN1", "GRIN2A", "GRIN2B",
    "SCN1A", "SCN2A", "SCN5A", "SCN9A",
    "KCNJ1", "KCNJ2", "KCNJ8", "KCNJ11",
    "KCNA1", "KCNA2", "KCNA5",
    "KCNH1", "KCNH2",
    "KCNQ1", "KCNQ2", "KCNQ3", "KCNQ4", "KCNQ5",
    "CACNA1A", "CACNA1C", "CACNA1D", "CACNA1G",
    "CACNA1H", "TRPA1", "TRPC1", "TRPC4", "TRPC5",
    "TRPC6", "TRPM8", "TRPV1", "TRPV2", "TRPV3", "TRPV4",
    "SLC6A1", "SLC6A2", "SLC6A3", "SLC6A4",
    "SLC2A1", "SLC2A3", "SLC2A4",
    "ABCB1", "ABCC1", "ABCG2",
    "ATP1A1", "ATP1A2", "ATP1A3", "ATP2A1", "ATP2A2",
    "TP53", "NFATC1", "NFKB1", "RELA",
    "STAT1", "STAT3", "STAT5A", "STAT6",
    "NFE2L2", "FOXO1", "FOXO3", "HIF1A",
    "NOTCH1", "NOTCH2", "NOTCH3",
    "SREBF1", "SREBF2", "EGR1", "JUN", "FOS",
    "ATF4", "DDIT3", "XBP1", "IRF1", "IRF3",
    "CEBPB", "CEBPD", "GATA3", "NR4A1", "NR4A2",
    "TNF", "IL1A", "IL1B", "IL2", "IL4", "IL5", "IL6", "IL8", "IL10",
    "CCL2", "CCL3", "CCL4", "CCL5",
    "CXCL8", "CXCL10", "CXCL12",
    "ICAM1", "VCAM1", "SELE", "SELP",
    "TLR2", "TLR3", "TLR4", "TLR7", "TLR9",
    "NLRP3", "ATG5", "ATG7", "BECN1",
    "APP", "MAPT", "PSEN1", "PSEN2",
    "BDNF", "NGF", "GDNF", "NT3",
    "SIGMAR1",
    "GFAP",
    "AIF1", "AHR", "ALDH1A1", "ARG1",
    "CCND1", "CCL2", "CDK4", "CXCR3",
    "HMOX2", "IFNG", "KEAP1",
    "MMP2", "MTOR", "NFE2L2", "NFKB1", "NOTCH1",
    "NOS2", "NOS3", "PDGFRA", "PPARG", "PTGS2",
    "RELA", "SIRT1", "SOD1", "SOD2", "SREBF1", "STAT3",
    "TGFB1", "TIMP1", "TNF", "TP53",
}

OXIDATIVE_STRESS_GENES = {
    "SOD1", "SOD2", "SOD3",
    "CAT", "GPX1", "GPX2", "GPX3", "GPX4",
    "GSR", "TXN", "TXNRD1", "TXNRD2",
    "NQO1", "NQO2", "HMOX1", "HMOX2",
    "GCLC", "GCLM", "GSS",
    "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6",
    "KEAP1", "NFE2L2",
    "SELENOP", "SELENOS",
}

APOPTOSIS_GENES = {
    "CASP3", "CASP6", "CASP7", "CASP8", "CASP9", "CASP10",
    "BAX", "BAK1", "BAD",
    "BCL2", "BCL2L1", "BCL2L11", "MCL1",
    "FAS", "FASLG", "TNFRSF1A", "TNFRSF1B",
    "FADD", "TRADD",
    "TP53", "PMAIP1", "NOXA", "PUMA",
    "BID", "BMF", "BIK",
    "DIABLO", "HTRA2", "CYCS", "APAF1",
    "XIAP", "BIRC2", "BIRC3",
}

CUPROPTOSIS_SPECIFIC_GENES = {
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHA1", "PDHB",
    "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7A", "ATP7B",
    "DLD", "DBT", "DLST", "PDHA2", "GCSH",
}

CIRI_PATHWAY_SUB = {
    "nfkb": {"NFKB1", "NFKB2", "RELA", "RELB", "REL", "IKBKB", "IKBKE",
             "TLR2", "TLR3", "TLR4", "TLR7", "TLR9", "TNF", "TNFRSF1A",
             "TRAF1", "TRAF2", "TRAF5", "TRAF6", "MYD88", "IRAK1", "IRAK4"},
    "nrf2": {"NFE2L2", "KEAP1", "HMOX1", "GCLC", "GCLM", "NQO1",
             "TXN", "TXNRD1", "TXNRD2"},
    "jak_stat": {"JAK1", "JAK2", "JAK3", "TYK2",
                 "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
                 "SOCS1", "SOCS2", "SOCS3"},
    "ferroptosis": {"GPX4", "ACSL4", "LPCAT3", "ALOX15", "ALOX12",
                    "TFRC", "SLC3A2", "SLC7A11", "FSP1", "AIFM2",
                    "DHODH", "GCH1"},
    "necroptosis": {"RIPK1", "RIPK3", "MLKL", "ZBP1",
                    "FADD", "CASP8", "CYLD", "TRADD"},
    "pyroptosis": {"NLRP1", "NLRP2", "NLRP3", "NLRC4", "AIM2",
                   "CASP1", "CASP4", "CASP5",
                   "GSDMD", "GSDME", "IL1B", "IL18", "NEK7", "P2RX7"},
    "apoptosis_ciri": {"CASP3", "CASP6", "CASP7", "CASP8", "CASP9",
                       "BAX", "BAK1", "BAD",
                       "BCL2", "BCL2L1", "MCL1",
                       "FAS", "FASLG", "TP53", "PMAIP1",
                       "DIABLO", "CYCS", "APAF1", "XIAP"},
    "autophagy_ciri": {"MTOR", "ULK1",
                       "ATG5", "ATG7", "ATG12", "ATG16L1",
                       "BECN1", "PIK3C3", "PIK3CA",
                       "SQSTM1", "OPTN", "TFEB", "FOXO1", "FOXO3"},
    "upr": {"ERN1", "EIF2AK3", "ATF6",
            "EIF2S1", "ATF4", "DDIT3", "XBP1", "HSPA5"},
    "cuproptosis": {"FDX1", "LIAS", "LIPT1", "DLAT", "PDHA1", "PDHB",
                    "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7A", "ATP7B",
                    "DLD", "DBT", "DLST", "PDHA2", "GCSH"},
}

CIRI_DEATH_RELATED_SUBS = {"ferroptosis", "necroptosis", "pyroptosis", "apoptosis_ciri", "cuproptosis"}

PROBE_BLACKLIST = {
    "1387883_A_AT", "1375219_A_AT", "1393852_AT",
    "1367651_AT", "1368430_AT", "1368518_AT",
}

UNIPROT_ALIASES = {
    "HLAE": "HLA-E",
    "PSD95": "DLG4",
    "NMDAR1": "GRIN1",
    "PSD-95": "DLG4",
    "NR1": "GRIN1",
}


class GeneFilterResult:
    __slots__ = [
        "gene", "tier", "score",
        "known_target", "brain_expressed", "brain_expressed_weight",
        "has_pdb", "network_hub", "druggable_family", "in_module",
        "cuproptosis", "ciri_pathway", "ciri_sub_pathways",
        "hard_blacklisted",
        "oxidative_stress", "apoptosis", "pathway_synergy",
        "metabolic_enzyme", "blacklisted", "blacklist_penalty",
        "intermediate_filament", "reason",
    ]

    def __init__(self, gene: str = ""):
        self.gene = gene
        self.tier = ""
        self.score = 0.0
        self.known_target = False
        self.brain_expressed = False
        self.brain_expressed_weight = 0.0
        self.has_pdb = False
        self.network_hub = False
        self.druggable_family = False
        self.in_module = False
        self.cuproptosis = False
        self.ciri_pathway = False
        self.ciri_sub_pathways: List[str] = []
        self.hard_blacklisted = False
        self.oxidative_stress = False
        self.apoptosis = False
        self.pathway_synergy = False
        self.metabolic_enzyme = False
        self.blacklisted = False
        self.blacklist_penalty = 0.0
        self.intermediate_filament = False
        self.reason = ""


def normalize_gene(gene: str):
    """Layer 0: UniProt标准化 + 探针ID过滤"""
    gene = gene.strip().upper()
    if gene in PROBE_BLACKLIST:
        return None
    alias = UNIPROT_ALIASES.get(gene)
    if alias is not None:
        return alias
    if gene == "":
        return None
    return gene


def check_known_target(gene: str) -> bool:
    """Layer 1: 已知靶点库锚定"""
    return gene in KNOWN_TARGET_SET


def check_brain_specificity(gene: str) -> Tuple[bool, float]:
    """Layer 2: 脑组织特异性分层加分"""
    if gene in CUPROPTOSIS_GENES:
        return True, 2.0
    if gene in NEURON_MARKERS:
        return True, 2.0
    if gene in GLIAL_MARKERS:
        return True, 1.5
    if gene in VASCULAR_MARKERS:
        return True, 1.5
    if gene in BRAIN_EXPRESSED_FLAT:
        return True, 1.0
    return False, 0.0


def check_blacklist(gene: str) -> Tuple[bool, float, bool]:
    """Layer 3: GO功能负向过滤 (梯度罚分 + Hard Blacklist)
    返回: (blacklisted, penalty, hard_blacklisted)
    """
    if gene in HARD_BLACKLIST:
        return True, 3.0, True
    if gene in BLACKLIST_PLASMA:
        return True, 2.0, False
    if gene in BLACKLIST_STRUCTURAL:
        return True, 2.0, False
    if gene in BLACKLIST_COLLAGEN:
        return True, 1.5, False
    return False, 0.0, False


def check_druggable_family(gene: str) -> bool:
    """Layer 4: 可成药家族检测 (仅startswith)"""
    for prefix in DRUGGABLE_PREFIXES:
        if gene.startswith(prefix):
            return True
    for prefix in DRUGGABLE_EXACT_PREFIXES:
        if gene.startswith(prefix) and len(gene) > len(prefix):
            return True
    return False


def check_has_pdb(gene: str) -> bool:
    """Layer 4b: PDB结构可解析性"""
    return gene in HAS_PDB_SET


def compute_all_ppi_degrees(ppi_df: pd.DataFrame, gene_set: Set[str]) -> Dict[str, int]:
    """计算候选基因在全局PPI网络中的真实degree (v4: 一端在候选集即可)

    v3缺陷: 仅计算子网络内部连接, 严重低估hub基因degree
    修复: 计算该基因在整个PPI网络中与任意基因的连接数
    """
    if ppi_df is None:
        return {}

    ppi_a = ppi_df["preferredName_A"].str.upper()
    ppi_b = ppi_df["preferredName_B"].str.upper()

    mask = ppi_a.isin(gene_set) | ppi_b.isin(gene_set)
    sub_ppi = ppi_df[mask]

    counts_a = sub_ppi["preferredName_A"].str.upper().value_counts()
    counts_b = sub_ppi["preferredName_B"].str.upper().value_counts()

    degree = counts_a.add(counts_b, fill_value=0).to_dict()

    for g in gene_set:
        degree.setdefault(g, 0)

    return degree


def compute_betweenness_approx(ppi_df: pd.DataFrame, gene_set: Set[str],
                                top_n: int = 30) -> Dict[str, float]:
    """近似计算betweenness centrality (使用PageRank作为代理, 避免networkx依赖)

    当degree阈值过高导致Hub判定为0时, 作为二级判定指标
    """
    if ppi_df is None:
        return {}

    try:
        import networkx as nx
    except ImportError:
        return {}

    ppi_a = ppi_df["preferredName_A"].str.upper()
    ppi_b = ppi_df["preferredName_B"].str.upper()

    mask = ppi_a.isin(gene_set) | ppi_b.isin(gene_set)
    sub_ppi = ppi_df[mask]

    gene_to_idx = {g: i for i, g in enumerate(gene_set)}
    idx_to_gene = {i: g for g, i in gene_to_idx.items()}

    G = nx.Graph()
    G.add_nodes_from(range(len(gene_set)))

    for _, row in sub_ppi.iterrows():
        a = row["preferredName_A"].upper()
        b = row["preferredName_B"].upper()
        if a in gene_to_idx and b in gene_to_idx:
            G.add_edge(gene_to_idx[a], gene_to_idx[b])

    if G.number_of_edges() == 0:
        return {}

    try:
        bc = nx.betweenness_centrality(G, k=min(top_n * 3, len(G)))
        return {idx_to_gene[idx]: val for idx, val in bc.items() if val > 0}
    except Exception:
        return {}


def check_ciri_sub_pathways(gene: str) -> List[str]:
    """Layer 6: CIRI子通路命中列表"""
    hit_sub = []
    for sub_name, sub_genes in CIRI_PATHWAY_SUB.items():
        if gene in sub_genes:
            hit_sub.append(sub_name)
    return hit_sub


def compute_score(result: GeneFilterResult) -> float:
    """计算confidence_score (v4: 子通路权重差异化)"""
    score = 0.0

    score += result.known_target * 3
    score += result.brain_expressed_weight
    score += result.network_hub * 2
    score += result.has_pdb * 2
    score += result.druggable_family * 1
    score += result.cuproptosis * 2
    score += result.in_module * 1

    if result.ciri_pathway:
        score += _ciri_pathway_weight(result.ciri_sub_pathways)

    score -= result.blacklist_penalty

    if result.metabolic_enzyme:
        score -= 1

    if result.pathway_synergy:
        score += 1

    return score


def _ciri_pathway_weight(sub_pathways: List[str]) -> float:
    """子通路命中权重差异化 (v4设计优化)

    死亡相关子通路(ferroptosis/necroptosis/pyroptosis/apoptosis_ciri): 1.5
    炎症/应激子通路(nfkb/nrf2/jak_stat): 1.0
    其他(autophagy_ciri/upr): 0.8
    """
    if not sub_pathways:
        return 0.0

    weights = []
    for sub in sub_pathways:
        if sub in CIRI_DEATH_RELATED_SUBS:
            weights.append(1.5)
        elif sub in {"nfkb", "nrf2", "jak_stat"}:
            weights.append(1.0)
        else:
            weights.append(0.8)

    return max(weights)


def assign_tier(result: GeneFilterResult) -> str:
    """分配Tier (v4: Hard Blacklist一票否决 + 提前跳过)"""
    if result.hard_blacklisted:
        return "Rejected"

    if result.blacklisted and not (result.cuproptosis or result.ciri_pathway):
        return "Rejected"

    if result.score >= 6:
        return "Tier1"
    elif result.score >= 3:
        return "Tier2"
    elif result.score >= 1:
        return "Tier3"
    else:
        return "Rejected"


def generate_reason(result: GeneFilterResult) -> str:
    """生成原因字符串 (v4: 修复空括号问题)"""
    reasons = []
    if result.known_target:
        reasons.append("known_target")
    if result.brain_expressed:
        reasons.append(f"brain(w={result.brain_expressed_weight:.1f})")
    if result.has_pdb:
        reasons.append("has_pdb")
    if result.network_hub:
        reasons.append("network_hub")
    if result.druggable_family:
        reasons.append("druggable_family")
    if result.cuproptosis:
        reasons.append("cuproptosis")
    if result.ciri_pathway:
        sub_str = ",".join(result.ciri_sub_pathways) if result.ciri_sub_pathways else "cuproptosis"
        reasons.append(f"ciri_pathway({sub_str})")
    if result.pathway_synergy:
        reasons.append("pathway_synergy")
    if result.in_module:
        reasons.append("in_module")
    if result.metabolic_enzyme:
        reasons.append("metabolic_enzyme")
    if result.hard_blacklisted:
        reasons.append("hard_blacklisted")
    elif result.blacklisted:
        reasons.append(f"blacklisted(penalty={result.blacklist_penalty:.1f})")
    if result.intermediate_filament:
        reasons.append("intermediate_filament")

    if not reasons:
        return "no_evidence"

    return ";".join(reasons)


def main():
    print("=" * 60)
    print("BCP靶点分级清洗 v5: 探针ID修复 + 通路协同修正 + 性能优化")
    print("=" * 60)

    try:
        from config import BCP_TARGETS
    except ImportError as e:
        print(f"错误: 无法从 config 导入 BCP_TARGETS ({e})")
        sys.exit(1)

    if not isinstance(BCP_TARGETS, (list, tuple, set)) or len(BCP_TARGETS) == 0:
        print("错误: BCP_TARGETS 为空或类型错误")
        sys.exit(1)

    try:
        raw_targets = list({str(g).strip().upper() for g in BCP_TARGETS if str(g).strip()})
    except Exception as e:
        print(f"错误: BCP_TARGETS 处理失败 ({e})")
        sys.exit(1)

    if len(raw_targets) == 0:
        print("错误: 过滤后无有效靶点")
        sys.exit(1)

    print(f"\n原始靶点: {len(raw_targets)} 个")

    ppi_df = None
    ppi_file = os.path.join(RESULTS_DIR, "stage5_ppi_mcode", "string_ppi.tsv")
    if os.path.exists(ppi_file):
        try:
            ppi_df = pd.read_csv(ppi_file, sep="\t")
            print(f"PPI网络: {len(ppi_df)} 条边")
        except Exception as e:
            print(f"警告: PPI文件读取失败 ({e}), 跳过PPI分析")
            ppi_df = None

    mcode_file = os.path.join(RESULTS_DIR, "stage5_ppi_mcode", "mcode_modules.json")
    mcode_modules = {}
    module_genes_set: Set[str] = set()
    if os.path.exists(mcode_file):
        try:
            with open(mcode_file, "r") as f:
                mcode_modules = json.load(f)
            for mod_info in mcode_modules.values():
                module_genes_set.update(g.upper() for g in mod_info.get("genes", []))
            print(f"MCODE模块: {len(mcode_modules)} 个, 模块内基因: {len(module_genes_set)}")
        except Exception as e:
            print(f"警告: MCODE文件读取失败 ({e})")

    gene_set = {normalize_gene(g) for g in raw_targets}
    gene_set.discard(None)

    degree_counts = compute_all_ppi_degrees(ppi_df, gene_set)

    betweenness_dict = compute_betweenness_approx(ppi_df, gene_set)

    if degree_counts:
        all_degrees = [d for d in degree_counts.values() if d > 0]
        p90_degree = float(np.percentile(all_degrees, 90)) if all_degrees else 0
        hub_threshold = max(p90_degree, 5)
        print(f"PPI degree (全局): 有效节点={len(all_degrees)}, "
              f"median={np.median(all_degrees) if all_degrees else 0:.1f}, "
              f"P90={p90_degree:.1f}, hub_threshold={hub_threshold}")

        n_above_thresh = sum(1 for d in all_degrees if d > hub_threshold)
        print(f"  degree>{hub_threshold}: {n_above_thresh} 个基因")

        if n_above_thresh == 0 and betweenness_dict:
            bc_values = [v for v in betweenness_dict.values() if v > 0]
            if bc_values:
                bc_p90 = float(np.percentile(bc_values, 90))
                print(f"  degree无Hub, 启用betweenness P90={bc_p90:.6f} 作为二级判定")
    else:
        hub_threshold = 5
        print("PPI网络无连接数据, hub_threshold=5 (默认)")

    bc_p90 = 0.0
    if betweenness_dict:
        bc_values = [v for v in betweenness_dict.values() if v > 0]
        if bc_values:
            bc_p90 = float(np.percentile(bc_values, 90))

    print(f"\n开始分层过滤...")

    results = []

    for gene in raw_targets:
        normalized = normalize_gene(gene)
        if normalized is None:
            continue

        result = GeneFilterResult(normalized)

        result.known_target = check_known_target(normalized)

        is_brain, brain_w = check_brain_specificity(normalized)
        result.brain_expressed = is_brain
        result.brain_expressed_weight = brain_w

        bl, penalty, hard_bl = check_blacklist(normalized)
        result.blacklisted = bl
        result.blacklist_penalty = penalty
        result.hard_blacklisted = hard_bl

        result.metabolic_enzyme = normalized in METABOLIC_ENZYMES
        result.intermediate_filament = normalized in INTERMEDIATE_FILAMENTS

        result.druggable_family = check_druggable_family(normalized)
        result.has_pdb = check_has_pdb(normalized)

        deg = degree_counts.get(normalized, 0)
        result.network_hub = deg > hub_threshold

        if not result.network_hub and betweenness_dict:
            bc_val = betweenness_dict.get(normalized, 0)
            if bc_val > bc_p90:
                result.network_hub = True

        result.in_module = normalized in module_genes_set

        result.cuproptosis = (normalized in CUPROPTOSIS_GENES or
                              normalized in CUPROPTOSIS_RELATED)

        ciri_subs = check_ciri_sub_pathways(normalized)
        result.ciri_sub_pathways = ciri_subs
        result.ciri_pathway = len(ciri_subs) > 0 or result.cuproptosis

        result.oxidative_stress = normalized in OXIDATIVE_STRESS_GENES
        result.apoptosis = normalized in APOPTOSIS_GENES

        is_cupro_specific = normalized in CUPROPTOSIS_SPECIFIC_GENES
        death_subs = set(ciri_subs) & CIRI_DEATH_RELATED_SUBS
        result.pathway_synergy = result.cuproptosis and bool(death_subs)

        result.score = compute_score(result)
        result.tier = assign_tier(result)
        result.reason = generate_reason(result)

        results.append(result)

    tier_order = {"Tier1": 0, "Tier2": 1, "Tier3": 2, "Rejected": 3}
    results.sort(key=lambda r: (tier_order.get(r.tier, 4), -r.score))

    df = pd.DataFrame([
        {
            "Gene": r.gene,
            "Tier": r.tier,
            "Score": round(r.score, 2),
            "known_target": r.known_target,
            "brain_expressed": r.brain_expressed,
            "brain_weight": round(r.brain_expressed_weight, 1),
            "has_pdb": r.has_pdb,
            "network_hub": r.network_hub,
            "druggable_family": r.druggable_family,
            "in_module": r.in_module,
            "cuproptosis": r.cuproptosis,
            "ciri_pathway": r.ciri_pathway,
            "ciri_sub_pathways": ",".join(r.ciri_sub_pathways) if r.ciri_sub_pathways else "",
            "pathway_synergy": r.pathway_synergy,
            "metabolic_enzyme": r.metabolic_enzyme,
            "blacklisted": r.blacklisted,
            "blacklist_penalty": round(r.blacklist_penalty, 1),
            "hard_blacklisted": r.hard_blacklisted,
            "intermediate_filament": r.intermediate_filament,
            "Reason": r.reason,
        }
        for r in results
    ])

    df = df.drop_duplicates(subset=["Gene"], keep="first")
    df.to_csv(os.path.join(OUTPUT_DIR, "tiered_targets.csv"), index=False)

    tier1_genes = df[df["Tier"] == "Tier1"]["Gene"].tolist()
    tier2_genes = df[df["Tier"] == "Tier2"]["Gene"].tolist()
    tier3_genes = df[df["Tier"] == "Tier3"]["Gene"].tolist()
    rejected_genes = df[df["Tier"] == "Rejected"]["Gene"].tolist()

    df_idx = df.set_index("Gene")

    with open(os.path.join(OUTPUT_DIR, "filter_report.txt"), "w", encoding="utf-8") as f:
        f.write("BCP靶点分级清洗报告 v4\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"原始靶点: {len(raw_targets)}\n")
        f.write(f"过滤探针ID后有效靶点: {len(results)}\n")
        f.write(f"去重后: {len(df)}\n\n")

        f.write("各层统计:\n")
        f.write(f"  Layer 1 - 已知靶点: {df['known_target'].sum()}/{len(df)}\n")
        f.write(f"  Layer 2 - 脑组织特异性: {df['brain_expressed'].sum()}/{len(df)}\n")
        f.write(f"  Layer 3 - 负向过滤(blacklisted): {df['blacklisted'].sum()}/{len(df)}\n")
        f.write(f"  Layer 3 - Hard Blacklist(绝对拒绝): {df['hard_blacklisted'].sum()}/{len(df)}\n")
        f.write(f"  Layer 3 - 代谢酶: {df['metabolic_enzyme'].sum()}/{len(df)}\n")
        f.write(f"  Layer 4 - PDB可解析: {df['has_pdb'].sum()}/{len(df)}\n")
        f.write(f"  Layer 4 - 可成药家族: {df['druggable_family'].sum()}/{len(df)}\n")
        f.write(f"  Layer 5 - PPI network_hub(thresh={hub_threshold}): {df['network_hub'].sum()}/{len(df)}\n")
        f.write(f"  Layer 5 - MCODE模块内: {df['in_module'].sum()}/{len(df)}\n")
        f.write(f"  Layer 6 - 铜死亡基因: {df['cuproptosis'].sum()}/{len(df)}\n")
        f.write(f"  Layer 6 - CIRI子通路: {df['ciri_pathway'].sum()}/{len(df)}\n")
        f.write(f"  Layer 6 - 通路协同(cupro+death_sub): {df['pathway_synergy'].sum()}/{len(df)}\n\n")

        for sub_name, sub_genes in CIRI_PATHWAY_SUB.items():
            n_in_targets = len(gene_set & sub_genes)
            f.write(f"  CIRI子通路-{sub_name}: {n_in_targets} 命中")
            if sub_name in CIRI_DEATH_RELATED_SUBS:
                f.write(" (死亡相关, 权重1.5)")
            elif sub_name in {"nfkb", "nrf2", "jak_stat"}:
                f.write(" (炎症/应激, 权重1.0)")
            else:
                f.write(" (权重0.8)")
            f.write("\n")

        f.write("\nTier分布:\n")
        f.write(f"  Tier1: {len(tier1_genes)}\n")
        f.write(f"  Tier2: {len(tier2_genes)}\n")
        f.write(f"  Tier3: {len(tier3_genes)}\n")
        f.write(f"  Rejected: {len(rejected_genes)}\n\n")

        f.write("Tier1基因 (score>=6, 高置信度):\n")
        for g in tier1_genes:
            row = df_idx.loc[g]
            f.write(f"  {g} (score={row['Score']}, reason={row['Reason']})\n")

        f.write("\nTier2基因 (score>=3, 中等置信度):\n")
        for g in tier2_genes:
            f.write(f"  {g}\n")

        f.write("\nTier3基因 (score>=1, 低置信度):\n")
        for g in tier3_genes:
            row = df_idx.loc[g]
            f.write(f"  {g} (score={row['Score']})\n")

        f.write("\nRejected基因及原因:\n")
        rej_df = df[df["Tier"] == "Rejected"]
        for _, row in rej_df.iterrows():
            f.write(f"  {row['Gene']}: {row['Reason']} (score={row['Score']})\n")

        f.write(f"\n\n建议: 使用Tier1+Tier2 ({len(tier1_genes) + len(tier2_genes)}个)替换config.BCP_TARGETS\n")

    print(f"\n分级完成:")
    print(f"  Tier1: {len(tier1_genes)} 个")
    print(f"  Tier2: {len(tier2_genes)} 个")
    print(f"  Tier3: {len(tier3_genes)} 个")
    print(f"  Rejected: {len(rejected_genes)} 个")
    print(f"\n输出: {OUTPUT_DIR}/tiered_targets.csv")
    print(f"输出: {OUTPUT_DIR}/filter_report.txt")

    filtered_targets = tier1_genes + tier2_genes
    print(f"\n建议替换config.BCP_TARGETS为: {len(filtered_targets)} 个(Tier1+Tier2)")

    return filtered_targets


if __name__ == "__main__":
    main()
