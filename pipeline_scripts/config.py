# ============================================================
# BCP x Cuproptosis x CIRI Target Screening System - Configuration
# ============================================================
# Version: v3.0 | Date: 2026-05-12
# GPU: RTX 5060 (8GB)
# ============================================================

import os
import torch

# ============================================================
# 1. Path Configuration
# ============================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")

for d in [DATA_DIR, RESULTS_DIR, SCRIPTS_DIR]:
    os.makedirs(d, exist_ok=True)

for sub in ["deg_analysis", "go_kegg", "ppi_network", "wgcna", 
            "single_cell", "gat_results", "ml_results", "shap_results",
            "docking", "figures"]:
    os.makedirs(os.path.join(RESULTS_DIR, sub), exist_ok=True)

# ============================================================
# 2. Data File Paths
# ============================================================
GSE61616_TOP_TABLE = r"C:\Users\Jy-Mentor-7\Downloads\GSE61616.top.table (2).tsv"
GSE61616_RAW_TAR = r"C:\Users\Jy-Mentor-7\Downloads\GSE61616_RAW (1).tar"
GSE61616_MATRIX = r"C:\Users\Jy-Mentor-7\Downloads\GSE61616_series_matrix (1).txt.gz"
GSE61616_PLATFORM1 = r"C:\Users\Jy-Mentor-7\Downloads\GPL1355-10794 (1).txt"
GSE61616_PLATFORM2 = r"C:\Users\Jy-Mentor-7\Downloads\GPL1261-56135.txt"

GSE97537_RAW_TAR = r"C:\Users\Jy-Mentor-7\Downloads\GSE97537_RAW.tar"
GSE97537_SOFT = r"C:\Users\Jy-Mentor-7\Downloads\GSE97537_family (1).soft.gz"
GSE97537_MATRIX = r"C:\Users\Jy-Mentor-7\Downloads\GSE97537_series_matrix (2).txt.gz"

GSE210986_RAW_TAR = r"C:\Users\Jy-Mentor-7\Downloads\GSE210986_RAW.tar"
GSE210986_MATRIX = r"C:\Users\Jy-Mentor-7\Downloads\GSE210986_series_matrix.txt.gz"

RAT_MOUSE_HUMAN_MAP = os.path.join(BASE_DIR, "大创", "大鼠 小鼠 人类映射库.txt")

# ============================================================
# 3. BCP Targets (FIX:[P1-7][cleaned of housekeeping noise])
# ============================================================
# Housekeeping genes to exclude (Eisenberg & Levanon, Trends Genet 2013, PMID:23213612)
HOUSEKEEPING_GENES = {
    "ALB", "TBP", "ACTA1", "GAPDHS", "GAPDH", "B2M", "ACTB", "RPLP0",
    "HPRT1", "PPIA", "YWHAZ", "GUSB", "TFRC", "POLR2A", "B2M",
    "PGK1", "TUBB", "RPL13A", "SDHA", "UBC"
}

# FIX:[P1-7][remove housekeeping genes from BCP targets]
_RAW_BCP_TARGETS = [
    "LYN", "FABP2", "HLA-E", "MB", "VDR", "AMY1A", "AMY1B", "AMY1C",
    "ESR1", "BCL9", "FNTB", "F3", "PTPRF", "ALDH9A1", "FABP3", "PARN",
    "GCH1", "RENBP", "HPD", "OAZ1", "PA2G4", "MAP2", "PRKCQ", "ATP7B",
    "POLR2D", "C3", "MVK", "CYP2R1", "MARF1", "ACP3", "LDHA",
    "ALK", "HPGDS", "S100A6", "ACADVL", "GALNT10", "SMAP1", "ABAT",
    "FABP1", "TPK1", "RBP2", "GAD2", "UQCRC1", "DLAT", "CLIP4", "CCNA2",
    "PNLIPRP2", "STK4", "RARA", "HSD17B4", "ALDH3A1", "GGPS1", "NMT1",
    "GAD1", "FAH", "HBA1", "HBA2", "PPP2R1A", "EPS8L2", "PDK4", "ECHS1",
    "PDXK", "SPSB2", "NDC80", "HASPIN", "MECR", "RBM39", "HSD17B10",
    "CHKA", "CEL", "OAT", "ITGA1", "IL5", "MYH11", "IGKC", "DCD", "PTPRJ",
    "FGA", "PITPNB", "REPS2", "PDHA1", "IRAK4", "DDC", "HARS1", "STARD13",
    "PTGR1", "RBM38", "STK16", "SMARCA4", "PTGR2", "SRC", "ZHX2", "HMGB3",
    "ERBB4", "ACAD11", "TDP1", "HBS1L", "CDC42", "MAPKAPK2",
    "ATG5", "TTPA", "PAH", "SEC24C", "NUDC", "HBEGF",
    "MSH2", "DPP4", "RGS9",
    "CITED2", "CASK", "CUL4B", "CTSD", "CD1B", "LEF1", "MKNK2", "PDCD6",
    "MAN2B1", "HNMT", "CHFR", "SEC13", "AKR1C3", "HIBADH", "CHAT",
    "ZNF292", "ACY1", "MGAT1", "MARK2", "OXSR1",
    "IYD", "BST1", "SLC9A1", "FDXR", "CROT",
    "CSAD", "RAB27B", "PDCD6IP", "SAT1", "LCK", "FOXP2", "EXOSC9",
    "SERPINB10", "TOP1", "SNTA1", "MMUT", "MYBPC3", "AP1G1", "ATP2A1",
    "PARP12", "IL10RA", "NUDCD2", "AOC3", "PAPOLA", "MTCYB", "PCTP",
    "SAT2", "IGFBP2", "CNDP2", "PTPRC", "XRCC6", "STAM2", "REV1", "PDE6G",
    "PTGES2", "ACADM", "SERPINB1", "CDC26", "PYCR1", "THOP1", "RHOC",
    "NOTCH1", "AOC1", "ACOX1", "DPYD", "LILRA2", "EIF4A1", "C5", "TCN2",
    "ECE1", "KCNE1", "FUT8", "CRYL1", "NOS1", "POLK", "ZEB2", "ZEB1",
    "XBP1", "VCAM1", "TSPO", "TP53", "TOP2A", "TNF", "TGFB1", "TBXAS1",
    "TBXA2R", "TACR2", "TACR1", "SYP", "SYN1", "STAT3", "STAT1", "SREBF1",
    "SRD5A1", "SQLE", "SOD2", "SOD1", "SNAI1", "SLC6A4", "SLC6A3",
    "SLC6A2", "SIGMAR1", "SCN9A", "S1PR1", "RELA", "PTGS2", "PTGS1",
    "PRSS1", "PPARG", "PPARA", "PNPLA2", "PGR", "PARP1", "PABPC1", "P2RX7",
    "OPRM1", "OPRK1", "OPRD1", "NR1H3", "NR1H2", "NOS3", "NOS2", "NLRP3",
    "NFKB1", "NFE2L2", "MTOR", "MTNR1B", "MTNR1A", "MMP9", "MMP13",
    "MLYCD", "MGLL", "MDM2", "MAPK9", "MAPK14", "MAPK10", "LSS", "KIF11",
    "KEAP1", "KDR", "KCNA5", "JAK2", "JAK1", "IRF1", "INS1", "IMPDH2",
    "IL6", "IL1B", "IL10", "IKBKB", "IGF1R", "IDO1", "ICMT", "HTR6",
    "HTR2C", "HTR2B", "HTR2A", "HSPA5", "HSD11B1", "HRH2", "HMOX1",
    "HIF1A", "HCRTR1", "HAVCR1", "GSR", "GPX1", "GPT", "GFAP", "GAP43",
    "GABRG2", "GABRB3", "GABRA5", "GABRA3", "GABRA2", "GABRA1", "FOXO1",
    "FLT4", "FASN", "FAS", "F2", "EPHX2", "EPHX1", "EGFR", "DRD4", "DRD3",
    "DRD2", "DRD1", "DGAT1", "DDIT3", "CYP3A4", "CYP2D6", "CYP2C9",
    "CYP2C19", "CYP1A2", "CYP1A1", "CYP19A1", "CTSV", "CTSS", "CTSL",
    "CTSK", "CTSF", "CTSC", "CTSB", "CTRC", "CTGF", "CSNK1D", "CPT2",
    "CPT1A", "CP", "COL1A1", "CNR2", "CHRM5", "CHRM4", "CHRM3", "CDK4",
    "CDC25B", "CDC25A", "CCR2", "CCND1", "CCL2", "CAT", "CASP9", "CASP8",
    "CASP3", "BRD4", "BRD3", "BRD2", "BCL2", "BAX", "ATF4", "ARG1",
    "ALOX5", "ALDH1A1", "AKT1", "AIF1", "AHR", "ADRB3", "ADRB1",
    "ADRA2B", "ADRA2A", "ADRA1D", "ADRA1A", "ADORA3", "ADORA2A",
    "ADORA1", "ACTA2", "ACLY", "ACACB", "ACACA", "FAAH", "NR1I3", "AR",
    "CHRM2", "ESR2", "BCHE", "CXCR3", "PTPN2", "MAOB", "UGT2B7", "GLI2",
    "GLI1", "RORC", "SREBF2", "NPC1L1", "FABP4", "TERT", "SHBG",
    "CYP17A1", "FABP5", "PPARD", "HMGCR", "CNR1", "CYP51A1", "SRD5A2",
    "NR3C2", "NR3C1", "CES1", "SERPINA6", "NR1I2", "HSD17B3", "SCD",
    "ACHE", "RBP4", "CCR5", "FDFT1", "PTGER1", "PTGER2", "CES2", "PREP",
    "FFAR1", "PTGES", "RARG", "RARB", "SPHK2", "PLA2G4A", "SPHK1", "RORA",
    "FNTA", "APP", "TIMP1", "PTPN6", "ICAM1", "MMP1", "EGR1", "TRPV1",
    "RAPGEF3", "PTPN1", "STAT5A", "XDH", "TRPA1", "ABHD6"
]

def _load_bcp_targets():
    """Load BCP targets: prefer cleaned Tier1+Tier2, fallback to raw (minus housekeeping)"""
    tier_file = os.path.join(RESULTS_DIR, "target_filter", "tiered_targets.csv")
    if os.path.exists(tier_file):
        import pandas as pd
        df = pd.read_csv(tier_file)
        tier12 = df[df["Tier"].isin(["Tier1", "Tier2"])]["Gene"].tolist()
        if tier12:
            return [g.upper() for g in tier12]
    cleaned = [g.upper() for g in _RAW_BCP_TARGETS if g.upper() not in HOUSEKEEPING_GENES]
    return list(set(cleaned))

BCP_TARGETS = _load_bcp_targets()
print(f"BCP targets (cleaned): {len(BCP_TARGETS)}")

# ============================================================
# 4. Cuproptosis Genes (FIX:[P1-7][trimmed to validated core])
# ============================================================
# Core 17 genes from Tsvetkov Science 2022 (PMID:35298263)
CUPROPTOSIS_GENES = [
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHA1", "PDHB",
    "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7A", "ATP7B",
    "DLD", "DBT", "DLST", "PDHA2", "GCSH"
]

# FIX:[P1-7][related genes restricted to direct copper metabolism, removing non-specific antioxidants]
CUPROPTOSIS_RELATED = [
    # Copper transport/chaperones (direct cuproptosis pathway)
    "ATOX1", "COX17", "CCS", "COX11", "SCO1", "SCO2",
    "STEAP1", "STEAP2", "STEAP3", "STEAP4",
    "CP", "COMMD1",
    # Metallothioneins (copper binding)
    "MT1A", "MT2A",
]

# Cuproptosis pathway hierarchy (Tsvetkov Science 2022, PMID:35298263)
# FIX:[P0-4][pathway coreness scoring for Stage8]
CUPROPTOSIS_PATHWAY_SCORES = {
    # Upstream regulators (FDX1-dependent pathway)
    "FDX1": 100,   # essential upstream regulator
    "LIAS": 95,    # lipoyltransferase substrate
    "LIPT1": 90,   # lipoyltransferase
    # Lipoylation effectors (TCA cycle enzymes)
    "DLAT": 85,    # dihydrolipoamide S-acetyltransferase
    "PDHA1": 85,   # pyruvate dehydrogenase E1 alpha
    "PDHB": 85,    # pyruvate dehydrogenase E1 beta
    "DLD": 80,     # dihydrolipoamide dehydrogenase
    "DBT": 80,     # branched-chain ketoacid dehydrogenase
    "DLST": 80,    # oxoglutarate dehydrogenase
    "PDHA2": 75,   # PDHA paralog
    "GCSH": 75,    # glycine cleavage system H
    # Copper transport/regulation
    "ATP7A": 65,   # copper transporter
    "ATP7B": 65,   # copper transporter
    "SLC31A1": 60, # copper importer (CTR1)
    "MTF1": 55,    # metal-regulatory transcription factor
    # Cell cycle/apoptosis downstream
    "CDKN2A": 50,  # cell cycle arrest
    "GLS": 45,     # glutaminase (metabolic reprogramming)
}

# ============================================================
# 5. Analysis Parameters
# ============================================================
DEG_LOG2FC_THRESHOLD = 1.0
DEG_ADJ_P_THRESHOLD = 0.05

WGCNA_SOFT_THRESHOLD = 6
WGCNA_MIN_MODULE_SIZE = 30
WGCNA_MEDISSIM = 0.25

ML_RANDOM_STATE = 42
ML_N_FOLDS = 5
ML_TOP_N_TARGETS = 20

SC_QC_MIN_GENES = 200
SC_QC_MAX_GENES = 5000
SC_QC_MIN_COUNTS = 500
SC_QC_MAX_COUNTS = 50000
SC_QC_MAX_MITO = 0.20
SC_RESOLUTION = 0.8
SC_N_HVGS = 2000
SC_N_PCS = 50
SC_KNN_NEIGHBORS = 15

LIMMA_ADJ_PVAL = 0.05
LIMMA_MIN_LOGFC = 1.0

GRN_MIN_CORR = 0.05
GRN_MIN_EFFECT_SIZE = 0.10
GRN_TOP_N_CORR = 30
GRN_GENE_POOL_SIZE = 350  # FIX:[P1-4][parameterized gene pool size]

# ============================================================
# 6. GPU Configuration
# ============================================================
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {DEVICE}")
if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"Memory: {torch.cuda.get_device_properties(0).total_mem / 1e9:.1f} GB")

torch.manual_seed(ML_RANDOM_STATE)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(ML_RANDOM_STATE)

# ============================================================
# 7. Output Format
# ============================================================
FIG_FORMAT = "pdf"
FIG_DPI = 300
FIG_SIZE = (10, 8)

print("=" * 60)
print("Configuration loaded (v3.0)")
print(f"BCP targets (cleaned): {len(BCP_TARGETS)}")
print(f"Cuproptosis core: {len(CUPROPTOSIS_GENES)}")
print(f"Cuproptosis related: {len(CUPROPTOSIS_RELATED)}")
print(f"Output format: {FIG_FORMAT}, DPI: {FIG_DPI}")
print("=" * 60)
