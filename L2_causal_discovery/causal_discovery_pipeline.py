#!/usr/bin/env python3
"""
L2 因果发现：PC + NOTEARS-MLP 双阶段管线
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
阶段1: PC算法 → 条件独立性检验定无向骨架 (CPDAG)
阶段2: NOTEARS-MLP → 可微分无环性约束定方向 → 完整DAG

参考文献:
  [1] Zheng, Aragam, Ravikumar, Xing. DAGs with NO TEARS. NeurIPS 2018.
  [2] Zheng, Dan, Aragam, Ravikumar, Xing. Learning Sparse Nonparametric DAGs. AISTATS 2020.
  [3] Spirtes, Glymour, Scheines. Causation, Prediction, and Search. 2000.
  [4] Tsvetkov et al. Copper induces cell death by targeting lipoylated TCA cycle proteins. Science 2022.
  [5] gCastle: A Python Toolbox for Causal Discovery. arXiv 2111.15155.

实现框架:
  - PC算法: gCastle
  - NOTEARS-MLP: 自实现 (gCastle封装)
  - 铜死亡参考通路: 从Tsvetkov 2022 + Cobine & Brady 2022 curated

自检标准:
  Acyclicity score < 0.01
  SHD vs curated cuproptosis pathway < 10
  5-fold CV prediction error < 0.1
"""

import numpy as np
import pandas as pd
import networkx as nx
import os
import json
import warnings
from copy import deepcopy
from scipy.linalg import expm

warnings.filterwarnings('ignore')
np.random.seed(42)

PROJECT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
CIRI_DIR = os.path.join(PROJECT_DIR, "ciri-cuproptosis-causal-discovery")
L2_DIR = os.path.join(CIRI_DIR, "l2_causal_discovery")
RESULTS_DIR = os.path.join(CIRI_DIR, "results", "L2_causal_discovery")
L1_RESULTS = os.path.join(CIRI_DIR, "results", "L1_phenotype_anchoring")

os.makedirs(L2_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# ===========================================================================
# 0. 铜死亡文献 curated 参考通路 (Ground Truth for SHD)
# ===========================================================================
# 来源: Tsvetkov 2022 (Science), Cobine & Brady 2022 (Mol Cell), 综述 (Frontiers 2025)
# 已知边: (上游 → 下游) 列表

LITERATURE_EDGES = [
    # ── Cu import / transport ──
    ("SLC31A1", "ATP7A"),       # CTR1 imports Cu, ATP7A exports
    ("ATP7A", "MT1A"),          # ATP7A delivers Cu to metallothioneins
    ("ATP7A", "MT2A"),
    ("SLC31A1", "ATOX1"),       # CTR1 → ATOX1 Cu chaperone
    ("ATOX1", "ATP7A"),         # ATOX1 delivers Cu to ATP7A

    # ── Cu chaperones → target enzymes ──
    ("CCS", "SOD1"),            # CCS delivers Cu to SOD1
    ("COX17", "MT-CO1"),        # COX17 delivers Cu to cytochrome c oxidase

    # ── Cuproptosis core pathway ──
    ("FDX1", "LIAS"),           # FDX1 reduces Fe-S cluster for LIAS
    ("LIAS", "DLAT"),           # LIAS lipoylates DLAT (lipoylation)
    ("LIAS", "DBT"),            # LIAS lipoylates DBT
    ("LIAS", "GCSH"),           # LIAS lipoylates GCSH (glycine cleavage system H)

    # ── TCA cycle / PDH complex ──
    ("DLAT", "PDHA1"),          # DLAT (E2) → PDHA1 (E1) in PDH complex
    ("DLAT", "PDHB"),
    ("DLD", "PDHA1"),           # DLD (E3) → E1
    ("DLD", "PDHB"),

    # ── Oxidative phosphorylation ──
    ("LIAS", "DLD"),            # LIAS → DLD (shared lipoylation domain)
    ("DLD", "DLST"),            # DLD ↔ DLST (E3 + E2 of α-KGDH)
    ("DLAT", "DLST"),           # Both E2 enzymes, shared mechanism

    # ── Fe-S cluster proteins ──
    ("FDX1", "FDXR"),           # FDX1 ↔ FDXR (ferredoxin reductase)
    ("FDX1", "CYP11A1"),        # FDX1 donates e- to CYP enzymes

    # ── Glutathione / Redox ──
    ("MTF1", "MT1A"),           # MTF1 (metal-responsive TF) → MT genes
    ("MTF1", "MT2A"),
    ("MTF1", "SLC31A1"),        # MTF1 regulates CTR1 transcription
    ("NFE2L2", "GCLC"),         # NFE2L2 (NRF2) → GCLC, GCLM for GSH synthesis
    ("NFE2L2", "GCLM"),
    ("NFE2L2", "SOD1"),         # NFE2L2 → antioxidant genes

    # ── Inflammatory / NLRP3 ──
    ("ATP7A", "NLRP3"),         # Cu efflux loss activates NLRP3 inflammasome
    ("NFE2L2", "NLRP3"),        # NRF2 inhibits NLRP3

    # ── Apoptosis intersection ──
    ("FDX1", "BAX"),            # FDX1 → BAX (mitochondrial apoptosis)
    ("SLC31A1", "BAX"),         # Cu overload → BAX activation

    # ── Negative feedback ──
    ("ATP7B", "SLC31A1"),       # ATP7B regulates Cu import
]

LITERATURE_NODES = sorted(set(
    [n for edge in LITERATURE_EDGES for n in edge]
))
GOLD_STANDARD = set(LITERATURE_EDGES)

print("=" * 70)
print("L2 因果发现: PC(NOTEARS-MLP)")
print("=" * 70)
print(f"\n文献通路: {len(GOLD_STANDARD)} 条已知边, {len(LITERATURE_NODES)} 个基因")
print()


# ===========================================================================
# 1. 数据准备: 从L1结果提取表达矩阵
# ===========================================================================
print("=" * 70)
print("[1/5] 数据准备\n")

# 1a. scRNA-seq细胞类型特异性log2FC矩阵
print("  1a. 加载 scRNA-seq celltype DEGs...")
ct_deg = pd.read_csv(os.path.join(L1_RESULTS, "celltype_cuproptosis_DEGs.csv"))
pivot_lfc = ct_deg.pivot_table(
    index="gene", columns="cell_type", values="log2FC"
)

cuproptosis_genes = sorted(ct_deg["gene"].unique())
print(f"  铜死亡基因数: {len(cuproptosis_genes)}")

# 排除 Unknown
if "Unknown" in pivot_lfc.columns:
    pivot_lfc = pivot_lfc.drop(columns=["Unknown"])

print(f"  细胞类型: {list(pivot_lfc.columns)}")
print(f"  log2FC矩阵: {pivot_lfc.shape}")

# 1b. 构建多条件伪观测数据 (bootstrap from celltype means)
#