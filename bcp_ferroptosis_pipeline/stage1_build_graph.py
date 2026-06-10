# -*- coding: utf-8 -*-
"""
阶段1 - 铁衰老基因数据库与异构图构建
==============================================
整合多源数据构建"石竹烯-铁衰老-脑缺血再灌注"异构图

数据源:
  - FerrDb V2: 铁死亡基因数据库 (http://www.zhounan.org/ferrdb/)
  - GenAge / CellAge: 衰老相关基因
  - GEO多数据集: CIRI差异表达基因
  - STRING: 蛋白互作网络
  - KEGG / Reactome: 信号通路
  - Cell Metabolism 2026: 灵长类铁衰老特征基因 (Vitamin C inhibits ACSL4...)

节点类型:
  - Drug: 石竹烯及其结构类似物
  - Gene/Protein: 铁死亡+衰老+CIRI相关靶点
  - Disease: 脑缺血再灌注
  - Pathway: 信号通路
  - Phenotype: 衰老/铁死亡表型

边类型:
  - drug_targets_gene: 药物-靶点关系
  - gene_interacts_gene: 蛋白互作 (PPI)
  - gene_belongs_to_pathway: 基因-通路
  - gene_associated_with_disease: 基因-疾病关联
  - pathway_related_to_phenotype: 通路-表型
"""

import os
import sys
import json
import csv
import gzip
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict

import numpy as np
import pandas as pd
import torch
from torch_geometric.data import HeteroData

# ============================================================
# 路径配置
# ============================================================
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data" / "ferroptosis_graph"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. 铁死亡核心基因集 (FerrDb V2 + 文献验证)
# ============================================================

# FerrDb V2 铁死亡驱动基因 (Driver genes)
FERROPTOSIS_DRIVERS = [
    "ACSL4",    # 核心: 长链脂酰辅酶A合成酶4 — 催化PUFA-PL合成，铁死亡执行者
    "LPCAT3",   # 溶血磷脂酰胆碱酰基转移酶3 — PUFA掺入膜磷脂
    "ALOX5",    # 花生四烯酸5-脂氧合酶 — 脂质过氧化
    "ALOX12",   # 花生四烯酸12-脂氧合酶
    "ALOX15",   # 花生四烯酸15-脂氧合酶
    "ALOXE3",   # 表皮型脂氧合酶3
    "TP53",     # p53 — 转录抑制SLC7A11，促进铁死亡
    "BAP1",     # BRCA1相关蛋白1 — 去泛素化H2A，抑制SLC7A11
    "CS",       # 柠檬酸合酶 — TCA循环
    "RPL8",     # 核糖体蛋白L8 — 铁死亡检查点
    "VDAC2",    # 电压依赖性阴离子通道2
    "VDAC3",    # 电压依赖性阴离子通道3
    "TFRC",     # 转铁蛋白受体 — 铁摄取
    "NCOA4",    # 核受体共激活因子4 — 铁蛋白自噬
    "HMOX1",    # 血红素加氧酶1 — 铁释放
    "SAT1",     # 精胺/精胺N1-乙酰转移酶1
    "GLS2",     # 谷氨酰胺酶2 — 谷氨酰胺分解
    "CARS",     # 半胱氨酰tRNA合成酶
    "ATL1",     # Atlastin GTPase 1
    "TTC35",    # Tetratricopeptide repeat domain 35
    "ATF4",     # 激活转录因子4 — 应激反应
    "CHAC1",    # ChaC谷胱甘肽特异性γ-谷氨酰环转移酶1
    "SLC1A5",   # 中性氨基酸转运蛋白
    "DPP4",     # 二肽基肽酶4
    "NFE2L2",   # NRF2 — 抗氧化主调控因子 (驱动取决于上下文)
    "HSPB1",    # 热休克蛋白B1 — 抑制铁摄取
    "CISD1",    # CDGSH铁硫结构域1 — 线粒体铁死亡
    "FANCD2",   # Fanconi贫血互补群D2
    "SLC7A11",  # xCT — 胱氨酸/谷氨酸反向转运蛋白 (核心抑制因子)
    "GPX4",     # 谷胱甘肽过氧化物酶4 — 铁死亡核心抑制因子
    "FTH1",     # 铁蛋白重链1 — 铁储存
    "AIFM2",    # FSP1 — 非GPX4铁死亡抑制通路
    "GCH1",     # GTP环化水解酶1 — BH4合成
    "SLC40A1",  # Ferroportin — 铁外排
    "PROM2",    # Prominin 2 — 铁外排囊泡
    "NF2",      # Merlin — Hippo通路
    "CAV1",     # Caveolin-1
    "PCBP1",    # Poly(rC)结合蛋白1 — 铁伴侣
    "PCBP2",    # Poly(rC)结合蛋白2
    "IREB2",    # 铁调节元件结合蛋白2
    "ACSF2",    # 酰基辅酶A合成酶家族成员2
    "GCLC",     # 谷氨酸-半胱氨酸连接酶催化亚基
    "GCLM",     # 谷氨酸-半胱氨酸连接酶调节亚基
    "GSS",      # 谷胱甘肽合成酶
    "SLC3A2",   # 溶质载体家族3成员2 (4F2hc)
    "KEAP1",    # Kelch样ECH相关蛋白1
    "SQSTM1",   # p62 — 自噬受体
    "BECN1",    # Beclin-1 — 自噬
    "ATG5",     # 自噬相关5
    "ATG7",     # 自噬相关7
    "MAP1LC3B", # LC3B
    "MTOR",     # 雷帕霉素靶蛋白
    "PRKAA1",   # AMPKα1
    "PRKAA2",   # AMPKα2
    "HIF1A",    # 缺氧诱导因子1α
    "EGLN1",    # PHD2 — 脯氨酸羟化酶
    "EPAS1",    # HIF2α
    "TF",       # 转铁蛋白
    "STEAP3",   # STEAP3金属还原酶
    "SLC11A2",  # DMT1 — 二价金属转运蛋白1
    "ACO1",     # IRP1 — 铁调节蛋白1
    "FBXL5",    # F-box和富含亮氨酸重复蛋白5
    "CDKN1A",   # p21
    "MDM2",     # MDM2原癌基因
    "MDM4",     # MDM4
    "AURKA",    # Aurora激酶A
    "BRD4",     # 溴结构域蛋白4
    "BRD7",     # 溴结构域蛋白7
    "ENPP2",    # 外核苷酸焦磷酸酶/磷酸二酯酶2
    "CRYAB",    # αB-晶状体蛋白
    "ARNTL",    # BMAL1 — 昼夜节律
    "KDM3B",    # 赖氨酸去甲基化酶3B
    "KDM4A",    # 赖氨酸去甲基化酶4A
    "KDM5C",    # 赖氨酸去甲基化酶5C
    "SIRT1",    # Sirtuin 1
    "NRF1",     # 核呼吸因子1
    "TFAM",     # 线粒体转录因子A
    "SIRT3",    # Sirtuin 3 — 线粒体
    "PPARG",    # PPARγ
    "PPARA",    # PPARα
    "PPARD",    # PPARδ
    "RXRA",     # RXRα
    "LXR",      # 肝X受体
    "PPARGC1A", # PGC-1α
    "NR1H3",    # LXRα
    "NR1H4",    # FXR
    "ABCA1",    # ATP结合盒亚家族A成员1
    "SCD",      # 硬脂酰辅酶A去饱和酶
    "FADS2",    # 脂肪酸去饱和酶2
    "ELOVL5",   # ELOVL脂肪酸延长酶5
    "ELOVL2",   # ELOVL脂肪酸延长酶2
    "ACACA",    # ACC1 — 乙酰辅酶A羧化酶
    "FASN",     # 脂肪酸合酶
    "ACSL3",    # 长链脂酰辅酶A合成酶3
    "ACSL5",    # 长链脂酰辅酶A合成酶5
    "ACSL6",    # 长链脂酰辅酶A合成酶6
    "AGPAT3",   # 1-酰基甘油-3-磷酸O-酰基转移酶3
    "LPCAT1",   # 溶血磷脂酰胆碱酰基转移酶1
    "PLA2G6",   # 磷脂酶A2 VI组
    "PLA2G4A",  # 胞质磷脂酶A2
    "PLA2G4C",  # 磷脂酶A2 IVC组
    "PNPLA2",   # ATGL — 脂肪甘油三酯脂肪酶
    "PNPLA8",   # 含patatin样磷脂酶结构域8
    "LIPE",     # HSL — 激素敏感脂肪酶
    "MGLL",     # 单酰甘油脂肪酶
    "CPT1A",    # 肉碱棕榈酰转移酶1A
    "CPT2",     # 肉碱棕榈酰转移酶2
    "ACADM",    # 中链酰基辅酶A脱氢酶
    "ACADVL",   # 极长链酰基辅酶A脱氢酶
    "ACOX1",    # 酰基辅酶A氧化酶1
    "HADHA",    # 羟酰辅酶A脱氢酶α亚基
]

# 铁死亡抑制基因 (Suppressor genes)
FERROPTOSIS_SUPPRESSORS = [
    "GPX4", "FSP1", "DHODH", "GCH1", "BH4", "SLC7A11",
    "NFE2L2", "SLC40A1", "PROM2", "FTH1", "FTL",
    "HSPB1", "CISD1", "CISD2", "FANCD2", "NF2",
    "CAV1", "PCBP1", "PCBP2", "GCLC", "GCLM", "GSS",
    "TXNRD1", "PRDX1", "PRDX6", "SOD1", "SOD2", "CAT",
    "GPX1", "GSR", "TXN", "SRXN1", "HMOX1",
]

# 铁死亡标志基因 (Marker genes)
FERROPTOSIS_MARKERS = [
    "PTGS2", "CHAC1", "SLC7A11", "GPX4", "ACSL4",
    "TFRC", "FTH1", "HMOX1", "SAT1", "ATF3",
    "DDIT3", "HSPA5", "TNFRSF10B", "TNFRSF10A",
]

ALL_FERROPTOSIS_GENES = list(set(
    FERROPTOSIS_DRIVERS + FERROPTOSIS_SUPPRESSORS + FERROPTOSIS_MARKERS
))

# ============================================================
# 2. 衰老相关基因 (GenAge + CellAge + 文献)
# ============================================================

# GenAge 人类衰老相关基因 (核心集)
AGING_GENES = [
    "TP53", "SIRT1", "SIRT3", "SIRT6", "FOXO1", "FOXO3",
    "MTOR", "IGF1R", "AKT1", "PIK3CA", "PTEN",
    "NFE2L2", "KEAP1", "SQSTM1", "ATG5", "ATG7",
    "BECN1", "MAP1LC3B", "ULK1", "ULK2",
    "CDKN1A", "CDKN2A", "CDKN2B", "RB1", "E2F1",
    "LMNA", "WRN", "BLM", "TERF1", "TERF2", "TERT",
    "PARP1", "SOD1", "SOD2", "CAT", "GPX1", "GPX4",
    "HMOX1", "NFKB1", "RELA", "IL6", "TNF", "IL1B",
    "TGFB1", "CCL2", "CXCL8", "HMGB1", "HMGB2",
    "PPARG", "PPARA", "PPARGC1A", "NR1H3",
    "HIF1A", "EPAS1", "VEGFA", "KDR",
    "GSK3B", "CTNNB1", "NOTCH1", "HES1",
    "SREBF1", "SREBF2", "FASN", "ACACA",
    "ATF4", "DDIT3", "ERN1", "XBP1", "HSPA5",
    "CD38", "NAMPT", "NNMT", "KAT2A", "HDAC1", "HDAC3",
    "DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2",
    "SASP", "CXCL1", "CXCL2", "MMP3", "MMP9",
    "ICAM1", "VCAM1", "SELP", "SELE",
    "CDKN1B", "CDK2", "CDK4", "CDK6", "CCND1",
    "AR", "ESR1", "ESR2", "PGR",
    "H2AFX", "ATM", "ATR", "CHEK1", "CHEK2",
    "BRCA1", "BRCA2", "RAD51", "XRCC6", "XRCC5",
    "MLH1", "MSH2", "MSH6", "PMS2",
    "GATA4", "MEF2A", "MEF2C", "MYC", "JUN", "FOS",
    "EGR1", "SP1", "STAT1", "STAT3", "STAT5A",
    "CREBBP", "EP300", "KMT2A", "KMT2D", "EZH2",
    "SUZ12", "BMI1", "RING1", "CBX4", "CBX8",
]

# Cell Metabolism 2026 灵长类铁衰老特征基因
# 基于论文 "Vitamin C inhibits ACSL4 to alleviate ferro-aging in primates"
# 基因来源于: C:\Users\Jy-Mentor-7\Desktop\铁衰老\铁衰老基因.txt

def _load_ferro_aging_genes(filepath: str = None) -> List[str]:
    """从文件加载铁衰老基因列表。"""
    if filepath is None:
        filepath = r"C:\Users\Jy-Mentor-7\Desktop\铁衰老\铁衰老基因.txt"
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            genes = [line.strip() for line in f if line.strip()]
        print(f"  从 {filepath} 加载了 {len(genes)} 个铁衰老基因")
        return genes
    except FileNotFoundError:
        print(f"  警告: 铁衰老基因文件未找到 {filepath}，使用内置基因集")
        return []

FERRO_AGING_GENES_FILE = _load_ferro_aging_genes()

FERRO_AGING_SIGNATURE = {
    "upregulated": [
        "ACSL4", "LPCAT3", "ALOX5", "ALOX15", "PTGS2",
        "TFRC", "HMOX1", "SAT1", "CHAC1", "ATF3",
        "DDIT3", "CDKN1A", "CDKN2A", "IL6", "TNF",
        "IL1B", "CCL2", "CXCL8", "MMP9", "ICAM1",
        "VCAM1", "NFKB1", "RELA", "HMGB1", "HIF1A",
        "NLRP3", "NOX4", "TLR4", "S100A8", "IFNG",
    ],
    "downregulated": [
        "GPX4", "FSP1", "SLC7A11", "FTH1", "FTL",
        "SLC40A1", "GCH1", "NFE2L2", "KEAP1",
        "SIRT1", "SIRT3", "SIRT6", "FOXO3",
        "SOD1", "SOD2", "CAT", "GSR",
        "NAMPT", "PPARGC1A",
    ],
}

ALL_AGING_GENES = list(set(
    AGING_GENES +
    FERRO_AGING_SIGNATURE["upregulated"] +
    FERRO_AGING_SIGNATURE["downregulated"] +
    FERRO_AGING_GENES_FILE  # 来自Cell Metabolism 2026论文的铁衰老基因
))


# ============================================================
# 3. ACSL4 中心调控网络 (文献挖掘)
# ============================================================

# ACSL4 一阶互作蛋白 (STRING PPI + BioGRID)
ACSL4_FIRST_NEIGHBORS = {
    # 直接互作蛋白
    "direct_interactors": [
        "LPCAT3",   # ACSL4合成PUFA-CoA → LPCAT3掺入膜磷脂 (协同)
        "GPX4",     # ACSL4促脂质过氧化 ↔ GPX4拮抗 (功能互作)
        "ALOX5",    # ACSL4提供底物 → ALOX5催化过氧化
        "ALOX15",   # 同上
        "SLC7A11",  # ACSL4依赖SLC7A11输入的半胱氨酸 (间接)
        "TFRC",     # 铁摄取 ↑ → ACSL4活性依赖铁
        "FTH1",     # 铁储存 ↓ → 游离铁 → ACSL4
        "PPARG",    # ACSL4受PPARγ转录调控
        "SP1",      # ACSL4启动子含SP1结合位点
        "NFE2L2",   # NRF2通过ARE调控ACSL4相关脂质代谢基因
    ],
    # 间接调控关系 (信号通路层次)
    "indirect_regulators": [
        "TP53",     # p53 → 转录抑制SLC7A11 → ACSL4依赖
        "NFKB1",    # NF-κB → 炎性正反馈 → 上调ACSL4
        "HIF1A",    # HIF-1α → 缺氧 → 铁代谢重编程
        "MTOR",     # mTORC1 → 脂质合成 → ACSL4底物可及性
        "PIK3CA",   # PI3K/AKT → 细胞存活 → 拮抗铁死亡
        "AKT1",     # AKT → 磷酸化调节
        "STAT3",    # STAT3 → ACSL4转录 (JAK/STAT)
        "SIRT1",    # SIRT1 → 去乙酰化 → 衰老调控
        "AMPK",     # PRKAA1/2 → 能量应激 → 脂质代谢
        "KEAP1",    # KEAP1-NRF2 → 氧化还原平衡
    ],
}


# ============================================================
# 4. 脑缺血再灌注 (CIRI) 相关基因
# ============================================================

# 从GEO多数据集 (GSE61616, GSE97537, GSE174574等) 整合的CIRI一致差异基因
# 基于前期分析结果
CIRI_DEGS = {
    "upregulated_in_ischemia": [
        "IL6", "TNF", "IL1B", "CCL2", "CXCL8", "ICAM1", "VCAM1",
        "MMP9", "MMP3", "TIMP1", "HMOX1", "HSPA5", "DDIT3",
        "ATF3", "ATF4", "FOS", "JUN", "EGR1", "NFKB1", "RELA",
        "STAT3", "HIF1A", "VEGFA", "PTGS2", "NOS2", "SERPINE1",
        "CDKN1A", "BAX", "CASP3", "CASP8", "CASP9", "BID",
        "CYCS", "APAF1", "AIFM1", "ENDOG", "PARP1",
        "TLR2", "TLR4", "MYD88", "HMGB1", "S100A8", "S100A9",
        "NLRP3", "IL18", "CASP1", "GSDMD",
        "TFRC", "ACSL4", "LPCAT3", "ALOX5", "ALOX15",
        "CHAC1", "SAT1", "PTGS2", "SLC7A11",
    ],
    "downregulated_in_ischemia": [
        "GPX4", "FSP1", "GCH1", "FTH1", "FTL", "SLC40A1",
        "NFE2L2", "KEAP1", "GCLC", "GCLM", "GSR", "SOD1", "SOD2",
        "CAT", "PRDX1", "TXN", "SRXN1",
        "BDNF", "NGF", "GDNF", "NTF3", "NTRK2",
        "SYN1", "SYP", "DLG4", "GRIN1", "GRIN2A", "GRIN2B",
        "GABRA1", "GABRB2", "GABRG2",
        "BCL2", "BCL2L1", "MCL1", "XIAP", "BIRC5",
        "SIRT1", "SIRT3", "PGC1A", "TFAM", "NRF1",
        "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP6",
        "CLDN5", "OCLN", "TJP1", "CDH5",
        "KDR", "FLT1", "TEK", "ANGPT1", "PDGFB",
        "CNR2", "OPRK1", "OPRD1", "OPRM1",
    ],
}

ALL_CIRI_GENES = list(set(
    CIRI_DEGS["upregulated_in_ischemia"] +
    CIRI_DEGS["downregulated_in_ischemia"]
))

# ============================================================
# 5. 石竹烯已知靶点 (来自config.py BCP_TARGETS)
# ============================================================
# 从工作区配置加载
sys.path.insert(0, str(BASE_DIR))
try:
    from config import BCP_TARGETS
except ImportError:
    BCP_TARGETS = [
        "CNR2", "TRPV1", "PPARG", "PPARA", "NFE2L2", "HMOX1",
        "PTGS2", "TNF", "IL6", "IL1B", "NFKB1", "RELA",
        "SLC7A11", "GPX4", "TP53", "AKT1", "MTOR", "STAT3",
        "MAPK14", "MAPK9", "MAPK10", "JAK2", "EGFR", "SRC",
        "CASP3", "BAX", "BCL2", "HIF1A", "VEGFA", "MMP9",
        "CAT", "SOD1", "SOD2", "GSR", "GCLC",
        "ACSL4",  # 假设石竹烯可通过CNR2→PPARγ→ACSL4轴调控
    ]

# ============================================================
# 6. 通路数据
# ============================================================

# 铁死亡核心通路
FERROPTOSIS_PATHWAYS = {
    "Ferroptosis (KEGG: hsa04216)": [
        "GPX4", "SLC7A11", "SLC3A2", "ACSL4", "LPCAT3",
        "TFRC", "STEAP3", "SLC11A2", "FTH1", "FTL",
        "HMOX1", "PCBP1", "PCBP2", "IREB2", "ACO1",
        "GCLC", "GCLM", "GSS", "GSR",
        "ALOX5", "ALOX12", "ALOX15", "PTGS2",
        "TP53", "BAP1", "SAT1", "CARS",
        "NFE2L2", "KEAP1", "ATF4", "DDIT3",
        "PRKAA1", "PRKAA2", "MTOR",
        "VDAC2", "VDAC3", "SLC40A1",
    ],
    "Glutathione metabolism (KEGG: hsa00480)": [
        "GCLC", "GCLM", "GSS", "GSR", "GPX4", "GPX1",
        "GGT1", "GGT5", "ANPEP", "OPLAH",
    ],
    "Lipid metabolism / PUFA biosynthesis": [
        "ACSL4", "ACSL3", "ACSL5", "ACSL6",
        "LPCAT3", "LPCAT1", "LPCAT2",
        "FADS1", "FADS2", "ELOVL2", "ELOVL5",
        "PLA2G4A", "PLA2G6", "PLA2G4C",
        "AGPAT3", "AGPAT4", "AGPAT5",
        "FASN", "ACACA", "SCD",
    ],
    "Iron homeostasis": [
        "TFRC", "TF", "STEAP3", "SLC11A2", "SLC40A1",
        "FTH1", "FTL", "HMOX1", "NCOA4", "IREB2",
        "ACO1", "FBXL5", "PCBP1", "PCBP2",
        "ISCU", "FXN", "BOLA3", "NFU1", "GLRX5",
    ],
    "Autophagy & ferroptosis": [
        "SQSTM1", "BECN1", "ATG5", "ATG7", "ATG13",
        "MAP1LC3B", "ULK1", "ULK2", "NCOA4", "RAB7A",
        "LAMP1", "LAMP2", "CTSD", "CTSL",
    ],
    "Necroptosis / Pyroptosis (交叉)": [
        "RIPK1", "RIPK3", "MLKL", "NLRP3", "CASP1",
        "GSDMD", "IL1B", "IL18", "HMGB1",
    ],
    "Aging-related pathways": [
        "TP53", "CDKN1A", "CDKN2A", "CDKN2B",
        "SIRT1", "SIRT3", "SIRT6", "FOXO1", "FOXO3",
        "MTOR", "IGF1R", "PIK3CA", "AKT1", "PTEN",
        "LMNA", "WRN", "TERT", "TERF1", "TERF2",
        "NFE2L2", "KEAP1", "PPARGC1A", "TFAM",
        "SOD1", "SOD2", "CAT",
        "NFKB1", "RELA", "IL6", "TNF",
    ],
    "CIRI-related pathways": [
        # 兴奋毒性
        "GRIN1", "GRIN2A", "GRIN2B", "GRIA1", "GRIA2",
        # 氧化应激
        "NOS1", "NOS2", "NOS3", "NOX1", "NOX2", "NOX4",
        "SOD1", "SOD2", "CAT", "PRDX1", "TXN",
        # 炎症
        "TLR2", "TLR4", "MYD88", "NFKB1", "RELA",
        "IL6", "TNF", "IL1B", "CCL2", "CXCL8",
        "PTGS2", "ALOX5",
        # 凋亡
        "BCL2", "BAX", "CASP3", "CASP8", "CASP9",
        "CYCS", "APAF1", "PARP1",
        # 血管
        "VEGFA", "KDR", "FLT1", "ANGPT1", "TEK",
        "CLDN5", "OCLN", "TJP1", "MMP9",
        # 神经营养
        "BDNF", "NGF", "NTRK2", "GDNF",
    ],
}

# ============================================================
# 7. 表型节点定义
# ============================================================

PHENOTYPES = {
    "ferroptosis": {
        "description": "铁死亡 — 铁依赖的脂质过氧化驱动的非凋亡细胞死亡",
        "hallmark_genes": ["GPX4", "ACSL4", "SLC7A11", "TFRC", "PTGS2"],
    },
    "cellular_senescence": {
        "description": "细胞衰老 — 不可逆的细胞周期停滞 + SASP",
        "hallmark_genes": ["CDKN1A", "CDKN2A", "TP53", "LMNA", "IL6"],
    },
    "ferro_aging": {
        "description": "铁衰老 — 铁死亡驱动的年龄相关退行性变",
        "hallmark_genes": ["ACSL4", "GPX4", "SLC7A11", "FTH1", "CDKN1A"],
    },
    "lipid_peroxidation": {
        "description": "脂质过氧化 — PUFA-PL过氧化级联",
        "hallmark_genes": ["ACSL4", "LPCAT3", "ALOX5", "ALOX15", "GPX4"],
    },
    "neuroinflammation": {
        "description": "神经炎症 — 缺血后胶质激活与炎症因子释放",
        "hallmark_genes": ["IL6", "TNF", "IL1B", "CCL2", "NFKB1"],
    },
    "oxidative_stress": {
        "description": "氧化应激 — ROS产生超过抗氧化防御能力",
        "hallmark_genes": ["NFE2L2", "SOD1", "SOD2", "CAT", "HMOX1"],
    },
    "mitochondrial_dysfunction": {
        "description": "线粒体功能障碍 — ATP耗竭与线粒体膜通透性转换",
        "hallmark_genes": ["TFAM", "PPARGC1A", "CYCS", "SIRT3", "MT-CO1"],
    },
    "blood_brain_barrier_disruption": {
        "description": "血脑屏障破坏 — 紧密连接蛋白降解",
        "hallmark_genes": ["CLDN5", "OCLN", "TJP1", "MMP9", "VEGFA"],
    },
}


# ============================================================
# 8. 异构图构建函数
# ============================================================

def build_ferroptosis_hetero_graph(
    include_string_ppi: bool = False,
    string_ppi_path: Optional[str] = None,
) -> HeteroData:
    """构建石竹烯-铁衰老-脑缺血异构图。

    节点类型:
      - drug (石竹烯)
      - gene (铁死亡+衰老+CIRI交集基因)
      - disease (脑缺血再灌注)
      - pathway (信号通路)
      - phenotype (表型)

    边类型:
      - drug_targets_gene: 石竹烯 → 靶点基因
      - gene_interacts_gene: PPI网络
      - gene_belongs_to_pathway: 基因 → 通路
      - gene_associated_with_disease: 基因 → 疾病
      - pathway_related_to_phenotype: 通路 → 表型
      - gene_regulates_phenotype: 基因 → 表型

    Returns:
        HeteroData: PyG异构图数据对象
    """
    data = HeteroData()

    # --- 收集所有唯一基因 ---
    all_genes = set()
    all_genes.update(ALL_FERROPTOSIS_GENES)
    all_genes.update(ALL_AGING_GENES)
    all_genes.update(ALL_CIRI_GENES)
    all_genes.update(BCP_TARGETS)
    all_genes.update(ACSL4_FIRST_NEIGHBORS["direct_interactors"])
    all_genes.update(ACSL4_FIRST_NEIGHBORS["indirect_regulators"])
    all_genes = sorted(all_genes)
    gene_to_idx = {g: i for i, g in enumerate(all_genes)}

    # --- 药物节点 (石竹烯) ---
    drug_names = ["Beta-caryophyllene"]
    drug_to_idx = {d: i for i, d in enumerate(drug_names)}

    # --- 疾病节点 (脑缺血再灌注) ---
    disease_names = ["Cerebral_Ischemia_Reperfusion_Injury"]
    disease_to_idx = {d: i for i, d in enumerate(disease_names)}

    # --- 通路节点 ---
    pathway_names = sorted(FERROPTOSIS_PATHWAYS.keys())
    pathway_to_idx = {p: i for i, p in enumerate(pathway_names)}

    # --- 表型节点 ---
    phenotype_names = sorted(PHENOTYPES.keys())
    phenotype_to_idx = {p: i for i, p in enumerate(phenotype_names)}

    # --- 存储节点数 ---
    data['drug'].num_nodes = len(drug_names)
    data['gene'].num_nodes = len(all_genes)
    data['disease'].num_nodes = len(disease_names)
    data['pathway'].num_nodes = len(pathway_names)
    data['phenotype'].num_nodes = len(phenotype_names)

    # --- 存储名称映射 ---
    data['gene'].names = all_genes
    data['drug'].names = drug_names
    data['disease'].names = disease_names
    data['pathway'].names = pathway_names
    data['phenotype'].names = phenotype_names

    # ============================================================
    # 构建边
    # ============================================================

    # --- drug_targets_gene: 石竹烯 → 靶点 ---
    drug_target_edges = []
    for target in BCP_TARGETS:
        if target in gene_to_idx:
            drug_target_edges.append((0, gene_to_idx[target]))  # drug=0 即石竹烯
    if drug_target_edges:
        data['drug', 'targets', 'gene'].edge_index = torch.tensor(
            drug_target_edges, dtype=torch.long
        ).t().contiguous()
        # 边权重：基于文献置信度 (CNR2最高，其余中等)
        edge_weights = []
        high_conf = {"CNR2", "TRPV1", "PPARG", "NFE2L2", "GPX4", "SLC7A11"}
        for _, g_idx in drug_target_edges:
            g_name = all_genes[g_idx]
            edge_weights.append(1.0 if g_name in high_conf else 0.7)
        data['drug', 'targets', 'gene'].edge_weight = torch.tensor(
            edge_weights, dtype=torch.float
        )

    # --- gene_interacts_gene: PPI网络 ---
    if include_string_ppi and string_ppi_path and os.path.exists(string_ppi_path):
        ppi_edges = _load_string_ppi(string_ppi_path, gene_to_idx)
    else:
        ppi_edges = _build_intrinsic_ppi(gene_to_idx, all_genes)
    if ppi_edges:
        data['gene', 'interacts', 'gene'].edge_index = torch.tensor(
            ppi_edges, dtype=torch.long
        ).t().contiguous()
        # 互作边权重 (ACSL4邻居边权重更高)
        acsl4_neighbors = set(ACSL4_FIRST_NEIGHBORS["direct_interactors"] + 
                             ACSL4_FIRST_NEIGHBORS["indirect_regulators"])
        ppi_weights = []
        for s, d in ppi_edges:
            sg, dg = all_genes[s], all_genes[d]
            if sg in acsl4_neighbors and dg == "ACSL4":
                ppi_weights.append(0.95)
            elif dg in acsl4_neighbors and sg == "ACSL4":
                ppi_weights.append(0.95)
            else:
                ppi_weights.append(0.6)
        data['gene', 'interacts', 'gene'].edge_weight = torch.tensor(
            ppi_weights, dtype=torch.float
        )

    # --- gene_belongs_to_pathway: 基因 → 通路 ---
    gp_edges = []
    gp_weights = []
    for pw_name, pw_genes in FERROPTOSIS_PATHWAYS.items():
        pw_idx = pathway_to_idx[pw_name]
        for g in pw_genes:
            if g in gene_to_idx:
                gp_edges.append((gene_to_idx[g], pw_idx))
                # ACSL4在铁死亡和脂质代谢通路中权值更高
                if g == "ACSL4" and ("Ferroptosis" in pw_name or "Lipid" in pw_name):
                    gp_weights.append(1.0)
                else:
                    gp_weights.append(0.8)
    if gp_edges:
        data['gene', 'belongs_to', 'pathway'].edge_index = torch.tensor(
            gp_edges, dtype=torch.long
        ).t().contiguous()
        data['gene', 'belongs_to', 'pathway'].edge_weight = torch.tensor(
            gp_weights, dtype=torch.float
        )

    # --- gene_associated_with_disease: 基因 → CIRI ---
    gd_edges = []
    gd_weights = []
    for g_name in ALL_CIRI_GENES:
        if g_name in gene_to_idx:
            gd_edges.append((gene_to_idx[g_name], 0))  # disease=0 即CIRI
            # 缺血中上调基因权值更高
            if g_name in CIRI_DEGS["upregulated_in_ischemia"]:
                gd_weights.append(0.9)
            elif g_name in CIRI_DEGS["downregulated_in_ischemia"]:
                gd_weights.append(0.7)
            else:
                gd_weights.append(0.5)
    if gd_edges:
        data['gene', 'associated_with', 'disease'].edge_index = torch.tensor(
            gd_edges, dtype=torch.long
        ).t().contiguous()
        data['gene', 'associated_with', 'disease'].edge_weight = torch.tensor(
            gd_weights, dtype=torch.float
        )

    # --- pathway_related_to_phenotype: 通路 → 表型 ---
    pp_edges = []
    pp_weights = []
    pathway_to_phenotype = {
        "Ferroptosis (KEGG: hsa04216)": ["ferroptosis", "ferro_aging", "lipid_peroxidation"],
        "Glutathione metabolism (KEGG: hsa00480)": ["oxidative_stress", "ferroptosis"],
        "Lipid metabolism / PUFA biosynthesis": ["lipid_peroxidation", "ferroptosis"],
        "Iron homeostasis": ["ferroptosis", "ferro_aging"],
        "Autophagy & ferroptosis": ["ferroptosis", "cellular_senescence"],
        "Necroptosis / Pyroptosis (交叉)": ["neuroinflammation"],
        "Aging-related pathways": ["cellular_senescence", "ferro_aging", "mitochondrial_dysfunction"],
        "CIRI-related pathways": [
            "neuroinflammation", "oxidative_stress",
            "mitochondrial_dysfunction", "blood_brain_barrier_disruption",
        ],
    }
    for pw_name, phenotypes in pathway_to_phenotype.items():
        if pw_name in pathway_to_idx:
            for pheno in phenotypes:
                if pheno in phenotype_to_idx:
                    pp_edges.append((pathway_to_idx[pw_name], phenotype_to_idx[pheno]))
                    pp_weights.append(0.85)
    if pp_edges:
        data['pathway', 'related_to', 'phenotype'].edge_index = torch.tensor(
            pp_edges, dtype=torch.long
        ).t().contiguous()
        data['pathway', 'related_to', 'phenotype'].edge_weight = torch.tensor(
            pp_weights, dtype=torch.float
        )

    # --- gene_regulates_phenotype: 基因 → 表型 ---
    gp_edges = []
    gp_weights = []
    # 使用hallmark基因关联
    for pheno_name, pheno_info in PHENOTYPES.items():
        pheno_idx = phenotype_to_idx[pheno_name]
        for g in pheno_info["hallmark_genes"]:
            if g in gene_to_idx:
                gp_edges.append((gene_to_idx[g], pheno_idx))
                gp_weights.append(0.9)
    if gp_edges:
        data['gene', 'regulates', 'phenotype'].edge_index = torch.tensor(
            gp_edges, dtype=torch.long
        ).t().contiguous()
        data['gene', 'regulates', 'phenotype'].edge_weight = torch.tensor(
            gp_weights, dtype=torch.float
        )

    # --- 存储元数据 ---
    data.metadata = {
        "node_types": ["drug", "gene", "disease", "pathway", "phenotype"],
        "edge_types": [
            ("drug", "targets", "gene"),
            ("gene", "interacts", "gene"),
            ("gene", "belongs_to", "pathway"),
            ("gene", "associated_with", "disease"),
            ("pathway", "related_to", "phenotype"),
            ("gene", "regulates", "phenotype"),
        ],
        "num_ferroptosis_genes": len(set(ALL_FERROPTOSIS_GENES) & set(all_genes)),
        "num_aging_genes": len(set(ALL_AGING_GENES) & set(all_genes)),
        "num_ciri_genes": len(set(ALL_CIRI_GENES) & set(all_genes)),
        "num_bcp_targets": len(set(BCP_TARGETS) & set(all_genes)),
        "acsl4_has_direct_neighbors": "ACSL4" in all_genes,
    }

    return data


def _build_intrinsic_ppi(
    gene_to_idx: Dict[str, int],
    all_genes: List[str],
) -> List[Tuple[int, int]]:
    """基于ACSL4中心网络 + 铁死亡通路共现构建内在PPI。

    在没有STRING外部数据时使用。
    """
    edges = []

    # ACSL4 一阶邻居
    acsl4_idx = gene_to_idx.get("ACSL4")
    if acsl4_idx is not None:
        for neighbor in (ACSL4_FIRST_NEIGHBORS["direct_interactors"] +
                         ACSL4_FIRST_NEIGHBORS["indirect_regulators"]):
            if neighbor in gene_to_idx:
                edges.append((acsl4_idx, gene_to_idx[neighbor]))

    # GPX4-SLC7A11 核心抑制轴
    if "GPX4" in gene_to_idx and "SLC7A11" in gene_to_idx:
        edges.append((gene_to_idx["GPX4"], gene_to_idx["SLC7A11"]))

    # FSP1 独立通路
    if "AIFM2" in gene_to_idx and "GPX4" in gene_to_idx:
        edges.append((gene_to_idx["AIFM2"], gene_to_idx["GPX4"]))

    # 脂质过氧化级联
    lipid_genes = ["ACSL4", "LPCAT3", "ALOX5", "ALOX15", "ALOX12", "PLA2G4A"]
    for i, g1 in enumerate(lipid_genes):
        for g2 in lipid_genes[i+1:]:
            if g1 in gene_to_idx and g2 in gene_to_idx:
                edges.append((gene_to_idx[g1], gene_to_idx[g2]))

    # 铁代谢核心
    iron_genes = ["TFRC", "FTH1", "FTL", "HMOX1", "SLC40A1", "IREB2"]
    for i, g1 in enumerate(iron_genes):
        for g2 in iron_genes[i+1:]:
            if g1 in gene_to_idx and g2 in gene_to_idx:
                edges.append((gene_to_idx[g1], gene_to_idx[g2]))

    # NRF2-ARE通路
    nrf2_genes = ["NFE2L2", "KEAP1", "GCLC", "GCLM", "HMOX1", "NQO1", "TXNRD1"]
    for i, g1 in enumerate(nrf2_genes):
        for g2 in nrf2_genes[i+1:]:
            if g1 in gene_to_idx and g2 in gene_to_idx:
                edges.append((gene_to_idx[g1], gene_to_idx[g2]))

    # 炎症核心
    inflam_genes = ["NFKB1", "RELA", "IL6", "TNF", "IL1B", "PTGS2"]
    for i, g1 in enumerate(inflam_genes):
        for g2 in inflam_genes[i+1:]:
            if g1 in gene_to_idx and g2 in gene_to_idx:
                edges.append((gene_to_idx[g1], gene_to_idx[g2]))

    # p53-p21衰老轴
    if "TP53" in gene_to_idx:
        for g in ["CDKN1A", "MDM2", "BAX"]:
            if g in gene_to_idx:
                edges.append((gene_to_idx["TP53"], gene_to_idx[g]))

    # 去重
    return list(set(edges))


def _load_string_ppi(
    ppi_path: str,
    gene_to_idx: Dict[str, int],
    min_score: int = 700,
) -> List[Tuple[int, int]]:
    """从STRING PPI文件加载并过滤互作关系。"""
    edges = []
    with open(ppi_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            p1, p2, score = parts[0], parts[1], int(parts[2])
            if score < min_score:
                continue
            if p1 in gene_to_idx and p2 in gene_to_idx:
                edges.append((gene_to_idx[p1], gene_to_idx[p2]))
    return edges


# ============================================================
# 9. 图统计分析
# ============================================================

def analyze_graph(data: HeteroData) -> Dict:
    """分析异构图的基本统计指标。"""
    stats = {}
    node_types, edge_types = data.metadata()
    stats["num_nodes"] = {nt: data[nt].num_nodes for nt in node_types}
    stats["num_edges"] = {}
    for et in edge_types:
        try:
            stats["num_edges"][str(et)] = data[et].edge_index.size(1)
        except Exception:
            stats["num_edges"][str(et)] = 0
    stats["total_edges"] = sum(stats["num_edges"].values())
    stats["total_nodes"] = sum(stats["num_nodes"].values())

    # ACSL4中心性分析
    acsl4_idx = data['gene'].names.index("ACSL4") if "ACSL4" in data['gene'].names else None
    if acsl4_idx is not None:
        acsl4_degree = 0
        for et in edge_types:
            try:
                ei = data[et].edge_index
                if et[0] == "gene" and et[2] == "gene":
                    acsl4_degree += (ei[0] == acsl4_idx).sum().item()
                    acsl4_degree += (ei[1] == acsl4_idx).sum().item()
                elif et[0] == "gene":
                    acsl4_degree += (ei[0] == acsl4_idx).sum().item()
                elif et[2] == "gene":
                    acsl4_degree += (ei[1] == acsl4_idx).sum().item()
            except Exception:
                pass
        stats["acsl4_degree"] = acsl4_degree

    # 交集统计
    ferro_set = set(ALL_FERROPTOSIS_GENES) & set(data['gene'].names)
    aging_set = set(ALL_AGING_GENES) & set(data['gene'].names)
    ciri_set = set(ALL_CIRI_GENES) & set(data['gene'].names)
    bcp_set = set(BCP_TARGETS) & set(data['gene'].names)

    stats["intersections"] = {
        "ferroptosis_only": len(ferro_set - aging_set - ciri_set),
        "ferro_aging_overlap": len(ferro_set & aging_set - ciri_set),
        "ferro_ciri_overlap": len(ferro_set & ciri_set - aging_set),
        "triple_overlap": len(ferro_set & aging_set & ciri_set),
        "bcp_ferro_overlap": len(bcp_set & ferro_set),
        "bcp_acsl4_neighbor_overlap": len(
            bcp_set & set(ACSL4_FIRST_NEIGHBORS["direct_interactors"] +
                         ACSL4_FIRST_NEIGHBORS["indirect_regulators"])
        ),
    }

    return stats


# ============================================================
# 10. 主函数
# ============================================================

if __name__ == "__main__":
    print("=" * 70)
    print("阶段1: 构建石竹烯-铁衰老-脑缺血异构图")
    print("=" * 70)

    # 构建异构图
    data = build_ferroptosis_hetero_graph(include_string_ppi=False)

    # 分析统计
    stats = analyze_graph(data)
    print(f"\n图结构统计:")
    print(f"  总节点数: {stats['total_nodes']}")
    print(f"  总边数: {stats['total_edges']}")
    print(f"\n  节点分布:")
    for nt, n in stats['num_nodes'].items():
        print(f"    {nt}: {n}")

    print(f"\n  边分布:")
    for et, n in stats['num_edges'].items():
        print(f"    {et[0]} → {et[2]} ({et[1]}): {n}")

    if 'acsl4_degree' in stats:
        print(f"\n  ACSL4 节点度数: {stats['acsl4_degree']}")

    print(f"\n  基因集交集:")
    for k, v in stats['intersections'].items():
        print(f"    {k}: {v}")

    # 保存图
    save_path = DATA_DIR / "ferroptosis_hetero_graph.pt"
    torch.save(data, save_path)
    print(f"\n✓ 异构图已保存至: {save_path}")

    # 保存统计JSON
    stats_path = DATA_DIR / "graph_statistics.json"
    with open(stats_path, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2, ensure_ascii=False, default=str)
    print(f"✓ 图统计已保存至: {stats_path}")

    # 导出节点列表
    gene_list_path = DATA_DIR / "gene_nodes.csv"
    pd.DataFrame({
        "gene": data['gene'].names,
        "in_ferroptosis": [g in ALL_FERROPTOSIS_GENES for g in data['gene'].names],
        "in_aging": [g in ALL_AGING_GENES for g in data['gene'].names],
        "in_ciri": [g in ALL_CIRI_GENES for g in data['gene'].names],
        "in_bcp_targets": [g in BCP_TARGETS for g in data['gene'].names],
        "is_acsl4_neighbor": [
            g in ACSL4_FIRST_NEIGHBORS["direct_interactors"] +
            ACSL4_FIRST_NEIGHBORS["indirect_regulators"]
            for g in data['gene'].names
        ],
    }).to_csv(gene_list_path, index=False)
    print(f"✓ 基因节点列表已保存至: {gene_list_path}")

    print("\n" + "=" * 70)
    print("阶段1 完成")
    print("=" * 70)