#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=====================================================================
IDSP (Iron-Driven Senescence Program) 三基因集构建
=====================================================================
基于权威数据库:
  - FerrDb V2 (Nucleic Acids Res 2022/2025) — 铁死亡基因
  - CellAge (Genome Biology 2020) — 衰老基因
  - SenMayo (Nature Communications 2022) — 体内衰老标志物
  - Liu et al. 2026 Cell Metabolism — 铁衰老基因集 (96基因)

设计原则:
  1. PURE_FERROPTOSIS: 只在铁死亡中富集，不与衰老共享的基因
  2. PURE_SENESCENCE: 只在衰老中富集，不与铁死亡共享的基因
  3. SHARED: 同时在铁死亡和衰老中起作用的基因 (含桥接基因)
  4. 每个PURE集 ≥ 50基因，确保GSVA评分统计学稳健
=====================================================================
"""

# ============================================================
# FerrDb V2 实验验证的铁死亡核心基因 (经实验证实的驱动/抑制基因)
# 来源: https://www.zhounan.org/ferrdb/current/
# 仅保留 "Validated" 级别基因
# ============================================================
FERROPTOSIS_VALIDATED = {
    # === System Xc- / GSH / GPX4 抗氧化轴 (9) ===
    "SLC7A11", "SLC3A2", "GPX4",
    "GCLC", "GCLM", "GSS",
    "CHAC1", "TXNRD1", "TXN",

    # === FSP1-CoQH2 / DHODH / GCH1 非GPX4防御 (3) ===
    "FSP1", "DHODH", "GCH1",

    # === 脂质代谢 & 脂质过氧化 (12) ===
    "ACSL4", "ACSL3", "LPCAT3", "LPCAT4",
    "ALOX5", "ALOX12", "ALOX15", "ALOX15B", "ALOXE3",
    "POR", "CYB5R1",
    "FADS1", "FADS2", "ELOVL5",

    # === 铁代谢 & 铁稳态 (12) ===
    "TFRC", "STEAP3", "SLC11A2",   # 铁摄取
    "FTH1", "FTL",                   # 铁储存
    "SLC40A1", "CP",                 # 铁外排
    "NCOA4",                         # 铁蛋白自噬
    "CISD1", "CISD2",                # 铁硫蛋白
    "IREB2",                         # 铁调控
    "PCBP1", "PCBP2",                # 铁伴侣
    "ABCB6", "ABCB7", "HEPH",

    # === 核心信号通路 & 转录调控 (15) ===
    "NFE2L2", "KEAP1",               # NRF2通路
    "TP53",                          # p53
    "NF2", "YAP1", "WWTR1",          # Hippo通路
    "MTOR",                          # mTOR
    "PRKAA1", "PRKAA2",             # AMPK
    "CDKN1A",                        # p21
    "RB1",                           # Rb
    "EIF2AK4",                       # GCN2
    "ATF3", "ATF4",                  # ATF

    # === 其他关键调控因子 (14) ===
    "HSPB1", "HSPA5", "HSP90AA1",   # 热休克蛋白
    "CAV1",                          # 小窝蛋白
    "RPL8", "ACSF2", "MYB",         # 其他驱动
    "SLC38A1",                       # 氨基酸转运
    "VDAC2", "VDAC3",               # 线粒体通道
    "MAP1LC3A", "GABARAPL1",        # 自噬
    "RRM2", "SQSTM1",               # 其他

    # === 铁死亡标志基因 (Marker) (5) ===
    "PTGS2", "HMOX1", "SOD1", "CHAC1", "HSPB1",
}

# 也包含经高通量筛选(Screened)但被广泛引用的基因
FERROPTOSIS_SCREENED = {
    "SELENOS", "DHFR", "ALOXE3",
    "ABCB6", "ABCB7", "HEPH", "CP",
    "FADS1", "FADS2", "ELOVL5",
    "HMOX2",
}

FERROPTOSIS_ALL = FERROPTOSIS_VALIDATED | FERROPTOSIS_SCREENED

# ============================================================
# 衰老基因集 (CellAge + SenMayo)
# 来源: 
#   CellAge: https://genomics.senescence.info/cells/ (279 genes)
#   SenMayo: Saul et al. Nat Commun 2022 (125 genes)
# ============================================================
SENESCENCE_CORE = {
    # === 核心细胞周期调控/衰老驱动 (10) ===
    "TP53",                          # p53, 核心衰老转录因子
    "CDKN1A",                        # p21, 经典衰老标志
    "CDKN2A",                        # p16, 体内衰老标志
    "CDKN1B",                        # p27
    "CDKN2B",                        # p15
    "RB1",                           # Rb
    "E2F1",                          # 增殖/衰老转录因子
    "MYC",                           # Myc
    "LMNB1",                         # 衰老时下降(反向标志)
    "GLB1",                          # SA-beta-Gal

    # === SASP: 细胞因子 (12) ===
    "IL6", "CXCL8", "IL1A", "IL1B", "TNF",
    "IL10", "IL13", "IL15", "IL18", "IL32",
    "CSF1", "CSF2",

    # === SASP: CCL趋化因子 (12) ===
    "CCL2", "CCL3", "CCL4", "CCL5",
    "CCL7", "CCL8", "CCL13", "CCL20",
    "CCL24", "CCL26", "CCL1", "CCL16",

    # === SASP: CXCL趋化因子 (9) ===
    "CXCL1", "CXCL2", "CXCL3",
    "CXCL5", "CXCL10", "CXCL12",
    "CXCL14", "CXCL16", "CX3CL1",

    # === SASP: 生长因子/调节因子 (12) ===
    "VEGFA", "HGF", "TGFB1",
    "EGF", "FGF1", "FGF2",
    "AREG", "EREG", "GDF15",
    "IGFBP1", "IGFBP3", "IGFBP7",
    "DKK1", "KITLG",

    # === SASP: MMP/蛋白酶 (9) ===
    "MMP1", "MMP2", "MMP3",
    "MMP9", "MMP10", "MMP12",
    "MMP13", "MMP14",
    "SERPINE1", "TIMP2",

    # === 表面标志/受体 (9) ===
    "ICAM1", "CXCR2",
    "FAS", "TNFRSF1A", "TNFRSF1B",
    "AXL", "CD9", "CD55",
    "NOTCH1", "NOTCH3", "JUN",

    # === 其他SASP/衰老相关 (8) ===
    "HMGB1", "SPP1", "ETS2", "CTSB",
    "CTNNB1", "PAPPA", "MIF", "PLAU",
}

# senmayo额外基因 (不在Core中但SenMayo收录)
SENESCENCE_SENMAYO = {
    "ACVR1B", "ADAMTS1", "ADAMTS4",
    "BMP2", "BMP6", "CALCA",
    "CCL11", "CCL13", "CCL22",
    "CTSK", "CXCL13", "DLL1",
    "EFNB1", "EFNB2", "ENG",
    "EPHA2", "EPHA3", "FGF7",
    "GEM", "ICAM2", "ICAM3",
    "IGFBP2", "IGFBP4", "IGFBP5", "IGFBP6",
    "IL11", "IL17A", "IL1R1",
    "INHBA", "ITGA1", "ITGA2",
    "ITGB3", "LAMA1", "LIF",
    "MMP11", "MMP14", "MMP16",
    "NRP1", "NRP2", "PDGFB",
    "PLAT", "PLAUR", "PTGFRN",
    "SCAMP1", "SCAMP4", "SEMA3F",
    "SERPINE2", "SERPINF1",
    "TIMP1", "TIMP3", "TNC",
    "TNFRSF11B", "TNFRSF12A",
    "VCAM1", "VEGFC", "WNT2",
    "WNT5A", "WNT7A",
}

SENESCENCE_ALL = SENESCENCE_CORE | SENESCENCE_SENMAYO

# ============================================================
# 交集分析：找出同时在铁死亡和衰老中起作用的共享基因
# ============================================================
SHARED_OVERLAP = FERROPTOSIS_ALL & SENESCENCE_ALL

print("=" * 60)
print(f"FerrDb铁死亡基因总数: {len(FERROPTOSIS_ALL)}")
print(f"  其中Validated: {len(FERROPTOSIS_VALIDATED)}")
print(f"衰老基因总数 (CellAge+SenMayo): {len(SENESCENCE_ALL)}")
print(f"  其中Core: {len(SENESCENCE_CORE)}")
print(f"  其中SenMayo追加: {len(SENESCENCE_SENMAYO)}")
print(f"铁死亡∩衰老交集: {len(SHARED_OVERLAP)} 基因")
print(f"交集基因: {sorted(SHARED_OVERLAP)}")
print("=" * 60)

# ============================================================
# 三基因集最终定义
# ============================================================

# --- 纯铁死亡基因 (FerrDb - 交集) ---
PURE_FERROPTOSIS = FERROPTOSIS_ALL - SHARED_OVERLAP

# --- 纯衰老基因 (Senescence - 交集) ---
PURE_SENESCENCE = SENESCENCE_ALL - SHARED_OVERLAP

# --- 共享基因 (交集 + 桥接基因) ---
# 桥接基因: 原论文铁衰老核心基因中连接铁代谢和衰老通路的
BRIDGE_GENES_SHARED = {
    # 原论文铁衰老核心基因 (Liu 2026), 不在FerrDb/CellAge中
    # 但在脑I/R中同时参与铁信号和衰老信号
    "CD74",    # 巨噬细胞迁移抑制因子受体, 铁-炎症桥接
    "S100A8",  # 钙卫蛋白, 铁螯合/炎症
    "IFNG",    # IFN-γ, 铁代谢调控/衰老相关
    "IRF1",    # 干扰素调控因子, 铁-免疫桥接
    "TLR4",    # TLR4, 铁/衰老相关炎症核心
    "NLRP3",   # 炎症小体, 铁过载/衰老共同激活
    "HIF1A",   # HIF-1α, 铁代谢/衰老氧感应
    "KEAP1",   # 既是铁死亡驱动, 也是衰老氧化应激感应
    "SOD1",    # SOD1, 氧化应激桥接
}

SHARED_GENES = SHARED_OVERLAP | BRIDGE_GENES_SHARED

# ============================================================
# 统计检验
# ============================================================
assert len(PURE_FERROPTOSIS) > 50, f"纯铁死亡基因集太小: {len(PURE_FERROPTOSIS)}"
assert len(PURE_SENESCENCE) > 50, f"纯衰老基因集太小: {len(PURE_SENESCENCE)}"
assert len(SHARED_GENES) > 10, f"共享基因集太小: {len(SHARED_GENES)}"
assert PURE_FERROPTOSIS.isdisjoint(PURE_SENESCENCE), "PURE_FERROPTOSIS和PURE_SENESCENCE不能重叠！"

print(f"\n{'='*60}")
print(f"最终三基因集统计:")
print(f"  PURE_FERROPTOSIS: {len(PURE_FERROPTOSIS)} 基因")
print(f"  PURE_SENESCENCE:  {len(PURE_SENESCENCE)} 基因")
print(f"  SHARED_GENES:     {len(SHARED_GENES)} 基因")
print(f"  总计:             {len(PURE_FERROPTOSIS|PURE_SENESCENCE|SHARED_GENES)} 基因")
print(f"{'='*60}")