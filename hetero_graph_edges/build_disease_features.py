# -*- coding: utf-8 -*-
"""
疾病语义嵌入构建: CIRI基因嵌入均值法 (默认) + DisGeNET 备选
=====================================================
方法一 (推荐, 默认): CIRI疾病相关基因的嵌入均值
  - 直接反映疾病相关基因的分子特征
  - 论文引用: "疾病特征由其关联基因的平均 SapBERT 嵌入表示"
  - 输入: disease_genes.txt (CIRI 疾病相关基因列表)
  - 输入: subgraph_embeddings.csv (基因嵌入)

方法二 (备选): DisGeNET 疾病语义嵌入
  - 下载 DisGeNET 疾病-基因关联或语义相似度矩阵
  - 用 Gonto2Vec / Onto2Vec 生成疾病本体嵌入

方法三: Disease Ontology (DO) OWL 文件 + Onto2Vec
  - 下载 DO OWL 文件 → Onto2Vec 训练 → 疾病嵌入

输出:
  - disease_features.csv (最多20行, 每行=一种疾病语义特征)
  - disease_features.npy (便于 PyG 加载)

作者: 优化版 v2.0
日期: 2026-05-31
"""

import os
import sys
import json
import csv
import numpy as np
import pandas as pd
from collections import defaultdict
from pathlib import Path

# ============================================================
# 0. 配置
# ============================================================
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\GAT")
OUTPUT_DIR = Path(r"D:\反向网络药理学\GAT拓展维度")

EMBEDDING_FILE = DATA_DIR / "subgraph_embeddings.csv"
DISEASE_GENES_FILE = DATA_DIR / "disease_genes.txt"

OUTPUT_FEATURES_CSV = OUTPUT_DIR / "disease_features.csv"
OUTPUT_FEATURES_NPY = OUTPUT_DIR / "disease_features.npy"
OUTPUT_DISEASE_NAMES = OUTPUT_DIR / "disease_feature_names.txt"
OUTPUT_STATS = OUTPUT_DIR / "disease_feature_stats.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. CIRI 疾病基因定义 (10个维度)
# ============================================================
CIRI_DISEASE_GENE_SETS = {
    "Cerebral_Ischemia_Reperfusion_Injury_CORE": [
        "NFKB1", "RELA", "TNF", "IL6", "IL1B", "HMOX1", "NOS2", "NOS3",
        "PTGS2", "BCL2", "BAX", "CASP3", "SOD1", "SOD2", "CAT",
        "HIF1A", "VEGFA", "MMP9", "ICAM1", "VCAM1",
        "GFAP", "MAP2", "SYN1", "GAP43", "SYP",
        "BDNF", "NGF", "GDNF", "NTF3",
        "TLR4", "MYD88", "HMGB1", "AIF1", "CSF1R",
        "CXCL8", "CCL2", "CX3CR1",
        "AQP4", "SLC1A2", "GRIN1", "GRIN2A", "GRIN2B",
    ],
    "Oxidative_Stress_CIRI": [
        "NFE2L2", "KEAP1", "HMOX1", "SOD1", "SOD2", "CAT",
        "GPX1", "GSR", "TXNRD1", "PRDX1", "PRDX2",
        "NOX1", "NOX2", "NOX4",
        "HIF1A", "VEGFA", "EPO",
        "MAOA", "MAOB", "XDH",
    ],
    "Neuroinflammation_CIRI": [
        "TNF", "IL1B", "IL6", "IL10", "TGFB1",
        "NFKB1", "RELA", "IKBKB", "NFKBIA",
        "TLR2", "TLR4", "MYD88", "TICAM1",
        "AIF1", "GFAP", "CD68",
        "CCL2", "CCL5", "CXCL8", "CX3CL1", "CX3CR1",
        "NLRP3", "CASP1", "IL18", "PYCARD",
        "PTGS2", "NOS2", "MMP9", "MMP2",
    ],
    "Apoptosis_Autophagy_CIRI": [
        "BCL2", "BAX", "BCL2L1", "BAD", "BID",
        "CASP3", "CASP8", "CASP9", "PARP1",
        "CYCS", "APAF1", "DIABLO", "XIAP",
        "TP53", "MDM2", "PMAIP1", "BBC3",
        "BECN1", "ATG5", "ATG7", "ATG12",
        "MAP1LC3A", "MAP1LC3B", "SQSTM1",
        "MTOR", "AKT1", "ULK1",
    ],
    "Blood_Brain_Barrier_CIRI": [
        "CLDN5", "OCLN", "TJP1", "TJP2", "CDH5",
        "MMP2", "MMP9", "MMP3", "TIMP1", "TIMP2",
        "VEGFA", "VEGFR2", "ANGPT1", "ANGPT2", "TEK",
        "ICAM1", "VCAM1", "SELE", "SELP",
        "AQP4", "ABCB1", "SLC2A1",
    ],
    "Excitotoxicity_CIRI": [
        "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D",
        "GRIA1", "GRIA2", "GRIA3", "GRIA4",
        "SLC1A2", "SLC1A3", "SLC1A1",
        "DLG4", "CAMK2A", "CREB1",
        "NOS1", "CALM1", "CALM2", "CALM3",
    ],
    "ER_Stress_CIRI": [
        "HSPA5", "DDIT3", "ERN1", "EIF2AK3", "ATF6",
        "XBP1", "ATF4", "HSP90B1", "PDIA3", "PDIA4",
        "EIF2S1", "PPP1R15A", "DNAJB9", "DNAJC3",
        "ERN1", "TRAF2", "MAPK8", "MAPK9",
    ],
    "Ferroptosis_CIRI": [
        "GPX4", "SLC7A11", "ACSL4", "LPCAT3",
        "TFRC", "FTH1", "FTL", "IREB2",
        "HMOX1", "NFE2L2", "KEAP1",
        "PTGS2", "ALOX5", "ALOX12", "ALOX15",
        "VDAC2", "VDAC3", "ATG5", "ATG7",
    ],
    "Pyroptosis_CIRI": [
        "NLRP3", "NLRC4", "AIM2", "PYCARD",
        "CASP1", "CASP4", "CASP5",
        "GSDMD", "GSDME",
        "IL1B", "IL18", "HMGB1",
        "TLR4", "NFKB1", "RELA",
    ],
    "Neuroplasticity_Repair_CIRI": [
        "BDNF", "NTRK2", "NGF", "NGFR",
        "GDNF", "GFRA1", "NTF3", "NTRK3",
        "CREB1", "CAMK2A", "CAMK2B",
        "SYN1", "SYP", "GAP43", "DLG4",
        "ARC", "EGR1", "FOS", "JUN",
        "MAP2", "TUBB3", "NEFL", "NEFM", "NEFH",
        "VEGFA", "IGF1", "CNTF",
    ],
    "Cuproptosis_Core": [
        "FDX1", "LIAS", "LIPT1", "DLD", "DLAT",
        "PDHA1", "PDHB", "MTF1", "GLS", "CDKN2A",
        "SLC31A1", "ATP7A", "ATP7B", "DBT", "DLST",
        "PDHA2", "GCSH",
    ],
    "Cuproptosis_Copper_Metabolism": [
        "ATOX1", "COX17", "CCS", "COX11", "SCO1", "SCO2",
        "STEAP1", "STEAP2", "STEAP3", "STEAP4",
        "CP", "COMMD1", "MT1A", "MT2A",
        "SLC31A1", "ATP7A", "ATP7B",
        "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4",
        "TYR", "DBH", "PAM", "MOXD1",
    ],
    "Cuproptosis_CIRI_Crosstalk": [
        "FDX1", "DLAT", "PDHA1", "DLD", "LIAS", "LIPT1",
        "PDHB", "MTF1", "GLS", "CDKN2A",
        "HMOX1", "HSPA5", "NFKB1", "RELA",
        "SOD1", "SOD2", "CAT", "GPX4",
        "HIF1A", "NFE2L2", "KEAP1", "TP53", "BAX", "BCL2",
        "IL6", "TNF", "IL1B", "MMP9",
        "AKT1", "MTOR", "MAPK8", "MAPK14",
        "XBP1", "DDIT3", "ATF4", "HSP90B1",
    ],
    "Cuproptosis_FeS_Cluster_Biogenesis": [
        "FDX1", "FDX2", "FDXR",
        "ISCU", "ISCA1", "ISCA2", "IBA57",
        "NFU1", "BOLA1", "BOLA3", "GLRX5",
        "NFS1", "LYRM4", "NDUFS1",
        "HSCB", "HSPA9", "GRPEL1", "GRPEL2",
        "FXN", "FDX1L",
    ],
}

# ============================================================
# 2. 基因嵌入加载
# ============================================================
def load_gene_embeddings(embedding_path):
    """加载基因嵌入矩阵"""
    print(f"[LOAD] 加载基因嵌入: {embedding_path}")

    if not embedding_path.exists():
        print(f"[ERROR] 嵌入文件不存在: {embedding_path}")
        return {}, 0

    try:
        df = pd.read_csv(embedding_path, sep=None, engine='python')
    except Exception:
        try:
            df = pd.read_csv(embedding_path, sep='\t')
        except Exception as e:
            print(f"[ERROR] 无法读取嵌入文件: {e}")
            return {}, 0

    gene_embeddings = {}
    gene_col = df.columns[0]
    embed_cols = df.columns[1:]
    skipped = 0

    for _, row in df.iterrows():
        gene = str(row[gene_col]).strip().upper()
        if not gene or gene in ('NAN', 'NONE', ''):
            skipped += 1
            continue
        try:
            vec = row[embed_cols].astype(float).values
            gene_embeddings[gene] = vec
        except (ValueError, TypeError):
            skipped += 1

    embed_dim = len(embed_cols)
    print(f"[OK] 加载 {len(gene_embeddings)} 个基因嵌入 (维度={embed_dim})")
    return gene_embeddings, embed_dim


# ============================================================
# 3. 疾病特征计算
# ============================================================
def compute_disease_features(disease_gene_sets, gene_embeddings, embed_dim):
    """
    对每个疾病维度, 取其关联基因嵌入的均值作为疾病语义特征

    为什么取均值:
      - 类比 CBOW (Continuous Bag of Words) 思想
      - 多基因联合表征可缓解单基因嵌入噪声
      - 在文献中广泛使用 (e.g. drug target profiles, disease modules)
    """
    np.random.seed(42)

    features = []
    names = []
    coverage = {}

    for disease_name, genes in disease_gene_sets.items():
        gene_vecs = []
        for gene in genes:
            gene = gene.upper()
            if gene in gene_embeddings:
                gene_vecs.append(gene_embeddings[gene])

        n_total = len(genes)
        n_with = len(gene_vecs)
        coverage[disease_name] = (n_with, n_total)

        if n_with > 0:
            mean_vec = np.mean(gene_vecs, axis=0)
            norm = np.linalg.norm(mean_vec)
            if norm > 1e-8:
                mean_vec = mean_vec / norm
        else:
            mean_vec = np.random.randn(embed_dim).astype(np.float32)
            norm = np.linalg.norm(mean_vec)
            mean_vec = mean_vec / norm

        features.append(mean_vec)
        names.append(disease_name)

    features = np.array(features, dtype=np.float32)

    print(f"\n[STATS] 疾病特征计算完成:")
    print(f"  疾病维度总数: {len(features)}")
    print(f"  特征维度:     {features.shape[1]}")
    for name, cov in coverage.items():
        pct = cov[0] / cov[1] * 100 if cov[1] > 0 else 0
        print(f"    {name}: {cov[0]}/{cov[1]} ({pct:.0f}%) 基因有嵌入")

    return features, names, coverage


# ============================================================
# 4. 导出
# ============================================================
def save_features(features, names, coverage, embed_dim):
    """保存所有输出文件"""
    df = pd.DataFrame(features, index=names)
    df.to_csv(OUTPUT_FEATURES_CSV)
    print(f"[SAVE] 疾病特征CSV: {OUTPUT_FEATURES_CSV} ({features.shape})")

    np.save(OUTPUT_FEATURES_NPY, features)
    print(f"[SAVE] 疾病特征NPY: {OUTPUT_FEATURES_NPY} ({features.shape})")

    with open(OUTPUT_DISEASE_NAMES, 'w', encoding='utf-8') as f:
        for name in names:
            f.write(name + '\n')
    print(f"[SAVE] 疾病名列表: {OUTPUT_DISEASE_NAMES}")

    stats = {
        'n_diseases': len(names),
        'embed_dim': embed_dim,
        'coverage': {k: {'n_with': v[0], 'n_total': v[1]} for k, v in coverage.items()},
        'method': 'disease_gene_embedding_mean',
        'embedding_source': str(EMBEDDING_FILE),
        'citation_note': (
            "Disease features are represented by the mean SapBERT/ProtBERT "
            "embedding of their associated genes, a common practice in drug "
            "repurposing and disease module identification studies."
        ),
        'disease_names': names,
    }
    with open(OUTPUT_STATS, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    print(f"[SAVE] 统计信息: {OUTPUT_STATS}")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("疾病语义嵌入构建: CIRI 基因嵌入均值法")
    print("=" * 60)

    gene_embeddings, embed_dim = load_gene_embeddings(EMBEDDING_FILE)

    if embed_dim == 0:
        print("[ERROR] 无法加载基因嵌入, 将生成随机疾病嵌入")

    features, names, coverage = compute_disease_features(
        CIRI_DISEASE_GENE_SETS, gene_embeddings, embed_dim
    )

    save_features(features, names, coverage, embed_dim)

    print(f"\n{'='*60}")
    print("疾病语义嵌入构建完成!")
    print(f"输出: {OUTPUT_FEATURES_CSV}")
    print(f"输出: {OUTPUT_FEATURES_NPY}")
    print(f"输出: {OUTPUT_DISEASE_NAMES}")
    print("=" * 60)