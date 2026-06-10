# -*- coding: utf-8 -*-
"""
基因特征增强脚本 v1.0
在不改变 GAT 架构下，追加额外特征到 SapBERT 嵌入后

特征类别:
  1. PPI 网络拓扑特征 (度、中介中心性、聚类系数、接近中心性)
  2. GTEx eQTL 特征 (Brain Cortex + Whole Blood eGenes)
  3. 蛋白质序列理化特征 (UniProt API + BioPython)

输出: enhanced_gene_features.csv / .npy
"""

import os
import sys
import csv
import json
import time
import pickle
import gzip
import io
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from collections import defaultdict, Counter
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx

# ============================================================
# 0. 配置
# ============================================================
GAT_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\GAT")
EQTL_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL")
OUTPUT_DIR = Path(r"D:\反向网络药理学\GAT拓展维度")
CACHE_DIR = OUTPUT_DIR / "cache"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

EMBEDDING_FILE = GAT_DIR / "subgraph_embeddings.csv"
PPI_FILE = GAT_DIR / "ppi_subgraph.csv"
GENES_FILE = GAT_DIR / "subgraph_genes.txt"
DRUG_TARGETS_FILE = GAT_DIR / "drug_targets.txt"
DISEASE_GENES_FILE = GAT_DIR / "disease_genes.txt"

BRAIN_EGENES = EQTL_DIR / "Brain_Cortex.v11.eGenes.txt"
BLOOD_EGENES = EQTL_DIR / "Whole_Blood.v11.eGenes.txt"

OUTPUT_CSV = OUTPUT_DIR / "enhanced_gene_features.csv"
OUTPUT_NPY = OUTPUT_DIR / "enhanced_gene_features.npy"
OUTPUT_DIM_MAP = OUTPUT_DIR / "feature_dimensions.json"

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)


# ============================================================
# 1. 加载现有数据
# ============================================================
def load_gene_list():
    """加载基因列表"""
    genes = []
    with open(GENES_FILE, 'r', encoding='utf-8') as f:
        next(f)
        for line in f:
            g = line.strip().upper()
            if g:
                genes.append(g)
    print(f"[LOAD] 基因列表: {len(genes)} 个基因")
    return genes


def load_sapbert_embeddings(genes):
    """加载 SapBERT 嵌入矩阵"""
    print(f"[LOAD] 加载 SapBERT 嵌入: {EMBEDDING_FILE}")
    df = pd.read_csv(EMBEDDING_FILE)
    gene_to_idx = {g.upper(): i for i, g in enumerate(genes)}

    n_genes = len(genes)
    embed_dim = df.shape[1] - 1

    emb_matrix = np.zeros((n_genes, embed_dim), dtype=np.float32)
    loaded = 0

    for _, row in df.iterrows():
        symbol = str(row.iloc[0]).strip().upper()
        if symbol in gene_to_idx:
            idx = gene_to_idx[symbol]
            emb_matrix[idx] = row.iloc[1:].values.astype(np.float32)
            loaded += 1

    print(f"[OK] 加载 {loaded}/{n_genes} 个基因嵌入 (维度={embed_dim})")
    return emb_matrix, embed_dim


def load_disease_and_drug_genes():
    """加载疾病基因和药物靶点集合"""
    disease = set()
    with open(DISEASE_GENES_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            g = line.strip().upper()
            if g:
                disease.add(g)

    drug = set()
    with open(DRUG_TARGETS_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            g = line.strip().upper()
            if g:
                drug.add(g)

    print(f"[LOAD] 疾病基因: {len(disease)}, 药物靶点: {len(drug)}")
    return disease, drug


# ============================================================
# 2. PPI 网络拓扑特征
# ============================================================
def compute_ppi_topology_features(genes):
    """
    从 PPI 图计算拓扑特征:
      - degree: 节点度
      - betweenness_approx: 中介中心性 (采样近似)
      - clustering_coeff: 聚类系数
      - closeness_approx: 接近中心性 (采样近似)
      - eigenvector_approx: 特征向量中心性 (子图近似)

    使用 NetworkX 分块计算以处理大图
    """
    cache_file = CACHE_DIR / "ppi_topology_features.npy"
    if cache_file.exists():
        cached = np.load(cache_file)
        if not np.allclose(cached, 0):
            print(f"[CACHE] 加载缓存的 PPI 拓扑特征")
            return cached
        else:
            print(f"[CACHE] 缓存无效(全零), 重新计算 PPI 拓扑特征")

    print(f"[PPI] 构建 PPI 图...")
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    gene_set = set(genes)

    G = nx.Graph()
    G.add_nodes_from(genes)

    edges_added = 0
    with open(PPI_FILE, 'r', encoding='utf-8') as f:
        first_line = f.readline()
        delim = '\t' if '\t' in first_line else ','
        for line in f:
            parts = line.strip().split(delim)
            if len(parts) < 2:
                continue
            a, b = parts[0].strip().upper(), parts[1].strip().upper()
            if a in gene_set and b in gene_set:
                score = float(parts[2]) if len(parts) >= 3 else 1.0
                G.add_edge(a, b, weight=score)
                edges_added += 1

    print(f"[PPI] 图构建完成: {G.number_of_nodes()} 节点, {G.number_of_edges()} 边")

    n = len(genes)
    feats = np.zeros((n, 12), dtype=np.float32)

    # 1. Degree (log1p 标准化)
    print("[PPI] 计算度...")
    degrees = dict(G.degree())
    for gene, deg in degrees.items():
        if gene in gene_to_idx:
            feats[gene_to_idx[gene], 0] = np.log1p(deg)

    # 2. Weighted degree (log1p)
    print("[PPI] 计算加权度...")
    w_degrees = dict(G.degree(weight='weight'))
    for gene, wdeg in w_degrees.items():
        if gene in gene_to_idx:
            feats[gene_to_idx[gene], 1] = np.log1p(wdeg)

    # 3. Clustering coefficient
    print("[PPI] 计算聚类系数...")
    clustering = nx.clustering(G)
    for gene, cc in clustering.items():
        if gene in gene_to_idx:
            feats[gene_to_idx[gene], 2] = cc

    # 4-5. Betweenness centrality (采样 k=2000)
    print("[PPI] 计算中介中心性 (k=2000 采样)...")
    try:
        bc = nx.betweenness_centrality(G, k=min(2000, G.number_of_nodes()), seed=RANDOM_SEED)
        for gene, val in bc.items():
            if gene in gene_to_idx:
                feats[gene_to_idx[gene], 3] = val
    except Exception as e:
        print(f"[WARN] 中介中心性计算失败: {e}")

    # 6. Closeness centrality (采样 k=500, 大规模图优化)
    print("[PPI] 计算接近中心性 (k=500 采样)...")
    try:
        sample_nodes = np.random.choice(list(G.nodes()), size=min(500, G.number_of_nodes()), replace=False)
        cc_dict = {}
        for idx, node in enumerate(sample_nodes):
            try:
                lengths = nx.single_source_shortest_path_length(G, node)
                total = sum(lengths.values())
                if total > 0:
                    cc_dict[node] = (len(lengths) - 1) / total
            except Exception:
                pass
            if (idx + 1) % 250 == 0:
                print(f"  [PPI] 接近中心性进度: {idx + 1}/{len(sample_nodes)}")
        for gene, val in cc_dict.items():
            if gene in gene_to_idx:
                feats[gene_to_idx[gene], 4] = val
    except Exception as e:
        print(f"[WARN] 接近中心性计算失败: {e}")

    # 7. Square clustering (跳过 - 大图太慢, 用聚类系数替代)
    print("[PPI] 跳过方形聚类系数 (用 triangles 替代)...")
    try:
        triangles = nx.triangles(G)
        for gene, val in triangles.items():
            if gene in gene_to_idx:
                feats[gene_to_idx[gene], 5] = np.log1p(val)
    except Exception as e:
        print(f"[WARN] triangles 计算失败: {e}")

    # 8-12. Degree statistics (neighbor stats)
    print("[PPI] 计算邻居统计特征...")
    for gene in genes:
        if gene in G:
            neighbors = list(G.neighbors(gene))
            if neighbors:
                n_degs = [G.degree(n) for n in neighbors]
                idx = gene_to_idx[gene]
                feats[idx, 6] = np.mean(n_degs)
                feats[idx, 7] = np.std(n_degs) if len(n_degs) > 1 else 0
                feats[idx, 8] = np.max(n_degs)
                feats[idx, 9] = np.min(n_degs)

                n_weights = [G.degree(n, weight='weight') for n in neighbors]
                feats[idx, 10] = np.mean(n_weights)
                feats[idx, 11] = np.std(n_weights) if len(n_weights) > 1 else 0

    node_to_gene = {v: k for k, v in gene_to_idx.items()}
    stats = {
        'degree': (feats[:, 0].min(), feats[:, 0].max(), feats[:, 0].mean()),
        'weighted_degree': (feats[:, 1].min(), feats[:, 1].max(), feats[:, 1].mean()),
        'clustering': (feats[:, 2].min(), feats[:, 2].max(), feats[:, 2].mean()),
        'betweenness': (feats[:, 3].min(), feats[:, 3].max(), feats[:, 3].mean()),
        'closeness': (feats[:, 4].min(), feats[:, 4].max(), feats[:, 4].mean()),
    }
    print(f"[PPI] 拓扑特征统计: {json.dumps(stats, indent=2, default=str)}")

    np.save(cache_file, feats)
    print(f"[PPI] 拓扑特征已保存: shape={feats.shape}")
    return feats


# ============================================================
# 3. GTEx eQTL 特征 (从 eGenes 文件)
# ============================================================
def load_egenes(filepath, tissue_label):
    """
    加载 GTEx eGenes 文件，提取每个基因的:
      - pval_beta: beta-approximated permutation p-value
      - slope: 效应大小
      - qval: Storey q-value
      - num_var: 检测的变异数

    返回 dict: {gene_symbol: [features]}
    """
    print(f"[eQTL] 加载 {tissue_label} eGenes: {filepath}")
    egenes = {}

    with open(filepath, 'r', encoding='utf-8') as f:
        header = next(f).strip().split('\t')
        name_idx = header.index('gene_name')
        pval_idx = header.index('pval_beta')
        slope_idx = header.index('slope')
        qval_idx = header.index('qval')
        num_var_idx = header.index('num_var') if 'num_var' in header else None

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < max(name_idx, pval_idx, slope_idx, qval_idx) + 1:
                continue
            gene_name = parts[name_idx].strip().upper()
            if not gene_name:
                continue

            try:
                pval = float(parts[pval_idx])
                slope = float(parts[slope_idx])
                qval = float(parts[qval_idx])
                num_var = float(parts[num_var_idx]) if num_var_idx else 1.0

                nlp = -np.log10(max(pval, 1e-300))
                abs_slope = abs(slope)
                nq = -np.log10(max(qval, 1e-300))

                egenes[gene_name] = [nlp, abs_slope, nq, num_var, 1.0]
            except (ValueError, IndexError):
                continue

    print(f"[eQTL] {tissue_label}: {len(egenes)} 个 eGenes")
    return egenes


def compute_eqtl_features(genes, brain_egenes, blood_egenes):
    """为每个基因构建 eQTL 特征向量"""
    cache_file = CACHE_DIR / "eqtl_features.npy"
    if cache_file.exists():
        print(f"[CACHE] 加载缓存的 eQTL 特征")
        return np.load(cache_file)

    n = len(genes)
    gene_to_idx = {g: i for i, g in enumerate(genes)}

    feats = np.zeros((n, 10), dtype=np.float32)

    brain_count = 0
    blood_count = 0

    for gene, idx in gene_to_idx.items():
        if gene in brain_egenes:
            vals = brain_egenes[gene]
            feats[idx, 0] = vals[0]
            feats[idx, 1] = vals[1]
            feats[idx, 2] = vals[2]
            feats[idx, 3] = vals[3]
            feats[idx, 4] = vals[4]
            brain_count += 1

        if gene in blood_egenes:
            vals = blood_egenes[gene]
            feats[idx, 5] = vals[0]
            feats[idx, 6] = vals[1]
            feats[idx, 7] = vals[2]
            feats[idx, 8] = vals[3]
            feats[idx, 9] = vals[4]
            blood_count += 1

    print(f"[eQTL] Brain 覆盖: {brain_count}/{n} ({100*brain_count/n:.1f}%)")
    print(f"[eQTL] Blood 覆盖: {blood_count}/{n} ({100*blood_count/n:.1f}%)")
    print(f"[eQTL] 任一组织覆盖: {np.sum((feats[:, 4] + feats[:, 9]) > 0)}/{n}")

    np.save(cache_file, feats)
    return feats


# ============================================================
# 4. 蛋白质序列理化特征 (UniProt API)
# ============================================================
def fetch_uniprot_mapping(gene_symbols, batch_size=10):
    """
    使用 UniProt REST API GET 批量查询蛋白质信息
    小批量避免 URL 过长
    """
    cache_file = CACHE_DIR / "uniprot_mapping.pkl"
    if cache_file.exists():
        print(f"[CACHE] 加载缓存的 UniProt 映射")
        with open(cache_file, 'rb') as f:
            mapping = pickle.load(f)
        if len(mapping) > 100:
            print(f"[CACHE] 有效映射: {len(mapping)} 个基因")
            return mapping
        else:
            print(f"[CACHE] 缓存映射数过少 ({len(mapping)}), 重新获取")

    print(f"[UniProt] 批量映射 {len(gene_symbols)} 个基因符号 (batch_size={batch_size})...")
    mapping = {}

    base_url = "https://rest.uniprot.org/uniprotkb/search"

    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i:i + batch_size]
        query_parts = [f"gene_exact:{g}" for g in batch]
        query = " OR ".join(query_parts)

        params = {
            'query': f"({query}) AND reviewed:true AND organism_id:9606",
            'format': 'tsv',
            'fields': 'accession,gene_names,length,sequence',
            'size': str(batch_size * 3),
        }
        query_string = urllib.parse.urlencode(params, doseq=True, quote_via=urllib.parse.quote_plus)
        url = f"{base_url}?{query_string}"

        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=60) as resp:
                content = resp.read().decode('utf-8')
                lines = content.strip().split('\n')
                if len(lines) < 2:
                    if (i // batch_size + 1) % 100 == 0:
                        print(f"[UniProt] 进度: {min(i + batch_size, len(gene_symbols))}/{len(gene_symbols)}, "
                              f"已映射: {len(mapping)}")
                    time.sleep(0.15)
                    continue
                for line in lines[1:]:
                    parts = line.split('\t')
                    if len(parts) < 4:
                        continue
                    uniprot_id = parts[0]
                    gene_names = parts[1].split()
                    length = int(parts[2]) if parts[2].isdigit() else 0
                    sequence = parts[3]
                    for gn in gene_names:
                        gn_upper = gn.upper()
                        if gn_upper in gene_symbols_set:
                            mapping[gn_upper] = {
                                'uniprot_id': uniprot_id,
                                'length': length,
                                'sequence': sequence,
                            }

            if (i // batch_size + 1) % 100 == 0:
                print(f"[UniProt] 进度: {min(i + batch_size, len(gene_symbols))}/{len(gene_symbols)}, "
                      f"已映射: {len(mapping)}")
            time.sleep(0.15)

        except Exception as e:
            if (i // batch_size + 1) % 100 == 0:
                print(f"[UniProt] 进度: {min(i + batch_size, len(gene_symbols))}/{len(gene_symbols)}, "
                      f"已映射: {len(mapping)}, 批次 {i} 错误: {e}")
            time.sleep(0.5)

    print(f"[UniProt] 映射完成: {len(mapping)}/{len(gene_symbols)} 个基因")

    with open(cache_file, 'wb') as f:
        pickle.dump(mapping, f)
    return mapping


def compute_protein_features_basic(seq):
    """本地计算基本蛋白质理化性质 (不依赖 BioPython)"""
    if not seq:
        return [0] * 20

    aa_weights = {
        'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
        'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
        'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
        'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
    }

    aa_pka = {
        'D': 3.9, 'E': 4.3, 'C': 8.3, 'Y': 10.1, 'H': 6.0,
        'K': 10.5, 'R': 12.5,
    }

    aa_hydropathy = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    }

    seq = seq.upper()
    length = len(seq)
    aa_counts = Counter(seq)

    # 1. 长度
    f_len = np.log1p(length)

    # 2. 分子量 (近似: 减去 (n-1)*18 水分子)
    mw = sum(aa_weights.get(aa, 0) for aa in seq) - (length - 1) * 18.015
    f_mw = mw / 1000.0

    # 3. 等电点 (简化 Henderson-Hasselbalch)
    n_asp = aa_counts.get('D', 0)
    n_glu = aa_counts.get('E', 0)
    n_cys = aa_counts.get('C', 0)
    n_tyr = aa_counts.get('Y', 0)
    n_his = aa_counts.get('H', 0)
    n_lys = aa_counts.get('K', 0)
    n_arg = aa_counts.get('R', 0)

    def _charge(ph):
        neg = n_asp * 10**(ph - aa_pka['D']) / (1 + 10**(ph - aa_pka['D'])) + \
              n_glu * 10**(ph - aa_pka['E']) / (1 + 10**(ph - aa_pka['E'])) + \
              n_cys * 10**(ph - aa_pka['C']) / (1 + 10**(ph - aa_pka['C'])) + \
              n_tyr * 10**(ph - aa_pka['Y']) / (1 + 10**(ph - aa_pka['Y']))
        pos = n_his * 10**(aa_pka['H'] - ph) / (1 + 10**(aa_pka['H'] - ph)) + \
              n_lys * 10**(aa_pka['K'] - ph) / (1 + 10**(aa_pka['K'] - ph)) + \
              n_arg * 10**(aa_pka['R'] - ph) / (1 + 10**(aa_pka['R'] - ph))
        return pos - neg

    lo, hi = 0.0, 14.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if _charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    f_pi = (lo + hi) / 2

    # 4. 疏水性 (GRAVY)
    gravy = sum(aa_hydropathy.get(aa, 0) for aa in seq) / length if length > 0 else 0
    f_gravy = gravy

    # 5. 芳香性
    aromatic = (aa_counts.get('F', 0) + aa_counts.get('Y', 0) + aa_counts.get('W', 0)) / length if length > 0 else 0

    # 6-25. 氨基酸组成 (20种)
    aa_list = 'ACDEFGHIKLMNPQRSTVWY'
    aa_comp = [aa_counts.get(aa, 0) / length if length > 0 else 0 for aa in aa_list]

    # 6. 不稳定指数 (简化)
    dipeptide_stability = {
        'W': 1.0, 'C': 1.0, 'M': 1.0, 'H': 1.0, 'Y': 1.0,
        'F': 1.0, 'Q': 1.0, 'N': 1.0, 'I': 1.0, 'V': 1.0,
        'P': 1.0, 'T': 1.0, 'K': 1.0, 'E': 1.0, 'D': 1.0,
        'A': 1.0, 'S': 1.0, 'R': 1.0, 'G': 1.0, 'L': 1.0,
    }
    instability = 0
    for j in range(length - 1):
        dipep = seq[j:j + 2]
        instability += dipeptide_stability.get(dipep[1], 1.0)
    f_instability = instability / max(length - 1, 1)

    return [f_len, f_mw, f_pi, f_gravy, aromatic, f_instability] + aa_comp


def compute_protein_features_gene_list(genes):
    """为基因列表计算蛋白质特征"""
    cache_file = CACHE_DIR / "protein_features.npy"
    if cache_file.exists():
        cached = np.load(cache_file)
        if not np.allclose(cached, 0):
            print(f"[CACHE] 加载缓存的蛋白质特征")
            return cached
        else:
            print(f"[CACHE] 蛋白质缓存无效(全零), 重新获取")

    global gene_symbols_set
    gene_symbols_set = set(genes)

    n = len(genes)
    feats = np.zeros((n, 26), dtype=np.float32)

    mapping = fetch_uniprot_mapping(genes)

    gene_to_idx = {g: i for i, g in enumerate(genes)}
    mapped = 0
    for gene, info in mapping.items():
        if gene in gene_to_idx:
            idx = gene_to_idx[gene]
            seq = info['sequence']
            seq_feats = compute_protein_features_basic(seq)
            feats[idx] = seq_feats
            mapped += 1

    print(f"[Protein] 蛋白质特征覆盖: {mapped}/{n} ({100*mapped/n:.1f}%)")

    np.save(cache_file, feats)
    return feats


# ============================================================
# 5. 特征组装与保存
# ============================================================
def normalize_features(feats_block, name):
    """对特征块做 Z-score 标准化 (处理 NaN/Inf)"""
    feats = feats_block.copy()
    for j in range(feats.shape[1]):
        col = feats[:, j]
        finite_mask = np.isfinite(col)
        if finite_mask.sum() < 2:
            continue
        mean = col[finite_mask].mean()
        std = col[finite_mask].std()
        if std < 1e-8:
            feats[:, j] = 0
        else:
            feats[:, j] = np.where(finite_mask, (col - mean) / std, 0)
    print(f"[NORM] {name}: shape={feats.shape}")
    return feats


def assemble_and_save():
    """主流程: 组装所有特征并保存"""
    print("=" * 60)
    print("基因特征增强 - 开始")
    print("=" * 60)

    genes = load_gene_list()
    disease_genes, drug_targets = load_disease_and_drug_genes()
    gene_to_idx = {g: i for i, g in enumerate(genes)}

    sapbert_emb, sapbert_dim = load_sapbert_embeddings(genes)

    # 2. PPI 拓扑特征
    t0 = time.time()
    ppi_feats = compute_ppi_topology_features(genes)
    ppi_feats = normalize_features(ppi_feats, "PPI topology")
    print(f"[TIME] PPI 拓扑特征: {time.time() - t0:.1f}s")

    # 3. eQTL 特征
    t0 = time.time()
    brain_egenes = load_egenes(BRAIN_EGENES, "Brain Cortex")
    blood_egenes = load_egenes(BLOOD_EGENES, "Whole Blood")
    eqtl_feats = compute_eqtl_features(genes, brain_egenes, blood_egenes)
    eqtl_feats = normalize_features(eqtl_feats, "eQTL")
    print(f"[TIME] eQTL 特征: {time.time() - t0:.1f}s")

    # 4. 蛋白质特征 (需要网络, 可能较慢)
    t0 = time.time()
    try:
        protein_feats = compute_protein_features_gene_list(genes)
        protein_feats = normalize_features(protein_feats, "Protein")
        print(f"[TIME] 蛋白质特征: {time.time() - t0:.1f}s")
        use_protein = True
    except Exception as e:
        print(f"[WARN] 蛋白质特征获取失败, 跳过: {e}")
        protein_feats = np.zeros((len(genes), 0), dtype=np.float32)
        use_protein = False

    # 5. 拼接所有特征
    all_feats_list = [sapbert_emb, ppi_feats, eqtl_feats]
    dim_names = [
        f"SapBERT_{sapbert_dim}",
        f"PPI_{ppi_feats.shape[1]}",
        f"eQTL_{eqtl_feats.shape[1]}",
    ]
    if use_protein:
        all_feats_list.append(protein_feats)
        dim_names.append(f"Protein_{protein_feats.shape[1]}")

    enhanced_matrix = np.hstack(all_feats_list)
    enhanced_matrix = np.nan_to_num(enhanced_matrix, nan=0.0, posinf=0.0, neginf=0.0)

    print(f"\n{'='*60}")
    print(f"增强特征矩阵: {enhanced_matrix.shape}")
    print(f"  原始 SapBERT: {sapbert_dim} 维")
    print(f"  + PPI 拓扑:   {ppi_feats.shape[1]} 维")
    print(f"  + eQTL:       {eqtl_feats.shape[1]} 维")
    if use_protein:
        print(f"  + 蛋白质:     {protein_feats.shape[1]} 维")
    print(f"  = 总计:       {enhanced_matrix.shape[1]} 维")

    # 保存为 CSV
    print(f"\n[SAVE] 保存增强特征 CSV: {OUTPUT_CSV}")
    feat_cols = [f"feat_{i:04d}" for i in range(enhanced_matrix.shape[1])]
    df_out = pd.DataFrame(enhanced_matrix, columns=feat_cols)
    df_out.insert(0, 'gene_symbol', genes)
    df_out.to_csv(OUTPUT_CSV, index=False, float_format='%.8g')
    print(f"[SAVE] CSV 已保存: {OUTPUT_CSV}")

    # 保存为 NPY
    np.save(OUTPUT_NPY, enhanced_matrix)
    print(f"[SAVE] NPY 已保存: {OUTPUT_NPY}")

    # 保存维度映射
    dim_map = {
        'total_dim': int(enhanced_matrix.shape[1]),
        'n_genes': int(enhanced_matrix.shape[0]),
        'blocks': dim_names,
        'sapbert_dim': sapbert_dim,
        'ppi_dim': int(ppi_feats.shape[1]),
        'eQTL_dim': int(eqtl_feats.shape[1]),
        'protein_dim': int(protein_feats.shape[1]) if use_protein else 0,
        'ppi_features': [
            'log1p_degree', 'log1p_weighted_degree', 'clustering_coeff',
            'betweenness_k2000', 'closeness_k500', 'log1p_triangles',
            'neighbor_mean_degree', 'neighbor_std_degree',
            'neighbor_max_degree', 'neighbor_min_degree',
            'neighbor_mean_weight', 'neighbor_std_weight',
        ],
        'eQTL_features': [
            'brain_nlog10_pval_beta', 'brain_abs_slope', 'brain_nlog10_qval',
            'brain_num_var', 'brain_is_eGene',
            'blood_nlog10_pval_beta', 'blood_abs_slope', 'blood_nlog10_qval',
            'blood_num_var', 'blood_is_eGene',
        ],
        'protein_features': [
            'log1p_length', 'mw_kda', 'pi', 'gravy', 'aromaticity',
            'instability_index',
            'aa_A', 'aa_C', 'aa_D', 'aa_E', 'aa_F',
            'aa_G', 'aa_H', 'aa_I', 'aa_K', 'aa_L',
            'aa_M', 'aa_N', 'aa_P', 'aa_Q', 'aa_R',
            'aa_S', 'aa_T', 'aa_V', 'aa_W', 'aa_Y',
        ] if use_protein else [],
    }
    with open(OUTPUT_DIM_MAP, 'w', encoding='utf-8') as f:
        json.dump(dim_map, f, indent=2, ensure_ascii=False)
    print(f"[SAVE] 维度映射: {OUTPUT_DIM_MAP}")

    # 防泄漏检查
    print(f"\n[LEAK] 防泄漏检查...")
    disease_set = disease_genes
    drug_set = drug_targets

    disease_idx = [gene_to_idx[g] for g in disease_set if g in gene_to_idx]
    drug_idx = [gene_to_idx[g] for g in drug_set if g in gene_to_idx]
    other_idx = [i for i in range(len(genes))
                 if i not in disease_idx and i not in drug_idx]

    print(f"  疾病基因: {len(disease_idx)}, 药物靶点: {len(drug_idx)}, "
          f"其他: {len(other_idx)}")

    # 检查 eQTL 特征是否与疾病状态独立
    if eqtl_feats.shape[1] > 0:
        for j in range(min(4, eqtl_feats.shape[1])):
            d_mean = eqtl_feats[disease_idx, j].mean()
            o_mean = eqtl_feats[other_idx, j].mean()
            print(f"  eQTL feat_{j}: 疾病均值={d_mean:.4f}, 其他均值={o_mean:.4f}")

    if ppi_feats.shape[1] > 0:
        for j in range(min(3, ppi_feats.shape[1])):
            d_mean = ppi_feats[disease_idx, j].mean()
            o_mean = ppi_feats[other_idx, j].mean()
            print(f"  PPI feat_{j}: 疾病均值={d_mean:.4f}, 其他均值={o_mean:.4f}")

    print(f"\n{'='*60}")
    print("基因特征增强 - 完成!")
    print(f"输出文件:")
    print(f"  {OUTPUT_CSV}")
    print(f"  {OUTPUT_NPY}")
    print(f"  {OUTPUT_DIM_MAP}")
    print("=" * 60)


if __name__ == "__main__":
    assemble_and_save()