# -*- coding: utf-8 -*-
"""
通路节点特征计算: 基因嵌入均值法
=====================================================
方法: 对每个通路取其所有关联基因的嵌入均值作为通路初始特征向量

输入:
  - pathway_nodes.csv (含 pathway_name, n_genes, 由 extract_gene_pathway_edges.py 生成)
  - subgraph_embeddings.csv (基因嵌入, 来自 SapBERT/ProtBERT)
  - gene_pathway_edges.txt (基因-通路关联, 由 extract_gene_pathway_edges.py 生成)

输出:
  - pathway_nodes.csv (更新, 追加嵌入向量列)
  - pathway_features.npy (NumPy 矩阵, shape=[n_pathways, embed_dim])
  - pathway_feature_stats.json (特征统计)

备选方案:
  - 方法二: Gonto2Vec / BioVec 本体嵌入 (较复杂)
  - 方法三: 通路内基因FPKM均值 (无预训练嵌入时)

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
PATHWAY_NODES_FILE = OUTPUT_DIR / "pathway_nodes.csv"
GENE_PATHWAY_EDGES = OUTPUT_DIR / "gene_pathway_edges.txt"

OUTPUT_PATHWAY_NODES = OUTPUT_DIR / "pathway_nodes.csv"
OUTPUT_FEATURES_NPY = OUTPUT_DIR / "pathway_features.npy"
OUTPUT_FEATURE_NAMES = OUTPUT_DIR / "pathway_feature_names.txt"
OUTPUT_STATS = OUTPUT_DIR / "pathway_feature_stats.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_gene_embeddings(embedding_path):
    """
    加载基因嵌入矩阵

    支持的格式:
      1. CSV: 首列=基因名, 其余列=嵌入向量 (有或没有表头)
      2. NPY: 纯数值矩阵 (需额外提供基因名文件)
    """
    print(f"[LOAD] 加载基因嵌入: {embedding_path}")

    if not embedding_path.exists():
        print(f"[ERROR] 嵌入文件不存在: {embedding_path}")
        print("[FALLBACK] 将使用随机初始化的通路特征")
        return None, None, 0

    try:
        df = pd.read_csv(embedding_path, sep=None, engine='python')
    except Exception:
        try:
            df = pd.read_csv(embedding_path, sep='\t')
        except Exception as e:
            print(f"[ERROR] 无法读取嵌入文件: {e}")
            return None, None, 0

    n_rows, n_cols = df.shape
    print(f"[INFO] 嵌入文件: {n_rows} 行 × {n_cols} 列")

    first_col_vals = df.iloc[:, 0]
    numeric_count = pd.to_numeric(first_col_vals, errors='coerce').notna().sum()

    if numeric_count > n_rows * 0.8:
        print("[INFO] 检测到纯数值矩阵 (无基因名列)")
        return None, None, 0

    if not isinstance(first_col_vals.iloc[0], str) or len(first_col_vals.iloc[0]) > 50:
        first_val = str(first_col_vals.iloc[0])
        has_header = any(c.isalpha() for c in first_val) and len(first_val) < 20
        if has_header:
            print(f"[INFO] 第一行可能是表头: '{first_val[:50]}'")
            df.columns = [first_val] + list(df.columns[1:])
            df = df.iloc[1:].reset_index(drop=True)

    gene_col = df.columns[0]
    embed_cols = df.columns[1:]

    gene_embeddings = {}
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
    print(f"[OK] 加载 {len(gene_embeddings)} 个基因的嵌入 (维度={embed_dim}), 跳过 {skipped} 行")

    return gene_embeddings, embed_dim


def load_gene_pathway_edges(edges_path):
    """加载基因-通路关联边"""
    print(f"[LOAD] 加载基因-通路关联: {edges_path}")

    pathway_genes = defaultdict(list)
    gene_pathways = defaultdict(list)

    if not edges_path.exists():
        print(f"[ERROR] 基因-通路关联文件不存在: {edges_path}")
        return pathway_genes, gene_pathways

    with open(edges_path, 'r', encoding='utf-8') as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            gene = parts[0].strip().upper()
            pathway = parts[1].strip()
            pathway_genes[pathway].append(gene)
            gene_pathways[gene].append(pathway)

    print(f"[OK] 加载 {len(pathway_genes)} 个通路, {sum(len(v) for v in pathway_genes.values())} 条关联")
    return pathway_genes, gene_pathways


def compute_pathway_features(pathway_genes, gene_embeddings, embed_dim):
    """
    计算通路特征: 通路内所有基因嵌入的均值

    策略:
      1. 存在嵌入的基因 → 取嵌入
      2. 不存在嵌入的基因 → 跳过
      3. 通路内无任何嵌入基因 → 随机单位向量 (避免零向量)

    返回:
      features:  (n_pathways, embed_dim) 浮点数矩阵
      names:     pathway 名称列表
      coverage:  {pathway: (n_with_embed, n_total)} 嵌入覆盖率
    """
    np.random.seed(42)

    if embed_dim == 0:
        print("[ERROR] 嵌入维度为0")
        return None, [], {}

    features = []
    names = []
    coverage = {}

    for pathway_name, genes in sorted(pathway_genes.items()):
        gene_vecs = []
        for gene in genes:
            if gene in gene_embeddings:
                gene_vecs.append(gene_embeddings[gene])

        n_total = len(genes)
        n_with = len(gene_vecs)
        coverage[pathway_name] = (n_with, n_total)

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
        names.append(pathway_name)

    features = np.array(features, dtype=np.float32)

    n_with_any = sum(1 for c in coverage.values() if c[0] > 0)
    print(f"\n[STATS] 通路特征计算完成:")
    print(f"  通路总数:     {len(features)}")
    print(f"  有嵌入基因:   {n_with_any} ({100*n_with_any/len(features):.1f}%)")
    print(f"  纯随机:       {len(features) - n_with_any}")
    print(f"  特征维度:     {features.shape[1]}")
    print(f"  平均覆盖率:   {np.mean([c[0]/c[1] for c in coverage.values()]):.1%}")

    return features, names, coverage


def save_features(features, names, coverage, embed_dim):
    """保存所有输出文件"""
    np.save(OUTPUT_FEATURES_NPY, features)
    print(f"[SAVE] 特征矩阵: {OUTPUT_FEATURES_NPY} ({features.shape})")

    with open(OUTPUT_FEATURE_NAMES, 'w', encoding='utf-8') as f:
        for name in names:
            f.write(name + '\n')
    print(f"[SAVE] 通路名列表: {OUTPUT_FEATURE_NAMES} ({len(names)} 个)")

    stats = {
        'n_pathways': len(names),
        'embed_dim': embed_dim,
        'n_with_embeddings': sum(1 for c in coverage.values() if c[0] > 0),
        'n_random': sum(1 for c in coverage.values() if c[0] == 0),
        'mean_coverage': np.mean([c[0] / c[1] for c in coverage.values()]),
        'method': 'gene_embedding_mean',
        'embedding_source': str(EMBEDDING_FILE),
        'top10_coverage': sorted(coverage.items(), key=lambda x: x[1][0] / x[1][1], reverse=True)[:10],
        'bottom10_coverage': sorted(coverage.items(), key=lambda x: x[1][0] / x[1][1])[:10],
    }
    with open(OUTPUT_STATS, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    print(f"[SAVE] 统计信息: {OUTPUT_STATS}")


def update_pathway_nodes_csv(names, coverage):
    """更新 pathway_nodes.csv, 追加特征信息列"""
    if not OUTPUT_PATHWAY_NODES.exists():
        print(f"[WARN] pathway_nodes.csv 不存在, 创建新文件")

    rows = []
    for i, name in enumerate(names):
        n_with, n_total = coverage[name]
        coverage_pct = n_with / n_total if n_total > 0 else 0.0
        rows.append({
            'pathway_name': name,
            'n_genes': n_total,
            'n_genes_with_embedding': n_with,
            'coverage': f"{coverage_pct:.2%}",
            'source': 'gene_mean',
        })

    with open(OUTPUT_PATHWAY_NODES, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'pathway_name', 'n_genes', 'n_genes_with_embedding',
            'coverage', 'source'
        ])
        writer.writeheader()
        writer.writerows(rows)

    print(f"[SAVE] 通路节点CSV: {OUTPUT_PATHWAY_NODES} ({len(rows)} 行)")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("通路节点特征计算: 基因嵌入均值法")
    print("=" * 60)

    gene_embeddings, embed_dim = load_gene_embeddings(EMBEDDING_FILE)

    pathway_genes, gene_pathways = load_gene_pathway_edges(GENE_PATHWAY_EDGES)

    if not pathway_genes:
        print("[ERROR] 没有基因-通路关联数据")
        print("[HINT] 请先运行 extract_gene_pathway_edges.py")
        sys.exit(1)

    features, names, coverage = compute_pathway_features(
        pathway_genes, gene_embeddings, embed_dim
    )

    if features is None:
        sys.exit(1)

    save_features(features, names, coverage, embed_dim)
    update_pathway_nodes_csv(names, coverage)

    print(f"\n{'='*60}")
    print("通路节点特征计算完成!")
    print(f"输出: {OUTPUT_FEATURES_NPY}")
    print(f"输出: {OUTPUT_PATHWAY_NODES}")
    print(f"输出: {OUTPUT_FEATURE_NAMES}")
    print("=" * 60)